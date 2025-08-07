# Streamlit application for air-cooled condenser design - DEBUGGING VERSION
import math
import streamlit as st
import pandas as pd
from CoolProp.CoolProp import PropsSI, get_global_param_string

st.set_page_config(layout="wide")
st.title("Air-Cooled Condenser Design Tool (Debugging Version)")

# ====================== INPUT SECTION ======================
with st.sidebar:
    st.header("🛠️ Geometry Inputs")
    tube_od_mm = st.number_input("Tube Outer Diameter (mm)", min_value=0.1, value=9.525)
    tube_thickness_mm = st.number_input("Tube Wall Thickness (mm)", min_value=0.01, value=0.35)
    row_pitch_mm = st.number_input("Row Pitch (mm)", min_value=0.1, value=25.4)
    tube_pitch_mm = st.number_input("Tube Pitch (mm)", min_value=0.1, value=25.4)
    fpi = st.number_input("Fins Per Inch", min_value=1, value=12)
    fin_thickness_mm = st.number_input("Fin Thickness (mm)", min_value=0.01, value=0.12)
    fin_material = st.selectbox("Fin Material", ["Aluminum", "Copper"])
    face_width_m = st.number_input("Coil Face Width (m)", min_value=0.1, value=1.0)
    face_height_m = st.number_input("Coil Face Height (m)", min_value=0.1, value=1.0)
    num_rows = st.number_input("Total Number of Rows", min_value=1, value=4)
    free_area_percent = st.slider("Free Flow Area %", 10, 100, 25)
    num_feeds = st.number_input("Number of Feeds", min_value=1, value=4)

    st.header("❄️ Refrigerant Inputs")
    fluid_list = get_global_param_string("FluidsList").split(',')
    refrigerants = sorted([f for f in fluid_list if f.startswith("R")])
    fluid = st.selectbox("Refrigerant", refrigerants, index=refrigerants.index("R134a") if "R134a" in refrigerants else 0)
    T1 = st.number_input("Inlet Superheat Temp (°C)", value=95.0)
    T3 = st.number_input("Outlet Subcooled Temp (°C)", value=52.0)
    T_cond = st.number_input("Condensing Temp (°C)", value=58.0)
    m_dot = st.number_input("Mass Flow Rate (kg/s)", min_value=0.01, value=0.6)
    air_temp = st.number_input("Air Inlet Temp (°C)", value=48.0)
    airflow_cmh = st.number_input("Air Flow (m³/hr)", min_value=1.0, value=10000.0)

# ====================== CALCULATIONS ======================
st.header("📊 Calculated Parameters")

# Derived Dimensions
tube_od_m = tube_od_mm / 1000
tube_thk_m = tube_thickness_mm / 1000
tube_id_m = max(0.001, tube_od_m - 2 * tube_thk_m)
tube_pitch_m = tube_pitch_mm / 1000
row_pitch_true_m = math.sqrt((tube_pitch_mm/2000)**2 + (row_pitch_mm/1000)**2)
fins_per_m = fpi * 39.37
fin_thk_m = fin_thickness_mm / 1000
fin_k = 235 if fin_material == "Aluminum" else 400
frontal_area_m2 = face_width_m * face_height_m
net_free_area = max(0.001, frontal_area_m2 * (free_area_percent/100))
airflow_m3s = max(0.001, airflow_cmh / 3600)
air_velocity_face = airflow_m3s / frontal_area_m2
air_velocity_fin = airflow_m3s / net_free_area

# Display geometry outputs
with st.expander("📐 Geometry Calculations", expanded=True):
    cols = st.columns(4)
    cols[0].metric("Tube OD", f"{tube_od_m*1000:.2f} mm")
    cols[1].metric("Tube ID", f"{tube_id_m*1000:.2f} mm")
    cols[2].metric("True Row Pitch", f"{row_pitch_true_m*1000:.1f} mm")
    cols[3].metric("Fins/m", f"{fins_per_m:.0f}")
    
    cols = st.columns(3)
    cols[0].metric("Frontal Area", f"{frontal_area_m2:.3f} m²")
    cols[1].metric("Free Flow Area", f"{net_free_area:.3f} m²")
    cols[2].metric("Face Velocity", f"{air_velocity_face:.2f} m/s")
    cols[0].metric("Fin Velocity", f"{air_velocity_fin:.2f} m/s")
    cols[1].metric("Fin Efficiency", f"{eta_fin:.3f}" if 'eta_fin' in locals() else "N/A")

# Tube and Fin Counts
tubes_per_row = max(1, math.floor(face_width_m / tube_pitch_m))
tube_length_per_tube = face_height_m
total_tubes = tubes_per_row * num_rows
total_tube_length = total_tubes * tube_length_per_tube
A_tube_ext = total_tube_length * math.pi * tube_od_m

# Fin Surface Area
num_fins = max(1, math.floor(face_height_m / 0.0254 * fpi))
L_fin = min(row_pitch_true_m, tube_pitch_m) / 2
A_fin_raw = row_pitch_true_m * num_rows * face_width_m * num_fins
A_hole = (math.pi/4 * tube_od_m**2) * total_tubes
A_fin = max(0, (A_fin_raw - A_hole) * 2)
m = math.sqrt((2*40)/(fin_k*fin_thk_m)) if fin_k*fin_thk_m > 0 else 0
eta_fin = math.tanh(m*L_fin)/(m*L_fin) if m*L_fin > 0 else 0
A_air_total = A_tube_ext + A_fin*eta_fin
A_air_per_m = A_air_total/total_tube_length if total_tube_length > 0 else 0

# Display surface areas
with st.expander("📏 Surface Areas", expanded=True):
    cols = st.columns(3)
    cols[0].metric("Tubes per Row", f"{tubes_per_row}")
    cols[1].metric("Total Tube Length", f"{total_tube_length:.1f} m")
    cols[2].metric("Tube Outer Area", f"{A_tube_ext:.2f} m²")
    
    cols = st.columns(3)
    cols[0].metric("Fin Area (raw)", f"{A_fin_raw:.2f} m²")
    cols[1].metric("Fin Area (effective)", f"{A_fin*eta_fin:.2f} m²")
    cols[2].metric("Total Air-Side Area", f"{A_air_total:.2f} m²")
    cols[0].metric("Area per Meter", f"{A_air_per_m:.2f} m²/m")

# Air Properties
try:
    T_K = air_temp + 273.15
    mu = PropsSI("V", "T", T_K, "P", 101325, "Air")
    k_air = PropsSI("L", "T", T_K, "P", 101325, "Air")
    rho_air = PropsSI("D", "T", T_K, "P", 101325, "Air")
    cp_air = PropsSI("C", "T", T_K, "P", 101325, "Air")/1000  # kJ/kg-K
    Pr = PropsSI("PRANDTL", "T", T_K, "P", 101325, "Air")
    Re = max(1, rho_air*air_velocity_fin*tube_od_m/mu)
    C = 0.27 if Re < 1000 else 0.021
    n = 0.63 if Re < 1000 else 0.84
    Nu = C * Re**n * Pr**0.36
    h_air = Nu*k_air/tube_od_m if tube_od_m > 0 else 0
    C_air = max(0.001, airflow_m3s*rho_air*cp_air)  # kW/K
    
    with st.expander("🌬️ Air Properties", expanded=True):
        cols = st.columns(4)
        cols[0].metric("Density", f"{rho_air:.3f} kg/m³")
        cols[1].metric("Viscosity", f"{mu*1e6:.2f} μPa·s")
        cols[2].metric("Thermal Cond.", f"{k_air:.4f} W/m·K")
        cols[3].metric("Specific Heat", f"{cp_air:.3f} kJ/kg·K")
        
        cols = st.columns(4)
        cols[0].metric("Reynolds No.", f"{Re:.0f}")
        cols[1].metric("Prandtl No.", f"{Pr:.3f}")
        cols[2].metric("Nusselt No.", f"{Nu:.1f}")
        cols[3].metric("h (air side)", f"{h_air:.1f} W/m²K")
        
except Exception as e:
    st.error(f"❌ Error calculating air properties: {str(e)}")
    st.stop()

# Refrigerant Properties
try:
    P_cond = PropsSI("P", "T", T_cond + 273.15, "Q", 0, fluid)
    h1 = PropsSI("H", "P", P_cond, "T", T1 + 273.15, fluid)
    h2 = PropsSI("H", "P", P_cond, "Q", 1, fluid)
    h3 = PropsSI("H", "P", P_cond, "Q", 0, fluid)
    h4 = PropsSI("H", "P", P_cond, "T", T3 + 273.15, fluid)
    cp_ref = PropsSI("C", "P", P_cond, "T", (T1+T3)/2 + 273.15, fluid)/1000  # kJ/kg-K
    
    with st.expander("🧊 Refrigerant Properties", expanded=True):
        cols = st.columns(3)
        cols[0].metric("Condensing Pressure", f"{P_cond/1000:.2f} kPa")
        cols[1].metric("Superheat Enthalpy", f"{h1/1000:.2f} kJ/kg")
        cols[2].metric("Vapor Enthalpy", f"{h2/1000:.2f} kJ/kg")
        
        cols = st.columns(3)
        cols[0].metric("Liquid Enthalpy", f"{h3/1000:.2f} kJ/kg")
        cols[1].metric("Subcool Enthalpy", f"{h4/1000:.2f} kJ/kg")
        cols[2].metric("Avg Cp", f"{cp_ref:.3f} kJ/kg·K")
        
except Exception as e:
    st.error(f"❌ Error calculating refrigerant properties: {str(e)}")
    st.stop()

# Heat Loads
Q_sens = m_dot * (h1-h2)/1000  # kW
Q_lat = m_dot * (h2-h3)/1000   # kW
Q_sub = m_dot * (h3-h4)/1000   # kW

# Heat Transfer Coefficients
h_ref_desuper = 1500  # W/m²-K
h_ref_cond = 3000     # W/m²-K
h_ref_subcool = 2000  # W/m²-K

def calculate_U(h_ref, h_air):
    """Calculate overall heat transfer coefficient"""
    R_total = (1/max(1, h_air)) + (1/max(1, h_ref)) + 0.0001  # Added fouling
    return 1/max(0.001, R_total)

U_desuper = calculate_U(h_ref_desuper, h_air)
U_cond = calculate_U(h_ref_cond, h_air)
U_subcool = calculate_U(h_ref_subcool, h_air)

with st.expander("🔥 Heat Transfer Parameters", expanded=True):
    cols = st.columns(3)
    cols[0].metric("Desuperheat U", f"{U_desuper:.1f} W/m²K")
    cols[1].metric("Condensing U", f"{U_cond:.1f} W/m²K")
    cols[2].metric("Subcooling U", f"{U_subcool:.1f} W/m²K")
    
    cols = st.columns(3)
    cols[0].metric("Desuperheat Q", f"{Q_sens:.2f} kW")
    cols[1].metric("Condensing Q", f"{Q_lat:.2f} kW")
    cols[2].metric("Subcooling Q", f"{Q_sub:.2f} kW")

# ====================== ZONE CALCULATIONS ======================
zone_data = [
    ("Subcooling", Q_sub, U_subcool, h_ref_subcool, T3),
    ("Condensing", Q_lat, U_cond, h_ref_cond, T_cond),
    ("Desuperheat", Q_sens, U_desuper, h_ref_desuper, T1)
]

zone_outputs = []
diagnostics = []
air_t = air_temp  # Initial air temperature

for label, Q_zone, U, h_ref, T_ref in zone_data:
    if Q_zone <= 0:
        zone_outputs.append((label, 0, 0, 0, 0, air_t))
        diagnostics.append((label, 0, 0, 0, 0, 0, 0, 0, 0))
        continue
        
    # Capacity rates
    if label == "Condensing":
        C_min = C_air
        C_star = 0
    else:
        C_r = m_dot * cp_ref
        C_min = min(C_air, C_r)
        C_star = C_min / max(0.001, max(C_air, C_r))
    
    # Effectiveness calculation
    delta_T = max(0.1, abs(T_ref - air_t))
    eps = min(0.999, Q_zone/(C_min*delta_T)) if (C_min*delta_T) > 0 else 0
    
    # NTU calculation
    try:
        if C_star == 0:  # Condensation
            NTU = -math.log(max(0.001, 1-eps))
        else:
            if eps >= 1:
                NTU = 100
            elif C_star == 1:
                NTU = eps/(1-eps + 1e-9)
            else:
                numerator = (eps-1)
                denominator = (eps*C_star-1)
                if denominator == 0 or numerator/denominator <= 0:
                    NTU = 100
                else:
                    NTU = (1/(C_star-1 + 1e-9)) * math.log(numerator/denominator)
    except:
        NTU = 100  # Fallback value
    
    # Area requirements
    A_required = NTU * C_min * 1000 / max(1, U)
    tube_len_zone = A_required / max(0.001, A_air_per_m)
    rows_zone = tube_len_zone / max(0.001, (tube_length_per_tube*tubes_per_row))
    
    # Update air temperature
    air_out = air_t + Q_zone/max(0.001, C_air)
    
    zone_outputs.append((label, Q_zone, A_required, tube_len_zone, rows_zone, air_out))
    diagnostics.append((label, C_min, C_star, delta_T, eps, NTU, U, A_air_per_m, tube_length_per_tube*tubes_per_row))
    air_t = air_out

# ====================== OUTPUT SECTION ======================
st.header("📋 Design Results")

# Main summary
cols = st.columns(4)
cols[0].metric("Total Heat Load", f"{Q_sens + Q_lat + Q_sub:.2f} kW", delta_color="off")
cols[1].metric("Required Area", f"{sum([z[2] for z in zone_outputs]):.2f} m²", 
               delta=f"{(sum([z[2] for z in zone_outputs])-A_air_total):.2f} vs provided")
cols[2].metric("Required Rows", f"{sum([z[4] for z in zone_outputs]):.1f}",
              delta=f"{(sum([z[4] for z in zone_outputs])-num_rows:.1f} vs provided")
cols[3].metric("Air Out Temp", f"{air_t:.1f} °C", delta_color="off")

# Detailed zone results
st.subheader("🔍 Zone-by-Zone Results")
df = pd.DataFrame(zone_outputs, columns=[
    "Zone", "Q (kW)", "Area (m²)", "Tube Length (m)", "Rows Needed", "Air Out Temp (°C)"
])
st.dataframe(df.style.format({
    "Q (kW)": "{:.2f}",
    "Area (m²)": "{:.2f}",
    "Tube Length (m)": "{:.2f}",
    "Rows Needed": "{:.2f}",
    "Air Out Temp (°C)": "{:.2f}"
}))

# Diagnostic data
st.subheader("🛠️ Diagnostic Parameters")
diag_df = pd.DataFrame(diagnostics, columns=[
    "Zone", "C_min (kW/K)", "C*", "ΔT (°C)", "ε", "NTU", "U (W/m²K)", 
    "A/m (m²/m)", "Row Area (m²/row)"
])
st.dataframe(diag_df.style.format({
    "C_min (kW/K)": "{:.3f}",
    "C*": "{:.3f}",
    "ΔT (°C)": "{:.1f}",
    "ε": "{:.3f}",
    "NTU": "{:.2f}",
    "U (W/m²K)": "{:.1f}",
    "A/m (m²/m)": "{:.2f}",
    "Row Area (m²/row)": "{:.2f}"
}))

# Air flow sequence
st.subheader("🌪️ Air Flow Sequence")
st.markdown("""
1. **Subcooling Section**: Ambient air first cools the subcooled liquid (coldest refrigerant)
2. **Condensing Section**: Air then flows through the condensing section  
3. **Desuperheating Section**: Finally passes through the desuperheating section (hottest refrigerant)

*Note: Refrigerant flows in the opposite direction for counterflow effectiveness*
""")

# Design verification
total_area = sum([z[2] for z in zone_outputs])
total_rows = sum([z[4] for z in zone_outputs])

if total_rows > num_rows * 1.1:
    st.error(f"⚠️ **Design Issue**: Requires {total_rows:.1f} rows but only {num_rows} provided")
    st.progress(min(1.0, num_rows/total_rows), text="Row Capacity")
    
    st.warning("""
    **Possible causes for high row requirement:**
    - Low air-side heat transfer coefficient (current: {h_air:.1f} W/m²K)
    - Small temperature differences (ΔT: {min([d[3] for d in diagnostics]):.1f}°C min)
    - High heat loads relative to airflow
    - Low fin efficiency (current: {eta_fin:.2f})
    """)
    
elif total_area > A_air_total * 1.1:
    st.warning(f"⚠️ **Design Note**: Requires {total_area:.1f} m² but only {A_air_total:.1f} provided")
else:
    st.success("✅ Design meets requirements")

# Download buttons
csv1 = df.to_csv(index=False)
csv2 = diag_df.to_csv(index=False)

cols = st.columns(2)
cols[0].download_button(
    label="Download Zone Results",
    data=csv1,
    file_name="condenser_zone_results.csv",
    mime="text/csv"
)
cols[1].download_button(
    label="Download Diagnostics",
    data=csv2,
    file_name="condenser_diagnostics.csv",
    mime="text/csv"
)