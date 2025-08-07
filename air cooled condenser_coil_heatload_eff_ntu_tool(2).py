# Streamlit application for air-cooled condenser design - ROBUST VERSION
import math
import streamlit as st
import pandas as pd
from CoolProp.CoolProp import PropsSI, get_global_param_string

st.set_page_config(layout="wide")
st.title("Air-Cooled Condenser Design Tool (Robust Version)")

# Geometry Inputs
st.sidebar.header("Geometry Inputs")
tube_od_mm = st.sidebar.number_input("Tube Outer Diameter (mm)", min_value=0.1, value=9.525)
tube_thickness_mm = st.sidebar.number_input("Tube Wall Thickness (mm)", min_value=0.01, value=0.35)
row_pitch_mm = st.sidebar.number_input("Row Pitch (mm)", min_value=0.1, value=25.4)
tube_pitch_mm = st.sidebar.number_input("Tube Pitch (mm)", min_value=0.1, value=25.4)
fpi = st.sidebar.number_input("Fins Per Inch", min_value=1, value=12)
fin_thickness_mm = st.sidebar.number_input("Fin Thickness (mm)", min_value=0.01, value=0.12)
fin_material = st.sidebar.selectbox("Fin Material", ["Aluminum", "Copper"])
face_width_m = st.sidebar.number_input("Coil Face Width (m)", min_value=0.1, value=1.0)
face_height_m = st.sidebar.number_input("Coil Face Height (m)", min_value=0.1, value=1.0)
num_rows = st.sidebar.number_input("Total Number of Rows", min_value=1, value=4)
free_area_percent = st.sidebar.slider("Free Flow Area %", 10, 100, 25)
num_feeds = st.sidebar.number_input("Number of Feeds", min_value=1, value=4)

# Refrigerant Inputs
st.sidebar.header("Refrigerant Inputs")
fluid_list = get_global_param_string("FluidsList").split(',')
refrigerants = sorted([f for f in fluid_list if f.startswith("R")])
fluid = st.sidebar.selectbox("Refrigerant", refrigerants, index=refrigerants.index("R134a") if "R134a" in refrigerants else 0)

T1 = st.sidebar.number_input("Inlet Superheat Temp (°C)", value=95.0)
T3 = st.sidebar.number_input("Outlet Subcooled Temp (°C)", value=52.0)
T_cond = st.sidebar.number_input("Condensing Temp (°C)", value=58.0)
m_dot = st.sidebar.number_input("Mass Flow Rate (kg/s)", min_value=0.01, value=0.6)
air_temp = st.sidebar.number_input("Air Inlet Temp (°C)", value=48.0)
airflow_cmh = st.sidebar.number_input("Air Flow (m³/hr)", min_value=1.0, value=10000.0)

# Derived Dimensions
tube_od_m = tube_od_mm / 1000
tube_thk_m = tube_thickness_mm / 1000
tube_id_m = max(0.001, tube_od_m - 2 * tube_thk_m)  # Ensure positive diameter
tube_pitch_m = tube_pitch_mm / 1000
row_pitch_true_m = math.sqrt((tube_pitch_mm / 2000) ** 2 + (row_pitch_mm / 1000) ** 2)
fins_per_m = fpi * 39.37
fin_thk_m = fin_thickness_mm / 1000
fin_k = 235 if fin_material == "Aluminum" else 400
frontal_area_m2 = face_width_m * face_height_m
net_free_area = max(0.001, frontal_area_m2 * (free_area_percent / 100))  # Prevent zero
airflow_m3s = max(0.001, airflow_cmh / 3600)  # Prevent zero
air_velocity_face = airflow_m3s / frontal_area_m2
air_velocity_fin = airflow_m3s / net_free_area

# Tube and Fin Counts
tubes_per_row = max(1, math.floor(face_width_m / tube_pitch_m))
tube_length_per_tube = face_height_m
total_tubes = tubes_per_row * num_rows
total_tube_length = total_tubes * tube_length_per_tube
A_tube_ext = total_tube_length * math.pi * tube_od_m

# Fin Surface Area with Correction
num_fins = max(1, math.floor(face_height_m / 0.0254 * fpi))
L_fin = min(row_pitch_true_m, tube_pitch_m) / 2
A_fin_raw = row_pitch_true_m * num_rows * face_width_m * num_fins
A_hole = (math.pi / 4 * tube_od_m ** 2) * total_tubes
A_fin = max(0, (A_fin_raw - A_hole) * 2)  # Prevent negative area
m = math.sqrt((2 * 40) / (fin_k * fin_thk_m)) if fin_k * fin_thk_m > 0 else 0
eta_fin = math.tanh(m * L_fin) / (m * L_fin) if m * L_fin > 0 else 0
A_air_total = A_tube_ext + A_fin * eta_fin
A_air_per_m = A_air_total / total_tube_length if total_tube_length > 0 else 0

# Air Properties
try:
    T_K = air_temp + 273.15
    mu = PropsSI("V", "T", T_K, "P", 101325, "Air")
    k_air = PropsSI("L", "T", T_K, "P", 101325, "Air")
    rho_air = PropsSI("D", "T", T_K, "P", 101325, "Air")
    cp_air = PropsSI("C", "T", T_K, "P", 101325, "Air") / 1000  # kJ/kg-K
    Pr = PropsSI("PRANDTL", "T", T_K, "P", 101325, "Air")
    Re = max(1, rho_air * air_velocity_fin * tube_od_m / mu)  # Prevent Re=0
    C = 0.27 if Re < 1000 else 0.021
    n = 0.63 if Re < 1000 else 0.84
    Nu = C * Re**n * Pr**0.36
    h_air = Nu * k_air / tube_od_m if tube_od_m > 0 else 0
except:
    st.error("Error calculating air properties. Check input temperatures.")
    st.stop()

# Refrigerant Properties
try:
    P_cond = PropsSI("P", "T", T_cond + 273.15, "Q", 0, fluid)
    h1 = PropsSI("H", "P", P_cond, "T", T1 + 273.15, fluid)
    h2 = PropsSI("H", "P", P_cond, "Q", 1, fluid)
    h3 = PropsSI("H", "P", P_cond, "Q", 0, fluid)
    h4 = PropsSI("H", "P", P_cond, "T", T3 + 273.15, fluid)
    cp_ref = PropsSI("C", "P", P_cond, "T", (T1 + T3)/2 + 273.15, fluid) / 1000  # kJ/kg-K
except:
    st.error("Error calculating refrigerant properties. Check refrigerant inputs.")
    st.stop()

Q_sens = m_dot * (h1 - h2) / 1000  # kW
Q_lat = m_dot * (h2 - h3) / 1000   # kW
Q_sub = m_dot * (h3 - h4) / 1000   # kW

# Refrigerant-side heat transfer coefficients (simplified)
h_ref_desuper = 1500  # W/m²-K for superheated vapor
h_ref_cond = 3000     # W/m²-K for condensation
h_ref_subcool = 2000  # W/m²-K for subcooled liquid

# Overall U calculations
def calculate_U(h_ref, h_air, A_ratio=1):
    """Calculate overall U including wall and fouling resistances"""
    try:
        R_total = (1/max(1, h_air)) + (1/max(1, h_ref)) + (0.0001)  # Added 0.0001 m²-K/W fouling
        return 1 / max(0.001, R_total)  # Prevent division by zero
    except:
        return 50.0  # Fallback U-value

U_desuper = calculate_U(h_ref_desuper, h_air)
U_cond = calculate_U(h_ref_cond, h_air)
U_subcool = calculate_U(h_ref_subcool, h_air)

# Effectiveness-NTU Calculations
C_air = max(0.001, airflow_m3s * rho_air * cp_air)  # kW/K (prevent zero)
zone_data = [
    ("Desuperheat", Q_sens, U_desuper, h_ref_desuper, T1),
    ("Condensing", Q_lat, U_cond, h_ref_cond, T_cond),
    ("Subcooling", Q_sub, U_subcool, h_ref_subcool, T3)
]

zone_outputs = []
air_t = air_temp
for label, Q_zone, U, h_ref, T_ref in zone_data:
    if Q_zone <= 0:
        zone_outputs.append((label, 0, 0, 0, 0, air_t))
        continue
        
    # Capacity rates
    if label == "Condensing":
        C_min = C_air
        C_star = 0
    else:
        C_r = m_dot * cp_ref
        C_min = min(C_air, C_r)
        C_star = C_min / max(0.001, max(C_air, C_r))  # Prevent division by zero
    
    # Effectiveness calculation with bounds checking
    delta_T = max(0.1, (T_ref - air_t))  # Ensure minimum ΔT of 0.1°C
    eps = min(0.999, Q_zone / (C_min * delta_T)) if C_min * delta_T > 0 else 0
    
    # NTU calculation with robust error handling
    try:
        if C_star == 0:  # Condensation
            NTU = -math.log(max(0.001, 1 - eps))  # Ensure positive argument for log
        else:  # Desuperheating/subcooling
            if eps >= 1:
                NTU = 100  # Large value for maximum heat transfer
            elif C_star == 1:
                NTU = eps/(1 - eps + 1e-9)  # Add small number to prevent division by zero
            else:
                numerator = (eps - 1)
                denominator = (eps * C_star - 1)
                if denominator == 0 or numerator/denominator <= 0:
                    NTU = 100
                else:
                    NTU = (1/(C_star - 1 + 1e-9)) * math.log(numerator/denominator)
    except:
        NTU = 100  # Fallback to large NTU if calculation fails
    
    # Area calculation with safeguards
    A_required = NTU * C_min * 1000 / max(1, U)  # Convert C_min from kW/K to W/K
    tube_len_zone = A_required / max(0.001, A_air_per_m)
    rows_zone = tube_len_zone / max(0.001, (tube_length_per_tube * tubes_per_row))
    
    # Air temperature update
    air_out = air_t + Q_zone / max(0.001, C_air)
    zone_outputs.append((label, Q_zone, A_required, tube_len_zone, rows_zone, air_out))
    air_t = air_out

# Outputs
st.subheader("Refrigerant Circuit Summary")
st.write(f"**Number of Circuits:** {max(1, int(tubes_per_row // num_feeds))}")
st.write(f"**Mass Flow Rate:** {m_dot:.3f} kg/s")
st.write(f"**Superheat Enthalpy (h1):** {h1/1000:.2f} kJ/kg")
st.write(f"**Saturated Vapor Enthalpy (h2):** {h2/1000:.2f} kJ/kg")
st.write(f"**Saturated Liquid Enthalpy (h3):** {h3/1000:.2f} kJ/kg")
st.write(f"**Subcooled Enthalpy (h4):** {h4/1000:.2f} kJ/kg")

st.subheader("Heat Load Summary")
st.write(f"**Desuperheating Load:** {Q_sens:.2f} kW")
st.write(f"**Condensing Load:** {Q_lat:.2f} kW")
st.write(f"**Subcooling Load:** {Q_sub:.2f} kW")
st.write(f"**Total Load:** {Q_sens + Q_lat + Q_sub:.2f} kW")

st.subheader("Air Side Summary")
st.write(f"**Air Flow Rate:** {airflow_m3s:.2f} m³/s")
st.write(f"**Air Capacity Rate:** {C_air:.3f} kW/K")
st.write(f"**Air-Side h:** {h_air:.2f} W/m²-K")
st.write(f"**Desuperheat U:** {U_desuper:.2f} W/m²-K")
st.write(f"**Condensing U:** {U_cond:.2f} W/m²-K")
st.write(f"**Subcooling U:** {U_subcool:.2f} W/m²-K")

st.subheader("Zone-by-Zone Results")
df = pd.DataFrame(zone_outputs, columns=["Zone", "Q (kW)", "Area (m²)", "Tube Length (m)", "Rows Required", "Air Temp Out (°C)"])
st.dataframe(df.style.format({
    "Q (kW)": "{:.2f}", 
    "Area (m²)": "{:.2f}", 
    "Tube Length (m)": "{:.2f}",
    "Rows Required": "{:.2f}",
    "Air Temp Out (°C)": "{:.2f}"
}))

st.subheader("Design Verification")
total_area = sum([z[2] for z in zone_outputs])
total_rows = sum([z[4] for z in zone_outputs])
st.write(f"**Total Required Area:** {total_area:.2f} m²")
st.write(f"**Total Required Rows:** {total_rows:.2f}")
st.write(f"**Provided Rows:** {num_rows}")
st.write(f"**Provided Area:** {A_air_total:.2f} m²")

if total_rows > num_rows * 1.1:
    st.error("Design Warning: Required rows exceed design rows by >10%")
elif total_area > A_air_total * 1.1:
    st.warning("Design Note: Required area exceeds provided area by >10%")
else:
    st.success("Design meets requirements")

# Add download button for results
csv = df.to_csv(index=False)
st.download_button(
    label="Download Results as CSV",
    data=csv,
    file_name="condenser_design_results.csv",
    mime="text/csv"
)