
import math
import streamlit as st
import pandas as pd
from CoolProp.CoolProp import PropsSI, get_global_param_string

st.title("Air-Cooled Condenser - NTU & Geometry Breakdown")

# Refrigerant inputs
st.sidebar.header("Refrigerant Inputs")
fluid_list = get_global_param_string("FluidsList").split(',')
refrigerants = sorted([f for f in fluid_list if f.startswith("R")])
fluid = st.sidebar.selectbox("Refrigerant", refrigerants, index=refrigerants.index("R134a") if "R134a" in refrigerants else 0)
m_dot = st.sidebar.number_input("Refrigerant Mass Flow Rate (kg/s)", value=0.6)
T_in_superheat = st.sidebar.number_input("Inlet Superheat Temp (°C)", value=95.0)
T_cond = st.sidebar.number_input("Condensation Temp (°C)", value=45.0)
T_out_subcool = st.sidebar.number_input("Subcooled Liquid Temp (°C)", value=50.0)

# Air inputs
st.sidebar.header("Air Inputs")
T_air_in = st.sidebar.number_input("Air Inlet Temp (°C)", value=35.0)

free_area_percent = st.slider("Free Flow Area (%)", min_value=10, max_value=100, value=25)

airflow_cmh = st.sidebar.number_input("Air Flow Rate (m³/h)", value=10000)

# Coil geometry
st.sidebar.header("Coil Geometry")
tube_od_mm = st.sidebar.number_input("Tube Outer Diameter (mm)", value=9.525)
tube_thickness_mm = st.sidebar.number_input("Tube Wall Thickness (mm)", value=0.35)
tube_pitch_mm = st.sidebar.number_input("Tube Pitch (horizontal) (mm)", value=25.4)
row_pitch_mm = st.sidebar.number_input("Row Pitch (depth) (mm)", value=25.4)
fin_thickness_mm = st.sidebar.number_input("Fin Thickness (mm)", value=0.12)
fpi = st.sidebar.number_input("Fins Per Inch (FPI)", value=12)
face_width_m = st.sidebar.number_input("Coil Face Width (m)", value=1.0)
face_height_m = st.sidebar.number_input("Coil Face Height (m)", value=1.0)
num_rows = st.sidebar.number_input("Number of Rows", value=4)

# U-values
st.sidebar.header("U Values (W/m²·K)")
U_subcool = st.sidebar.number_input("U - Subcooling", value=40.0)
U_cond = st.sidebar.number_input("U - Condensing", value=60.0)
U_desuper = st.sidebar.number_input("U - Desuperheating", value=40.0)

# Geometry calculations
tube_od = tube_od_mm / 1000
tube_pitch = tube_pitch_mm / 1000
row_pitch = row_pitch_mm / 1000
fin_thickness = fin_thickness_mm / 1000
fins_per_m_vertical = fpi / 0.0254
fin_depth = row_pitch * num_rows

tubes_per_row = int(face_width_m / tube_pitch)
tube_length_per_tube = face_height_m
total_tubes = tubes_per_row * num_rows
total_tube_length = tube_length_per_tube * total_tubes
tube_ext_area = total_tube_length * math.pi * tube_od

number_of_fins = int(face_height_m * fins_per_m_vertical)
area_per_fin = 2 * face_width_m * fin_depth
total_fin_area = number_of_fins * area_per_fin

total_air_side_area = tube_ext_area + total_fin_area

# Refrigerant enthalpies
P_cond = PropsSI("P", "T", T_cond + 273.15, "Q", 0, fluid)
h1 = PropsSI("H", "T", T_in_superheat + 273.15, "P", P_cond, fluid)
h2 = PropsSI("H", "T", T_cond + 273.15, "Q", 1, fluid)
h3 = PropsSI("H", "T", T_cond + 273.15, "Q", 0, fluid)
h4 = PropsSI("H", "T", T_out_subcool + 273.15, "P", P_cond, fluid)

# Air properties
T_air_K = T_air_in + 273.15
cp_air = PropsSI("C", "T", T_air_K, "P", 101325, "Air")
rho_air = PropsSI("D", "T", T_air_K, "P", 101325, "Air")
m_dot_air = airflow_cmh / 3600 * rho_air

def zone_calc(name, q_kw, T_ref, T_air_in, U):
    Q = q_kw * 1000
    deltaT = (T_ref + 273.15) - (T_air_in + 273.15)
    C_air = m_dot_air * cp_air
    eps = Q / (C_air * deltaT) if deltaT > 0 else 0
    eps = min(eps, 0.99)
    NTU = -math.log(1 - eps)
    UA = NTU * C_air
    A = UA / U
    T_air_out = T_air_in + Q / C_air
    return {
        "Zone": name,
        "Q (kW)": q_kw,
        "Effectiveness": eps,
        "NTU": NTU,
        "UA (W/K)": UA,
        "U (W/m²·K)": U,
        "Area Required (m²)": A,
        "T_air_out (°C)": T_air_out
    }

# Heat loads
q_desuper = m_dot * (h1 - h2) / 1000
q_cond = m_dot * (h2 - h3) / 1000
q_subcool = m_dot * (h3 - h4) / 1000

# NTU calculation
res_subcool = zone_calc("Subcooling", q_subcool, T_out_subcool, T_air_in, U_subcool)
res_cond = zone_calc("Condensing", q_cond, T_cond, res_subcool["T_air_out (°C)"], U_cond)
res_desuper = zone_calc("Desuperheating", q_desuper, T_in_superheat, res_cond["T_air_out (°C)"], U_desuper)

# Area/tube length ratio and tube rows per zone
area_per_m_tube = total_air_side_area / total_tube_length
for res in [res_subcool, res_cond, res_desuper]:
    res["Tube Length (m)"] = res["Area Required (m²)"] / area_per_m_tube
    res["Tube Rows Required"] = res["Tube Length (m)"] / (tube_length_per_tube * tubes_per_row)

df = pd.DataFrame([res_subcool, res_cond, res_desuper])

st.header("Zone-wise NTU + Area + Tube Geometry")
st.dataframe(df)

st.header("Geometry Breakdown")
st.write(f"**Tubes per Row:** {tubes_per_row}")
st.write(f"**Tube Length per Tube:** {tube_length_per_tube:.2f} m")
st.write(f"**Total Tubes:** {total_tubes}")
st.write(f"**Total Tube Length:** {total_tube_length:.2f} m")
st.write(f"**Number of Fins:** {number_of_fins}")
st.write(f"**Area per Fin:** {area_per_fin:.4f} m²")
st.write(f"**Total Fin Area:** {total_fin_area:.2f} m²")
st.write(f"**Tube External Area:** {tube_ext_area:.2f} m²")
st.write(f"**Total Air-Side Area:** {total_air_side_area:.2f} m²")
st.write(f"**Area per meter of tube:** {area_per_m_tube:.4f} m²/m")

st.header("Refrigerant Enthalpies")
st.write(f"h1 (Inlet Superheated): {h1:.2f} J/kg")
st.write(f"h2 (Saturated Vapor): {h2:.2f} J/kg")
st.write(f"h3 (Saturated Liquid): {h3:.2f} J/kg")
st.write(f"h4 (Subcooled): {h4:.2f} J/kg")


# Air velocity and Reynolds number section
face_area_m2 = face_width_m * face_height_m
airflow_m3s = airflow_cmh / 3600
air_velocity_face = airflow_m3s / face_area_m2
air_velocity_fin = airflow_m3s / (face_area_m2 * (free_area_percent / 100))

T_air_K = T_air_in + 273.15
rho_air = PropsSI("D", "T", T_air_K, "P", 101325, "Air")
mu_air = PropsSI("V", "T", T_air_K, "P", 101325, "Air")

Re_fin = rho_air * air_velocity_fin * tube_od / mu_air

st.header("Air Flow Velocities and Reynolds Number")

# --------------------------------

# --------------------------------
# Zukauskas Correlation for Tube Banks
# --------------------------------
# Compute Reynolds number based on tube OD
# --- Air properties for Re, Nu, h ---
T_K = air_temp_C + 273.15
P = 101325  # Pa
rho = PropsSI('D', 'T', T_K, 'P', P, 'Air')
mu = PropsSI('V', 'T', T_K, 'P', P, 'Air')
k = PropsSI('L', 'T', T_K, 'P', P, 'Air')
cp = PropsSI('C', 'T', T_K, 'P', P, 'Air')
Re = rho * air_velocity_fin * tube_od_m / mu
Pr = cp * mu / k  # Prandtl number

# Zukauskas for staggered tube bank, valid 100 < Re < 10000
C_z = 0.193
m_z = 0.618
Nu = C_z * (Re ** m_z) * (Pr ** (1/3))
h_air = Nu * k / tube_od_m
U = h_air  # assuming air-side dominates

st.subheader("Zukauskas Heat Transfer Estimation")
st.write(f"**Reynolds Number (Re):** {Re:.2f}")
st.write(f"**Prandtl Number (Pr):** {Pr:.3f}")
st.write(f"**Nusselt Number (Nu):** {Nu:.2f}")
st.write(f"**Air-side Heat Transfer Coefficient (h):** {h_air:.2f} W/m²-K")
st.write(f"**Estimated Overall U-value:** {U:.2f} W/m²-K")
st.write(f"**Air Face Velocity:** {air_velocity_face:.2f} m/s")
st.write(f"**Air Velocity in Fin Passage:** {air_velocity_fin:.2f} m/s")
st.write(f"**Reynolds Number (based on tube OD):** {Re_fin:.0f}")
