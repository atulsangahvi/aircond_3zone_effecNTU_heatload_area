# Adding full output breakdown as per user's requirement
import math
import streamlit as st
import pandas as pd
from CoolProp.CoolProp import PropsSI, get_global_param_string

st.title("Comprehensive Air-Cooled Condenser Design Tool")

# ----------------------
# Inputs: Geometry
# ----------------------
st.sidebar.header("Geometry Inputs")
tube_od_mm = st.sidebar.number_input("Tube Outer Diameter (mm)", value=9.525)
tube_thickness_mm = st.sidebar.number_input("Tube Wall Thickness (mm)", value=0.35)
row_pitch_mm = st.sidebar.number_input("Row Pitch (mm)", value=25.4)
tube_pitch_mm = st.sidebar.number_input("Tube Pitch (mm)", value=25.4)
fpi = st.sidebar.number_input("Fins Per Inch", value=12)
fin_thickness_mm = st.sidebar.number_input("Fin Thickness (mm)", value=0.12)
fin_material = st.sidebar.selectbox("Fin Material", ["Aluminum", "Copper"])
face_width_m = st.sidebar.number_input("Coil Face Width (m)", value=1.0)
face_height_m = st.sidebar.number_input("Coil Face Height (m)", value=1.0)
num_rows = st.sidebar.number_input("Total Number of Rows", value=4, step=1)
free_area_percent = st.sidebar.slider("Free Flow Area %", 10, 100, 25)
num_feeds = st.sidebar.number_input("Number of Feeds", value=9)

# ----------------------
# Inputs: Refrigerant
# ----------------------
st.sidebar.header("Refrigerant Inputs")
fluid_list = get_global_param_string("FluidsList").split(',')
refrigerants = sorted([f for f in fluid_list if f.startswith("R")])
fluid = st.sidebar.selectbox("Refrigerant", refrigerants, index=refrigerants.index("R134a") if "R134a" in refrigerants else 0)

T1 = st.sidebar.number_input("Inlet Superheat Temp (°C)", value=95.0)
T3 = st.sidebar.number_input("Outlet Subcooled Temp (°C)", value=52.7)
T_cond = st.sidebar.number_input("Condensing Temp (°C)", value=55.0)
m_dot = st.sidebar.number_input("Mass Flow Rate (kg/s)", value=0.6)
air_inlet_temp = st.sidebar.number_input("Air Inlet Temp (°C)", value=35.0)
airflow_cmh = st.sidebar.number_input("Air Flow (m³/hr)", value=10000)

# ----------------------
# U-value (Optional)
# ----------------------
st.sidebar.header("U-Value Inputs (Optional)")
U_user = st.sidebar.checkbox("Use custom U-values")
U_desuper = st.sidebar.number_input("U-value Desuperheating (W/m²K)", value=45.0)
U_cond = st.sidebar.number_input("U-value Condensation (W/m²K)", value=55.0)
U_subcool = st.sidebar.number_input("U-value Subcooling (W/m²K)", value=65.0)

# ----------------------
# Geometry Calculations
# ----------------------
tube_od_m = tube_od_mm / 1000
tube_thk_m = tube_thickness_mm / 1000
tube_id_m = tube_od_m - 2 * tube_thk_m
row_pitch_m = row_pitch_mm / 1000
tube_pitch_m = tube_pitch_mm / 1000
fin_thk_m = fin_thickness_mm / 1000
fins_per_m = fpi * 39.37
fin_k = 235 if fin_material == "Aluminum" else 400
frontal_area_m2 = face_width_m * face_height_m
net_free_area = frontal_area_m2 * (free_area_percent / 100)
airflow_m3s = airflow_cmh / 3600

# Tube layout
corrected_row_pitch_m = math.sqrt((tube_pitch_m / 2)**2 + row_pitch_m**2)
tubes_per_row = math.floor(face_width_m / tube_pitch_m)
tube_length_per_tube = face_height_m
total_tubes = tubes_per_row * num_rows
total_tube_length = total_tubes * tube_length_per_tube
A_tube_ext = total_tube_length * math.pi * tube_od_m

# Fin Area
num_fins = math.floor(face_height_m / 0.0254 * fpi)
area_per_fin = row_pitch_m * face_width_m
A_fin_raw = area_per_fin * num_fins
tube_hole_area = total_tubes * (math.pi / 4) * tube_od_m**2
A_fin = A_fin_raw - tube_hole_area

# Fin Efficiency
L_fin = min(row_pitch_m, tube_pitch_m) / 2
m = math.sqrt((2 * 40) / (fin_k * fin_thk_m))
eta_fin = math.tanh(m * L_fin) / (m * L_fin)
A_air_total = A_tube_ext + A_fin * eta_fin
A_air_per_m = A_air_total / total_tube_length

# Air properties
T_K = air_inlet_temp + 273.15
mu = PropsSI("V", "T", T_K, "P", 101325, "Air")
k_air = PropsSI("L", "T", T_K, "P", 101325, "Air")
rho_air = PropsSI("D", "T", T_K, "P", 101325, "Air")
cp_air = PropsSI("C", "T", T_K, "P", 101325, "Air")
Pr = mu * cp_air / k_air

# Velocities and Nu
air_velocity_face = airflow_m3s / frontal_area_m2
air_velocity_fin = airflow_m3s / net_free_area
Re = rho_air * air_velocity_fin * tube_od_m / mu
C = 0.27 if Re < 1000 else 0.021
n = 0.63 if Re < 1000 else 0.84
Nu = C * Re**n * Pr**0.36
h_air = Nu * k_air / tube_od_m
U_final = lambda U_fixed: U_fixed if U_user else h_air

# Refrigerant calculations
P_cond = PropsSI("P", "T", T_cond + 273.15, "Q", 0, fluid)
h1 = PropsSI("H", "P", P_cond, "T", T1 + 273.15, fluid)
h2 = PropsSI("H", "P", P_cond, "Q", 1, fluid)
h3 = PropsSI("H", "P", P_cond, "Q", 0, fluid)
h4 = PropsSI("H", "P", P_cond, "T", T3 + 273.15, fluid)
cp_ref = PropsSI("C", "P", P_cond, "T", (T1 + T3)/2 + 273.15, fluid) / 1000
Q_sens = m_dot * (h1 - h2) / 1000
Q_lat = m_dot * (h2 - h3) / 1000
Q_sub = m_dot * (h3 - h4) / 1000

# Circuits and refrigerant velocity
mass_flow_per_circuit = m_dot / num_feeds
A_flow_refrig = (math.pi / 4) * tube_id_m**2
v_refrig = mass_flow_per_circuit / (A_flow_refrig * rho_air)

# Display
st.subheader("Refrigerant Side Summary")
st.write(f"Number of Circuits: {num_feeds}")
st.write(f"Mass Flow per Circuit: {mass_flow_per_circuit:.3f} kg/s")
st.write(f"Tube Inner Diameter: {tube_id_m*1000:.2f} mm")
st.write(f"Refrigerant Velocity: {v_refrig:.2f} m/s")

st.subheader("Air Side Summary")
st.write(f"Tubes per Row: {tubes_per_row}")
st.write(f"Total Tube Length: {total_tube_length:.2f} m")
st.write(f"Total Tube Area: {A_tube_ext:.2f} m²")
st.write(f"Fin Area (corrected): {A_fin:.2f} m²")
st.write(f"Total Air-Side Area: {A_air_total:.2f} m²")
st.write(f"Face Velocity: {air_velocity_face:.2f} m/s")
st.write(f"Fin Passage Velocity: {air_velocity_fin:.2f} m/s")
st.write(f"Reynolds Number: {Re:.0f}")
st.write(f"Nusselt Number: {Nu:.1f}")
st.write(f"Air-Side h (W/m²K): {h_air:.2f}")
st.write(f"U Estimated: {h_air:.2f} W/m²K")
