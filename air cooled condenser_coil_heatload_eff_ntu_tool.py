
import math
import streamlit as st
import pandas as pd
from CoolProp.CoolProp import PropsSI, get_global_param_string

st.title("Air-Cooled Condenser: Geometry, Heat Load, Pressure Drop")

# -------------------------------
# Coil Geometry Input Parameters
# -------------------------------
st.header("Coil Geometry")

tube_od_mm = st.number_input("Tube Outer Diameter (mm)", value=9.525)
tube_thickness_mm = st.number_input("Tube Wall Thickness (mm)", value=0.35)
row_pitch_mm = st.number_input("Row Pitch (mm)", value=25.4)
tube_pitch_mm = st.number_input("Tube Pitch (mm)", value=25.4)
fin_thickness_mm = st.number_input("Fin Thickness (mm)", value=0.12)
fpi = st.number_input("Fins per Inch (FPI)", value=12)
num_rows = st.number_input("Number of Tube Rows", value=4)
face_width_m = st.number_input("Coil Face Width (m)", value=1.2)
face_height_m = st.number_input("Coil Face Height (m)", value=1.0)
airflow_cmh = st.number_input("Airflow Rate (m³/h)", value=10000)
air_temp_C = st.number_input("Air Inlet Temperature (°C)", value=35.0)
free_area_percent = st.slider("Free Flow Area (%)", 10, 100, 30)
fin_material = st.selectbox("Fin Material", ["Aluminum", "Copper"])

# Unit conversions
tube_od_m = tube_od_mm / 1000
tube_thickness_m = tube_thickness_mm / 1000
row_pitch_m = row_pitch_mm / 1000
tube_pitch_m = tube_pitch_mm / 1000
fin_thickness_m = fin_thickness_mm / 1000
fins_per_m = fpi * 39.37
airflow_m3s = airflow_cmh / 3600
face_area_m2 = face_width_m * face_height_m
net_free_flow_area = face_area_m2 * (free_area_percent / 100)

# Tube layout
tubes_per_row = int(face_width_m / tube_pitch_m)
total_tubes = tubes_per_row * int(num_rows)
tube_length_per_tube = face_height_m
total_tube_length = total_tubes * tube_length_per_tube

# Fin geometry
num_fins = int(face_height_m * fpi / 0.0254)
area_per_fin = face_width_m * (num_rows * row_pitch_m)
total_fin_area = num_fins * area_per_fin
tube_external_area = math.pi * tube_od_m * total_tube_length
total_air_side_area = tube_external_area + total_fin_area
area_per_meter_tube = total_air_side_area / total_tube_length

# Air velocities
air_velocity_face = airflow_m3s / face_area_m2
air_velocity_fin = airflow_m3s / net_free_flow_area

# Air properties and Re
T_K = air_temp_C + 273.15
rho = PropsSI("D", "T", T_K, "P", 101325, "Air")
mu = PropsSI("V", "T", T_K, "P", 101325, "Air")
k_air = PropsSI("L", "T", T_K, "P", 101325, "Air")
cp_air = PropsSI("C", "T", T_K, "P", 101325, "Air")
Pr = cp_air * mu / k_air
Re = rho * air_velocity_fin * tube_od_m / mu

# Zukauskas Nu
C, m, n = (0.9, 0.4, 0.36) if Re < 1000 else (0.52, 0.5, 0.36)
Nu = C * (Re**m) * (Pr**n)
h_air = Nu * k_air / tube_od_m

# Fin efficiency
k_fin = 205 if fin_material == "Aluminum" else 385
L_fin = 0.5 * (math.sqrt(row_pitch_m**2 + tube_pitch_m**2) - tube_od_m)
m_fin = math.sqrt(2 * h_air / (k_fin * fin_thickness_m))
fin_eff = math.tanh(m_fin * L_fin) / (m_fin * L_fin)

# -------------------------------
# Refrigerant Heat Load
# -------------------------------
st.header("Refrigerant Heat Load Calculation")

fluid_list = get_global_param_string("FluidsList").split(',')
refrigerants = sorted([f for f in fluid_list if f.startswith("R")])
fluid = st.selectbox("Select Refrigerant", refrigerants, index=refrigerants.index("R134a") if "R134a" in refrigerants else 0)

m_dot = st.number_input("Refrigerant Mass Flow Rate (kg/s)", value=0.6)
T_superheat = st.number_input("Superheated Inlet Temp (°C)", value=95.0)
T_subcool = st.number_input("Subcooled Outlet Temp (°C)", value=52.7)
T_cond = st.number_input("Condensation Temp (°C)", value=55.0)
P_cond = PropsSI("P", "T", T_cond + 273.15, "Q", 0, fluid)

T1 = T_superheat + 273.15
T3 = T_subcool + 273.15
T_sat = T_cond + 273.15

h1 = PropsSI("H", "P", P_cond, "T", T1, fluid)
h2 = PropsSI("H", "P", P_cond, "Q", 1, fluid)
h3 = PropsSI("H", "P", P_cond, "Q", 0, fluid)
h4 = PropsSI("H", "P", P_cond, "T", T3, fluid)

Q_desuper = m_dot * (h1 - h2) / 1000
Q_cond = m_dot * (h2 - h3) / 1000
Q_subcool = m_dot * (h3 - h4) / 1000
Q_total = Q_desuper + Q_cond + Q_subcool

# -------------------------------
# Refrigerant Velocity in Tubes
# -------------------------------
st.subheader("Refrigerant Circuit Velocity")
num_feeds = st.number_input("Number of Feeds/Circuits", value=4, min_value=1)
circuits = max(1, tubes_per_row // num_feeds)
mass_flow_per_circuit = m_dot / circuits
tube_id_m = tube_od_m - 2 * tube_thickness_m
A_tube_internal = (math.pi / 4) * tube_id_m**2
rho_vap = PropsSI("D", "P", P_cond, "Q", 1, fluid)
rho_liq = PropsSI("D", "P", P_cond, "Q", 0, fluid)
rho_avg = (rho_vap + rho_liq) / 2
v_refrigerant = mass_flow_per_circuit / (rho_avg * A_tube_internal)

# -------------------------------
# Output
# -------------------------------
st.header("Output Summary")

st.subheader("Geometry Breakdown")
st.write(f"Tubes per Row: {tubes_per_row}")
st.write(f"Tube Length per Tube: {tube_length_per_tube:.2f} m")
st.write(f"Total Tubes: {total_tubes}")
st.write(f"Total Tube Length: {total_tube_length:.2f} m")
st.write(f"Number of Fins: {num_fins}")
st.write(f"Area per Fin: {area_per_fin:.4f} m²")
st.write(f"Total Fin Area: {total_fin_area:.2f} m²")
st.write(f"Tube External Area: {tube_external_area:.2f} m²")
st.write(f"Total Air-Side Area: {total_air_side_area:.2f} m²")
st.write(f"Area per meter of tube: {area_per_meter_tube:.4f} m²/m")

st.subheader("Air Flow and Heat Transfer")
st.write(f"Air Face Velocity: {air_velocity_face:.2f} m/s")
st.write(f"Air Velocity in Fin Passage: {air_velocity_fin:.2f} m/s")
st.write(f"Reynolds Number (Re): {Re:.2f}")
st.write(f"Nusselt Number (Nu): {Nu:.2f}")
st.write(f"Air-side h (W/m²-K): {h_air:.2f}")
st.write(f"Fin Efficiency: {fin_eff:.4f}")

st.subheader("Refrigerant Enthalpies")
st.write(f"h1 (Inlet Superheated): {h1:.2f} J/kg")
st.write(f"h2 (Saturated Vapor): {h2:.2f} J/kg")
st.write(f"h3 (Saturated Liquid): {h3:.2f} J/kg")
st.write(f"h4 (Subcooled): {h4:.2f} J/kg")

st.subheader("Heat Load Breakdown")
st.write(f"Sensible (Desuperheating): {Q_desuper:.2f} kW")
st.write(f"Latent (Condensing): {Q_cond:.2f} kW")
st.write(f"Subcooling: {Q_subcool:.2f} kW")
st.write(f"Total Heat Removed: {Q_total:.2f} kW")

st.subheader("Refrigerant Circuit Velocity")
st.write(f"Number of Circuits: {circuits}")
st.write(f"Mass Flow per Circuit: {mass_flow_per_circuit:.4f} kg/s")
st.write(f"Tube Inner Diameter: {tube_id_m*1000:.2f} mm")
st.write(f"Refrigerant Velocity in Tube: {v_refrigerant:.2f} m/s")
