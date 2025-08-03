
import math
import streamlit as st
from CoolProp.CoolProp import PropsSI, get_global_param_string

st.title("Full Condenser Coil + Refrigerant Heat Load + NTU Effectiveness Tool")

# ------------------------------
# Coil Geometry Inputs
# ------------------------------
st.header("Air-Side Coil Geometry")
tube_od_mm = st.number_input("Tube Outer Diameter (mm)", value=9.525)
tube_thickness_mm = st.number_input("Tube Wall Thickness (mm)", value=0.35)
row_pitch_mm = st.number_input("Row Pitch (mm)", value=25.4)
tube_pitch_mm = st.number_input("Tube to Tube Pitch (mm)", value=25.4)
fin_thickness_mm = st.number_input("Fin Thickness (mm)", value=0.12)
fpi = st.number_input("Fins per Inch (FPI)", value=12, step=1)
num_rows = st.number_input("Number of Tube Rows", value=4, step=1)
face_width_m = st.number_input("Coil Face Width (m)", value=1.0)
face_height_m = st.number_input("Coil Face Height (m)", value=1.0)
free_area_percent = st.slider("Free Flow Area (%)", min_value=10, max_value=100, value=25)

# ------------------------------
# Air Flow Inputs
# ------------------------------
st.header("Air-Side Operating Conditions")
air_flow_cmh = st.number_input("Air Flow Rate (m³/h)", value=10000)
T_air_in = st.number_input("Air Inlet Temperature (°C)", value=35.0)

# Convert air-side geometry
tube_od = tube_od_mm / 1000
row_pitch = row_pitch_mm / 1000
tube_pitch = tube_pitch_mm / 1000
fin_thickness = fin_thickness_mm / 1000
fins_per_m = fpi * 39.3701
frontal_area = face_width_m * face_height_m
fin_depth = row_pitch * num_rows
tubes_per_row = math.floor(face_width_m / tube_pitch)
total_tubes = tubes_per_row * num_rows
tube_ext_area = total_tubes * (math.pi * tube_od)
fin_area_per_fin = 2 * face_width_m * fin_depth
gross_fin_area = fin_area_per_fin * fins_per_m
hole_area_per_tube = (math.pi / 4) * tube_od**2
total_hole_area = hole_area_per_tube * total_tubes * fins_per_m
net_fin_area = gross_fin_area - total_hole_area
total_air_side_area = (tube_ext_area + net_fin_area) * face_height_m
net_free_flow_area = frontal_area * (free_area_percent / 100)

# ------------------------------
# Air Properties & Velocity
# ------------------------------
air_flow_m3s = air_flow_cmh / 3600
T_air_K = T_air_in + 273.15
rho_air = PropsSI("D", "T", T_air_K, "P", 101325, "Air")
mu_air = PropsSI("V", "T", T_air_K, "P", 101325, "Air")
cp_air = PropsSI("C", "T", T_air_K, "P", 101325, "Air")
m_dot_air = air_flow_m3s * rho_air
air_velocity = air_flow_m3s / net_free_flow_area if net_free_flow_area > 0 else 0
Re_air = rho_air * air_velocity * tube_od / mu_air
f_air = 0.25 * Re_air**-0.25 if Re_air > 0 else 0
dp_air = (f_air * num_rows * rho_air * air_velocity**2) / 2

st.subheader("Air-Side Calculations")
st.write(f"Tubes per Row: {tubes_per_row}")
st.write(f"Total Tubes: {total_tubes}")
st.write(f"Tube External Area: {tube_ext_area:.4f} m²")
st.write(f"Net Fin Area: {net_fin_area:.4f} m²")
st.write(f"Total Air-Side Area: {total_air_side_area:.4f} m²")
st.write(f"Free Flow Area A_min: {net_free_flow_area:.4f} m²")
st.write(f"Air Velocity: {air_velocity:.2f} m/s")
st.write(f"Air Density: {rho_air:.3f} kg/m³")
st.write(f"Air Viscosity: {mu_air:.7f} Pa·s")
st.write(f"Reynolds Number: {Re_air:.1f}")
st.write(f"Friction Factor: {f_air:.4f}")
st.write(f"Air-side Pressure Drop: {dp_air:.2f} Pa")

# ------------------------------
# Refrigerant & NTU Section
# ------------------------------
st.header("Refrigerant Heat Load + NTU Analysis")

fluid_list = get_global_param_string("FluidsList").split(',')
refrigerants = sorted([f for f in fluid_list if f.startswith("R")])
fluid = st.selectbox("Select Refrigerant", refrigerants, index=refrigerants.index("R134a") if "R134a" in refrigerants else 0)
m_dot_ref = st.number_input("Refrigerant Mass Flow Rate (kg/s)", value=0.6)

# Helper function
def effectiveness_crossflow(NTU, R):
    if R == 0:
        return 1 - math.exp(-NTU)
    return 1 - math.exp((1/R)*(math.exp(-R*NTU) - 1))

def NTU_from_eps(eps, R):
    if R == 0:
        return -math.log(1 - eps)
    return -(1/R) * math.log(1 - R * math.log(1 - eps))

def zone(name, q_kW, T_ref_in_C, T_ref_out_C, T_air_in_C):
    T_ref_in_K = T_ref_in_C + 273.15
    T_ref_out_K = T_ref_out_C + 273.15
    T_air_in_K = T_air_in_C + 273.15
    h_in = PropsSI("H", "T", T_ref_in_K, "P", 101325, fluid)
    h_out = PropsSI("H", "T", T_ref_out_K, "P", 101325, fluid)
    cp_ref = (h_in - h_out)/(T_ref_in_K - T_ref_out_K)
    C_ref = m_dot_ref * cp_ref
    C_air = m_dot_air * cp_air
    C_min = min(C_ref, C_air)
    R = C_min / max(C_ref, C_air)
    Q = q_kW * 1000
    eps = Q / (C_min * (T_ref_in_K - T_air_in_K))
    NTU = NTU_from_eps(eps, R)
    UA = NTU * C_min
    T_air_out_K = T_air_in_K + Q / (m_dot_air * cp_air)
    st.subheader(name)
    st.write(f"Air In: {T_air_in_C:.2f} °C → Out: {T_air_out_K - 273.15:.2f} °C")
    st.write(f"ε = {eps:.3f} | NTU = {NTU:.3f} | R = {R:.3f}")
    st.write(f"UA Required: {UA:.1f} W/K")
    return T_air_out_K - 273.15

# Zone 1: Desuperheating
q1 = st.number_input("Desuperheating Load (kW)", value=5.0)
T1_ref_in = st.number_input("Ref Inlet Temp (°C) - Desuperheat", value=95.0)
T1_ref_out = st.number_input("Ref Outlet Temp (°C) - Desuperheat", value=60.0)
T_air_out1 = zone("Desuperheating Zone", q1, T1_ref_in, T1_ref_out, T_air_in)

# Zone 2: Condensing
q2 = st.number_input("Condensing Load (kW)", value=90.0)
T_cond = st.number_input("Condensing Temp (°C)", value=54.0)
C_air2 = m_dot_air * cp_air
eps2 = (q2 * 1000) / (C_air2 * (T_cond + 273.15 - (T_air_out1 + 273.15)))
NTU2 = -math.log(1 - eps2)
UA2 = NTU2 * C_air2
T_air_out2 = T_air_out1 + q2 * 1000 / (m_dot_air * cp_air)
st.subheader("Condensing Zone")
st.write(f"Air In: {T_air_out1:.2f} °C → Out: {T_air_out2:.2f} °C")
st.write(f"ε = {eps2:.3f} | NTU = {NTU2:.3f}")
st.write(f"UA Required: {UA2:.1f} W/K")

# Zone 3: Subcooling
q3 = st.number_input("Subcooling Load (kW)", value=15.0)
T3_ref_in = st.number_input("Ref Inlet Temp (°C) - Subcooling", value=52.0)
T3_ref_out = st.number_input("Ref Outlet Temp (°C) - Subcooling", value=45.0)
T_air_out3 = zone("Subcooling Zone", q3, T3_ref_in, T3_ref_out, T_air_out2)
