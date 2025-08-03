
import math
import streamlit as st
import pandas as pd
from CoolProp.CoolProp import PropsSI, get_global_param_string

st.title("Air-Cooled Condenser Design Tool (NTU + Geometry Breakdown)")

# Refrigerant setup
st.sidebar.header("Refrigerant Inputs")
fluid_list = get_global_param_string("FluidsList").split(',')
refrigerants = sorted([f for f in fluid_list if f.startswith("R")])
fluid = st.sidebar.selectbox("Refrigerant", refrigerants, index=refrigerants.index("R134a") if "R134a" in refrigerants else 0)
m_dot = st.sidebar.number_input("Refrigerant Mass Flow Rate (kg/s)", value=0.6)
T_in_superheat = st.sidebar.number_input("Inlet Superheat Temp (°C)", value=95.0)
T_cond = st.sidebar.number_input("Condensation Temp (°C)", value=45.0)
T_out_subcool = st.sidebar.number_input("Subcooled Liquid Temp (°C)", value=50.0)

# Air side
st.sidebar.header("Air Inputs")
T_air_in = st.sidebar.number_input("Air Inlet Temp (°C)", value=35.0)
airflow_cmh = st.sidebar.number_input("Air Flow Rate (m³/h)", value=10000)

# Coil Geometry
st.sidebar.header("Coil Geometry")
tube_od_mm = st.sidebar.number_input("Tube Outer Diameter (mm)", value=9.525)
tube_thickness_mm = st.sidebar.number_input("Tube Wall Thickness (mm)", value=0.35)
tube_pitch_mm = st.sidebar.number_input("Tube Pitch (mm)", value=25.4)
row_pitch_mm = st.sidebar.number_input("Row Pitch (mm)", value=22.0)
fin_thickness_mm = st.sidebar.number_input("Fin Thickness (mm)", value=0.12)
fpi = st.sidebar.number_input("Fins Per Inch", value=12)
face_width_m = st.sidebar.number_input("Coil Face Width (m)", value=1.0)
face_height_m = st.sidebar.number_input("Coil Face Height (m)", value=1.0)

# U values
st.sidebar.header("U Values (W/m²·K)")
U_subcool = st.sidebar.number_input("U - Subcooling", value=40.0)
U_cond = st.sidebar.number_input("U - Condensing", value=60.0)
U_desuper = st.sidebar.number_input("U - Desuperheating", value=40.0)

# Enthalpies and pressures
P_cond = PropsSI("P", "T", T_cond + 273.15, "Q", 0, fluid)
h1 = PropsSI("H", "T", T_in_superheat + 273.15, "P", P_cond, fluid)
h2 = PropsSI("H", "T", T_cond + 273.15, "Q", 1, fluid)
h3 = PropsSI("H", "T", T_cond + 273.15, "Q", 0, fluid)
h4 = PropsSI("H", "T", T_out_subcool + 273.15, "P", P_cond, fluid)

# Zone-wise heat loads
q_desuper = m_dot * (h1 - h2) / 1000
q_cond = m_dot * (h2 - h3) / 1000
q_subcool = m_dot * (h3 - h4) / 1000

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

# NTU calculation zone by zone
res_subcool = zone_calc("Subcooling", q_subcool, T_out_subcool, T_air_in, U_subcool)
res_cond = zone_calc("Condensing", q_cond, T_cond, res_subcool["T_air_out (°C)"], U_cond)
res_desuper = zone_calc("Desuperheating", q_desuper, T_in_superheat, res_cond["T_air_out (°C)"], U_desuper)

# Geometry calculations
tube_od = tube_od_mm / 1000
tube_pitch = tube_pitch_mm / 1000
row_pitch = row_pitch_mm / 1000
fins_per_m = fpi * 39.37
fin_thickness = fin_thickness_mm / 1000

tubes_per_row = int(face_width_m / tube_pitch)
total_tubes = tubes_per_row * int(face_height_m / row_pitch)
tube_ext_area = total_tubes * (math.pi * tube_od * face_height_m)
fin_area_per_fin = 2 * face_width_m * face_height_m
total_fins = fins_per_m * face_width_m
total_fin_area = total_fins * fin_area_per_fin
total_air_side_area = total_fin_area + tube_ext_area

total_tube_length = total_tubes * face_height_m
area_per_m_tube = total_air_side_area / total_tube_length

for res in [res_subcool, res_cond, res_desuper]:
    res["Tube Length (m)"] = res["Area Required (m²)"] / area_per_m_tube
    res["Tube Rows Required"] = res["Tube Length (m)"] / (face_height_m * tubes_per_row)

df = pd.DataFrame([res_subcool, res_cond, res_desuper])

st.header("Zone-wise NTU + Area + Tube Geometry")
st.dataframe(df)

st.header("Air Side Geometry Summary")
st.write(f"**Tube External Area:** {tube_ext_area:.2f} m²")
st.write(f"**Fin Area:** {total_fin_area:.2f} m²")
st.write(f"**Total Air Side Area:** {total_air_side_area:.2f} m²")
st.write(f"**Total Tube Length:** {total_tube_length:.2f} m")
st.write(f"**Area per meter of tube:** {area_per_m_tube:.4f} m²/m")

st.header("Refrigerant Enthalpies")
st.write(f"h1 (Inlet Superheated): {h1:.2f} J/kg")
st.write(f"h2 (Saturated Vapor): {h2:.2f} J/kg")
st.write(f"h3 (Saturated Liquid): {h3:.2f} J/kg")
st.write(f"h4 (Subcooled): {h4:.2f} J/kg")
