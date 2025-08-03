
import math
import streamlit as st
import pandas as pd
from CoolProp.CoolProp import PropsSI

st.title("Air-Cooled Condenser (Subcooling → Condensing → Desuperheating) NTU Tool")

# ----------------- USER INPUT -----------------
st.sidebar.header("Air Side Inputs")
face_area = st.sidebar.number_input("Frontal Coil Area (m²)", value=1.0)
free_flow_percent = st.sidebar.slider("Free Flow Area (%)", 10, 100, 25)
airflow_cmh = st.sidebar.number_input("Air Flow Rate (m³/h)", value=10000)
T_air_in = st.sidebar.number_input("Inlet Air Temperature (°C)", value=35.0)

st.sidebar.header("Coil Geometry")
num_rows = st.sidebar.number_input("Number of Tube Rows", min_value=1, max_value=6, value=4)
tube_od_mm = st.sidebar.number_input("Tube Outer Diameter (mm)", value=9.525)
fin_spacing_mm = st.sidebar.number_input("Fin Spacing (mm)", value=1.7)

st.sidebar.header("Refrigerant & Zone Loads")
fluid = st.sidebar.selectbox("Refrigerant", ["R134a", "R407C"])
q_subcool = st.sidebar.number_input("Subcooling Load (kW)", value=10.0)
q_cond = st.sidebar.number_input("Condensation Load (kW)", value=90.0)
q_desuper = st.sidebar.number_input("Desuperheating Load (kW)", value=10.0)
T_ref_subcool = st.sidebar.number_input("Refrigerant Temp (Subcooling) (°C)", value=50.0)
T_ref_cond = st.sidebar.number_input("Refrigerant Temp (Condensation) (°C)", value=45.0)
T_ref_desuper = st.sidebar.number_input("Refrigerant Temp (Desuperheating) (°C)", value=60.0)

# ----------------- CALCULATIONS -----------------
# Geometry
tube_od = tube_od_mm / 1000
A_min = face_area * (free_flow_percent / 100)
V_face = airflow_cmh / 3600 / A_min
rho_air = PropsSI("D", "T", T_air_in + 273.15, "P", 101325, "Air")
mu_air = PropsSI("V", "T", T_air_in + 273.15, "P", 101325, "Air")
cp_air = PropsSI("C", "T", T_air_in + 273.15, "P", 101325, "Air")
m_dot_air = airflow_cmh / 3600 * rho_air

# Hydraulic diameter
s_f = fin_spacing_mm / 1000
Dh = 4 * s_f * tube_od / (2 * (s_f + tube_od))
Re = rho_air * V_face * Dh / mu_air

# Rich empirical friction factor (approximate)
f_rich = 0.38 * Re**-0.25
dP_per_row = 0.5 * f_rich * rho_air * V_face**2
total_dP = dP_per_row * num_rows

# Effectiveness-NTU method
def calc_zone(q_kw, T_ref_C, T_air_C, m_dot_air, cp_air, zone_name):
    Q = q_kw * 1000
    C_air = m_dot_air * cp_air
    deltaT = (T_ref_C + 273.15) - (T_air_C + 273.15)
    eps = Q / (C_air * deltaT) if deltaT > 0 else 0
    eps = min(eps, 0.99)
    NTU = -math.log(1 - eps)
    UA = NTU * C_air
    T_air_out = T_air_C + Q / C_air
    return {
        "Zone": zone_name,
        "Effectiveness": eps,
        "NTU": NTU,
        "UA (W/K)": UA,
        "T_air_out (°C)": T_air_out
    }

results = []

# Zone 1: Subcooling
res1 = calc_zone(q_subcool, T_ref_subcool, T_air_in, m_dot_air, cp_air, "Subcooling")
results.append(res1)

# Zone 2: Condensing
res2 = calc_zone(q_cond, T_ref_cond, res1["T_air_out (°C)"], m_dot_air, cp_air, "Condensing")
results.append(res2)

# Zone 3: Desuperheating
res3 = calc_zone(q_desuper, T_ref_desuper, res2["T_air_out (°C)"], m_dot_air, cp_air, "Desuperheating")
results.append(res3)

# ----------------- OUTPUT -----------------
st.header("Results: Effectiveness NTU Per Zone")
st.dataframe(pd.DataFrame(results))

st.header("Air Side Performance")
st.write(f"**Free Flow Area A_min:** {A_min:.4f} m²")
st.write(f"**Face Velocity:** {V_face:.2f} m/s")
st.write(f"**Reynolds Number (air):** {Re:.0f}")
st.write(f"**Friction Factor (Rich):** {f_rich:.4f}")
st.write(f"**Pressure Drop per Row:** {dP_per_row:.2f} Pa")
st.write(f"**Total Air-side Pressure Drop:** {total_dP:.2f} Pa")
