
import math
import streamlit as st
import pandas as pd
from CoolProp.CoolProp import PropsSI, get_global_param_string

st.title("Air-Cooled Condenser NTU & Area Tool (Subcool → Condense → Desuperheat)")

# Refrigerant setup
st.sidebar.header("Refrigerant Inputs")
fluid_list = get_global_param_string("FluidsList").split(',')
refrigerants = sorted([f for f in fluid_list if f.startswith("R")])
fluid = st.sidebar.selectbox("Refrigerant", refrigerants, index=refrigerants.index("R134a") if "R134a" in refrigerants else 0)
m_dot = st.sidebar.number_input("Mass Flow Rate of Refrigerant (kg/s)", value=0.6)
T_in_superheat = st.sidebar.number_input("Refrigerant Inlet Temp (°C)", value=95.0)
T_cond = st.sidebar.number_input("Condensation Temp (°C)", value=45.0)
T_out_subcool = st.sidebar.number_input("Subcooled Liquid Temp (°C)", value=50.0)

# Pressure at condensation
P_cond = PropsSI("P", "T", T_cond + 273.15, "Q", 0, fluid)

# Enthalpies
h1 = PropsSI("H", "T", T_in_superheat + 273.15, "P", P_cond, fluid)
h2 = PropsSI("H", "T", T_cond + 273.15, "Q", 1, fluid)
h3 = PropsSI("H", "T", T_cond + 273.15, "Q", 0, fluid)
h4 = PropsSI("H", "T", T_out_subcool + 273.15, "P", P_cond, fluid)

# Zone heat loads (kW)
q_desuper = m_dot * (h1 - h2) / 1000
q_cond = m_dot * (h2 - h3) / 1000
q_subcool = m_dot * (h3 - h4) / 1000

# Air side
st.sidebar.header("Air Side Inputs")
T_air_in = st.sidebar.number_input("Air Inlet Temp (°C)", value=35.0)
airflow_cmh = st.sidebar.number_input("Air Flow Rate (m³/h)", value=10000)
cp_air = PropsSI("C", "T", T_air_in + 273.15, "P", 101325, "Air")
rho_air = PropsSI("D", "T", T_air_in + 273.15, "P", 101325, "Air")
m_dot_air = airflow_cmh / 3600 * rho_air

# U values
st.sidebar.header("Heat Transfer Coefficient Inputs")
U_subcool = st.sidebar.number_input("U - Subcooling (W/m²·K)", value=40.0)
U_cond = st.sidebar.number_input("U - Condensation (W/m²·K)", value=60.0)
U_desuper = st.sidebar.number_input("U - Desuperheating (W/m²·K)", value=40.0)

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

# Zone calculations in reverse order of air flow
res_subcool = zone_calc("Subcooling", q_subcool, T_out_subcool, T_air_in, U_subcool)
res_cond = zone_calc("Condensation", q_cond, T_cond, res_subcool["T_air_out (°C)"], U_cond)
res_desuper = zone_calc("Desuperheating", q_desuper, T_in_superheat, res_cond["T_air_out (°C)"], U_desuper)

df = pd.DataFrame([res_subcool, res_cond, res_desuper])

st.header("Zone-wise NTU & Area Results")
st.dataframe(df)

st.header("Refrigerant State Enthalpies")
st.write(f"**h1 (Superheated Inlet):** {h1:.2f} J/kg")
st.write(f"**h2 (Saturated Vapor):** {h2:.2f} J/kg")
st.write(f"**h3 (Saturated Liquid):** {h3:.2f} J/kg")
st.write(f"**h4 (Subcooled Liquid):** {h4:.2f} J/kg")
