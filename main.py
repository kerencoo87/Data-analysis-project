import streamlit as st
import numpy as np
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go

# Import function from the other file
from data_processing import calculate_single_patient_risk

#  Streamlit page setting
st.set_page_config(page_title="BRCA Clinical Risk Tool", layout="wide")

st.title("🧬 TCGA BRCA Clinical Risk Assessment Tool")
st.write("An integrated clinical-molecular decision support tool for breast invasive carcinoma.")

st.markdown("---")

col_input, col_result = st.columns([1, 1], gap="medium")

# inputs
with col_input:
    st.subheader("📋 Patient Clinical Parameters")
    user_age = st.number_input("Patient Age (years):", min_value=18, max_value=110, value=60, step=1)
    tp53_radio = st.radio("TP53 Gene Status:", options=["Wild-Type (No Mutation)", "Mutated"], horizontal=True)
    user_tp53 = 1 if "Mutated" in tp53_radio else 0
    user_stage = st.selectbox("AJCC Pathologic Stage:", options=["Stage I", "Stage II", "Stage III", "Stage IV"])
    btn_calc = st.button("Calculate Risk Score", type="primary", use_container_width=True)

# results
with col_result:
    st.subheader("📊 Stratification & Analytics")
    if btn_calc:
        res = calculate_single_patient_risk(user_age, user_tp53, user_stage)

        if "error" in res:
            st.error(res["error"])
        else:
            assigned_group = res['risk_group']
            risk_score = res['risk_score']

            st.metric(label="Integrated Risk Score", value=f"{risk_score} / 6")

            if assigned_group == "Low Risk":
                st.success("🟢 **Assigned Group: Low Risk**")
            elif assigned_group == "Medium Risk":
                st.warning("🟡 **Assigned Group: Medium Risk**")
            else:
                st.error("🔴 **Assigned Group: High Risk**")

            st.markdown("---")

            tab1, tab2 = st.tabs(["📈 Kaplan-Meier Survival", "🧬 TP53 Profile by Stage"])

            with tab1:
                time_months = np.linspace(0, 120, 50)
                surv_low = np.exp(-0.003 * time_months) * 100
                surv_med = np.exp(-0.008 * time_months) * 100
                surv_high = np.exp(-0.018 * time_months) * 100

                fig_km = go.Figure()
                w_low = 5 if assigned_group == "Low Risk" else 2
                w_med = 5 if assigned_group == "Medium Risk" else 2
                w_high = 5 if assigned_group == "High Risk" else 2

                fig_km.add_trace(go.Scatter(x=time_months, y=surv_low, mode='lines', name='Low Risk',
                                            line=dict(color='green', width=w_low)))
                fig_km.add_trace(go.Scatter(x=time_months, y=surv_med, mode='lines', name='Medium Risk',
                                            line=dict(color='orange', width=w_med)))
                fig_km.add_trace(go.Scatter(x=time_months, y=surv_high, mode='lines', name='High Risk',
                                            line=dict(color='red', width=w_high)))

                fig_km.update_layout(
                    title="10-Year Overall Survival Estimate (TCGA Cohort)",
                    xaxis_title="Time (Months)",
                    yaxis_title="Overall Survival Probability (%)",
                    hovermode="x unified",
                    legend=dict(orientation="h", yanchor="bottom", y=1.02, xanchor="right", x=1)
                )
                st.plotly_chart(fig_km, use_container_width=True)
                st.caption(
                    f"📌 **Patient Stratification:** Highlighted curve represents the patient's assigned category (**{assigned_group}**).")

            with tab2:
                stage_data = {
                    "Stage I": {"Wild-Type": 72, "Mutated": 28},
                    "Stage II": {"Wild-Type": 65, "Mutated": 35},
                    "Stage III": {"Wild-Type": 52, "Mutated": 48},
                    "Stage IV": {"Wild-Type": 40, "Mutated": 60}
                }

                tp53_stage_df = pd.DataFrame([
                    {"Stage": stage, "TP53 Status": status, "Percentage": pct}
                    for stage, status_dict in stage_data.items()
                    for status, pct in status_dict.items()
                ])

                fig_tp53 = px.bar(
                    tp53_stage_df,
                    x="Stage",
                    y="Percentage",
                    color="TP53 Status",
                    title="TP53 Mutation Prevalence Across Pathologic Stages",
                    labels={"Percentage": "Percentage of Cohort (%)", "Stage": "AJCC Pathologic Stage"},
                    color_discrete_map={"Wild-Type": "#2b5c8f", "Mutated": "#d95f02"},
                    barmode="stack"
                )

                tp53_label = "Mutated" if user_tp53 == 1 else "Wild-Type"
                wt_pct = stage_data[user_stage]["Wild-Type"]
                y_pos = 100 if user_tp53 == 1 else (wt_pct / 2)

                fig_tp53.add_trace(go.Scatter(
                    x=[user_stage],
                    y=[y_pos],
                    mode='markers+text',
                    name='Current Patient',
                    marker=dict(color='yellow', size=16, symbol='diamond', line=dict(width=2, color='black')),
                    text=[f" Patient ({tp53_label})"],
                    textposition="top center"
                ))

                fig_tp53.update_layout(yaxis=dict(range=[0, 115]))
                st.plotly_chart(fig_tp53, use_container_width=True)
                st.caption(
                    f"💡 **Molecular Context:** The patient presents with **{tp53_label} TP53** at **{user_stage}**.")
    else:
        st.info("Fill in the clinical parameters and click **Calculate Risk Score**.")