import streamlit as st
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
from lifelines import KaplanMeierFitter

from api_data import get_data
from data_processing import calculate_single_patient_risk, calculate_risk_groups

# 1. Load data and assign stage column
df_final = get_data()
df_final['STAGE_CATEGORY'] = df_final['AJCC_PATHOLOGIC_TUMOR_STAGE']

# 2. Calculate risk groups for the full cohort
df_k_m = calculate_risk_groups(df_final)

# Streamlit Configuration
st.set_page_config(page_title="BRCA Clinical Risk Tool", layout="wide")
st.title("🧬 TCGA BRCA Clinical Risk Assessment Tool")
st.write("An integrated clinical-molecular decision support tool for breast invasive carcinoma.")
st.markdown("---")

# Initialize session state to track calculation status
if 'calculated' not in st.session_state:
    st.session_state['calculated'] = False

def reset_calculation():
    """Resets the calculation state when any input parameter changes."""
    st.session_state['calculated'] = False

col_input, col_result = st.columns([1, 1], gap="medium")

with col_input:
    st.subheader("📋 Patient Clinical Parameters")
    user_age = st.number_input("Patient Age (years):", min_value=0, value=60, step=1, on_change=reset_calculation)
    tp53_radio = st.radio("TP53 Gene Status:", options=["Wild-Type (No Mutation)", "Mutated"], horizontal=True, on_change=reset_calculation)
    user_tp53 = 1 if "Mutated" in tp53_radio else 0
    user_stage = st.selectbox("AJCC Pathologic Stage:", options=["Stage I", "Stage II", "Stage III", "Stage IV"], on_change=reset_calculation)
    user_survival = st.number_input("Patient OS (months):", min_value=0, value=0, step=1, on_change=reset_calculation)
    btn_calc = st.button("Calculate Risk Score", type="primary", use_container_width=True)
    if btn_calc:
            st.session_state['calculated'] = True
with col_result:
    st.subheader("📊 Stratification & Analytics")
    is_age_valid = (user_age is not None) and (18 <= user_age <= 110)
    is_survival_valid = (user_survival is not None) and (0 <= user_survival <= 300)
    # 2. Render warning if invalid
    if not is_age_valid:
        st.error("❌ **Invalid Input:** Patient Age must be between 18 and 110 years. Previous results cleared.")
        st.info("Please correct the age parameter and click **Calculate Risk Score**.")
    elif not is_survival_valid:
        st.error(
            "❌ **Invalid Input:** Overall Survival time must be between 0 and 360 months. Previous results cleared.")
        st.info("Please correct the survival time parameter and click **Calculate Risk Score**.")

    # 3. Render results ONLY if inputs are valid AND button was clicked
    elif st.session_state.get('calculated', False):
    #if btn_calc:
        res = calculate_single_patient_risk(user_age, user_tp53, user_stage)

        if "error" in res:
            st.error(res["error"])
            assigned_group = None
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


    # Swapped tab order
        tab1, tab2 = st.tabs([
                "📊 Data-Driven Kaplan-Meier",
            "🧬 TP53 Profile by Stage"
        ])


        # ---------------------------------------------------------
        # TAB 1: Empirical TCGA Data-Based Kaplan-Meier
        # ---------------------------------------------------------
        with tab1:
            kmf = KaplanMeierFitter()
            fig_km = go.Figure()

            color_map = {
                "Low Risk": "green",
                "Medium Risk": "orange",
                "High Risk": "red"
            }

            patient_survival_prob = None

            for group in ["Low Risk", "Medium Risk", "High Risk"]:
                group_df = df_k_m[df_k_m["RISK_GROUP"] == group]

                if group_df.empty:
                    continue

                kmf.fit(
                    durations=group_df["OS_MONTHS"],
                    event_observed=group_df["IS_DECEASED"],
                    label=group
                )

                km_timeline = kmf.survival_function_.index.values
                km_survival = kmf.survival_function_.iloc[:, 0].values * 100

                line_width = 5 if group == assigned_group else 2
                line_color = color_map.get(group, "blue")

                fig_km.add_trace(go.Scatter(
                    x=km_timeline,
                    y=km_survival,
                    mode='lines',
                    line=dict(color=line_color, width=line_width, shape='hv'),
                    name=f"{group} (n={len(group_df)})"
                ))

                if group == assigned_group:
                    patient_survival_prob = kmf.predict(user_survival) * 100

            if patient_survival_prob is not None:
                fig_km.add_trace(go.Scatter(
                    x=[user_survival],
                    y=[patient_survival_prob],
                    mode='markers+text',
                    name='Current Patient',
                    marker=dict(color='purple', size=14, symbol='diamond', line=dict(width=2, color='white')),
                    text=[f" Patient ({user_survival}m, {patient_survival_prob:.1f}%)"],
                    textposition="top right"
                ))

            fig_km.update_layout(
                title="10-Year Overall Survival Estimate (TCGA Empirical Data)",
                xaxis_title="Time (Months)",
                yaxis_title="Overall Survival Probability (%)",
                hovermode="x unified",
                legend=dict(orientation="h", yanchor="bottom", y=1.02, xanchor="right", x=1)
            )

            st.plotly_chart(fig_km, use_container_width=True)
            if assigned_group:
                st.caption(
                    f"📌 **Empirical Cohort:** Highlighted curve represents **{assigned_group}** based on direct TCGA-BRCA observations."
                )

        # ---------------------------------------------------------
        # TAB 2: TP53 Profile by Stage
        # ---------------------------------------------------------
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

            # if btn_calc and assigned_group:
            #     tp53_label = "Mutated" if user_tp53 == 1 else "Wild-Type"
            #     wt_pct = stage_data[user_stage]["Wild-Type"]
            #     y_pos = 100 if user_tp53 == 1 else (wt_pct / 2)
            #
            #     fig_tp53.add_trace(go.Scatter(
            #         x=[user_stage],
            #         y=[y_pos],
            #         mode='markers+text',
            #         name='Current Patient',
            #         # marker=dict(color='yellow', size=16, symbol='diamond', line=dict(width=2, color='black')),
            #         text=[f" Patient ({tp53_label})"],
            #         textposition="top center"
            #     ))

            fig_tp53.update_layout(yaxis=dict(range=[0, 115]))
            st.plotly_chart(fig_tp53, use_container_width=True)
            if btn_calc and assigned_group:
                st.caption(
                    f"💡 **Molecular Context:** The patient presents with **{tp53_label} TP53** at **{user_stage}**."
                )
    else:
        assigned_group = None
        st.info("Fill in the clinical parameters and click **Calculate Risk Score**.")