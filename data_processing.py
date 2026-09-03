import pandas as pd
import numpy as np

# ---------------------------------------------------------
# 1. Categorize tumor stages
# ---------------------------------------------------------

# Function that categorizes the stages
def categorize_stage_fixed(val):

    if pd.isna(val):
        return "Unknown"

    val_str = str(val).upper().strip()

    # 1. Stage IV
    if "STAGE IV" in val_str or "T4" in val_str:
        return "Advanced (Stage IV/T4)"

    # 2. Stage III
    elif "STAGE III" in val_str or "T3" in val_str:
        return "Advanced (Stage III/T3)"

    # 3. Stage II
    elif "STAGE II" in val_str or "T2" in val_str:
        return "Intermediate (Stage II/T2)"

    # 4. Stage I
    elif "STAGE I" in val_str or "T1" in val_str:
        return "Early (Stage I/T1)"

    else:
        return "Other/Unknown"

# ---------------------------------------------------------
# 2. Add stage category to the DataFrame
# ---------------------------------------------------------

def add_stage_category(df):

    # Using the function for the categories
    df = df.copy()

    df["STAGE_CATEGORY"] = (
        df["AJCC_PATHOLOGIC_TUMOR_STAGE"]
        .map(categorize_stage_fixed)
    )

    return df


# ---------------------------------------------------------
# 3. Create age category
# ---------------------------------------------------------

def add_age_category(df):

    df = df.copy()

    # Patients aged 60 or above receive the value 1.
    # Patients below 60 receive the value 0.
    df["AGE_60_PLUS"] = (
        df["AGE"] >= 60
    ).astype(int)

    return df


# ---------------------------------------------------------
# 4. Calculate risk score and risk groups
# ---------------------------------------------------------

def calculate_risk_groups(df):

    # A copy of the DataFrame for risk analysis
    df_risk = df.copy()

    # Calculation of the risk score by age:
    # age above or equal to 60 gets 1 point
    # age below 60 gets 0 points
    age_points = np.where(
        df_risk["AGE"] >= 60,
        1,
        0
    )

    # Calculation of the risk score by mutation:
    # presence of TP53 mutation gets 2 points
    tp53_points = np.where(
        df_risk["TP53_MUTATION"] == 1,
        2,
        0
    )

    # Calculation of the risk score by STAGE_CATEGORY
    #
    # Stage I = 0
    # Stage II = 1
    # Stage III/IV = 3

    stage_conditions = [

        df_risk["STAGE_CATEGORY"].str.contains(
            "Stage I/T1",
            na=False
        ),

        df_risk["STAGE_CATEGORY"].str.contains(
            "Stage II/T2",
            na=False
        ),

        df_risk["STAGE_CATEGORY"].str.contains(
            "Stage III|Stage IV",
            na=False
        ),
    ]

    stage_choices = [0, 1, 3]

    stage_points = np.select(
        stage_conditions,
        stage_choices,
        default=0
    )

    # Calculation of the total risk score
    df_risk["RISK_SCORE"] = (
        age_points
        + tp53_points
        + stage_points
    )

    # Split groups according to risk score:
    #
    # 0-1 = Low risk
    # 2-3 = Medium risk
    # 4-6 = High risk

    group_conditions = [

        df_risk["RISK_SCORE"] <= 1,

        (df_risk["RISK_SCORE"] >= 2)
        & (df_risk["RISK_SCORE"] <= 3),

        df_risk["RISK_SCORE"] >= 4,
    ]

    group_choices = [
        "Low Risk",
        "Medium Risk",
        "High Risk"
    ]

    df_risk["RISK_GROUP"] = np.select(
        group_conditions,
        group_choices,
        default="Unknown"
    )

    return df_risk


def calculate_single_patient_risk(age, tp53, stage):
    # Column names aligned with calculate_risk_groups expectations
    patient_df = pd.DataFrame([{
        "AGE": float(age),
        "TP53_MUTATION": int(tp53),
        "AJCC_PATHOLOGIC_TUMOR_STAGE": str(stage)
    }])
    patient_df.columns = patient_df.columns.str.strip()
    # Calculate Stage Category and Risk
    patient_df = add_stage_category(patient_df)
    result_df = calculate_risk_groups(patient_df)

    # Return result dictionary for the single patient
    return {
        "risk_score": int(result_df["RISK_SCORE"].iloc[0]),
        "risk_group": str(result_df["RISK_GROUP"].iloc[0])
    }