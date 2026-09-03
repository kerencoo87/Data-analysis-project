import pandas as pd


# ---------------------------------------------------------
# 1. Mortality according to age
# ---------------------------------------------------------

def mortality_by_age(df):

    # Analysis of the deceased percentage according to age.
    # This will help us understand whether age should be
    # included as a risk factor.

    age_stats = (
        df.groupby("AGE_60_PLUS")["IS_DECEASED"]
        .agg(["count", "mean"])
        .round(2)
    )

    # Convert the mean mortality value to a percentage
    age_stats["mortality %"] = (
        age_stats["mean"] * 100
    ).round(2)

    return age_stats


# ---------------------------------------------------------
# 2. Mortality according to TP53 mutation
# ---------------------------------------------------------

def mortality_by_tp53(df):

    # Analysis of the deceased percentage according to
    # TP53 mutation status.

    tp53_stats = (
        df.groupby("TP53_MUTATION")["IS_DECEASED"]
        .agg(["count", "mean"])
        .round(2)
    )

    # Convert the mean mortality value to a percentage
    tp53_stats["mortality %"] = (
        tp53_stats["mean"] * 100
    ).round(2)

    return tp53_stats


# ---------------------------------------------------------
# 3. Mortality according to tumor stage
# ---------------------------------------------------------

def mortality_by_stage(df):

    # Analysis of the deceased percentage according to
    # stage status.

    stage_stats = (
        df.groupby("STAGE_CATEGORY")["IS_DECEASED"]
        .agg(["count", "mean"])
        .round(2)
    )

    # Convert the mean mortality value to a percentage
    stage_stats["mortality %"] = (
        stage_stats["mean"] * 100
    ).round(2)

    return stage_stats


# ---------------------------------------------------------
# 4. Mortality according to risk group
# ---------------------------------------------------------

def mortality_by_risk_group(df):

    # Analysis of the mortality percentage according
    # to the risk groups.

    risk_summary = (
        df.groupby("RISK_GROUP")["IS_DECEASED"]
        .agg(["count", "mean"])
    )

    # Convert the mean mortality value to a percentage
    risk_summary["mortality_pct"] = (
        risk_summary["mean"] * 100
    ).round(2)

    risk_summary["mean"] = (
        risk_summary["mean"].round(2)
    )

    return risk_summary