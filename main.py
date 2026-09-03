from api_data import get_data

from data_processing import (
    add_stage_category,
    add_age_category,
    calculate_risk_groups
)

from analysis import (
    mortality_by_age,
    mortality_by_tp53,
    mortality_by_stage,
    mortality_by_risk_group
)


def main():

    # ---------------------------------------------------------
    # 1. Get the data from the API
    # ---------------------------------------------------------

    df_final = get_data()

    print("Initial df_final:")
    print(df_final.head())


    # ---------------------------------------------------------
    # 2. Data processing
    # ---------------------------------------------------------

    # Categorize patients according to tumor stage
    df_final = add_stage_category(df_final)

    # Create age category
    df_final = add_age_category(df_final)


    # ---------------------------------------------------------
    # 3. Analyze mortality according to individual
    #    risk factors
    # ---------------------------------------------------------

    age_stats = mortality_by_age(df_final)

    print("\nDeceased percentage according to age:")
    print(age_stats)


    tp53_stats = mortality_by_tp53(df_final)

    print("\nDeceased percentage according to TP53 mutation:")
    print(tp53_stats)


    stage_stats = mortality_by_stage(df_final)

    print("\nDeceased percentage according to stage:")
    print(stage_stats)


    # ---------------------------------------------------------
    # 4. Calculate risk score and create risk groups
    # ---------------------------------------------------------

    df_risk = calculate_risk_groups(df_final)


    # ---------------------------------------------------------
    # 5. Analyze mortality according to risk group
    # ---------------------------------------------------------

    risk_summary = mortality_by_risk_group(df_risk)

    print("\nSummary of mortality % for each risk group:")
    print(risk_summary)


if __name__ == "__main__":
    main()

