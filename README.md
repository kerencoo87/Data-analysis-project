
# TCGA BRCA Clinical Risk Assessment Tool

An integrated clinical-molecular decision support tool for Breast Invasive Carcinoma (BRCA) based on Empirical The Cancer Genome Atlas (TCGA) dataset analysis.

## Project Overview
This application provides clinicians and researchers with an interactive interface to stratify breast cancer patients into risk categories based on key clinical parameters and genetic alterations (specifically `TP53` mutation status). 
The tool integrates data processing, survival probability estimation via Kaplan-Meier curves, and molecular profiling across different AJCC pathologic stages to support data-driven insights.

## Live App
Access the deployed Streamlit web application here:
https://data-analysis-project-bfwrlacsqq6amv64hozgz3.streamlit.app/

## How to Use the Project

### 1. Patient Parameter Input
* **Patient Age (years):** Enter patient age (Valid range: 18–110 years).
* **TP53 Gene Status:** Select whether the `TP53` gene is Wild-Type or Mutated.
* **AJCC Pathologic Stage:** Choose the pathologic stage (Stage I-Stage IV).
* **Patient OS (months):** Input overall survival tracking duration in months (Valid range: 0–300 months).

### 2. Risk Calculation & Analytics
* Click **Calculate Risk Score** to perform risk stratification.
* View the assigned Risk Group (Low, Medium, or High Risk) and calculated Integrated Risk Score.
* Explore interactive tabs:
  * **Data-Driven Kaplan-Meier:** Visualize the 10-year overall survival probability curve relative to the TCGA cohort with the current patient highlighted.
  * **TP53 Profile by Stage:** Inspect molecular mutation prevalence across AJCC pathologic stages.

### Execution
Execute the main application file using Streamlit:
```bash
streamlit run app.py

