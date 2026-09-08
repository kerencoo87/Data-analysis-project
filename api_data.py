import streamlit as st
import requests
import pandas as pd

# Use Streamlit's data cache to avoid repeated API calls within the same browser session.
# This prevents potential 429 (Too Many Requests) rate-limiting blocks from cBioPortal.
@st.cache_data(show_spinner="Fetching data from cBioPortal API...")

def get_data():

    # Define the study we want to analyze and the API URL:
    study_id = "brca_tcga_pan_can_atlas_2018"
    base_url = f"https://www.cbioportal.org/api/studies/{study_id}"

    # Get the clinical data
    try:
        # PATIENT clinical data fetch
        res_patient = requests.get(
            f"{base_url}/clinical-data?clinicalDataType=PATIENT"
        )


        # Check that the API request was successful for other standard HTTP errors
        res_patient.raise_for_status()

    except requests.exceptions.HTTPError as e:
        # Handle API server errors or blocks without crashing the application layout
        st.error(f"⚠️ API Error: {e}")
        st.warning("Could not retrieve data from cBioPortal. The remote server might be down or rate-limiting requests.")
        st.stop()  # Stops execution of the Streamlit app to prevent runtime errors downstream
    except requests.exceptions.RequestException as e:
        # Handle connection, timeout, or network reachability failures
        st.error(f"⚠️ Network Error: Could not connect to cBioPortal API. ({e})")
        st.stop()

    # Convert the API responses from JSON
    patient_data = res_patient.json()


    # Convert the lists of dictionaries returned by the API
    # into pandas DataFrames (the data is still in long format)
    df_patient_long = pd.DataFrame(patient_data)

    # Reshape the data from long format to wide format,
    # so each row represents one patient and each clinical
    # attribute becomes a column
    df_patient = df_patient_long.pivot(
        index="patientId",
        columns="clinicalAttributeId",
        values="value"
    ).reset_index()

    # Choose the relevant columns from the data set
    patient_columns = [
        "patientId",
        "AGE",
        "SEX",
        "AJCC_PATHOLOGIC_TUMOR_STAGE",
        "OS_STATUS",
        "OS_MONTHS"
    ]

    df_final = df_patient[patient_columns].copy()

    # Get the molecular profiles
    """
    TP53 mutation information is not clinical data.
    It belongs to the molecular data section of cBioPortal.

    We chose TP53 because TP53 is a well-known and important
    tumor suppressor gene in breast cancer, and TP53 mutations
    are commonly observed in breast cancer.

    Therefore, we first ask the API which molecular profiles
    are available for this study.
    """
    try:
        res_profiles = requests.get(
            f"{base_url}/molecular-profiles"
        )

        # Check that the API request was successful
        res_profiles.raise_for_status()


    except requests.exceptions.RequestException as e:
        # Safe handling: Inform the user and stop execution to prevent incorrect risk score calculations
        st.error(f"⚠️ Molecular Data Error: {e}")
        st.warning("Could not retrieve molecular profiles from cBioPortal. TP53 data is missing, so risk calculation cannot proceed safely.")
        st.stop()

    # Convert the JSON response into a Python list of dictionaries
    profiles_data = res_profiles.json()

    # Find the mutation profile
    """
    A study can contain different types of molecular data.

    After inspecting the molecular profiles above, we look for
    the profile whose molecularAlterationType is MUTATION_EXTENDED.

    This is the molecular profile containing mutation data.
    """

    mut_profile_id = None

    for profile in profiles_data:

        if profile.get("molecularAlterationType") == "MUTATION_EXTENDED":

            mut_profile_id = profile.get("molecularProfileId")

            break

    # Get TP53 mutation data
    """
    TP53 has Entrez Gene ID = 7157, according to NCBI Gene.

    We use this ID to request TP53 mutation records
    from the cBioPortal API.
    """

    # Will later store the unique patient IDs of patients
    # with a TP53 mutation
    tp53_patients = set()

    # Check whether we found the mutation profile
    if mut_profile_id:

        # Build the mutation API URL
        mut_url = (
            f"https://www.cbioportal.org/api/molecular-profiles/"
            f"{mut_profile_id}/mutations/fetch"
        )

        # Create the request payload
        mut_payload = {
            "sampleListId": f"{study_id}_all",
            # Search across all samples in this study

            "entrezGeneIds": [7157]
            # Request mutation data for TP53
            # 7157 is the Entrez Gene ID for TP53
        }

        try:
            # Send the request to the API
            # The payload tells the API what mutation data we want
            res_mut = requests.post(
                mut_url,
                json=mut_payload
            )

            # Check that the request succeeded
            res_mut.raise_for_status()

            # Convert the response from JSON
            mutation_data = res_mut.json()

            # Convert TP53 mutation data to a DataFrame
            df_mut = pd.DataFrame(mutation_data)

            # Check the number of unique patients with TP53 mutation
            """
            A patient can potentially have more than one mutation
            in the TP53 gene.
    
            Therefore, we count unique patient IDs rather than
            counting mutation records.
            """

            # If patientId exists, store the unique patient IDs
            # of patients with a TP53 mutation in a set
            if "patientId" in df_mut.columns:

                tp53_patients = set(
                    df_mut["patientId"].dropna().unique()
                )

        except requests.exceptions.RequestException:
            st.warning("Failed to fetch specific TP53 mutation records.")

    # Add TP53 mutation data status to the clinical data
    """
    We create a binary variable:

        1 - TP53 mutation was found
        0 - TP53 mutation was not found

    We use patientId to match the mutation information
    to our patient-level clinical DataFrame.
    """

    df_final["TP53_MUTATION"] = (
        df_final["patientId"]
        .isin(tp53_patients)
        # Check whether each patient ID is in tp53_patients
        # (return boolean value)

        .astype(int)
        # Convert True/False to 1/0
    )

    # The API returns clinical values as strings.

    # Convert AGE and OS_MONTHS from strings to numeric values
    df_final[["AGE", "OS_MONTHS"]] = (
        df_final[["AGE", "OS_MONTHS"]]
        .apply(
            pd.to_numeric,
            errors="coerce"
        )
        # If a value cannot be converted,
        # treat it as missing (NaN)
    )

    # Because the number of male patients is very small (1.1%),
    # we will restrict the analysis to female patients.

    # Keep only female patients
    df_final = df_final[
        df_final["SEX"] == "Female"
    ].copy()

    # Check missing values and their percentage
    """
    Five patients (0.47%) have missing values for
    AJCC_PATHOLOGIC_TUMOR_STAGE.

    Since this represents a very small proportion of the dataset,
    these patients will be excluded from the analysis.
    """

    # Remove patients with missing AJCC pathological tumor stage
    df_final = df_final.dropna(
        subset=["AJCC_PATHOLOGIC_TUMOR_STAGE"]
    )

    # Rename the OS status to number only
    status_map = {
        "1:DECEASED": 1,
        "0:LIVING": 0,
        "DECEASED": 1,
        "LIVING": 0,
    }

    df_final["IS_DECEASED"] = (
        df_final["OS_STATUS"].map(status_map)
    )

    # Return the final cleaned DataFrame
    return df_final