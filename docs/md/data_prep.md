# Input data preparation

## Requirements for input data

***DiNetxify*** enables 3D disease network analysis on cohort data from electronic health records (EHR) and offers three study designs: the **standard cohort**, which compares individuals with a specific disease or exposure (e.g., depression or smoking) against the general population; the **matched cohort**, which pairs subjects on key characteristics to reduce confounding; and the **exposed-only cohort**, which examines disease networks within a defined group (e.g., older adults) without an unexposed comparison group.

To begin using ***DiNetxify***, two datasets are required: a **phenotype data** file containing each participant’s baseline information, and one or more **medical record data** files extracted from an EHR database that list diagnoses (codes and dates) for all cohort individuals over the study period. The specific requirements for these datasets are:

- **Phenotype data:** A CSV (or TSV) file with a header row, where each row represents a participant. Required columns are:

  - **Participant ID** – Unique identifier for each individual.

  - **Index date** – Start of follow-up (e.g., date of exposure or baseline).

  - **End date** – End of follow-up (e.g., last visit, death, or study completion).

  - **Exposure** – Binary indicator (1 = exposed, 0 = unexposed) for standard and matched cohorts (omit for exposed-only cohorts).

  - **Match ID** – Identifier for matched sets (only for matched cohort designs).

  - **Sex** – Biological sex (1 = female, 0 = male).

  - **Additional covariates (optional)** – Any number of extra variables for adjustment or matching (e.g., age, BMI, education).

For all required columns, missing values are not permitted and dates must follow a consistent format (default `YYYY-MM-DD`). The **Sex** and **Exposure** fields must use the specified 1/0 coding. You may include unlimited additional covariates; their types (binary, categorical, or continuous) will be auto-detected and processed accordingly (e.g., one-hot encoding for categorical variables). Missing values in continuous covariates are dropped, whereas missing categorical values are treated as a separate “NA” category.

- **Medical records data**: One or more CSV/TSV files (each with a header row), listing diagnosis events for the participants. Each record (row) must include:

  - **Participant ID** – The same unique ID used in the phenotype data, linking each record to an individual.
  - **Diagnosis code** – A diagnosis code (e.g., ICD-10 or ICD-9 code).
  - **Date of diagnosis** – The date of that diagnosis/event (format consistent with the phenotype dates, e.g., `YYYY-MM-DD`).

  The **medical record data** should be in a “long” format (multiple rows per participant if they have multiple diagnoses). ***DiNetxify*** will automatically filter these records to include only those within each individual’s follow-up period (from index date up to end date). **Do not pre-filter** the medical record by date or by first occurrence — provide the complete set of diagnoses for each participant, and let the software handle the filtering and mapping. Each medical record file should use a single coding system for diagnoses. Currently supported code versions are ICD-9 (WHO and CM) and ICD-10 (WHO and CM). If your data uses a different coding system, you will need to map it to one of the supported formats beforehand.

## Dummy dataset overview

A [dummy dataset](https://github.com/HZcohort/DiNetxify/tree/main/tests/data) is provided to help you become familiar with the input format and to allow you to run through the full analysis workflow before using your own data. It simulates a matched-cohort study of 10,000 exposed individuals and 50,000 matched unexposed individuals, along with their entire follow-up EHR records.

> **Note:** All participant characteristics and diagnoses in this dummy dataset are randomly generated. The ICD-9 and ICD-10 codes correspond to real classifications, and the analysis may yield seemingly significant associations, but these results **do not** reflect true medical findings! They are for instructional purposes only.

- The dummy dataset consists of three CSV files:
  - **`dummy_phenotype.csv`** – Simulated baseline characteristics for 60,000 individuals, containing:
    
    - **ID** – Unique participant identifier.
    - **date_start**, **date_end** – Follow-up start and end dates.
    - **exposure** – Exposure status (0 = unexposed, 1 = exposed).
    - **group_id** – Matching group identifier (each exposed is matched with unexposed in groups).
    - **sex** – Biological sex (1 = female, 0 = male).
    - **age** – Baseline age (years).
    - **BMI** – Body mass index category.
  - **`dummy_EHR_ICD9.csv`** – Simulated EHR diagnoses coded in ICD-9 (10,188 records), containing:
  
      - **ID** – Participant ID (matches the phenotype file).
  
      - **dia_date** – Diagnosis date.
  
      - **diag_icd9** – ICD-9 diagnosis code.
  
  
  - **`dummy_EHR_ICD10.csv`** – Simulated EHR diagnoses coded in ICD-10 (1,048,576 records), containing:
  
      - **ID** – Participant ID.
  
      - **dia_date** – Diagnosis date.
  
      - **diag_icd10** – ICD-10 diagnosis code.
  

Using this dummy dataset, you can practice the workflow and verify that the tool runs correctly. In the following sections, we will demonstrate the analysis steps using the dummy data.