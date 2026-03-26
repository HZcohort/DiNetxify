# Input data preparation

## Overview

`DiNetxify` requires two types of input:

1. A phenotype file describing the cohort
2. One or more medical-record files containing diagnosis codes and diagnosis dates

The package supports three study designs:

- `cohort`
- `matched cohort`
- `exposed-only cohort`

Before loading any data, create a `DiseaseNetworkData` object.

```python
import DiNetxify as dnt

data = dnt.DiseaseNetworkData(
    study_design="cohort",
    phecode_level=1,
    min_required_icd_codes=1,
    date_fmt="%Y-%m-%d",
    phecode_version="1.2",
)
```

Key initialization arguments:

- `study_design`: One of `'cohort'`, `'matched cohort'`, or `'exposed-only cohort'`.
- `phecode_level`: Use `1` or `2`.
- `min_required_icd_codes`: Minimum number of mapped ICD codes required for a phecode to count as valid.
- `date_fmt`: Date format used in the phenotype file. This is also the default for medical-record files unless you override it later.
- `phecode_version`: Version `1.2` is the recommended general-purpose option.

## Phenotype data

Phenotype data must be provided as a CSV or TSV file with a header row and one row per participant.

### Required columns

For `cohort`:

- `Participant ID`
- `Exposure`
- `Sex`
- `Index date`
- `End date`

For `matched cohort`:

- `Participant ID`
- `Exposure`
- `Sex`
- `Index date`
- `End date`
- `Match ID`

For `exposed-only cohort`:

- `Participant ID`
- `Sex`
- `Index date`
- `End date`

You may also provide any number of additional covariates, such as age, BMI, smoking, or education.

### Input rules

- Required columns cannot contain missing values.
- `Exposure` must be coded as `1` for exposed and `0` for unexposed when the study design includes an exposure group.
- `Sex` must be coded as `1` for female and `0` for male.
- Dates must use a consistent format.
- Covariate types are detected automatically and converted internally.
- Continuous covariates with missing values lead to participant removal during loading.
- Covariate names must not conflict with reserved internal variable names.

### Load phenotype data

```python
col_dict = {
    "Participant ID": "ID",
    "Exposure": "exposure",
    "Sex": "sex",
    "Index date": "date_start",
    "End date": "date_end",
}

covariates = ["age", "BMI"]

data.phenotype_data(
    phenotype_data_path="tests/data/dummy_phenotype.csv",
    column_names=col_dict,
    covariates=covariates,
    is_single_sex=False,
    force=False,
)
```

Notes:

- For a matched cohort, include `"Match ID"` in `column_names`.
- For an exposed-only cohort, omit `"Exposure"` from `column_names`.
- If your cohort contains only one sex, set `is_single_sex=True`.
- If you want to overwrite already loaded data, set `force=True`.

## Medical records data

Medical records must be provided in long format, with one diagnosis event per row.

### Required columns

Each medical-record file must contain:

- `Participant ID`
- `Diagnosis code`
- `Date of diagnosis`

### Supported diagnosis code systems

Each file must use exactly one supported coding system:

- `ICD-9-CM`
- `ICD-9-WHO`
- `ICD-10-CM`
- `ICD-10-WHO`

If you have multiple coding systems, load them as separate files by calling `merge_medical_records()` multiple times.

### Input rules

- Use the same participant IDs as in the phenotype file.
- Do not mix ICD-9 and ICD-10 codes in one file.
- Do not mix CM and WHO code systems in one file.
- Do not pre-filter diagnoses to first occurrence only.
- Do not pre-filter diagnoses to the follow-up period; `DiNetxify` handles follow-up filtering internally.
- ICD-10 codes may be provided with or without decimal points.
- ICD-9 codes may be provided in decimal or short format.

### Load medical records

```python
data.merge_medical_records(
    medical_records_data_path="tests/data/dummy_EHR_ICD9.csv",
    diagnosis_code="ICD-9-WHO",
    column_names={
        "Participant ID": "ID",
        "Diagnosis code": "diag_icd9",
        "Date of diagnosis": "dia_date",
    },
    date_fmt=None,
    chunksize=1000000,
    diagnosis_code_exclusion=[],
)

data.merge_medical_records(
    medical_records_data_path="tests/data/dummy_EHR_ICD10.csv",
    diagnosis_code="ICD-10-WHO",
    column_names={
        "Participant ID": "ID",
        "Diagnosis code": "diag_icd10",
        "Date of diagnosis": "dia_date",
    },
)
```

Notes:

- If `date_fmt=None`, the medical-record file uses the same date format as the phenotype file.
- `chunksize` controls how many rows are processed at a time and is useful for large files.
- `diagnosis_code_exclusion` lets you exclude specific diagnosis codes before mapping.
- After each merge, the package updates diagnosis, diagnosis-count, and history information inside the `DiseaseNetworkData` object.

## Dummy dataset

A dummy dataset is included under [tests/data](https://github.com/HZcohort/DiNetxify/tree/main/tests/data) so you can test the full workflow before using your own data.

The dataset contains:

- `dummy_phenotype.csv`: 60,000 participants in a matched-cohort-style example
- `dummy_EHR_ICD9.csv`: 10,188 ICD-9 diagnosis records
- `dummy_EHR_ICD10.csv`: 1,668,795 ICD-10 diagnosis records

Important columns in the dummy phenotype file:

- `ID`: Participant identifier
- `group_id`: Matching group identifier
- `exposure`: Exposure status
- `date_start`: Follow-up start date
- `date_end`: Follow-up end date
- `age`: Baseline age
- `sex`: Biological sex
- `BMI`: BMI category

Important columns in the dummy medical-record files:

- `ID`: Participant identifier
- `dia_date`: Diagnosis date
- `diag_icd9` or `diag_icd10`: Diagnosis code

> **Note:** The dummy data are simulated for demonstration only. They are useful for learning the workflow and checking that the software runs, but they should not be interpreted as real clinical findings.
