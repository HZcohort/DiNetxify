# Data harmonization

Data harmonization in ***DiNetxify*** means loading phenotype data and medical-record data into a single `DiseaseNetworkData` object. During this process, diagnosis codes are mapped to phecodes, dates are standardized, follow-up windows are applied, and phenotype covariates are converted into analysis-ready variables.

## Initializing the data object

First, import ***DiNetxify*** and create a `DiseaseNetworkData` object with your study design and phecode settings.

```python
import DiNetxify as dnt

# Matched cohort
data = dnt.DiseaseNetworkData(
    study_design="matched cohort",
    phecode_level=1,
)

# Standard cohort
data = dnt.DiseaseNetworkData(
    study_design="cohort",
    phecode_level=1,
)

# Exposed-only cohort
data = dnt.DiseaseNetworkData(
    study_design="exposed-only cohort",
    phecode_level=1,
)
```

- **study_design**: One of `'cohort'`, `'matched cohort'`, or `'exposed-only cohort'`. Default is `'cohort'`.
- **phecode_level**: Use `1` or `2`. Level 1 is broader and usually more stable for smaller datasets; level 2 is more granular and often more suitable for larger cohorts.

**Optional parameters:**

- **min_required_icd_codes**: Minimum number of ICD records mapping to the same phecode before that phecode is counted as present for an individual. Default is `1`.
- **date_fmt**: Date format used in the phenotype file. This is also the default format for medical-record files unless overridden later. Default is `'%Y-%m-%d'`.
- **phecode_version**: Phecode mapping version. `'1.2'` is the recommended general-purpose choice. `'1.3a'` is also available for special use cases.

## Load phenotype data

After creating the data object, use `phenotype_data()` to load the cohort phenotype data.

```python
# Matched cohort
col_dict = {
    "Participant ID": "ID",
    "Exposure": "exposure",
    "Sex": "sex",
    "Index date": "date_start",
    "End date": "date_end",
    "Match ID": "group_id",
}
covariates = ["age", "BMI"]
data.phenotype_data(
    phenotype_data_path="tests/data/dummy_phenotype.csv",
    column_names=col_dict,
    covariates=covariates,
)

# Standard cohort
col_dict = {
    "Participant ID": "ID",
    "Exposure": "exposure",
    "Sex": "sex",
    "Index date": "date_start",
    "End date": "date_end",
}
data.phenotype_data(
    phenotype_data_path="tests/data/dummy_phenotype.csv",
    column_names=col_dict,
    covariates=["age", "BMI"],
)

# Exposed-only cohort
col_dict = {
    "Participant ID": "ID",
    "Sex": "sex",
    "Index date": "date_start",
    "End date": "date_end",
}
data.phenotype_data(
    phenotype_data_path="tests/data/dummy_phenotype.csv",
    column_names=col_dict,
    covariates=["age", "BMI"],
)
```

- **phenotype_data_path**: Path to a CSV or TSV phenotype file.
- **column_names**: Dictionary mapping required DiNetxify field names to the corresponding column names in your dataset.
- **covariates**: List of additional phenotype variables to load.

**Optional parameters:**

- **is_single_sex**: Set to `True` if the cohort contains only one sex. Default is `False`.
- **force**: If `True`, overwrite phenotype and medical-record information already stored in the object. Default is `False`.

**Important phenotype input rules:**

- For `cohort` and `matched cohort`, `Exposure` must be coded as `1` for exposed and `0` for unexposed.
- `Sex` must be coded as `1` for female and `0` for male.
- Required columns cannot contain missing values.
- Continuous covariates with missing values will remove those participants during loading.
- Categorical covariates with missing values are retained and treated as an `"NA"` category.
- Column names used for covariates must not conflict with DiNetxify reserved variables.

**After loading phenotype data:**

You can inspect the cohort summary by printing the object:

```python
print(data)
```

You can also generate a phenotype summary table with:

```python
table1_df = data.Table1()
print(table1_df.head())
```

`Table1()` summarizes sex, covariates included and follow-up time. For `cohort` and `matched cohort`, the table compares exposed and unexposed groups; for `exposed-only cohort`, it produces a single-group summary.

## Load medical record data

After the phenotype data are loaded, merge one or more diagnosis files with `merge_medical_records()`. Call the method once per file if your records use different ICD systems.

```python
data.merge_medical_records(
    medical_records_data_path="tests/data/dummy_EHR_ICD10.csv",
    diagnosis_code="ICD-10-WHO",
    column_names={
        "Participant ID": "ID",
        "Diagnosis code": "diag_icd10",
        "Date of diagnosis": "dia_date",
    },
)

data.merge_medical_records(
    medical_records_data_path="tests/data/dummy_EHR_ICD9.csv",
    diagnosis_code="ICD-9-WHO",
    column_names={
        "Participant ID": "ID",
        "Diagnosis code": "diag_icd9",
        "Date of diagnosis": "dia_date",
    },
    diagnosis_code_exclusion=[],
)
```

- **medical_records_data_path**: Path to a CSV or TSV diagnosis file.
- **diagnosis_code**: One ICD system per file. Supported values are `'ICD-9-CM'`, `'ICD-9-WHO'`, `'ICD-10-CM'`, and `'ICD-10-WHO'`.
- **column_names**: Dictionary mapping `'Participant ID'`, `'Diagnosis code'`, and `'Date of diagnosis'` to the actual columns in your file.

**Optional parameters:**

- **date_fmt**: Date format for the diagnosis file. If `None`, DiNetxify uses the `date_fmt` defined in the `DiseaseNetworkData` object.
- **chunksize**: Number of rows processed per chunk. Default is `1000000`.
- **diagnosis_code_exclusion**: List of diagnosis codes to exclude before phecode mapping.

**Important medical-record input rules:**

- Use the same participant IDs as the phenotype file.
- Do not mix ICD-9 and ICD-10 codes in one file.
- Do not mix CM and WHO coding systems in one file.
- Do not restrict records to first occurrences only.
- Do not pre-filter records to the follow-up period; DiNetxify handles follow-up filtering internally.

**During loading, the package reports:**

- how many rows were read
- how many rows were excluded because of missing values or ID/code filtering
- how many diagnosis codes mapped directly or after truncation
- how many codes were not mapped to any phecode
- how many mapped records were invalid because the diagnosis date fell outside the usable follow-up window

**After loading medical record data:**

Print the object again to see the merged-data summary:

```python
print(data)
```

This summary includes the number of merged files, the total number of processed diagnosis records, mean numbers of recorded phecodes during and before follow-up, and any mapping or follow-up warnings collected during harmonization.

## Save DiseaseNetworkData object

Once phenotype data and medical records have been harmonized, you can save the `DiseaseNetworkData` object for later reuse.

```python
# Save as gzip-compressed pickle
data.save("results/cohort_data")

# Save as compressed NumPy archive
data.save_npz("results/cohort_data")
```

You do not need to add the extension yourself. `save()` appends `.pkl.gz`, and `save_npz()` writes a `.npz` file.

## Reload DiseaseNetworkData object

To reload a saved object, create a new `DiseaseNetworkData` instance and call `load()` or `load_npz()`.

```python
import DiNetxify as dnt

data = dnt.DiseaseNetworkData()

# Load from gzip-compressed pickle
data.load("results/cohort_data")

# Or load from compressed NumPy archive
data.load_npz("results/cohort_data")
```

Use `force=True` only if you are loading into an object that already contains data and you want to overwrite it.
