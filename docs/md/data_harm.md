# Data harmonization

Data harmonization involves loading and merging the **phenotype data** and **medical record data** into a single `DiseaseNetworkData` object for analysis. During this process, the software ensures consistent coding (e.g., converting diagnosis codes to phecodes) and standardized formatting (e.g., datetime parsing for diagnosis and follow-up periods).

## Initializing the data object

First, import the ***DiNetxify*** package and instantiate a `DiseaseNetworkData` object with your chosen study design, phecode level, and any optional parameters. For example:

```python
import DiNetxify as dnt  

# For a matched cohort study  
data = dnt.DiseaseNetworkData(  
    study_design='matched cohort',  
    phecode_level=1,  
)  

# For a standard cohort study  
data = dnt.DiseaseNetworkData(  
    study_design='cohort',  
    phecode_level=1,  
)  

# For an exposed-only cohort study  
data = dnt.DiseaseNetworkData(  
    study_design='exposed-only cohort',  
    phecode_level=1,  
)  
```

- **study_design** – Type of study design. Options: `'cohort'`, `'matched cohort'`, or `'exposed-only cohort'`. *(Default: 'cohort')*.
- **phecode_level** – Level of phecode to use for grouping diagnoses. Level 1 provides broader categories (~585 conditions) while level 2 offers more detailed categories (~1257 conditions). For smaller datasets, level 1 is recommended to maintain statistical power; for larger datasets, level 2 can provide finer granularity. *(Options: 1 or 2; Default: 1)*.

**Optional parameters:**

- **min_required_icd_codes** – Minimum number of ICD diagnosis records mapping to the same phecode for that phecode to be considered “present” in an individual. For example, `min_required_icd_codes=2` means a single occurrence of a code isn’t enough to count the person as having that phecode; at least two occurrences are required. Ensure your medical record are comprehensive (not limited to first occurrences) if using this parameter. *(Default: 1)*.
- **date_fmt** – Date format of the **Index date** and **End date** columns in your phenotype data. *(Default: '%Y-%m-%d', i.e. YYYY-MM-DD)*.
- **phecode_version** – Phecode version for mapping diagnosis codes. Currently, version `'1.2'` is the recommended official version (with mapping files for ICD-9-CM/WHO and ICD-10-CM/WHO). An unofficial `'1.3a'` is available in the package for special use cases but is **not** recommended for general use. *(Default: '1.2')*.

## Load phenotype data

After initializing the `DiseaseNetworkData` object, use the `phenotype_data()` method to load your phenotype data file. You need to provide the file path, a dictionary mapping the required column names to your file’s column headers, and a list of any additional covariate column names.

Below are examples demonstrating how to load the dummy phenotype dataset under different study designs. The dummy file is structured for a matched cohort study, but it can be adapted for other designs by dropping certain columns when mapping: if you omit the **Match ID** column in the mapping, the data will be treated as a standard cohort (ignoring the matching groups); if you omit both **Match ID** and **Exposure**, it will be treated as an exposed-only cohort (all individuals considered “exposed”).

```python
# Load phenotype data for a matched cohort study  
col_dict = {  
    'Participant ID': 'ID',  
    'Exposure': 'exposure',  
    'Sex': 'sex',  
    'Index date': 'date_start',  
    'End date': 'date_end',  
    'Match ID': 'group_id'  
}  
covariates_list = ['age', 'BMI']  
data.phenotype_data(  
    phenotype_data_path=r"/test/data/dummy_phenotype.csv",  
    column_names=col_dict,  
    covariates=covariates_list  
)  

# Load phenotype data for a standard cohort study (no matching)  
col_dict = {  
    'Participant ID': 'ID',  
    'Exposure': 'exposure',  
    'Sex': 'sex',  
    'Index date': 'date_start',  
    'End date': 'date_end'  
}  
covariates_list = ['age', 'BMI']  
data.phenotype_data(  
    phenotype_data_path=r"/test/data/dummy_phenotype.csv",  
    column_names=col_dict,  
    covariates=covariates_list  
)  

# Load phenotype data for an exposed-only cohort study (only exposed group, no comparator)  
col_dict = {  
    'Participant ID': 'ID',  
    'Sex': 'sex',  
    'Index date': 'date_start',  
    'End date': 'date_end'  
}  
covariates_list = ['age', 'BMI']  
data.phenotype_data(  
    phenotype_data_path=r"/test/data/dummy_phenotype.csv",  
    column_names=col_dict,  
    covariates=covariates_list  
)  
```

- **phenotype_data_path** – Path to your phenotype data file (CSV or TSV).
- **column_names** – Dictionary mapping the required column names (`'Participant ID'`, `'Index date'`, `'End date'`, `'Sex'`, and depending on design `'Exposure'` and `'Match ID'`) to the corresponding column headers in your file. Include `'Exposure'` for cohort and matched cohort designs, and include `'Match ID'` only for matched cohorts.
- **covariates** – List of additional covariate column names to load (if any). Use an empty list `[]` if there are none. The function will automatically detect each covariate’s type and process it appropriately (e.g., encode categorical variables). For continuous covariates, any rows with missing values will be dropped; for categorical covariates, missing values will be categorized as "NA".

**Optional parameters:**

- **is_single_sex** – If your cohort contains only one sex (all male or all female), set this to `True` so the software knows to treat the Sex column accordingly. *(Default: False)*.
- **force** – If `False`, the method will raise an error if phenotype data has already been loaded into this `DiseaseNetworkData` object (to prevent accidental overwrite). Setting `force=True` will overwrite any existing data in the object with the new data. *(Default: False)*.

**After loading phenotype data:**

Once the phenotype data is loaded, you can inspect the basic characteristics by printing the `data` object:

```python
print(data)  
# Example output (for a matched cohort study):  
"""  
DiNetxify.DiseaseNetworkData  

Study design: matched cohort  

Phenotype data  
Total number of individuals: 60,000 (10,000 exposed and 50,000 unexposed)  
The average group size is: 6.00  
Average follow-up years: 10.44 (exposed) and 10.46 (unexposed)  

Warning: 102 exposed individuals and 440 unexposed individuals have negative or zero follow-up time.  
Consider removing them before merge.  
"""  

```

The printed summary confirms the number of individuals, breakdown by exposure, average matching group size (for matched cohorts), and average follow-up times. Warnings are provided if any participants have non-positive follow-up lengths, which you may want to address (e.g., by removing those individuals) before proceeding.

Additionally, you can generate a basic descriptive table (Table 1) of the phenotype data using the `Table1()` method. This returns a pandas DataFrame summarizing each variable (e.g., medians/IQRs for continuous variables, counts/percentages for categorical variables) and performing simple statistical comparisons between exposed and unexposed groups:

```python
table1_df = data.Table1()  
print(table1_df)  
# Example (truncated) output:  
"""  
                   Variable exposure=1 (n=10,000) exposure=0 (n=50,000)  Test and p-value
0        _age (median, IQR)   57.08 (48.91-65.32)   57.05 (48.87-65.35)  Mann-Whitney U test p-value=9.824e-01
1   follow_up (median, IQR)     9.18 (5.77-13.70)     9.22 (5.80-13.75)  Mann-Whitney U test p-value=6.806e-01
2                sex (n, %)                                               
3                sex=Female        5,045 (50.45%)       25,225 (50.45%)   
4                  sex=Male        4,955 (49.55%)       24,775 (49.55%)  Chi-squared test p-value=1.000e+00 
5                BMI (n, %)                                              ...
6                    BMI=c2        1,945 (19.45%)       10,170 (20.34%)   
7                    BMI=c4        2,022 (20.22%)       10,022 (20.04%)   
8                    BMI=c5        2,002 (20.02%)       10,031 (20.06%)   
9                    BMI=c1        1,999 (19.99%)        9,952 (19.90%)   
10                   BMI=c3        2,032 (20.32%)        9,825 (19.65%)  Chi-squared test p-value=2.552e-01     
"""  
```

This Table 1 gives a quick overview of how the exposed and unexposed groups compare on key variables. You can save this `DataFrame` to a CSV/TSV/Excel file using pandas if needed.

## Load medical record data

After loading the phenotype data, use the `merge_medical_records()` method to load and merge each medical record file. You will call this method for each separate file (e.g., one for ICD-10 and one for ICD-9 in our dummy data). Provide the file path, specify the ICD coding standard used in that file, and a dictionary mapping required columns. The following example code shows how to load the dummy EHR ICD-10 and ICD-9 files:

```python
# Merge the first medical record file (dummy_EHR_ICD10.csv)  
data.merge_medical_records(  
    medical_records_data_path=r"/test/data/dummy_EHR_ICD10.csv",  
    diagnosis_code='ICD-10-WHO',  
    column_names={  
        'Participant ID': 'ID',  
        'Diagnosis code': 'diag_icd10',  
        'Date of diagnosis': 'dia_date'  
    }  
)  

# Merge the second medical record file (dummy_EHR_ICD9.csv)  
data.merge_medical_records(  
    medical_records_data_path=r"/test/data/dummy_EHR_ICD9.csv",  
    diagnosis_code="ICD-9-WHO",  
    column_names={  
        'Participant ID': 'ID',  
        'Diagnosis code': 'diag_icd9',  
        'Date of diagnosis': 'dia_date'  
    }  
)  

```

- **medical_records_data_path** – Path to a medical record data file (CSV or TSV).
- **diagnosis_code** – The diagnosis coding system used in that file. Options include `'ICD-9-CM'`, `'ICD-9-WHO'`, `'ICD-10-CM'`, `'ICD-10-WHO'` (case-sensitive).
- **column_names** – Dictionary mapping the required column names (`'Participant ID'`, `'Diagnosis code'`, `'Date of diagnosis'`) to your file’s column headers.

**Optional parameters:**

- **date_fmt** – Date format of the **Date of diagnosis** column in this file. If not provided, it defaults to the same format used for phenotype dates (`date_fmt` specified in the `DiseaseNetworkData` initialization).
- **chunksize** – If the file is very large, you can specify a number of rows to read per chunk (the function will stream through the file in chunks to manage memory usage). *(Default: 1,000,000 rows per chunk.)*

**During data loading:**

As each medical record file is processed, ***DiNetxify*** will output progress messages and basic stats. For example:

```python
"""
1,000,000 records read, 1,000,000 left after filltering on participant ID/exclusion list of diagnosis codes, 0 records with missing values excluded.
1,668,795 records read, 1,668,795 left after filltering on participant ID/exclusion list of diagnosis codes, 0 records with missing values excluded.
Total: 1,668,795 diagnosis records processed, 0 records with missing values were excluded.
1,286,386 diagnosis records mapped to phecode without truncating.
0 diagnosis records mapped to phecode after truncating to 4 digits.
72,073 diagnosis records mapped to phecode after truncating to 3 digits.
302,908 diagnosis records not mapped to any phecode.
Phecode diagnosis records successfully merged (18,486 invalid records were not merged, typically with diagnosis date later than date of follow-up end)

1 medical records data already merged, merging with a new one.
10,188 records read, 10,188 left after filltering on participant ID/exclusion list of diagnosis codes, 0 records with missing values excluded.
Total: 10,188 diagnosis records processed, 0 records with missing values were excluded.
9,711 diagnosis records mapped to phecode without truncating.
0 diagnosis records mapped to phecode after truncating to 4 digits.
266 diagnosis records mapped to phecode after truncating to 3 digits.
211 diagnosis records not mapped to any phecode.
Phecode diagnosis records successfully merged (0 invalid records were not merged, typically with diagnosis date later than date of follow-up end)
"""
```

From these logs, you can see how many records were read and included, how many were excluded (e.g., missing values or out-of-follow-up-range dates), and how many diagnosis codes were successfully mapped to phecodes versus not mapped. The logs also indicate when multiple files are being merged sequentially.

**After loading medical record data:**

After merging all medical record files, you can print the `data` object again to see a summary of the combined dataset:

```python
print(data)  
# Example output (matched cohort study):  
"""  
Merged Medical records
1,678,983 diagnosis records from 2 medical records file were merged (0 with missing values).
Average number of disease diagnosis during follow-up: 18.99 (exposed) and 7.31 (unexposed)
Average number of disease diagnosis before follow-up: 8.40 (exposed) and 3.46 (unexposed)

Warning: 102 exposed individuals and 440 unexposed individuals have negative or zero follow-up time.  
Consider removing them before merge.  
Warning: 18.15% of ICD-10-WHO codes were not mapped to phecodes for file /test/data/dummy_EHR_ICD10.csv.  
Warning: 2.07% of ICD-9-WHO codes were not mapped to phecodes for file /test/data/dummy_EHR_ICD9.csv.  
"""  
```

This output confirms the number of diagnosis records merged and provides average counts of diagnoses per person (during and before follow-up, by exposure group). Warnings indicate the percentage of codes that could not be mapped to a phecode for each file, so you’re aware of any unmapped codes.

## Save DiseaseNetworkData object

At this stage, after loading phenotype and medical record data, you may want to save the `DiseaseNetworkData` object for later use. Saving allows you to reuse the prepared data without re-reading and processing raw files each time, facilitating reproducibility and easy sharing of the processed data. ***DiNetxify*** provides two methods: `save()` (which uses Python’s pickle serialization, saving to a compressed `.pkl.gz` file) and `save_npz()` (which saves to a compressed NumPy `.npz` file). You can use either or both depending on your needs. For example:

```python
# Save the data object to a gzipped pickle file  
data.save('/your/project/path/cohort_data')  
# (This will produce a file named "cohort_data.pkl.gz")  

# Save the data object to a NumPy .npz file  
data.save_npz('/your/project/path/cohort_data')  
# (This will produce a file named "cohort_data.npz")
```

You do not need to add the file extension in the path; the functions will append `.pkl.gz` or `.npz` automatically. Make sure to choose a directory where you have write permissions and enough storage space (the files can be large if your dataset is large).

## Reload DiseaseNetworkData object

If you have previously saved a `DiseaseNetworkData` object, you can reload it instead of re-reading all input files. This is especially useful for large datasets or when sharing the processed object with collaborators. To reload, first instantiate a new `DiseaseNetworkData` object with the same `study_design` and `phecode_level` that the data was created with, then call the corresponding load function (`load()` or `load_npz()`). For example:

```python
import DiNetxify as dnt  

# Create a new DiseaseNetworkData object with the same design/parameters  
data = dnt.DiseaseNetworkData(  
    study_design='matched cohort',  
    phecode_level=1,  
)  

# Load from a .pkl.gz file  
data.load('/your/project/path/cohort_data', force=True)  

# Or load from a .npz file  
data.load_npz('/your/project/path/cohort_data', force=True)  
```
