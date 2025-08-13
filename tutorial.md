# Three-dimensional disease network analysis using DiNetxify
## Table of contents

- [1. Input data preparation](#1-input-data-preparation)  
  - [1.1 Requirements for input data](#11-requirements-for-input-data)  
  - [1.2 Dummy dataset overview](#12-dummy-dataset-overview)  
- [2. Data Harmonization](#2-data-harmonization)  
  - [2.1 Initializing the data object](#21-initializing-the-data-object)  
  - [2.2 Load phenotype data](#22-load-phenotype-data)  
  - [2.3 Load medical records data](#23-load-medical-records-data)  
  - [2.4 Save DiseaseNetworkData object](#24-save-diseasenetworkdata-object)
  - [2.5 Reload DiseaseNetworkData object](#25-reload-diseasenetworkdata-object)
- [3. Data Analysis](#3-data-analysis)  
  - [3.1 One-step analysis](#31-one-step-analysis)
  - [3.2 Step-by-step analysis](#32-step-by-step-analysis)
    - [3.2.1 PheWAS Analysis](#321-phewas-analysis)
    - [3.2.2 Disease pair generation](#322-disease-pair-generation)
    - [3.2.3 Comorbidity strength estimation](#323-comorbidity-strength-estimation)
    - [3.2.4 Binomial test](#324-binomial-test)
    - [3.2.5 Comorbidity network analysis](#325-comorbidity-network-analysis)
    - [3.2.6 Disease trajectory analysis](#326-disease-trajectory-analysis)
- [4. Visualization](#4-visualization)  
  - [4.1 Initializing the plot object](#41-initializing-the-plot-object)  
  - [4.2 PheWAS plot](#42-phewas-plot)  
  - [4.3 Comorbidity network plot](#43-comorbidity-network-plot)  
  - [4.4 Disease trajectory plot](#44-disease-trajectory-plot)  
  - [4.5 Three dimension plot](#45-three-dimension-plot)  
- [API Reference](#api-reference)  
  - [Class DiseaseNetworkData](#class-diseasenetworkdata)  
    - [phenotype_data](#phenotype_data)  
    - [Table1](#table1)  
    - [merge_medical_records](#merge_medical_records)  
    - [get_attribute](#get_attribute)  
    - [medical_records_to_dataframe](#medical_records_to_dataframe)  
    - [modify_phecode_level](#modify_phecode_level)  
    - [disease_pair](#disease_pair) 
    - [save](#save)
    - [load](#load)
    - [save_npz](#save_npz)
    - [load_npz](#load_npz)
  - [Analysis functions](#analysis-functions)  
    - [disease_network_pipeline](#disease_network_pipeline)  
    - [phewas](#phewas)  
    - [phewas_multipletests](#phewas_multipletests)  
    - [comorbidity_strength](#comorbidity_strength)  
    - [comorbidity_strength_multipletests](#comorbidity_strength_multipletests)  
    - [binomial_test](#binomial_test)  
    - [binomial_multipletests](#binomial_multipletests)  
    - [comorbidity_network](#comorbidity_network)  
    - [comorbidity_multipletests](#comorbidity_multipletests)  
    - [disease_trajectory](#disease_trajectory)  
    - [trajectory_multipletests](#trajectory_multipletests)  
  - [Class Plot](#class-plot)  
    - [three_dimension_plot](#three_dimension_plot)  
    - [comorbidity_network_plot](#comorbidity_network_plot)  
    - [trajectory_plot](#trajectory_plot)  
    - [phewas_plot](#phewas_plot)

## 1. Input data preparation

### 1.1 Requirements for input data

***DiNetxify*** enables 3D disease network analysis on cohort data from electronic health records (EHR) and offers three study designs: the **standard cohort**, which compares individuals with a specific disease or exposure (e.g., depression or smoking) against the general population; the **matched cohort**, which pairs subjects on key characteristics to reduce confounding; and the **exposed-only cohort**, which examines disease networks within a defined group (e.g., older adults) without an unexposed comparison group.

To begin using ***DiNetxify***, two datasets are required: a **phenotype data** file containing each participant’s baseline information, and one or more **medical records data** files extracted from an EHR database that list diagnoses (codes and dates) for all cohort individuals over the study period. The specific requirements for these datasets are:

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

  - **Participant ID **– The same unique ID used in the phenotype data, linking each record to an individual.
  - **Diagnosis code** – A diagnosis code (e.g., ICD-10 or ICD-9 code).
  - **Date of diagnosis** – The date of that diagnosis/event (format consistent with the phenotype dates, e.g., `YYYY-MM-DD`).

  The **medical records data** should be in a “long” format (multiple rows per participant if they have multiple diagnoses). ***DiNetxify*** will automatically filter these records to include only those within each individual’s follow-up period (from index date up to end date). **Do not pre-filter** the medical records by date or by first occurrence — provide the complete set of diagnoses for each participant, and let the software handle the filtering and mapping. Each medical records file should use a single coding system for diagnoses. Currently supported code versions are ICD-9 (WHO and CM) and ICD-10 (WHO and CM). If your data uses a different coding system, you will need to map it to one of the supported formats beforehand.

### 1.2 Dummy dataset overview

A [dummy dataset](https://github.com/HZcohort/DiNetxify/tree/main/tests/data) is provided to help you become familiar with the input format and to allow you to run through the full analysis workflow before using your own data. It simulates a matched-cohort study of 10,000 exposed individuals and 50,000 matched unexposed individuals, along with their entire follow-up EHR records.

> **Note:** All participant characteristics and diagnoses in this dummy dataset are randomly generated. The ICD-9 and ICD-10 codes correspond to real classifications, and the analysis may yield seemingly significant associations, but these results do **not** reflect true medical findings. They are for instructional purposes only.

- The dummy dataset consists of three CSV files:
  - **`dummy_phenotype.csv`** – Simulated baseline characteristics for 60,000 individuals, containing:
    - **ID** – Unique participant identifier.
    - **date_start**, **date_end** – Follow-up start and end dates.
    - **exposure** – Exposure status (0 = unexposed, 1 = exposed).
    - **group_id** – Matching group identifier (each exposed is matched with unexposed in groups).
    - **sex** – Biological sex (1 = female, 0 = male).
    - **age** – Baseline age (years).
    - **BMI** – Body mass index category.
  - **`dummy_EHR_ICD9.csv`** – Simulated EHR diagnoses coded in ICD-9 (10,188 records). Columns:
    - **ID** – Participant ID (matches the phenotype file).
    - **dia_date** – Diagnosis date.
    - **diag_icd9** – ICD-9 diagnosis code.
  - **`dummy_EHR_ICD10.csv`** – Simulated EHR diagnoses coded in ICD-10 (1,048,576 records). Columns:
    - **ID** – Participant ID.
    - **dia_date** – Diagnosis date.
    - **diag_icd10** – ICD-10 diagnosis code.

Using this dummy dataset, you can practice the workflow and verify that the tool runs correctly. In the following sections, we will demonstrate the analysis steps using the dummy data.

## 2. Data harmonization

Data harmonization involves loading and merging the **phenotype data** and **medical records data** into a single `DiseaseNetworkData` object for analysis. During this process, the software ensures consistent coding (e.g., converting diagnosis codes to phecodes) and standardized formatting (e.g., datetime parsing for diagnosis and follow-up periods).

### 2.1 Initializing the data object

First, import the ***DiNetxify*** package and instantiate a `DiseaseNetworkData` object with your chosen study design, phecode level, and any optional parameters. For example:

```python
import DiNetxify as dnt  

# For a standard cohort study  
data = dnt.DiseaseNetworkData(  
    study_design='cohort',  
    phecode_level=1,  
)  

# For a matched cohort study  
data = dnt.DiseaseNetworkData(  
    study_design='matched cohort',  
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

### 2.2 Load phenotype data

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
                   Variable        exposure=1 (n=10,000)   exposure=0 (n=50,000)            Test and p-value  
0        age (median, IQR)       57.08 (48.91–65.32)       57.05 (48.87–65.35)    Mann-Whitney U p=0.9824  
1   follow_up (median, IQR)       9.18 (5.77–13.70)         9.22 (5.80–13.75)    Mann-Whitney U p=0.6806  
2                sex (n, %)  
3                sex=Female          5,045 (50.45%)           25,225 (50.45%)   …  
...  
"""  
```

This Table 1 gives a quick overview of how the exposed and unexposed groups compare on key variables. You can save this `DataFrame` to a CSV/TSV/Excel file using pandas if needed.

### 2.3 Load medical record data

After loading the phenotype data, use the `merge_medical_records()` method to load and merge each medical records file. You will call this method for each separate file (e.g., one for ICD-10 and one for ICD-9 in our dummy data). Provide the file path, specify the ICD coding standard used in that file, and a dictionary mapping required columns. The following example code shows how to load the dummy EHR ICD-10 and ICD-9 files:

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

# Merge the second medical records file (dummy_EHR_ICD9.csv)  
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

- **medical_records_data_path** – Path to a medical records data file (CSV or TSV).
- **diagnosis_code** – The diagnosis coding system used in that file. Options include `'ICD-9-CM'`, `'ICD-9-WHO'`, `'ICD-10-CM'`, `'ICD-10-WHO'` (case-sensitive).
- **column_names** – Dictionary mapping the required column names (`'Participant ID'`, `'Diagnosis code'`, `'Date of diagnosis'`) to your file’s column headers.

**Optional parameters:**

- **date_fmt** – Date format of the **Date of diagnosis** column in this file. If not provided, it defaults to the same format used for phenotype dates (`date_fmt` specified in the `DiseaseNetworkData` initialization).
- **chunksize** – If the file is very large, you can specify a number of rows to read per chunk (the function will stream through the file in chunks to manage memory usage). *(Default: 1,000,000 rows per chunk.)*

**During data loading:**

As each medical records file is processed, ***DiNetxify*** will output progress messages and basic stats. For example:

```python
"""
1,000,000 records read (1,000,000 included after filtering on participant ID), 0 records with missing values excluded.  
1,668,795 records read (1,668,795 included after filtering on participant ID), 0 records with missing values excluded.  
Total: 1,668,795 diagnosis records processed, 0 records with missing values were excluded.  
1,286,386 diagnosis records mapped to phecode without truncating.  
0 diagnosis records mapped to phecode after truncating to 4 digits.  
72,073 diagnosis records mapped to phecode after truncating to 3 digits.  
302,908 diagnosis records not mapped to any phecode.  
Phecode diagnosis records successfully merged (18,486 invalid records were not merged, typically due to diagnosis date beyond follow-up end).  

1 medical record file already merged, merging with a new one.  
10,188 records read (10,188 included after filtering on participant ID), 0 records with missing values excluded.  
Total: 10,188 diagnosis records processed, 0 records with missing values were excluded.  
9,711 diagnosis records mapped to phecode without truncating.  
0 diagnosis records mapped to phecode after truncating to 4 digits.  
266 diagnosis records mapped to phecode after truncating to 3 digits.  
211 diagnosis records not mapped to any phecode.  
Phecode diagnosis records successfully merged (0 invalid records were not merged).  
"""
```

From these logs, you can see how many records were read and included, how many were excluded (e.g., missing values or out-of-follow-up-range dates), and how many diagnosis codes were successfully mapped to phecodes versus not mapped. The logs also indicate when multiple files are being merged sequentially.

**After loading medical record data:**

After merging all medical records files, you can print the `data` object again to see a summary of the combined dataset:

```python
print(data)  
# Example output (matched cohort study):  
"""  
Merged Medical records  
2 medical record files with 1,678,983 diagnosis records were merged (0 with missing values).  
Average number of disease diagnoses during follow-up: 18.99 (exposed) and 7.31 (unexposed)  
Average number of disease diagnoses before follow-up: 8.40 (exposed) and 3.46 (unexposed)  

Warning: 102 exposed individuals and 440 unexposed individuals have negative or zero follow-up time.  
Consider removing them before merge.  
Warning: 18.15% of ICD-10-WHO codes were not mapped to phecodes for file /test/data/dummy_EHR_ICD10.csv.  
Warning: 2.07% of ICD-9-WHO codes were not mapped to phecodes for file /test/data/dummy_EHR_ICD9.csv.  
"""  
```

This output confirms the number of diagnosis records merged and provides average counts of diagnoses per person (during and before follow-up, by exposure group). Warnings indicate the percentage of codes that could not be mapped to a phecode for each file, so you’re aware of any unmapped codes.

### 2.4 Save DiseaseNetworkData object

At this stage, after loading phenotype and medical records data, you may want to save the `DiseaseNetworkData` object for later use. Saving allows you to reuse the prepared data without re-reading and processing raw files each time, facilitating reproducibility and easy sharing of the processed data. ***DiNetxify*** provides two methods: `save()` (which uses Python’s pickle serialization, saving to a compressed `.pkl.gz` file) and `save_npz()` (which saves to a compressed NumPy `.npz` file). You can use either or both depending on your needs. For example:

```python
# Save the data object to a gzipped pickle file  
data.save('/your/project/path/cohort_data')  
# (This will produce a file named "cohort_data.pkl.gz")  

# Save the data object to a NumPy .npz file  
data.save_npz('/your/project/path/cohort_data')  
# (This will produce a file named "cohort_data.npz")
```

You do not need to add the file extension in the path; the functions will append `.pkl.gz` or `.npz` automatically. Make sure to choose a directory where you have write permissions and enough storage space (the files can be large if your dataset is large).

### 2.5 Reload DiseaseNetworkData object

If you have previously saved a `DiseaseNetworkData` object, you can reload it instead of re-reading all input files. This is especially useful for large datasets or when sharing the processed object with collaborators. To reload, first instantiate a new `DiseaseNetworkData` object with the same `study_design` and `phecode_level` that the data was created with, then call the corresponding load function (`load()` or `load_npz()`). For example:

```python
import DiNetxify as dnt  

# Create a new DiseaseNetworkData object with the same design/parameters  
data = dnt.DiseaseNetworkData(  
    study_design='cohort',  
    phecode_level=1,  
)  

# Load from a .pkl.gz file  
data.load('/your/project/path/cohort_data')  

# Or load from a .npz file  
data.load_npz('/your/project/path/cohort_data')  
```

## 3. Data Analysis

Once the data is prepared and stored in a `DiseaseNetworkData` object, DiNetxify offers two approaches to perform the disease network analysis:

1. **One-step analysis:** a comprehensive pipeline that automates the entire sequence of analyses (PheWAS → disease pair generation → comorbidity strength estimation → binomial test → comorbidity network analysis → disease trajectory analysis) with one function call. This is convenient and ensures all steps are performed in the correct order with default or specified parameters.
2. **Step-by-step analysis:** individual functions for each analysis component, allowing you to run and inspect each step separately. This approach offers more control and flexibility (e.g., to tweak parameters at each step or to examine intermediate results), at the expense of writing a bit more code.

Both approaches output their results as pandas DataFrames, which you can further analyze or export (to CSV/Excel, etc.) using pandas. The one-step pipeline minimizes redundant computations and code, but does not allow modifying certain internal parameters beyond what its arguments expose. The step-by-step approach is more verbose but lets you adjust and understand each phase of the analysis in detail. We’ll demonstrate both.

### 3.1 One-step analysis

With your `DiseaseNetworkData` ready (let’s call it `data`), you can perform the entire analysis in one go using the `disease_network_pipeline()` function. This function returns all major result DataFrames and can also save output files. The example below illustrates using `disease_network_pipeline()`; it is applicable to any of the three study designs (the function internally adapts to the design specified in `data`):

```python
# Reminder:
# When using multiprocessing, ensure that the code is enclosed within the following block.
# This prevents entering a never ending loop of new process creation.
from DiNetxify import disease_network_pipeline  

if __name__ == "__main__":  # Required when using multiprocessing on Windows/Mac  
    phewas_result, com_strength_result, com_network_result, binomial_result, trajectory_result = disease_network_pipeline(
        data=data, n_process=2,
        n_threshold_phewas=100,
        n_threshold_comorbidity=100, 
        comorbidity analysis,
        output_dir="/your/project/path/results/",
        project_prefix="disease_network",
        keep_positive_associations=False,
        save_intermediate_data=False,
        system_exl=['symptoms', 'others', 'injuries & poisonings'],
        pipeline_mode="v1",
        method="RPCN",
        covariates=['BMI', 'age'],
        matching_var_dict={'sex': 'exact'},
        matching_n=2,
        min_interval_days=0,
        max_interval_days=float('inf'),
        enforce_temporal_order=False,
        correction='bonferroni',
        cutoff=0.05) 
```

- **Parameters (key arguments in `disease_network_pipeline`):**

    - **data** – The `DiseaseNetworkData` object containing your loaded cohort data.
    - **n_process** – Number of parallel processes for computation. Use `1` for single-threaded execution, or higher to speed up analysis with multiprocessing (especially beneficial for large datasets). *(No default; you must specify this.)*
    - **n_threshold_phewas** – Minimum number of exposed cases required for a disease (phecode) to be included in the PheWAS analysis. This filters out very rare outcomes. (This value is passed to the internal `phewas()` function.)
    - **n_threshold_comorbidity** – Minimum number of exposed individuals in whom a given disease pair co-occurs (considering both temporal and non-temporal occurrences) to include that pair in the comorbidity strength analysis. (Passed to `comorbidity_strength()`.)
    - **output_dir** – Directory path for saving output files. The pipeline will create result files here (e.g., CSVs of significant results, log files, etc.). Use an absolute path or a path relative to your working directory.
    - **project_prefix** – A string prefix for naming output files. For example, if `project_prefix="disease_network"`, output files might be named like `disease_network_phewas_results.csv`, etc.
    - **keep_positive_associations** – If set to `True`, the pipeline will filter results to retain only “positive” associations: diseases with hazard ratio (HR) > 1 in the PheWAS and disease pairs with positive correlation in comorbidity analysis. *(Default: False – retains all significant associations regardless of direction.)*
    - **save_intermediate_data** – If `True`, intermediate `DiseaseNetworkData` objects (specifically those created during disease pair generation) will be saved to disk. This can be useful for debugging or inspecting intermediate steps, but will consume additional disk space. *(Default: False)*.
    - **system_exl** – A list of phecode disease *systems* to exclude from all analyses. If certain categories of diseases (systems) are not of interest or should be filtered out (e.g., ‘symptoms’ or ‘injuries & poisonings’), list them here. If set to None or an empty list, no system is excluded. *(Default: None)*. Valid system names include: *circulatory system, congenital anomalies, dermatologic, digestive, endocrine/metabolic, genitourinary, hematopoietic, infectious diseases, injuries & poisonings, mental disorders, musculoskeletal, neoplasms, neurological, pregnancy complications, respiratory, sense organs, symptoms, others*.
    - **pipeline_mode** – Specifies the order of analyses. Two modes are available:
      - **'v1'**: PheWAS → comorbidity strength → binomial test → *then parallel/complementary*: comorbidity network analysis and disease trajectory analysis. (In this mode, the binomial test is run on all eligible disease pairs without considering network results, so trajectory and network analyses can be done independently.)
      - **'v2'**: PheWAS → comorbidity strength → comorbidity network analysis → binomial test → disease trajectory analysis. (In this mode, only the disease pairs deemed significant in the network analysis are subjected to the binomial test and subsequently used for trajectory analysis. This means the trajectory analysis focuses on a subset defined by network results.)
         *(Default: 'v1')*. Choose 'v2' if you want a more stringent approach where trajectory analysis is conditional on network significance; otherwise, 'v1' covers all pairs passing earlier filters.
    - **method** – The method used for comorbidity network and disease trajectory analyses. Options are:
      - **'RPCN'** – *Regularized Partial Correlation Network*. This method uses a regularized logistic regression framework including all other diseases as covariates (with L1 penalty) to evaluate direct disease-disease associations, adjusting for covariates. *(This is the default and our recommended approach.)*
      - **'PCN_PCA'** – *Partial Correlation Network with PCA*. Similar to RPCN but applies principal component analysis to reduce dimensionality of the “other diseases” covariates before computing the network. This can simplify the model when there are many diseases.
      - **'CN'** – *Correlation Network*. A simpler approach using standard logistic regression for each disease pair (plus covariates) without partialling out other diseases. Essentially assesses correlation of each pair independently.
         *(Default: 'RPCN')*. The choice of method will affect how comorbidity networks and trajectories are inferred.
    - **covariates** – List of covariate names to adjust for in the PheWAS, network, and trajectory analyses (in addition to the required *sex* variable, which is always adjusted for). These should match the covariate names you provided in `phenotype_data()`. For example, `['BMI', 'age']` as shown above. *(Default: None, meaning adjust for sex only along with any default adjustments like matching factors in matched cohorts.)*
    - **matching_var_dict** – A dictionary specifying how controls are matched to cases for the trajectory analysis (which uses an incidence density sampling approach). Keys are variable names to match on, and values specify matching criteria: for categorical or binary variables use `'exact'`; for continuous variables, provide a numeric tolerance. **Important:** use `'sex'` (literally) to match on sex (even if your original column name was different) because the data object uses a standardized 'sex' field. Other covariates should be referred to by their original names. *(Default: {'sex': 'exact'}, meaning match controls to cases by sex.)*.
    - **matching_n** – Maximum number of matched controls to select for each case in the trajectory analysis. For example, `matching_n=2` tries to find up to 2 controls per case. *(Default: 2)*.
    - **min_interval_days** – Minimum time interval in days required between two diagnoses to consider one occurring *before* the other for temporal (trajectory) analysis. If the time between D1 and D2 diagnoses in an individual is less than or equal to this threshold, the pair is treated as effectively simultaneous (and thus not counted as a temporal sequence). *(Default: 0 days)*. Setting a positive number here can exclude very closely timed diagnoses from being considered as one preceding the other.
    - **max_interval_days** – Maximum time interval in days to consider for disease pairs. If the gap between two diagnoses is larger than this, that pair occurrence might be ignored for certain analyses. This applies to both temporal (ordered) and non-temporal pair considerations. *(Default: infinity, i.e., no maximum gap applied.)*.
    - **enforce_temporal_order** – If `True`, the pipeline will enforce strict temporal ordering in the trajectory analyses and related significance testing: any individual who has a disease pair in the opposite order (D2 before D1) may be excluded from certain calculations, and the binomial test will only consider pairs where a clear ordering can be established. In practice, setting this to True means the binomial test will ignore individuals who have the diseases in both orders, and the trajectory analysis will also respect the specified `min_interval_days` and `max_interval_days` strictly. *(Default: False)*.
    - **correction** – The multiple hypothesis testing correction method to apply to p-values (where applicable). This uses methods from `statsmodels.stats.multitest.multipletests`. Options include `'bonferroni'`, `'holm'`, `'fdr_bh'`, etc. *(Default: 'bonferroni')*.
    - **cutoff** – Significance cutoff for adjusted p-values. *(Default: 0.05)*. Any results with adjusted p-value above this threshold will be considered non-significant and typically filtered out from the final results.

    The `disease_network_pipeline` will return a tuple of DataFrames: in order, these correspond to **phewas_result**, **com_strength_result** (comorbidity strength), **com_network_result** (comorbidity network), **binomial_result**, and **trajectory_result**. In addition to returning these, the function writes out certain results and logs to files in `output_dir` for your records. You can adjust which results to focus on based on your research question (for instance, the comorbidity network and trajectory analyses might produce quite a lot of output; you could choose pipeline_mode ‘v1’ to consider them separately).

#### After one-step analysis:

The **one-step analysis** generates comprehensive results presented in sequential order: PheWAS analysis, comorbidity strength estimation, comorbidity network analysis, binomial testing, and disease trajectory analysis. Each result includes detailed variable descriptions to facilitate interpretation.

**Result of PheWAS analysis**

```python
# print the result to show some details of PheWAS analysis result
print(phewas_result)

# This will output something like (e.g., for a matched cohort study):
     phecode                                            disease               system  ...       phewas_p  phewas_p_significance phewas_p_adjusted
0        8.0                               Intestinal infection  infectious diseases  ...  4.141110e-109                   True     1.192640e-106 
1       10.0                                       Tuberculosis  infectious diseases  ...   3.515722e-55                   True      1.012528e-52 
3       38.0                                         Septicemia  infectious diseases  ...   2.963084e-40                   True      8.533681e-38 
4       41.0                            Bacterial infection NOS  infectious diseases  ...  1.601040e-139                   True     4.610994e-137 
7       70.0                                    Viral hepatitis  infectious diseases  ...  4.783200e-129                   True     1.377562e-126 
..       ...                                                ...                  ...  ...            ...                    ...               ... 
419    792.0  Abnormal Papanicolaou smear of cervix and cerv...        genitourinary  ...            NaN                  False               NaN 
420    796.0           Elevated prostate specific antigen [PSA]        genitourinary  ...            NaN                  False               NaN 
421    860.0                Bone marrow or stem cell transplant            neoplasms  ...            NaN                  False               NaN 
422    931.0  Contact dermatitis and other eczema due to pla...         dermatologic  ...            NaN                  False               NaN 
426    980.0  Encounter for long-term (current) use of antib...  infectious diseases  ...            NaN                  False               NaN
```

| Variable Name           | Type      | Description |
|------------------------|-----------|-------------|
| `phecode`              | String    | Disease code (Phecode) used in PheWAS analysis |
| `disease`              | String    | Disease name corresponding to the Phecode |
| `system`               | String    | Phecode disease system corresponding to the Phecode (e.g., infectious diseases) |
| `sex`                  | String    | Sex-specific for the disease (e.g., Both, Male, Female) |
| `N_cases_exposed`      | Integer   | Number of individuals diagnosed with the disease in the exposed group |
| `describe`             | String    | Descriptions of the model fitting state and removed covariates with reasons |
| `exposed_group`        | String    | Incidence rate (unit: per 1,000 person-years) in the exposed group (i.e., 372/96.25 (3.86), which means number of individuals diagnosed with the disease divided by time at risk) |
| `unexposed_group`      | String    | Incidence rate (unit: per 1,000 person-years) in the unexposed group (i.e., 372/96.25 (3.86), which means number of individuals diagnosed with the disease divided by time at risk) |
| `phewas_coef`          | Float     | Estimated coefficient from the model |
| `phewas_se`            | Float     | Standard error of the estimated coefficient |
| `phewas_p`             | Float     | P-value indicating statistical significance of the coefficient |
| `phewas_p_significance`| Boolean   | Indicates whether the result is statistically significant based on adjusted p-value (True/False) |
| `phewas_p_adjusted`    | Float     | Adjusted p-value accounting for multiple comparisons |

**Result of comorbidity strength estimation**

```python
print(com_strength_result)
# This will output something like (e.g., for a matched cohort study):
       phecode_d1  phecode_d2 name_disease_pair  N_exposed  n_total  ...  sex_d2  phi_p_significance  phi_p_adjusted  RR_p_significance  RR_p_adjusted
0            71.0       172.0        71.0-172.0      10000     9048  ...    Both               False    1.000000e+00              False       1.0
1           117.0       759.0       117.0-759.0      10000     9333  ...    Both               False    1.000000e+00              False       1.0
2           327.0       371.0       327.0-371.0      10000     9183  ...    Both                True    1.398268e-08              False       1.0
3           333.0       345.0       333.0-345.0      10000     7918  ...    Both                True    9.509075e-03              False       1.0
4           473.0       522.0       473.0-522.0      10000     8802  ...    Both               False    1.000000e+00              False       1.0
...           ...         ...               ...        ...      ...  ...     ...                 ...             ...                ...       ...
41275       383.0       565.0       383.0-565.0      10000     9391  ...    Both               False             NaN              False       NaN
41285       339.0       418.0       339.0-418.0      10000     9724  ...    Both               False             NaN              False       NaN
41286       372.0       754.0       372.0-754.0      10000     8746  ...    Both               False             NaN              False       NaN
41302       212.0       312.0       212.0-312.0      10000     9777  ...    Both               False             NaN              False       NaN
41314       242.0       739.0       242.0-739.0      10000     8627  ...    Both               False             NaN              False       NaN
```

| Variable Name            | Type    | Description |
|-------------------------|---------|-------------|
| `phecode_d1`            | Integer | Phecode for disease 1 in the disease pair |
| `phecode_d2`            | Integer | Phecode for disease 2 in the disease pair |
| `name_disease_pair`     | String  | Name of the disease pair (format: "D1-D2") |
| `N_exposed`             | Integer | Total number of individuals in exposed group |
| `n_total`               | Integer | Number of individuals in the sub-cohort (Individuals who had no history of comorbidities with disease 1 and no history of comorbidities associated with disease 2, and met the sex-specific eligibility criteria for both diseases) |
| `n_d1d2_diagnosis`      | Integer | Number of individuals diagnosed with both diseases |
| `n_d1_diagnosis`        | Integer | Number of individuals diagnosed with disease 1 |
| `n_d2_diagnosis`        | Integer | Number of individuals diagnosed with disease 2 |
| `n_d1d2_nontemporal`    | Integer | Number of individuals diagnosed with both disease 1 and disease 2 with a defined non-temporal order (i.e., the time between the two diagnosis is less than and equal to the `min_interval_days` parameter or more than the `max_interval_days` parameter) |
| `n_d1d2_temporal`       | Integer | Number of individuals diagnosed with disease 1 followed by disease 2 in a defined temporal order (i.e., the interval time between the two diagnosis is more than the `min_interval_days` parameter, and less than and equal to the `max_interval_days` parameter) |
| `n_d2d1_temporal`       | Integer | Number of individuals diagnosed with disease 2 followed by disease 1 in a defined temporal order (i.e., the interval time between the two diagnosis is more than the `min_interval_days` parameter, and less than and equal to the `max_interval_days` parameter) |
| `phi_coef`              | Float   | Phi coefficient (φ), the association between disease 1 and disease 2 |
| `phi_p`                 | Float   | P-value for Phi coefficient significance |
| `RR`                    | Float   | Relative risk estimate for the disease pair |
| `RR_p`                  | Float   | P-value for relative risk |
| `phi_p_adjusted`        | Float   | Adjusted P-value for Phi coefficient (multiple comparisons) |
| `RR_p_adjusted`         | Float   | Adjusted P-value for relative risk (multiple comparisons) |
| `phi_p_significance`    | Boolean | Whether the Phi is statistically significant based on adjusted p-value |
| `RR_p_significance`     | Boolean | Whether the RR is statistically significant based on adjusted p-value |
| `disease_d1`            | String  | Name of disease 1 |
| `system_d1`             | String  | Phecode disease system related to disease 1 |
| `sex_d1`                | String  | Sex-specific for disease 1 |
| `disease_d2`            | String  | Name of disease 2 |
| `system_d2`             | String  | Phecode disease system related to disease 2 |
| `sex_d2`                | String  | Sex-specific for disease 2 |

**Result of binomial test**

```python
print(binomial_result)
# This will output something like (e.g., for a matched cohort study):
     phecode_d1  phecode_d2 name_disease_pair  n_d1d2_nontemporal  ...             system_d2  sex_d2  binomial_p_significance  binomial_p_adjusted
0         743.0       573.0      743.0->573.0                  39  ...             digestive    Both                    False                  1.0
1         752.0       627.0      752.0->627.0                  30  ...         genitourinary  Female                    False                  1.0
2         287.0       496.0      287.0->496.0                  19  ...           respiratory    Both                    False                  1.0
3         557.0       365.0      557.0->365.0                   9  ...          sense organs    Both                    False                  1.0
4         733.0       528.0      733.0->528.0                  30  ...             digestive    Both                    False                  1.0
..          ...         ...               ...                 ...  ...                   ...     ...                      ...                  ...
160       364.0       626.0      364.0->626.0                   8  ...         genitourinary  Female                    False                  1.0
161       287.0       369.0      287.0->369.0                  24  ...          sense organs    Both                    False                  1.0
162       269.0       481.0      269.0->481.0                  19  ...           respiratory    Both                    False                  1.0
163       287.0       749.0      287.0->749.0                  20  ...  congenital anomalies    Both                    False                  1.0
164       287.0       528.0      287.0->528.0                  28  ...             digestive    Both                    False                  1.0
```

| Variable Name               | Type    | Description |
|----------------------------|---------|-------------|
| `phecode_d1`               | Float   | Phecode for disease 1 in the temporal disease pair |
| `phecode_d2`               | Float   | Phecode for disease 2 in the temporal disease pair |
| `name_disease_pair`        | String  | Name of the temporal disease pair (e.g., D1->D2) |
| `n_d1d2_nontemporal`       | Float   | Number of individuals diagnosed with both disease 1 and disease 2 with a defined non-temporal order (i.e., the time between the two diagnosis is less than and equal to the **min_interval_days** parameter or more than the **max_interval_days** parameter) |
| `n_d1d2_temporal`          | Float   | Number of individuals diagnosed (D1 followed by D2) with a defined temporal order (i.e., the time between the two diagnosis is more than the **min_interval_days** parameter, and less than and equal to the **max_interval_days** parameter) |
| `n_d2d1_temporal`          | Float   | Number of individuals diagnosed (D2 followed by D1) with a defined temporal order (i.e., the time between the two diagnosis is more than the **min_interval_days** parameter, and less than and equal to the **max_interval_days** parameter) |
| `binomial_p`               | Float   | P-value from the binomial test for directionality |
| `binomial_proportion`      | Float   | Proportion of D1 followed by D2 among all temporal disease pairs (**n_d1d2_temporal**/(**n_d1d2_temporal**+**n_d2d1_temporal**)) |
| `binomial_proportion_ci`   | String  | Confidence interval for the binomial proportion |
| `disease_d1`               | String  | Name of disease 1 |
| `system_d1`                | String  | Phecode disease system for disease 1 |
| `sex_d1`                   | String  | Sex-specific for disease 1 |
| `disease_d2`               | String  | Name of disease 2 |
| `system_d2`                | String  | Phecode disease system for disease 2 |
| `sex_d2`                   | String  | Sex-specific for disease 2 |
| `binomial_p_significance` | Boolean | Indicates whether the result is statistically significant based on adjusted p-value |
| `binomial_p_adjusted`     | Float   | Adjusted p-value for multiple comparisons |

**Result of comorbidity network analysis**

```python
print(com_network_result)
# This will output something like (e.g., for a matched cohort study):
     phecode_d1  phecode_d2 name_disease_pair  N_exposed  ...             system_d2  sex_d2 comorbidity_p_significance comorbidity_p_adjusted
0         573.0       743.0       573.0-743.0      10000  ...       musculoskeletal    Both                       True           4.444539e-76     
1         627.0       752.0       627.0-752.0      10000  ...  congenital anomalies    Both                       True           3.258098e-17     
2         287.0       496.0       287.0-496.0      10000  ...           respiratory    Both                       True           3.293334e-36     
3         365.0       557.0       365.0-557.0      10000  ...             digestive    Both                       True          6.826683e-111
4         528.0       733.0       528.0-733.0      10000  ...       musculoskeletal    Both                       True           1.819824e-56
..          ...         ...               ...        ...  ...                   ...     ...                        ...                    ...
160       364.0       626.0       364.0-626.0      10000  ...         genitourinary  Female                       True           5.244532e-12
161       287.0       369.0       287.0-369.0      10000  ...          sense organs    Both                       True           1.010592e-57
162       269.0       481.0       269.0-481.0      10000  ...           respiratory    Both                       True           1.133811e-62
163       287.0       749.0       287.0-749.0      10000  ...  congenital anomalies    Both                       True           4.266938e-61
164       287.0       528.0       287.0-528.0      10000  ...             digestive    Both                       True           1.914133e-82
```

| Variable Name                | Type    | Description |
|-----------------------------|---------|-------------|
| `phecode_d1`                | Float   | Phecode for disease 1 in the non-temporal disease pair |
| `phecode_d2`                | Float   | Phecode for disease 2 in the non-temporal disease pair |
| `name_disease_pair`         | String  | Name of the non-temporal disease pair (e.g., "D1-D2") |
| `N_exposed`                 | Integer | Total number of individuals in exposed group |
| `n_total`                   | Integer | Number of individuals in the sub-cohort (Individuals who had no history of comorbidities with disease 1 and no history of comorbidities associated with disease 2, and met the sex-specific eligibility criteria for both diseases) |
| `n_exposed/n_cases`         | String  | Number of exposed (the number of individuals diagnosed with D1 followed by D2 in exposed group) divided by number of case (the number of individuals diagnosed with D2 in exposed group) |
| `n_exposed/n_controls`      | String  | Number of exposed (the number of individuals diagnosed with D1 not followed by D2 in exposed group) divided by number of case (the number of individuals diagnosed without D2 in exposed group) |
| `comorbidity_network_method`| String  | Method used for comorbidity network analysis |
| `describe`                  | String  | Description of the model fitting, removed covariates in the model, and reasons for removal of covariates in the model |
| `co_vars_list`              | String  | List of covariates used in the model |
| `co_vars_zvalues`           | String  | Z-values for each covariate in the model |
| `comorbidity_beta`          | Float   | Estimated coefficient from the comorbidity model |
| `comorbidity_se`            | Float   | Standard error of the estimated coefficient |
| `comorbidity_p`             | Float   | P-value for the comorbidity coefficient |
| `comorbidity_aic`           | Float   | Akaike information criterion for the model |
| `disease_d1`                | String  | Name of the disease 1 |
| `system_d1`                 | String  | Phecode disease system for the disease 1 |
| `sex_d1`                    | String  | Sex-specific for the disease 1 |
| `disease_d2`                | String  | Name of the disease 2 |
| `system_d2`                 | String  | Phecode disease system for the disease 2 |
| `sex_d2`                    | String  | Sex-specific for the disease 2 |
| `comorbidity_p_significance`| Boolean | Whether the result is statistically significant based on adjusted p-value |
| `comorbidity_p_adjusted`    | Float   | Adjusted p-value accounting for multiple comparisons |
| Additions for RPCN method |
| `alpha`     | Float   | Hyperparameter used for l1-norm (Weight multiplying the l1 penalty term) |
| Additions for PCN_PCA method |
| `pc_sum_variance_explained` | Float   | The cumulative proportion of variance in a dataset that is accounted for by a selected number of principal components in a Principal Component Analysis (sum of explained variance for principal components) |

**Result of disease trajectory analysis**

```python
print(trajectory_result)
# This will output something like (e.g., for a matched cohort study):
   phecode_d1  phecode_d2 name_disease_pair  N_exposed  n_total n_exposed/n_cases n_exposed/n_controls  ...            system_d1 sex_d1                                         disease_d2             system_d2  sex_d2  trajectory_p_significance  trajectory_p_adjusted
0       287.0       155.0       287.0-155.0      10000     1662            74/282              53/1380  ...        hematopoietic   Both         Cancer of liver and intrahepatic bile duct             neoplasms    Both                       True           1.862810e-25
1       244.0        90.0        244.0-90.0      10000     2459           117/428              59/2031  ...  endocrine/metabolic   Both  Sexually transmitted infections (not HIV or he...   infectious diseases    Both                       True           2.099407e-41
2       198.0       155.0       198.0-155.0      10000     1521            71/270              38/1251  ...            neoplasms   Both         Cancer of liver and intrahepatic bile duct             neoplasms    Both                       True           4.350023e-25
3       145.0        90.0        145.0-90.0      10000     2373           139/419             101/1954  ...            neoplasms   Both  Sexually transmitted infections (not HIV or he...   infectious diseases    Both                       True           3.468080e-44
4       722.0       557.0       722.0-557.0      10000     2514           134/446             144/2068  ...      musculoskeletal   Both              Intestinal malabsorption (non-celiac)             digestive    Both                       True           1.579863e-35
5       702.0       155.0       702.0-155.0      10000     1573            68/275              67/1298  ...         dermatologic   Both         Cancer of liver and intrahepatic bile duct             neoplasms    Both                       True           2.593236e-19
6       705.0       752.0       705.0-752.0      10000     9013          458/1749            1222/7264  ...         dermatologic   Both                Nervous system congenital anomalies  congenital anomalies    Both                       True           1.160539e-18
7       153.0       560.0       153.0-560.0      10000     2028           141/359              64/1669  ...            neoplasms   Both   Intestinal obstruction without mention of hernia             digestive    Both                       True           8.010304e-50
8       383.0       596.0       383.0-596.0      10000     5556            97/955              48/4601  ...         sense organs   Both                         Other disorders of bladder         genitourinary    Both                       True           2.532090e-36
```

| Variable Name               | Type    | Description |
|-----------------------------|---------|-------------|
| `phecode_d1`                | Float   | Phecode for disease 1 in the temporal disease pair |
| `phecode_d2`                | Float   | Phecode for disease 2 in the temporal disease pair |
| `name_disease_pair`         | String  | Name of the temporal disease pair (e.g., "D1 → D2") |
| `N_exposed`                 | Integer | Total number of individuals in exposed group |
| `n_total`                   | Integer | Number of individuals in the sub-cohort (Individuals who had no history of comorbidities with disease 1 and no history of comorbidities associated with disease 2, and met the sex-specific eligibility criteria for both diseases) |
| `n_exposed/n_cases`         | String  | Number of exposed (the number of individuals diagnosed with D1 followed by D2 in exposed group) divided by number of case (the number of individuals diagnosed with D2 in exposed group) |
| `n_exposed/n_controls`      | String  | Number of exposed (the number of individuals diagnosed with D1 not followed by D2 in exposed group) divided by number of case (the number of individuals diagnosed without D2 in exposed group) |
| `trajectory_method`         | String  | Method used for disease trajectory analysis |
| `describe`                  | String  | Description of the model fitting, removed covariates in the model, and reasons for removal of covariates in the model |
| `co_vars_list`              | String  | List of covariates included in the model |
| `co_vars_zvalues`           | String  | Z-values for each covariate in the model |
| `trajectory_beta`           | Float   | Estimated coefficient from the model |
| `trajectory_se`             | Float   | Standard error of the estimated coefficient |
| `trajectory_p`              | Float   | P-value for the coefficient |
| `trajectory_aic`            | Float   | Akaike information criterion for the model |
| `disease_d1`                | String  | Name of the disease 1 |
| `system_d1`                 | String  | Phecode disease system for the disease 1 |
| `sex_d1`                    | String  | Sex-specific for the disease 1 |
| `disease_d2`                | String  | Name of the disease 2 |
| `system_d2`                 | String  | Phecode disease system for the disease 2 |
| `sex_d2`                    | String  | Sex-specific for the disease 2 |
| `trajectory_p_significance` | Boolean | Whether the result is statistically significant based on adjusted p-value |
| `trajectory_p_adjusted`     | Float   | Adjusted p-value accounting for multiple comparisons |
| Additions for RPCN method |
| `alpha`     | Float   | Hyperparameter used for l1-norm (Weight multiplying the l1 penalty term) |
| Additions for PCN_PCA method |
| `pc_sum_variance_explained` | Float   | The cumulative proportion of variance in a dataset that is accounted for by a selected number of principal components in a Principal Component Analysis (sum of explained variance for principal components) |

#### Save the results:
All analysis outputs are standardized `pandas.DataFrame` objects, these can be exported to multiple formats (i.e., .csv, .xlsx, .feather). The following example code show how to save results to CSV files.

```python
# For example: save the results to .csv file
phewas_result.to_csv("/your/project/path/phewas_result.csv")
com_strength_result.to_csv("/your/project/path/com_strength_result.csv")
com_network_result.to_csv("/your/project/path/com_network_result.csv")
binomial_result.to_csv("/your/project/path/binomial_result.csv")
trajectory_result.to_csv("/your/project/path/trajectory_result.csv")
```

### 3.2 Step-by-step analysis

#### 3.2.1 PheWAS Analysis

The **step-by-step analysis** starts with a Phenome-wide Association Study (PheWAS) based on `DiseaseNetworkData` object, which aims to identify outcome phecode diseases significantly and positively associated with the exposure. In a standard cohort study, an unconditional cox regression model adjusted for covariates is employed to estimate the effects of exposure on each outcome phecode disease. In a matched-cohort study, a stratified cox regression model adjusting for covariates, is utilized to determine the association between exposure and each outcome phecode disease. For a exposed-only cohort study, the incidence rate of each outcome phecode disease is calculated within the overall study population.

The following example code show how to use `DiNetxify.phewas()` to PheWAS analysis in a matched-cohort study.

```python
# Reminder:
# When using multiprocessing, ensure that the code is enclosed within the following block.
# This prevents entering a never ending loop of new process creation.
if __name__ == "__main__":
    phewas_result = dnt.phewas(
        data=data,                                            
        proportion_threshold=0.01,                            
        n_process=2,                                          
        system_exl=['symptoms', 'others', 'injuries & poisonings', 'pregnancy complications'],                                                   
        covariates=['age', 'sex', 'BMI'],                                                    
        correction='bonferroni',                               
        lifelines_disable=True,                               
        log_file='/your/project/path/dep.log'
    )
```

- **data** – the `DiseaseNetworkData` object.
- **proportion_threshold** – the minimum proportion of cases within the exposed group required for a phecode to be included in the PheWAS analysis. If the proportion of cases is below this threshold, the phecode is excluded from the analysis. **proportion_threshold** and **n_threshold** are mutually exclusive. Default is `None`.
- **n_process** - number of parallel processes to use for analysis. Multiprocessing is enabled when set to greater than one. Default is `1`.
- **system_exl** - phecode systems to exclude from analysis. *Note:* Mutually exclusive with **system_inc**. Options: Same as **system_inc**. Default is `None`.
- **covariates** – list of phenotypic covariates to include in the model. Default is `None`.
- **correction** - method for p-value correction (from `statsmodels.stats.multitest.multipletests`).  
  - Options:
    - none: No correction  
    - bonferroni: One-step correction  
    - sidak: One-step correction  
    - holm-sidak: Step-down method using Sidak adjustments  
    - holm: Step-down method using Bonferroni adjustments  
    - simes-hochberg: Step-up method (independent)  
    - hommel: Closed method based on Simes tests (non-negative)  
    - fdr_bh: Benjamini/Hochberg (non-negative)  
    - fdr_by: Benjamini/Yekutieli (negative)  
    - fdr_tsbh: Two stage FDR correction (non-negative)  
    - fdr_tsbky: Two stage FDR correction (non-negative)  
  - *Reference:* [statsmodels documentation](https://www.statsmodels.org/dev/generated/statsmodels.stats.multitest.multipletests.html). Default is `bonferroni`.
- **lifelines_disable** - whether to disable lifelines. Lifelines provide more robust fitting but require longer computation time. Default is `False`.
- **log_file** - path and prefix for log file. If None, logs are written to temporary directory with prefix `DiseaseNet_phewas_`. Default is `None`.

**Optional parameters**:

- **cutoff** - significance threshold for adjusted PheWAS p-values. Default is `0.05`.
- **system_inc** - phecode systems to include in analysis. *Note:* Mutually exclusive with **system_exl**. Options: `circulatory`, `congenital anomalies`, `dermatologic`, `digestive`, `endocrine/metabolic`, `genitourinary`, `hematopoietic`, `infectious diseases`, `injuries & poisonings`, `mental disorders`, `musculoskeletal`, `neoplasms`, `neurological`, `pregnancy complications`, `respiratory`, `sense organs`, `symptoms`, `others`. Default is `None`.
- **phecode_inc** - Specific phecodes to include in analysis. *Note:* Mutually exclusive with **phecode_exl**. Default is `None`.
- **phecode_exl** - Specific phecodes to exclude from analysis. *Note:* Mutually exclusive with **phecode_inc**. Default is `None`.
- **n_threshold** - The minimum number of cases within the exposed group required for a phecode to be included in the PheWAS analysis. If the number of cases is below this threshold, the phecode is excluded. *Note:* **n_threshold** and **proportion_threshold** are mutually exclusive. Default is `None`.

**After PheWAS analysis**:

By running the `DiNetxify.phewas()` function, you will obtain a result in the format of a `pandas.DataFrame`. This result is identical to the `phewas_result` generated by the **one-step analysis** described above, and thus will not be introduced again here. In standard cohort or matched-cohort studies, outcome phecode diseases with a hazard ratio (HR) greater than 1 are typically selected for further analysis. For exposed-only cohort studies, it is usually necessary to define an empirical incidence threshold to identify significantly occurring phecode diseases for subsequent analysis. Additionally, you have the option to apply various statistical methods to perform further multiple hypothesis testing corrections.

The following example code demonstrates how to filter outcome phecode diseases with a hazard ratio (HR) greater than 1 from the PheWAS result and how to apply different statistical methods for multiple hypothesis testing corrections.

```python
# Filter results for diseases with Hazard Ratio (HR) > 1 if necessary
phewas_result = phewas_result[phewas_result['phewas_coef'] > 0]

# Redo p-value adjustment using False Discovery Rate (FDR) with other methods
# For example: Benjamini-Hochberg method
phewas_result = dnt.phewas_multipletests(
    df=phewas_result,                                     
    correction='fdr_bh',                                  
    cutoff=0.05                                           
)
```

**Save the result of PheWAS analysis**:

Same as the **one-step analysis**, the result of PheWAS analysis can be exported to multiple formats (i.e., .csv, .xlsx, .feather). The following example code show how to save the result to CSV files.

```python
# For example: save the entire PheWAS results to a CSV file
phewas_result.to_csv('/your/project/path/dep_phewas.csv')  
```

#### 3.2.2 Disease pair generation

In the step of generating disease pairs, the `DiseaseNetworkData.disease_pair()` function utilizes only two variables from the PheWAS results: the phecode disease identifier and a Boolean value (`True` or `False`) indicating whether the analysis results are significantly associated. Therefore, if you prefer to study specific phecode diseases only, you can independently create a `pandas.DataFrame` containing these two variables by any other method, thereby bypassing the initial **PheWAS analysis step entirely.

Additionally, to distinguish whether disease pairs have temporal associations, the `DiseaseNetworkData.disease_pair()` function introduces two optional parameters: **min_interval_days** and **max_interval_days**. For example, if **min_interval_days** is set to 30 days and **max_interval_days** is set to 365.25 × 5 days, temporal disease pairs (e.g., D1 → D2) must fulfill the condition that the time interval between occurrences of diseases D1 and D2 exceeds **min_interval_days** and is less than or equal to **max_interval_days**.

The following example code show how to use `DiseaseNetworkData.disease_pair()` to disease pair construction.

```python
# Disease pair construction of the `DiseaseNetworkData` object
data.disease_pair(
    phewas_result=phewas_result,               
    min_interval_days=30,                      
    max_interval_days=365.25*5,                
    force=True                              
)
```

- **phewas_result** - `pd.DataFrame` containing PheWAS analysis results produced by the `DiNetxify.phewas()` function.
- **min_interval_days** - minimum required time interval (in days) between diagnosis dates when constructing temporal D1 → D2 disease pairs. Individuals with D1 and D2 diagnoses interval ≤ this value are considered to have non-temporal pairs. Default is `0`.
- **max_interval_days** - maximum allowed time interval (in days) between diagnosis dates when constructing disease pairs. Individuals with interval more than this value are excluded from temporal analysis. Default is `np.inf`.
- **force** - if `True`, overwrites existing data attributes. If `False`, raises error when data exists. Default is `False`.

**Optional Parameters**:

- **n_process** - number of processes for parallel processing. Values more than 1 enable multiprocessing. Default is `1`.
- **phecode_col** - column name for phecode identifiers in `phewas_result`. Default is `'phecode'`.  
- **significance_col** - column name for PheWAS significance values. Default is `'phewas_p_significance'`.

**After disease pair generation**

After running the `DiseaseNetworkData.disease_pair()` function, you will obtain a `DiseaseNetworkData` object that has been updated in place, specifically generating the `DiseaseNetworkData.trajectory` attribute for subsequent analyses. Typically, once this step is completed, the resulting `DiseaseNetworkData` object should be saved in a compressed file format (e.g., .npz for NumPy-based storage or .pkl.gz for gzipped Python object persistence) to facilitate cross-platform data transfer and reproducibility of experimental results. The following example code show how to save the object to .npz file.

```python
# For example: save the updated data object with disease pairs
data.save('/your/project/path/dep_withtra')
```

#### 3.2.3 Comorbidity strength estimation

After constructing disease pairs, the updated `DiseaseNetworkData` object is used with the `DiNetxify.comorbidity_strength()` function to calculate the relative risk (RR) and phi-correlation for each disease pair, which serves as the comorbidity strength between diseases. Additionally, this function performs certain statistical computations (e.g., counting the number of individuals diagnosed with disease D1 and the number of individuals exhibiting temporal D1 → D2 disease pairs) that help determine the temporal directionality of the disease pairs.

The following example code show how to use `DiNetxify.comorbidity_strength()` function.

```python
# Reminder:
# When using multiprocessing, ensure that the code is enclosed within the following block.
# This prevents entering a never ending loop of new process creation.

if __name__ == "__main__":
    com_strength_result = dnt.comorbidity_strength(
        data=data,                                     
        proportion_threshold=0.001,                    
        n_process=2,                                   
        log_file='/your/project/path/dep.log'          
    )
```

- **data** - `DiseaseNetworkData` object containing the processed disease pairs data.
- **proportion_threshold** - minimum proportion of exposed individuals required for disease pair co-occurrence to be included in analysis. Disease pairs below this threshold are excluded. *Note:* Mutually exclusive with `n_threshold`. Default is `None`.
- **n_process** - number of parallel processes for analysis. Values >1 enable multiprocessing. Default is `1`.
- **log_file** - path/prefix for log file. If `None`, uses temporary directory with prefix `DiseaseNet_com_strength_`. Default is `None`.

**Optional Parameters**:

- **n_threshold** - minimum number of exposed individuals required for disease pair co-occurrence to be included in analysis. Disease pairs below this threshold are excluded. *Note:* Mutually exclusive with proportion_threshold. Default is `None`.
- **correction_phi** - P-value correction method for phi-correlation (from `statsmodels.stats.multitest`). 
  - Options:
    - none: No correction  
    - bonferroni: One-step correction  
    - sidak: One-step correction  
    - holm-sidak: Step-down method using Sidak adjustments  
    - holm: Step-down method using Bonferroni adjustments  
    - simes-hochberg: Step-up method (independent)  
    - hommel: Closed method based on Simes tests (non-negative)  
    - fdr_bh: Benjamini/Hochberg (non-negative)  
    - fdr_by: Benjamini/Yekutieli (negative)  
    - fdr_tsbh: Two stage FDR correction (non-negative)  
    - fdr_tsbky: Two stage FDR correction (non-negative) 
- **cutoff_phi** - significance threshold for adjusted phi-correlation p-values. Default is `0.05`.
- **correction_RR** - P-value correction method for relative risk (same methods as `correction_phi`). Default is `bonferroni`.
- **cutoff_RR** - significance threshold for adjusted RR p-values. Default is `0.05`.

**After comorbidity strength estimation**:

After running the `DiNetxify.comorbidity_strength()` function, the result will be provided in a pandas.DataFrame format, which is identical to the `com_strength_result` obtained from the **one-step analysis** described earlier. Typically, disease pairs with a RR greater than 1 and a phi-correlation greater than 0 are selected for subsequent construction of disease trajectories and comorbidity networks.

The following example code demonstrates how to filter outcome phecode diseases with a RR greater than 1 and a phi-correlation greater than 0 from the comorbidity strength estimation result and how to apply different statistical methods for multiple hypothesis testing corrections.

```python
# Further filter based on Relative Risk (RR) > 1 and phi-correlation > 0 if necessary
com_strength_result = com_strength_result[
    (com_strength_result['phi'] > 0) & (com_strength_result['RR'] > 1)
]

# Redo p-value adjustment using False Discovery Rate (FDR) with other methods
# For example: Benjamini-Hochberg method
com_strength_result = dnt.comorbidity_strength_multipletests(
    df=com_strength_result,                         
    correction_phi='fdr_bh',                       
    correction_RR='fdr_bh',                         
    cutoff_phi=0.05,                                
    cutoff_RR=0.05       
)
```

**Save the result of comorbidity strength estimation**:

Same as the **one-step analysis**, the result of comorbidity strength estimation can be exported to multiple formats (i.e., .csv, .xlsx, .feather). The following example code show how to save the result to CSV file.

```python
# For example: save the comorbidity strength estimation results to a CSV file
com_strength_result.to_csv('/your/project/path/dep_com_strength.csv')  
```

#### 3.2.4 Binomial test

After calculating the comorbidity strength for each disease pair, a **binomial test** is performed to assess the temporal directionality of the disease pairs. The **binomial test** is based on some variables of `com_strength_result` (Including `phecode_d1_col`, `phecode_d2_col`, `significance_phi_col`, `significance_RR_col`, `n_nontemporal_col`, `n_temporal_d1d2_col`, and `n_temporal_d2d1_col`) and the updated `DiseaseNetworkData` object. For example, for a D1–D2 disease pair, the test evaluates whether a temporal pattern exists in the direction of D1 → D2 or D2 → D1. The following example code show how to use `DiNetxify.binomial_test()` function.

```python
binomial_result = dnt.binomial_test(
    data=data,                                        
    comorbidity_strength_result=com_strength_result,
    n_process=1,                                      
    enforce_temporal_order=True,                      
    log_file='/your/project/path/dep.log'            
)
```

- **data** - `DiseaseNetworkData` object containing processed disease network data.
- **comorbidity_strength_result** - `DataFrame` containing comorbidity strength analysis results from `DiNetxify.comorbidity_strength()`.
- **n_process** - number of parallel processes. Note: Multiprocessing is disabled for this analysis. Default is `1`.
- **enforce_temporal_order** - if `True`, excludes individuals with non-temporal D1-D2 pairs. If `False`, includes all individuals. Default is `False`.
- **log_file** - path/prefix for log file. If `None`, uses temporary directory with prefix `DiseaseNet_binomial_test_`. Default is `None`.

**Optional Parameters**:

- **comorbidity_network_result** - `DataFrame` containing comorbidity network analysis results from `DiNetxify.comorbidity_network()`. When provided, limits binomial test to significant disease pairs. Default is `None`.
- **correction** - P-value correction method for binomial tests. 
  - Options:  
    - none: no correction  
    - bonferroni: one-step correction  
    - sidak: one-step correction  
    - holm-sidak: Step-down method using Sidak adjustments  
    - holm: Step-down method using Bonferroni adjustments  
    - simes-hochberg: Step-up method (independent)  
    - hommel: Closed method based on Simes tests (non-negative)  
    - fdr_bh: Benjamini/Hochberg (non-negative)  
    - fdr_by: Benjamini/Yekutieli (negative)  
    - fdr_tsbh: Two stage FDR correction (non-negative)  
    - fdr_tsbky: Two stage FDR correction (non-negative)  
- **cutoff** - significance threshold for adjusted binomial p-values. Default is `0.05`.
- **phecode_d1_col** - column for disease 1 phecode. Default is `'phecode_d1'`.  
- **phecode_d2_col** - column for disease 2 phecode. Default is `'phecode_d2'`.
- **n_nontemporal_col** - column for non-temporal disease pair counts. Default is `'n_d1d2_nontemporal'`.  
- **n_temporal_d1d2_col** - column for D1 → D2 temporal disease pair counts. Default is `'n_d1d2_temporal'`.  
- **n_temporal_d2d1_col** - column for D2 → D1 temporal disease pair counts. Default is `'n_d2d1_temporal'`.  
- **significance_phi_col** - column for phi-correlation significance. Default is `'phi_p_significance'`.  
- **significance_RR_col** - column for RR significance. Default is `'RR_p_significance'`.
- **significance_coef_col** - column for comorbidity significance. Default is `'comorbidity_p_significance'`.

**After binomial test**:

After running the `DiNetxify.binomial_result()` function, the results will be provided in a `pandas.DataFrame` format, which is identical to the `binomial_result` obtained from the **one-step analysis** described earlier. The disease pairs in the `binomial_result` that exhibit significant temporal directionality will be used for constructing disease trajectories.

The following example code demonstrates how to apply different statistical methods for multiple hypothesis testing corrections.

```python
# Redo p-value adjustment using False Discovery Rate (FDR) with other methods
# For example: Benjamini-Hochberg method
binomial_result = dnt.binomial_multipletests(
    df=binomial_result,                              
    correction='fdr_bh',                             
    cutoff=0.05                                      
)
```

**Save the result of binomial test**:

Same as the **one-step analysis**, the result of binomial test can be exported to multiple formats (i.e., .csv, .xlsx, .feather). The following example code show how to save the result to CSV files.

```python
# For example: save the binomial test results to a CSV file
binomial_result.to_csv('/your/project/path/dep_binomial.csv')
```

#### 3.2.5 Comorbidity network analysis

In comorbidity network analysis, disease pairs filtered using the `DiNetxify.comorbidity_strength()` function are subjected to logistic regression model to identify significantly associated disease pairs. The **comorbidity network analysis** is based on some variables of `com_strength_result` (Including `phecode_d1_col`, `phecode_d2_col`, `significance_phi_col`, `significance_RR_col`) and some variables of `binomial_result` (Including `phecode_d1_col`, `phecode_d2_col`, and `significance_binomial_col`). You can choose any of the methods (`CN`, `RPCN`, and `PCN_PCA`) provided by `DiNetxify.comorbidity_network()` function. 

The first method is the Correlation Network ("CN"), which uses an unconditional logistic regression model to assess the association of each disease pair and identify significantly related pairs. The second method is the Regularized Partial Correlation Network ("RPCN"), which also relies on an unconditional logistic regression model combined with L1 regularization. In this method, other diseases involved in non-temporal D1–D2 pairs are included as covariates (for example, when evaluating the association between D1 and D2, additional diseases such as D3 and D4 are considered as covariates) and L1 regularization is applied to reduce the influence of these additional covariates. The third method is the Partial Correlation Network with PCA ("PCN_PCA"), where principal component analysis (PCA) is applied to the additional diseases, and the top principal components are included as covariates in the unconditional logistic regression model. Empirical studies have shown that the results obtained using the "PCN_PCA" method are comparable to those from the "RPCN" method, while offering faster computation.

The following example code show how to use `DiNetxify.comorbidity_network()` function.

```python
# Reminder:
# When using multiprocessing, ensure that the code is enclosed within the following block.
# This prevents entering a never ending loop of new process creation.

if __name__ == "__main__":
    comorbidity_result = dnt.comorbidity_network(
        data=data,                                       
        comorbidity_strength_result=com_strength_result, 
        binomial_test_result=binomial_result,            
        n_process=2,                                     
        covariates=['age', 'BMI', 'sex'],                  
        method='CN',                                    
        log_file='/your/project/path/dep.log'            
    )
```

- **data** - `DiseaseNetworkData` object containing processed disease network data.
- **comorbidity_strength_result** - `DataFrame` containing comorbidity strength analysis results from `DiNetxify.comorbidity_strength()`.
- **binomial_test_result** - `DataFrame` containing binomial test analysis results from `DiNetxify.binomial_test()`. Default is `None`.
- **n_process** - number of parallel processes. Values more than 1 enable multiprocessing. Default is `1`.
- **covariates** - list of phenotypic covariates to include. If `None`, `sex` includes in **covariates**. Default is `None`.
- **method** - comorbidity network analysis method to use. 
  - Options: 
    - `RPCN`: regularized partial correlation Network. 
        - **alpha** - L1 penalty weight. Default is `None`  
        - **auto_penalty** - auto-determine alpha. Default is `True`  
        - **alpha_range** - alpha search range. Default is `(1,15)`  
        - **scaling_factor** - alpha scaling factor. Default is `1`
    - `PCN_PCA`: partial correlation network with PCA. 
        - **n_PC** - number of principal components. Default is `5`  
        - **explained_variance** - variance threshold. Default is `None`
    - `CN`: correlation network. Default is `RPCN`.
- **log_file** - path/prefix for log file. If `None`, uses temporary directory with prefix `DiseaseNet_comorbidity_network_`. Default is `None`.

**Optional parameters**:

- **correction** - P-value correction method. 
  - Options:  
    - none: no correction  
    - bonferroni: one-step correction  
    - sidak: one-step correction  
    - holm-sidak: step-down method using Sidak adjustments  
    - holm: step-down method using Bonferroni adjustments  
    - simes-hochberg: step-up method (independent)  
    - hommel: closed method based on Simes tests (non-negative)  
    - fdr_bh: benjamini/Hochberg (non-negative)  
    - fdr_by: benjamini/Yekutieli (negative)  
    - fdr_tsbh: two stage FDR correction (non-negative)  
    - fdr_tsbky: two stage FDR correction (non-negative)  
- **cutoff** - significance threshold for adjusted p-values. Default is `0.05`.
- **phecode_d1_col** - column for disease 1 phecode. Default is `phecode_d1`.  
- **phecode_d2_col** - column for disease 2 phecode. Default is `phecode_d2`.  
- **significance_phi_col** - column for phi-correlation. Default is `phi_p_significance`.  
- **significance_RR_col** - column for RR. Default is `RR_p_significance`.  
- **significance_binomial_col** - column for binomial test. Default is `binomial_p_significance`.

**After comorbidity network analysis**:

After running the `DiNetxify.comorbidity_network()` function, the result will be provided in a `pandas.DataFrame` format, which is identical to the `com_network_result` obtained from the **one-step analysis** described earlier. At this stage, the significantly associated non-temporal disease pairs have been identified, collectively forming a comorbidity network.

The following example code demonstrates how to apply different statistical methods for multiple hypothesis testing corrections.

```python
# Redo p-value adjustment using False Discovery Rate (FDR) with other methods
# For example: Benjamini-Hochberg method
comorbidity_result = dnt.comorbidity_multipletests(
    df=comorbidity_result,                            
    correction='fdr_bh',                              
    cutoff=0.05                                      
)
```

**Save the result of comorbidity network analysis**:

Same as the **one-step analysis**, the result of comorbidity network analysis can be exported to multiple formats (i.e., .csv, .xlsx, .feather). The following example code show how to save the result to CSV files.

```python
# For example: save the comorbidity network analysis results to a CSV file
comorbidity_result.to_csv('/your/project/path/dep_comorbidity.csv')
```

#### 3.2.6 Disease trajectory analysis

In disease trajectory analysis, the significantly temporal disease pairs filtered using the `DiNetxify.binomial_result()` function are subjected to logistic regression validation to identify significantly associated temporal disease pairs. The **disease trajectory analysis** is based on some variables of `com_strength_result` (Including `phecode_d1_col`, `phecode_d2_col`, `significance_phi_col`, `significance_RR_col`) and some variables of `binomial_result` (Including `phecode_d1_col`, `phecode_d2_col`, and `significance_binomial_col`). You can select any of the validation methods provided by the `DiNetxify.disease_trajectory()` function (`CN`, `RPCN`, and `PCN_PCA`).

All methods of disease trajectory analysis are based on conditional logistic regression models, and the procedures for each method are consistent with the corresponding approaches used in **comorbidity network analysis**.

```python
# Reminder:
# When using multiprocessing, ensure that the code is enclosed within the following block.
# This prevents entering a never ending loop of new process creation.

if __name__ == "__main__":
    trajectory_result = dnt.disease_trajectory(
        data=data,                                       
        comorbidity_strength_result=com_strength_result,
        binomial_test_result=binomial_result,          
        method='RPCN',                                   
        n_process=2,                                     
        matching_var_dict={'age': 2, 'sex': 'exact'},   
        matching_n=5,                                    
        enforce_time_interval=False,                   
        covariates=['age', 'BMI'],                   
        log_file='/your/project/path/dep.log'            
    )
```

- **data** - `DiseaseNetworkData` object containing processed disease network data.
- **comorbidity_strength_result** - `DataFrame` containing comorbidity strength analysis results from `DiNetxify.comorbidity_strength()`.
- **binomial_test_result** - `DataFrame` containing binomial test analysis results from `DiNetxify.binomial_test()`.
- **method** - disease trajectory method: 
 - `RPCN`: regularized partial correlation network. 
    - **alpha**: L1 penalty weight (ignored if auto_penalty). Default is `None`
    - **auto_penalty**: auto-determine optimal alpha. Default is `True`
    - **alpha_range**: alpha search range. Default is `(1,15)`
    - **scaling_factor**: alpha scaling factor. Default is `1`
 - `PCN_PCA`: partial correlation network with PCA. 
    - **n_PC**: principal components count. Default is `5`
    - **explained_variance**: variance threshold (overrides n_PC). Default is `None`
 - `CN`: correlation network. Default is `RPCN`.
- **matching_var_dict** - dictionary specifying matching variables and criteria: Categorical/binary: `{'var':'exact'}`. Continuous: `{'var':max_diff}` (scalar > 0). Always use `'sex'` for sex matching. Default is `{'sex':'exact'}`.
- **matching_n** - maximum matched controls per case. Default is `2`.
- **covariates** - list of phenotypic covariates (use `sex` for sex). Exclude matching variables. Default is `None`.
- **n_process** - number of parallel processes (>1 enables multiprocessing). Default is `1`.
- **log_file** - log file path/prefix. If `None`, uses temp dir with `DiseaseNet_trajectory_` prefix. Default is `None`.
- **enforce_time_interval** - apply min/max time intervals for D2 outcome determination. Default is `True`.

**Optional Parameters**:

- **max_n_cases** - Maximum D2 cases to include (random sampling if exceeded). Default is `np.inf`.
- **global_sampling** - `True` for single sampling across all pairs, `False` for per-pair sampling. Default is `False`.
- **correction** - P-value correction method. 
  - Options:  
    - none: no correction  
    - bonferroni: one-step correction  
    - sidak: one-step correction  
    - holm-sidak: step-down method using Sidak adjustments  
    - holm: step-down method using Bonferroni adjustments  
    - simes-hochberg: step-up method (independent)  
    - hommel: closed method based on Simes tests (non-negative)  
    - fdr_bh: benjamini/Hochberg (non-negative)  
    - fdr_by: benjamini/Yekutieli (negative)  
    - fdr_tsbh: two stage FDR correction (non-negative)  
    - fdr_tsbky: two stage FDR correction (non-negative)
- **cutoff** - significance threshold for adjusted p-values. Default is `0.05`.
- **phecode_d1_col**: disease 1 phecode column. Default is `'phecode_d1'`
- **phecode_d2_col**: disease 2 phecode column. Default is `'phecode_d2'`
- **significance_phi_col**: Phi-correlation column. Default is `'phi_p_significance'`
- **significance_RR_col**: RR column. Default is `'RR_p_significance'`
- **significance_binomial_col**: binomial test column. Default is `'binomial_p_significance'`

**After disease trajectory analysis**:

After running the `DiNetxify.disease_trajectory()` function, the result will be provided in a `pandas.DataFrame` format, which is identical to the `trajectory_result` obtained from the **one-step analysis** described earlier and will not be reintroduced here. At this stage, the significantly associated temporal disease pairs have been identified, collectively forming multiple disease trajectories.

The following example code demonstrates how to apply different statistical methods for multiple hypothesis testing corrections.

```python
# Redo p-value adjustment using False Discovery Rate (FDR) with other methods
# For example: Benjamini-Hochberg method
trajectory_result = dnt.trajectory_multipletests(
    df=trajectory_result,                             
    correction='fdr_bh',                            
    cutoff=0.05                                       
)
```

**Save the result of disease trajectory analysis**:

Same as the **one-step analysis**, the result of disease trajectory analysis can be exported to multiple formats (i.e., .csv, .xlsx, .feather). The following example code show how to save the result to CSV files.

```python
# For example: save the disease trajectory analysis analysis results to a CSV file
trajectory_result.to_csv('/your/project/path/dep_trajectory.csv')
```

## 4. Visualization

The visualization part provides four methods: PheWAS plot, comorbidity network plot, disease trajectory plot, and three dimension plot.

### 4.1 Initializing the plot object

To get started with visualization, first import the `Plot` class from `DiNetxify.visualization`. You'll then create a `Plot` object by passing in your analysis results (**PheWAS analysis**, **comorbidity network analysis**, and **disease trajectory analysis**). The following example code demonstrates how to initialize a `Plot` object for a standard cohort/matched-cohort study.

```python
# Create Plot object
from DiNetxify.visualization import Plot

# cohort/matched cohort
result_plot = Plot(
    phewas_result=phewas_result,                         
    comorbidity_network_result=comorbidity_result,       
    disease_trajectory_result=trajectory_result,         
    exposure_name='Short LTL',                                      
    exposure_size=15,                                    
    exposure_location=(0,0,0),                           
)

# or exposed-only cohort
result_plot = Plot(
    phewas_result=phewas_result,                         
    comorbidity_network_result=comorbidity_result,       
    disease_trajectory_result=trajectory_result,         
    exposure_name=None,                                       
    exposure_size=None,                                 
    exposure_location=None,                              
)
```

- **comorbidity_result** - comorbidity network analysis results (`pandas.DataFrame`) that must be included D1 disease column, D2 disease column, disease pair identifier column, beta value column of the fitted model and significance identifier of adjusted p-value column.
- **trajectory_result** - disease trajectory analysis results (`pandas.DataFrame`) that must be included D1 disease column, D2 disease column, beta value column of the fitted model and significance identifier of adjusted p-value column.
- **phewas_result** - PheWAS analysis results (`pandas.DataFrame`) that must be included phecode disease column, case count of phecode disease column, disease system column, and significance identifier of adjusted p-value column.
- **exposure_name** - name of exposure. Default is `None`. If `None`, it means that this is an exposed-only cohort study.
- **exposure_location** - custom 3D coordinates (x, y, z) for exposure node positioning. Default is `None`. If `None`, it will be auto-positioned at (0,0,0).
- **exposure_size** - relative size scaling factor for exposure node. Default is `None`. If `None`, it means that this is an exposed-only cohort study.

#### Optional Parameters:

- **source** - D1 disease column. Default is `'phecode_d1'`.
- **target** - D2 disease column. Default is `'phecode_d2'`.
- **phewas_phecode** - phecode disease column. Default `'phecode'`.
- **phewas_number** - case count of disease column. Default `'N_cases_exposed'`.
- **system_col** - disease system column. Default `'system'`.
- **col_disease_pair**: disease pair identifier column. Default `'name_disease_pair'`.
- **filter_phewas_col**: significance identifier column of PheWAS analysis results. Default `'phewas_p_significance'`.
- **filter_comorbidity_col**: significance identifier column of comorbidity network analysis results. Default `'comorbidity_p_significance'`.
- **filter_trajectory_col**: significance identifier column of disease trajectory analysis results. Default `'trajectory_p_significance'`.
- **SYSTEM** - list of phecode systems to visualize. Available systems (17 total):
  ['neoplasms', 'genitourinary', 'digestive', 'respiratory', 'infectious diseases', 'mental disorders', 'musculoskeletal',
   'hematopoietic', 'dermatologic', 'circulatory system', 'neurological', 'endocrine/metabolic', 'sense organs',
   'injuries & poisonings', 'congenital anomalies', 'symptoms', 'others']
- **COLOR** - custom colors for systems (hex, RGB, or named colors). Must match SYSTEM order. Default is 
    ['#F46D5A', '#5DA5DA', '#5EBCD1', '#C1D37F', '#CE5A57', '#A5C5D9', '#F5B36D', '#7FCDBB', '#ED9A8D',
    '#94B447', '#8C564B', '#E7CB94', '#8C9EB2', '#E0E0E0', '#F1C40F', '#9B59B6', '#4ECDC4', '#6A5ACD']

#### After initializing the Plot object:

After initializing the `Plot` object, you can use the `result_plot` object to visualization.

### 4.2 PheWAS plot

The `phewas_plot()` function provides a visualization of the phecode diseases included in the disease network analysis (**comorbidity network analysis** and **disease trajectory analysis**). For a standard cohort study or a matched-cohort study, this method displays the effect of exposure on outcome phecode diseases, represented by hazard ratios (HR). For an exposed-only cohort study, the method presents the number of individuals affected by each phecode disease. Additionally, this function supports three image formats (.png, .svg, and .jpg). The following example code demonstrates how to use `phewas_plot()` function.

```python
# phewas plot
result_plot.phewas_plot(
    path="/your/project/path/phewas_plot.png",
    is_exposure_only=False,                         
)
```
- **path** - output file path for saving the plot (including filename and extension).
- **is_exposure_only** - boolean flag indicating whether the plot is for an exposure-only cohort study. Default is `False`. If `False`, it means that this is not an exposed-only cohort study.

#### Optional parameters:

- **col_coef** - column name containing effect size coefficients (e.g., hazard ratios). Default is `'phewas_coef'`.
- **col_system** - column name containing disease system/category classifications. Default is `'system'`.
- **col_se** - column name containing standard errors for effect sizes. Default is `'phewas_se'`.
- **col_disease** - column name containing disease names/descriptions. Default is `'disease'`.
- **col_exposure** - column name containing case counts for exposed individuals. Default is `'N_cases_exposed'`.
- **disese_font_size** - column name containing font size of diseases. Default is `10`.
- **system_font_size** - column name containing font size of disease system. Default is `17`.
- **dpi** - image resolution in dots per inch for output files. Default is `200`.

### 4.3 Comorbidity network plot

The `comorbidity_network_plot()` function visualizes the associations of phecode diseases within the same comorbidity community and across different comorbidity communities. In this method, the Louvain algorithm is first applied to group each outcome phecode disease into distinct comorbidity communities based on their associations with other phecode diseases. Each community is then arranged within a separate sector of the plot, with clear spacing between different communities. In the visualization, the phecode disease systems to which each disease belongs are represented by distinct colors, allowing for a clear identification of which types of phecode diseases are predominant within each comorbidity community. 

The following example code demonstrates how to use `comorbidity_network_plot()` function.

```python
# comorbidity network visualization
result_plot.comorbidity_network_plot(
    path="/your/project/path/comorbidity_network.html"
)
```
- **path** - output file path for saving the HTML visualization (including filename and extension).  

#### Optional parameters:

- **max_radius** - maximum radial position for nodes (in pixels). Controls outer boundary of the network. Default is `90.0`.
- **min_radius** - minimum radial position for nodes (in pixels). Controls inner boundary. Default is `35.0`
- **layer_distance** - spacing between concentric circles (in pixels). Affects radial grouping. Default is `40.0`.
- **size_reduction** - scaling factor for node diameters (0-1 range). Adjusts visual prominence. Default is `0.5`
- **line_width** - stroke width for comorbidity connections (in pixels). Default is `1.0`
- **line_color** - color specification for comorbidity lines. Accepts: Named colors (e.g., `steelblue`). HEX codes (e.g., `#4682B4`). RGB tuples (e.g., `(70, 130, 180)`). Default is `'black'`
- **cluster_reduction_ratio** - compression factor (0-1) for cluster tightness. Lower values create denser groupings. Default is `0.4`.
- **cluster_weight** - edge attribute used for clustering calculations. Typically the association strength metric. Default is `'comorbidity_beta'`
- **font_style** - font family for all text elements. Use web-safe fonts or loaded font families. Default is `Times New Roman`.

### 4.4 Disease trajectory plot

The `trajectory_plot()` function provides a visualization of all temporal disease pairs within the same comorbidity community. In this visualization, multiple temporal disease pairs are sequentially connected (e.g., D1 → D2 and D2 → D3 are combined to form D1 → D2 → D3), creating a complete disease trajectory.

The following example code demonstrates how to use `trajectory_plot()` function.

```python
# Disease trajectory visualization
result_plot.trajectory_plot(
    path="/your/project/path/"
)
```

- **path** - directory path where output visualization images will be saved. *Note:* Include trailing slash for proper path resolution (e.g., `/output/plots/`)

#### Optional parameters:

- **source** - D1 disease column. Default is `'phecode_d1'`.
- **target** - D2 disease column. Default is `'phecode_d2'`.
- **dpi** - image resolution in dots per inch for output files. Default is `500`.
- **cluster_weight** specifies the edge weight metric used for network clustering calculations. Default is `'comorbidity_beta'`.

#### After disease trajectory plot:

After running the `trajectory_plot()` function, multiple PNG images will be generated, with each image displaying the disease trajectories within a specific comorbidity community.

### 4.5 Three dimension plot

The `three_dimension_plot()` function combines the visualization approaches of both the `comorbidity_network_plot()` and `trajectory_plot()` functions, displaying the two types of plots on two perpendicular planes (the x–y plane and the y–z plane). From a top-down perspective, the plot displays the comorbidity network, while from a horizontal perspective, it presents the disease trajectory.

The following example code demonstrates how to use `three_dimension_plot()` function.

```python
# three dimension visualization
result_plot.three_dimension_plot(
    path="/your/project/path/"
)
```

- **path** - absolute or relative file path to save the interactive HTML visualization.

#### Optional Parameters:

- **max_radius** - maximum radial distance (in pixels) from center for node placement. Default is `180.0`.
- **min_radius** - minimum radial distance (in pixels) from center for node placement. Default is `35.0`.
- **layer_distance** - vertical spacing (in pixels) between concentric layers. Default is `40.0`.
- **layout_width** - total figure width in pixels. Default is `900.0`.
- **layout_height** - total figure height in pixels. Default is `900.0`.
- **line_color** - color specification for trajectory pathways. Accepts: CSS color names (e.g., `steelblue`). HEX codes (e.g., `#4682B4`). RGB tuples (e.g., `(70, 130, 180)`). Default is `'black'`.
- **line_width** - stroke width (in pixels) for trajectory lines. Default is `1.0`.
- **size_reduction** - multiplicative scaling factor for node diameters (0.1-1.0). Default is `0.5`.
- **cluster_reduction_ratio** - compression factor (0.1-1.0) for cluster density. Default is `0.4`.
- **cluster_weight** - edge metric used for clustering calculations. Default is `'comorbidity_beta'`.
- **font_style** - font family for all text elements. Use web-safe fonts. Default is `'Times New Roman'`
- **font_size** - base font size in points for all text elements. Default is `15.0`.

## API Reference
### Class `DiseaseNetworkData`

```python
class DiseaseNetworkData(
    study_design: str = 'cohort',
    phecode_level: int = 1,
    min_required_icd_codes: int = 1,
    date_fmt: str = '%Y-%m-%d',
    phecode_version: str = '1.2'
)
```

A class for handling disease network data creation and operations, for use in DiNetxify module.

**Parameters:**
- `study_design` (`str`): Specify the type of study design, either "cohort", "matched cohort", or "exposed-only cohort". Defaults to `'cohort'`.
- `phecode_level` (`int`): The level of phecode to use for analysis, where level 1 (with a total of 585 medical conditions) corresponds to 3-digit ICD-10 codes and level 2 (with a total of 1257 medical conditions) to 4-digit ICD-10 codes. Level 2 phecodes offer a more granular analysis with potentially smaller sample sizes per disease category. For larger studies, level 2 phecodes may enhance result interpretation. For smaller studies, level 1 is recommended to maintain statistical power. Defaults to `1`.
- `min_required_icd_codes` (`int`): The minimum number of ICD codes mapping to a specific phecode required for the phecode to be considered valid. For example, if set to 2, a single diagnosis record will not be sufficient to count as an occurrence. Ensure that your medical records are complete (i.e., not limited to only the first occurrence for each code) when using this parameter. Defaults to `1`.
- `date_fmt` (`str`): The format of the date fields in your phenotype and medical records data. Defaults to `'%Y-%m-%d'`.
- `phecode_version` (`str`): The version of the phecode system used for converting diagnosis codes. Version 1.2 is the official version of the phecode system, with mapping files available for ICD-9-CM, ICD-9-WHO, ICD-10-CM, and ICD-10-WHO codes. While option 1.3a is provided, it’s an unofficial version and not recommended for general use. Defaults to `'1.2'`.

---

#### Instance Methods

##### `phenotype_data`
```python
phenotype_data(
    self,
    phenotype_data_path: str,
    column_names: dict,
    covariates: list,
    is_single_sex: bool = False,
    force: bool = False
) -> None
```
Load phenotype data into the object.

**Parameters:**
- `phenotype_data_path` (`str`): Path to CSV/TSV phenotype file with header row.
- `column_names` (`dict`): Mapping of dataset column names. Required keys: `'Participant ID'`, `'Index date'`, `'End date'`, `'Exposure'`, `'Sex'`, `'Match ID'`.
- `covariates` (`list`): List of additional covariate names (e.g., `['age', 'BMI']`).
- `is_single_sex` (`bool`): True if dataset contains only one sex. Defaults to `False`.
- `force` (`bool`): If True, overwrite existing data attributes. Defaults to `False`.

**Returns:**
- `None`

---

##### `Table1`
```python
Table1(
    self,
    continuous_stat_mode: str = 'auto'
) -> pd.DataFrame
```
Generate a descriptive summary table of phenotype data.

**Parameters:**
- `continuous_stat_mode` (`str`): Method for continuous variable statistics. Choices:
  - `auto`: Automatic normality-based choice.
  - `normal`: Mean and standard deviation.
  - `nonnormal`: Median and interquartile range.
  Defaults to `'auto'`.

**Returns:**
- `pd.DataFrame`

---

##### `merge_medical_records`
```python
merge_medical_records(
    self,
    medical_records_data_path: str,
    diagnosis_code: str,
    column_names: dict,
    date_fmt: str = None,
    chunksize: int = 1000000
) -> None
```
Load one or more medical records datasets.

**Parameters:**
- `medical_records_data_path` (`str`): Path to CSV/TSV medical records file.
- `diagnosis_code` (`str`): Code type: `'ICD-9-CM'`, `'ICD-9-WHO'`, `'ICD-10-CM'`, or `'ICD-10-WHO'`.
- `column_names` (`dict`): Mapping for dataset columns. Required keys: `'Participant ID'`, `'Diagnosis code'`, `'Date of diagnosis'`.
- `date_fmt` (`str`): Date format (defaults to phenotype data format). Defaults to `None`.
- `chunksize` (`int`): Rows per chunk for large files. Defaults to `1000000`.

**Returns:**
- `None`

---

##### `get_attribute`
```python
get_attribute(
    self,
    attr_name: str
) -> any
```
Retrieve the value of a private or protected attribute.

**Parameters:**
- `attr_name` (`str`): Name of the attribute to retrieve.

**Returns:**
- Attribute value (`any`)

---

##### `medical_records_to_dataframe`
```python
concat(
    self, 
    phecode_list: list,
    medical_history: bool=False
) -> DiseaseNetworkData
```
Convert stored medical records into a tidy pandas DataFrame.

**Parameters:**
- `phecode_list` (`list`): List of phecodes to extract from the medical records. Only phecodes valid for the current phecode_level are accepted.
- `medical_history` (`bool`): Include a binary history column for each phecode if set to True. Default to `False`
      

**Returns:**
- `pd.DataFrame`

---

##### `modify_phecode_level`
```python
modify_phecode_level(
    self,
    phecode_level: int
) -> None
```
Update the phecode level setting.

**Parameters:**
- `phecode_level` (`int`): New phecode level (1 or 2).

**Returns:**
- `None`

---

##### `disease_pair`
```python
disease_pair(
    self,
    phewas_result: pd.DataFrame,
    min_interval_days: int = 0,
    max_interval_days: float = float('inf'),
    force: bool = False,
    n_process: int = 1,
    **kwargs
) -> None
```
Construct temporal and non-temporal disease pairs.

**Parameters:**
- `phewas_result` (`pd.DataFrame`): DataFrame from `phewas()`.
- `min_interval_days` (`int`): Minimum days between diagnoses. Defaults to `0`.
- `max_interval_days` (`float`): Maximum days between diagnoses. Defaults to `inf`.
- `force` (`bool`): Overwrite existing data. Defaults to `False`.
- `n_process` (`int`): Number of parallel processes. Defaults to `1`.
- `**kwargs`: Additional mappings:
  - `phecode_col` (`str`): Column for phecode. Defaults to `'phecode'`.
  - `significance_col` (`str`): Column for significance. Defaults to `'phewas_p_significance'`.

**Returns:**
- `None`

---

##### `save`
```python
save(
    self,
    file: str
) -> None
```
Save object state to a gzip-compressed pickle file (`.pkl.gz`).

**Parameters:**
- `file` (`str`): Filename or prefix (adds `.pkl.gz`).

**Returns:**
- `None`

---

##### `load`
```python
load(
    self,
    file: str,
    force: bool = False
) -> None
```
Load object state from a gzip-compressed pickle file.

**Parameters:**
- `file` (`str`): Filename or prefix (adds `.pkl.gz`).
- `force` (`bool`): Overwrite if True. Defaults to `False`.

**Returns:**
- `None`

---

##### `save_npz`
```python
save_npz(
    self,
    file: str
) -> None
```
Save object state to a NumPy `.npz` file.

**Parameters:**
- `file` (`str`): Filename or prefix (adds `.npz`).

**Returns:**
- `None`

---

##### `load_npz`
```python
load_npz(
    self,
    file: str,
    force: bool = False
) -> None
```
Load object state from a NumPy `.npz` file.

**Parameters:**
- `file` (`str`): Filename or prefix (adds `.npz`).
- `force` (`bool`): Overwrite if True. Defaults to `False`.

**Returns:**
- `None`

### Analysis Functions

#### Function: `disease_network_pipeline`

```python
disease_network_pipeline(
    data: DiseaseNetworkData,
    n_process: int,
    n_threshold_phewas: int,
    n_threshold_comorbidity: int,
    output_dir: str,
    project_prefix: str,
    keep_positive_associations: bool = False,
    save_intermediate_data: bool = False,
    system_exl: list = None,
    pipeline_mode: str = 'v1',
    method: str = 'RPCN',
    covariates: list = None,
    matching_var_dict: dict = {'sex':'exact'},
    matching_n: int = 2,
    min_interval_days: int = 0,
    max_interval_days: float = float('inf'),
    enforce_temporal_order: bool = False,
    correction: str = 'bonferroni',
    cutoff: float = 0.05,
    **kwargs
) -> dict
```

**Parameters:**
- `data` (`DiseaseNetworkData`): The DiseaseNetworkData object.
- `n_process` (`int`): Specifies the number of parallel processes to use. Defaults to required.
- `n_threshold_phewas` (`int`): Minimum cases in exposed group for PheWAS inclusion. Passed to `phewas()`.
- `n_threshold_comorbidity` (`int`): Minimum co-occurrences for comorbidity strength. Passed to `comorbidity_strength()`.
- `output_dir` (`str`): Directory path for pipeline outputs.
- `project_prefix` (`str`): Prefix for naming outputs.
- `keep_positive_associations` (`bool`): Retain only positive associations. Defaults to `False`.
- `save_intermediate_data` (`bool`): Save intermediate data objects. Defaults to `False`.
- `system_exl` (`list`): Phecode systems to exclude. Defaults to `None`.
- `pipeline_mode` (`str`): Analysis order mode (`'v1'` or `'v2'`). Defaults to `'v1'`.
- `method` (`str`): Comorbidity network / trajectory method (`'RPCN'`, `'PCN_PCA'`, `'CN'`). Defaults to `'RPCN'`.
- `covariates` (`list`): Covariates for models. Defaults to `None`.
- `matching_var_dict` (`dict`): Matching variables and criteria. Defaults to `{'sex':'exact'}`.
- `matching_n` (`int`): Number of matched controls per case. Defaults to `2`.
- `min_interval_days` (`int`): Minimum days between diagnoses. Defaults to `0`.
- `max_interval_days` (`float`): Maximum days between diagnoses. Defaults to `inf`.
- `enforce_temporal_order` (`bool`): Enforce temporal order in testing. Defaults to `False`.
- `correction` (`str`): p-value correction method. Defaults to `'bonferroni'`.
- `cutoff` (`float`): Significance threshold. Defaults to `0.05`.
- `**kwargs`:
  - `alpha` (`float`): L1 penalty weight. Defaults per method.
  - `auto_penalty` (`bool`): Auto-select alpha. Defaults to `True`.
  - `alpha_range` (`tuple`): Search range for alpha. Defaults to `(1,15)`.
  - `scaling_factor` (`float`): Scaling factor for alpha. Defaults to `1`.
  - `n_PC` (`int`): Number of principal components. Defaults to `5`.
  - `explained_variance` (`float`): Variance threshold for PCs.

**Returns:**  
- `dict`: Summary of significant results count.

---

#### Function: `phewas`

```python
phewas(
    data: DiseaseNetworkData,
    covariates: list = None,
    proportion_threshold: float = None,
    n_threshold: int = None,
    n_process: int = 1,
    correction: str = 'bonferroni',
    cutoff: float = 0.05,
    system_inc: list = None,
    system_exl: list = None,
    phecode_inc: list = None,
    phecode_exl: list = None,
    log_file: str = None,
    lifelines_disable: bool = False
) -> pd.DataFrame
```

**Parameters:**
- `data` (`DiseaseNetworkData`): Input data object.
- `covariates` (`list`): Phenotypic covariates. Defaults to `None`.
- `proportion_threshold` (`float`): Minimum proportion of cases. Mutually exclusive with `n_threshold`. Defaults to `None`.
- `n_threshold` (`int`): Minimum case count. Mutually exclusive with `proportion_threshold`. Defaults to `None`.
- `n_process` (`int`): Parallel processes. Defaults to `1`.
- `correction` (`str`): p-value correction method. Defaults to `'bonferroni'`.
- `cutoff` (`float`): Significance threshold. Defaults to `0.05`.
- `system_inc` (`list`): Systems to include. Defaults to `None`.
- `system_exl` (`list`): Systems to exclude. Defaults to `None`.
- `phecode_inc` (`list`): Specific phecodes to include. Defaults to `None`.
- `phecode_exl` (`list`): Specific phecodes to exclude. Defaults to `None`.
- `log_file` (`str`): Log file prefix. Defaults to `None`.
- `lifelines_disable` (`bool`): Disable lifelines. Defaults to `False`.

**Returns:**  
- `pd.DataFrame`: PheWAS results.

---

#### Function: `phewas_multipletests`

```python
phewas_multipletests(
    df: pd.DataFrame,
    correction: str = 'bonferroni',
    cutoff: float = 0.05
) -> pd.DataFrame
```

**Parameters:**
- `df` (`pd.DataFrame`): Input results from `phewas()`.
- `correction` (`str`): p-value correction method. Defaults to `'bonferroni'`.
- `cutoff` (`float`): Significance threshold. Defaults to `0.05`.

**Returns:**  
- `pd.DataFrame`: Adjusted results.

#### Function: `comorbidity_strength`

```python
comorbidity_strength(
    data: DiseaseNetworkData,
    proportion_threshold: float = None,
    n_threshold: int = None,
    n_process: int = 1,
    log_file: str = None,
    correction_phi: str = 'bonferroni',
    cutoff_phi: float = 0.05,
    correction_RR: str = 'bonferroni',
    cutoff_RR: float = 0.05
) -> pd.DataFrame
```

**Parameters:**

- `data` (*DiseaseNetworkData*): DiseaseNetworkData object.
- `proportion_threshold` (*float*): The minimum proportion of individuals in the exposed group in which a disease pair must co-occur (temporal or non-temporal) to be included in the comorbidity strength estimation. If the proportion of co-occurrence is below this threshold, the disease pair is excluded from the analysis. proportion_threshold and n_threshold are mutually exclusive.
- `n_threshold` (*int*): The minimum number of individuals in the exposed group in which a disease pair must co-occur (temporal or non-temporal) to be included in the comorbidity strength estimation. If the number of co-occurrences is below this threshold, the disease pair is excluded from the analysis. n_threshold and proportion_threshold are mutually exclusive.
- `n_process` (*int, default=1*): Specifies the number of parallel processes to use for the analysis. Multiprocessing is enabled when `n_process` is set to a value greater than one.
- `correction_phi` (*str, default='bonferroni'*): Method for phi-correlation p-value correction from the statsmodels.stats.multitest.multipletests. 
  - Available methods are:
    - none : no correction 
    - bonferroni : one-step correction 
    - sidak : one-step correction 
    - holm-sidak : step down method using Sidak adjustments 
    - holm : step-down method using Bonferroni adjustments 
    - simes-hochberg : step-up method (independent) 
    - hommel : closed method based on Simes tests (non-negative) 
    - fdr_bh : Benjamini/Hochberg (non-negative) 
    - fdr_by : Benjamini/Yekutieli (negative) fdr_tsbh : two stage fdr correction (non-negative) 
    - fdr_tsbky : two stage fdr correction (non-negative) 
  - See https://www.statsmodels.org/dev/generated/statsmodels.stats.multitest.multipletests.html for more details.
- `cutoff_phi` (*float, default=0.05*): The significance threshold for adjusted phi-correlatio p-values.
- `correction_RR` (*str, default='bonferroni'*): Method for RR p-value correction from the statsmodels.stats.multitest.multipletests. 
  - Available methods are:
    - none : no correction 
    - bonferroni : one-step correction 
    - sidak : one-step correction 
    - holm-sidak : step down method using Sidak adjustments 
    - holm : step-down method using Bonferroni adjustments 
    - simes-hochberg : step-up method (independent) 
    - hommel : closed method based on Simes tests (non-negative) 
    - fdr_bh : Benjamini/Hochberg (non-negative) 
    - fdr_by : Benjamini/Yekutieli (negative) fdr_tsbh : two stage fdr correction (non-negative) 
    - fdr_tsbky : two stage fdr correction (non-negative) 
  - See https://www.statsmodels.org/dev/generated/statsmodels.stats.multitest.multipletests.html for more details.
- `cutoff_RR` (*float, default=0.05*): The significance threshold for adjusted RR p-values.
- `log_file` (*str, default=None*): Path and prefix for the text file where log will be recorded. If None, the log will be written to the temporary files directory with file prefix of DiseaseNet_com_strength_.


#### Function: `comorbidity_strength_multipletests`

```python
comorbidity_strength_multipletests(
    df: pd.DataFrame,
    correction_phi: str = 'bonferroni',
    cutoff_phi: float = 0.05,
    correction_RR: str = 'bonferroni',
    cutoff_RR: float = 0.05
) -> pd.DataFrame
```

**Parameters:**

- `df` (*pd.DataFrame*): DataFrame containing the results from the comorbidity_strength function.
- `correction_phi` (*str, default='bonferroni'*): Method for phi-correlation p-value correction from the statsmodels.stats.multitest.multipletests.   
  - Available methods are:
    - none : no correction 
    - bonferroni : one-step correction 
    - sidak : one-step correction 
    - holm-sidak : step down method using Sidak adjustments 
    - holm : step-down method using Bonferroni adjustments 
    - simes-hochberg : step-up method (independent) 
    - hommel : closed method based on Simes tests (non-negative) 
    - fdr_bh : Benjamini/Hochberg (non-negative) 
    - fdr_by : Benjamini/Yekutieli (negative) fdr_tsbh : two stage fdr correction (non-negative) 
    - fdr_tsbky : two stage fdr correction (non-negative) 
  - See https://www.statsmodels.org/dev/generated/statsmodels.stats.multitest.multipletests.html for more details.
- `cutoff_phi` (*float, default=0.05*): The significance threshold for adjusted phi-correlatio p-values.
- `correction_RR` (*str, default='bonferroni'*): Method for RR p-value correction from the statsmodels.stats.multitest.multipletests. 
  - Available methods are:
    - none : no correction 
    - bonferroni : one-step correction 
    - sidak : one-step correction 
    - holm-sidak : step down method using Sidak adjustments 
    - holm : step-down method using Bonferroni adjustments 
    - simes-hochberg : step-up method (independent) 
    - hommel : closed method based on Simes tests (non-negative) 
    - fdr_bh : Benjamini/Hochberg (non-negative) 
    - fdr_by : Benjamini/Yekutieli (negative) fdr_tsbh : two stage fdr correction (non-negative) 
    - fdr_tsbky : two stage fdr correction (non-negative) 
  - See https://www.statsmodels.org/dev/generated/statsmodels.stats.multitest.multipletests.html for more details.
- `cutoff_RR` (*float, default=0.05*): The significance threshold for adjusted RR p-values.

#### Function: `binomial_test`

```python
binomial_test(
    data: DiseaseNetworkData,
    comorbidity_strength_result: pd.DataFrame,
    comorbidity_network_result: pd.DataFrame = None,
    n_process: int = 1,
    log_file: str = None,
    correction: str = 'bonferroni',
    cutoff: float = 0.05,
    enforce_temporal_order: bool = False,
    **kwargs
) -> pd.DataFrame
```

**Parameters:**

- `data` (*DiseaseNetworkData*): DiseaseNetworkData object.
- `comorbidity_strength_result` (*pd.DataFrame*): DataFrame containing comorbidity strength analysis results produced by the 'DiNetxify.comorbidity_strength' function.
- `comorbidity_network_result` (*pd.DataFrame, default=None*): DataFrame containing comorbidity network analysis results produced by the 'DiNetxify.comorbidity_network' function. When provided, the binomial test is limited to disease pairs deemed significant in the comorbidity network analysis.
- `n_process` (*int, default=1*): Multiprocessing is disabled for this analysis.
- `correction` (*str, default='bonferroni'*): Method for binomial p-value correction from the statsmodels.stats.multitest.multipletests. 
  - Available methods are:
    - none : no correction 
    - bonferroni : one-step correction 
    - sidak : one-step correction 
    - holm-sidak : step down method using Sidak adjustments 
    - holm : step-down method using Bonferroni adjustments 
    - simes-hochberg : step-up method (independent) 
    - hommel : closed method based on Simes tests (non-negative) 
    - fdr_bh : Benjamini/Hochberg (non-negative) 
    - fdr_by : Benjamini/Yekutieli (negative) fdr_tsbh : two stage fdr correction (non-negative) 
    - fdr_tsbky : two stage fdr correction (non-negative) 
  - See https://www.statsmodels.org/dev/generated/statsmodels.stats.multitest.multipletests.html for more details.
- `cutoff` (*float, default=0.05*): The significance threshold for adjusted binomial p-values.
- `log_file` (*str, default=None*): Path and prefix for the text file where log will be recorded. If None, the log will be written to the temporary files directory with file prefix of DiseaseNet_binomial_test_.
- `enforce_temporal_order` (*bool, default=False*): If True, exclude individuals with non-temporal D1-D2 pair when performing the test. If False, include all individuals, including those with non-temporal D1-D2 pair. 
- `**kwargs` 
  - `phecode_d1_col` : str, default='phecode_d1' Name of the column in 'comorbidity_strength_result' and 'comorbidity_network_result' that specifies the phecode identifiers for disease 1 of the disease pair. 
  - `phecode_d2_col` : str, default='phecode_d2' Name of the column in 'comorbidity_strength_result' and 'comorbidity_network_result' that specifies the phecode identifiers for disease 2 of the disease pair. 
  - `n_nontemporal_col` : str, default='n_d1d2_nontemporal' Name of the column in 'comorbidity_strength_result' that specifies the number of individuals with non-temporal d1-d2 disease pair 
  - `n_temporal_d1d2_col` : str, default='n_d1d2_temporal' Name of the column in 'comorbidity_strength_result' that specifies the number of individuals with temporal d1->d2 disease pair.
  - `n_temporal_d2d1_col` : str, default='n_d2d1_temporal' Name of the column in 'comorbidity_strength_result' that specifies the number of individuals with temporal d2->d1 disease pair. 
  - `significance_phi_col` : str, default='phi_p_significance' Name of the column in 'comorbidity_strength_result' that indicates the significance of phi-correlation for each disease pair. 
  - `significance_RR_col` : str, default='RR_p_significance' Name of the column in 'comorbidity_strength_result' that indicates the significance of RR for each disease pair. 
  - `significance_coef_col` : str, default='comorbidity_p_significance' Name of the column in 'comorbidity_network_result' that indicates the significance of comorbidity network analysis for each disease pair.


#### Function: `binomial_multipletests`

```python
binomial_multipletests(
    df: pd.DataFrame,
    correction: str = 'bonferroni',
    cutoff: float = 0.05
) -> pd.DataFrame
```

**Parameters:**

- `df` (*pd.DataFrame*): DataFrame containing the results from the comorbidity_strength function.
- `correction` (*str, default='bonferroni'*): Method for binomial p-value correction from the statsmodels.stats.multitest.multipletests. 
  - Available methods are:
    - none : no correction 
    - bonferroni : one-step correction 
    - sidak : one-step correction 
    - holm-sidak : step down method using Sidak adjustments 
    - holm : step-down method using Bonferroni adjustments 
    - simes-hochberg : step-up method (independent) 
    - hommel : closed method based on Simes tests (non-negative) 
    - fdr_bh : Benjamini/Hochberg (non-negative) 
    - fdr_by : Benjamini/Yekutieli (negative) fdr_tsbh : two stage fdr correction (non-negative) 
    - fdr_tsbky : two stage fdr correction (non-negative) 
  - See https://www.statsmodels.org/dev/generated/statsmodels.stats.multitest.multipletests.html for more details.
- `cutoff` (*float, default=0.05*): The significance threshold for adjusted binomial p-values.


#### Function: `comorbidity_network`
```python
comorbidity_network(
    data: DiseaseNetworkData,
    comorbidity_strength_result: pd.DataFrame,
    binomial_test_result: pd.DataFrame = None,
    method: str = 'RPCN',
    covariates: list = None,
    n_process: int = 1,
    log_file: str = None,
    correction: str = 'bonferroni',
    cutoff: float = 0.05,
    **kwargs
) -> pd.DataFrame
```

**Parameters:**

- `data` (*DiseaseNetworkData*): DiseaseNetworkData object.
- `comorbidity_strength_result` (*pd.DataFrame*): DataFrame containing comorbidity strength analysis results produced by the 'DiNetxify.comorbidity_strength' function.
- `binomial_test_result` (*pd.DataFrame, default=None*): DataFrame containing binomial test analysis results produced by the 'DiNetxify.binomial_test' function.
- `method` (*str, default='RPCN'*): Specifies the comorbidity network analysis method to use. Choices are: - 'RPCN: Regularized Partial Correlation Network. - 'PCN_PCA: Partial Correlation Network with PCA. - 'CN': Correlation Network. **Additional Options for RPCN:** - 'alpha' : non-negative scalar The weight multiplying the l1 penalty term for other diseases covariates. Ignored if 'auto_penalty' is enabled. - 'auto_penalty' : bool, default=True If 'True', automatically determine the optimal 'alpha' based on model AIC value. - 'alpha_range' : tuple, default=(1,15) When 'auto_penalty' is True, search the optimal 'alpha' in this range. - 'scaling_factor' : positive scalar, default=1 The scaling factor for the alpha when 'auto_penalty' is True. **Additional Options for PCN_PCA:** - 'n_PC' : int, default=5 Fixed number of principal components to include in each model. - 'explained_variance' : float Determines the number of principal components based on the cumulative explained variance. Overrides 'n_PC' if specified.
- `covariates` (*list, default=None*): List of phenotypic covariates to include in the model. By default, includes ['sex'] and all covariates specified in the 'DiNetxify.DiseaseNetworkData.phenotype_data()' function. To include the required variable sex as a covariate, always use 'sex' instead of its original column name. For other covariates specified in the 'DiNetxify.DiseaseNetworkData.phenotype_data()' function, use their original column names.
- `n_process` (*int, default=1*): Specifies the number of parallel processes to use for the analysis. Multiprocessing is enabled when `n_process` is set to a value greater than one.
- `correction` (*str, default='bonferroni'*): Method for binomial p-value correction from the statsmodels.stats.multitest.multipletests. 
  - Available methods are:
    - none : no correction 
    - bonferroni : one-step correction 
    - sidak : one-step correction 
    - holm-sidak : step down method using Sidak adjustments 
    - holm : step-down method using Bonferroni adjustments 
    - simes-hochberg : step-up method (independent) 
    - hommel : closed method based on Simes tests (non-negative) 
    - fdr_bh : Benjamini/Hochberg (non-negative) 
    - fdr_by : Benjamini/Yekutieli (negative) fdr_tsbh : two stage fdr correction (non-negative) 
    - fdr_tsbky : two stage fdr correction (non-negative) 
  - See https://www.statsmodels.org/dev/generated/statsmodels.stats.multitest.multipletests.html for more details.
- `cutoff` (*float, default=0.05*): The significance threshold for adjusted comorbidity network analysis p-values.
- `log_file` (*str, default=None*): Path and prefix for the text file where log will be recorded. If None, the log will be written to the temporary files directory with file prefix of DiseaseNet_comorbidity_network_. 
- `**kwargs`  
  - `phecode_d1_col` : str, default='phecode_d1' Name of the column in 'comorbidity_strength_result' and 'binomial_test_result' that specifies the phecode identifiers for disease 1 of the disease pair. 
  - `phecode_d2_col` : str, default='phecode_d2' Name of the column in 'comorbidity_strength_result' and 'binomial_test_result' that specifies the phecode identifiers for disease 2 of the disease pair. 
  - `significance_phi_col` : str, default='phi_p_significance' Name of the column in 'comorbidity_strength_result' that indicates the significance of phi-correlation for each disease pair. 
  - `significance_RR_col` : str, default='RR_p_significance' Name of the column in 'comorbidity_strength_result' that indicates the significance of RR for each disease pair. 
  - `significance_binomial_col` : str default='binomial_p_significance' Name of the column in 'binomial_test_result' that indicates the significance of binomial test for each disease pair.
  - `alpha` : non-negative scalar The weight multiplying the l1 penalty term for other diseases covariates. Ignored if 'auto_penalty' is enabled. 
  - `auto_penalty` : bool, default=True If 'True', automatically determines the best 'alpha' based on model AIC value. 
  - `alpha_range` : tuple, default=(1,15) When 'auto_penalty' is True, search the optimal 'alpha' in this range. 
  - `scaling_factor` : positive scalar, default=1 The scaling factor for the alpha when 'auto_penalty' is True.
  - `n_PC` : int, default=5 Fixed number of principal components to include in each model. 
  - `explained_variance` : float Cumulative explained variance threshold to determine the number of principal components. Overrides `'n_PC'` if specified.


#### Function: `comorbidity_multipletests`

```python
comorbidity_multipletests(
    df: pd.DataFrame,
    correction: str = 'bonferroni',
    cutoff: float = 0.05
) -> pd.DataFrame
```

**Parameters:**

- `df` (*pd.DataFrame*): DataFrame containing the results from the 'comorbidity_network' function.
- `correction` (*str, default='bonferroni'*): Method for binomial p-value correction from the statsmodels.stats.multitest.multipletests. 
  - Available methods are:
    - none : no correction 
    - bonferroni : one-step correction 
    - sidak : one-step correction 
    - holm-sidak : step down method using Sidak adjustments 
    - holm : step-down method using Bonferroni adjustments 
    - simes-hochberg : step-up method (independent) 
    - hommel : closed method based on Simes tests (non-negative) 
    - fdr_bh : Benjamini/Hochberg (non-negative) 
    - fdr_by : Benjamini/Yekutieli (negative) fdr_tsbh : two stage fdr correction (non-negative) 
    - fdr_tsbky : two stage fdr correction (non-negative) 
  - See https://www.statsmodels.org/dev/generated/statsmodels.stats.multitest.multipletests.html for more details.
- `cutoff` (*float, default=0.05*): The significance threshold for adjusted binomial p-values.


#### Function: `disease_trajectory`

```python
disease_trajectory(
    data: DiseaseNetworkData,
    comorbidity_strength_result: pd.DataFrame,
    binomial_test_result: pd.DataFrame,
    method: str = 'RPCN',
    matching_var_dict: dict = {'sex':'exact'},
    matching_n: int = 2,
    max_n_cases: float = np.inf,
    global_sampling: bool = False,
    covariates: list = None,
    n_process: int = 1,
    log_file: str = None,
    correction: str = 'bonferroni',
    cutoff: float = 0.05,
    **kwargs
) -> pd.DataFrame
```

**Parameters:**

- `data` (*DiseaseNetworkData*): DESCRIPTION.
- `comorbidity_strength_result` (*pd.DataFrame*): DataFrame containing comorbidity strength analysis results produced by the 'DiNetxify.comorbidity_strength' function.
- `binomial_test_result` (*pd.DataFrame*): DataFrame containing binomial test analysis results produced by the 'DiNetxify.binomial_test' function.
- `method` (*str, default='RPCN'*): Specifies the comorbidity network analysis method to use. Choices are: 
  - `'RPCN'`: Regularized Partial Correlation Network.
  - `'PCN_PCA'`: Partial Correlation Network with PCA. 
  - `'CN'`: Correlation Network. 

- `matching_var_dict` (*dict, default={'sex':'exact'}*): Specifies the matching variables and the criteria used for incidence density sampling. For categorical and binary variables, the matching criteria should always be `'exact'`. For continuous variables, provide a scalar greater than 0 as the matching criterion, indicating the maximum allowed difference when matching. To include the required variable sex as a matching variable, always use `'sex'` instead of its original column name. For other covariates specified in the `DiNetxify.DiseaseNetworkData.phenotype_data()` function, use their original column names.
- `matching_n` (*int, default=2*): Specifies the maximum number of matched controls for each case.
- `max_n_cases` (*int, default=np.inf*): Specifies the maximum number of D2 cases to include in the analysis. If the number of D2 cases exceeds this value, a random sample of cases will be selected.
- `global_sampling` (*bool, default=False*): Indicates whether to perform independent incidence density sampling for each D1→D2 pair (if False), or to perform a single incidence density sampling for all Dx→D2 pairs with separate regression models for each D1→D2 pair (if True). Global sampling is recommended when processing large datasets, though it might reduce result heterogeneity.
- `covariates` (*list, default=None*): List of phenotypic covariates to include in the model. By default, includes all covariates specified in the `DiNetxify.DiseaseNetworkData.phenotype_data()` function. Categorical and binary variables used for matching should not be included as covariates. Continuous variables used for matching can be included as covariates, but caution is advised. To include the required variable sex as a covariate, always use `sex` instead of its original column name. For other covariates specified in the `DiNetxify.DiseaseNetworkData.phenotype_data()` function, use their original column names.
- `n_process` (*int, default=1*): Specifies the number of parallel processes to use for the analysis. Multiprocessing is enabled when `n_process` is set to a value greater than one.
- `correction` (*str, default='bonferroni'*): Method for binomial p-value correction from the statsmodels.stats.multitest.multipletests. 
  - Available methods are:
    - none : no correction 
    - bonferroni : one-step correction 
    - sidak : one-step correction 
    - holm-sidak : step down method using Sidak adjustments 
    - holm : step-down method using Bonferroni adjustments 
    - simes-hochberg : step-up method (independent) 
    - hommel : closed method based on Simes tests (non-negative) 
    - fdr_bh : Benjamini/Hochberg (non-negative) 
    - fdr_by : Benjamini/Yekutieli (negative) fdr_tsbh : two stage fdr correction (non-negative) 
    - fdr_tsbky : two stage fdr correction (non-negative) 
  - See https://www.statsmodels.org/dev/generated/statsmodels.stats.multitest.multipletests.html for more details.
- `cutoff` (*float, default=0.05*): The significance threshold for adjusted comorbidity network analysis p-values.
- `log_file` (*str, default=None*): Path and prefix for the text file where log will be recorded. If None, the log will be written to the temporary files directory with file prefix of DiseaseNet_trajectory_. 
- `**kwargs` Analysis option 
  - `enforce_time_interval` : bool, default=True If set to True, applies the specified minimum and maximum time intervals when determining the D2 outcome among individuals diagnosed with D1. These time interval requirements should be defined using the `DiNetxify.DiseaseNetworkData.disease_pair()` function.
  - `phecode_d1_col` : str, default='phecode_d1' Name of the column in `comorbidity_strength_result` and `binomial_test_result` that specifies the phecode identifiers for disease 1 of the disease pair.
  - `phecode_d2_col` : str, default='phecode_d2' Name of the column in `comorbidity_strength_result` and `binomial_test_result` that specifies the phecode identifiers for disease 2 of the disease pair. 
  - `significance_phi_col` : str, default='phi_p_significance' Name of the column in `comorbidity_strength_result` that indicates the significance of phi-correlation for each disease pair. 
  - `significance_RR_col` : str, default='RR_p_significance' Name of the column in `comorbidity_strength_result` that indicates the significance of RR for each disease pair. 
  - `significance_binomial_col` : str default='binomial_p_significance' Name of the column in `binomial_test_result` that indicates the significance of binomial test for each disease pair. 
  - `alpha` : non-negative scalar The weight multiplying the l1 penalty term for other diseases covariates. Ignored if `auto_penalty` is enabled. 
  - `auto_penalty` : bool, default=True If `True`, automatically determines the best `alpha` based on model AIC value. 
  - `alpha_range` : tuple, default=(1,15) When `auto_penalty` is True, search the optimal `alpha` in this range. 
  - `scaling_factor` : positive scalar, default=1 The scaling factor for the `alpha` when 'auto_penalty' is True. 
  - `n_PC` : int, default=5 Fixed number of principal components to include in each model. 
  - `explained_variance` : float Cumulative explained variance threshold to determine the number of principal components. Overrides `'n_PC'` if specified.

#### Function: `trajectory_multipletests`

```python
trajectory_multipletests(
    df: pd.DataFrame,
    correction: str = 'bonferroni',
    cutoff: float = 0.05
) -> pd.DataFrame
```

**Parameters:**

- `df` (*pd.DataFrame*): DataFrame containing the results from the 'disease_trajectory' function.
- `correction` (*str, default='bonferroni'*): Method for binomial p-value correction from the statsmodels.stats.multitest.multipletests. 
  - Available methods are:
    - none : no correction 
    - bonferroni : one-step correction 
    - sidak : one-step correction 
    - holm-sidak : step down method using Sidak adjustments 
    - holm : step-down method using Bonferroni adjustments 
    - simes-hochberg : step-up method (independent) 
    - hommel : closed method based on Simes tests (non-negative) 
    - fdr_bh : Benjamini/Hochberg (non-negative) 
    - fdr_by : Benjamini/Yekutieli (negative) fdr_tsbh : two stage fdr correction (non-negative) 
    - fdr_tsbky : two stage fdr correction (non-negative) 
  - See https://www.statsmodels.org/dev/generated/statsmodels.stats.multitest.multipletests.html for more details.
- `cutoff` (*float, default=0.05*): The significance threshold for adjusted binomial p-values.

### Class `Plot`

```python
class Plot(
    comorbidity_result: pd.DataFrame,
    trajectory_result: pd.DataFrame,
    phewas_result: pd.DataFrame,
    exposure_name: str = None,
    exposure_location: Tuple[float, float, float] = None,
    exposure_size: float = None,
    source: str = 'phecode_d1',
    target: str = 'phecode_d2',
    phewas_phecode: str = 'phecode',
    phewas_number: str = 'N_cases_exposed',
    system_col: str = 'system',
    col_disease_pair: str = 'name_disease_pair',
    filter_phewas_col: str = 'phewas_p_significance',
    filter_comorbidity_col: str = 'comorbidity_p_significance',
    filter_trajectory_col: str = 'trajectory_p_significance',
)
```

A class for integrating and visualizing disease relationships from PHEWAS, comorbidity network, and trajectory analyses.

**Constructor Parameters:**
- `comorbidity_result` (`pd.DataFrame`): Non-temporal disease pairs with association metrics and significance flag.
- `trajectory_result` (`pd.DataFrame`): Temporal disease pairs (source→target) with metrics and significance flag.
- `phewas_result` (`pd.DataFrame`): PheWAS results including phecode, effect sizes, case counts, and system classifications.
- `exposure_name` (`float`, optional): Name of exposure. Default is `None`. If `None`, it means that this is an exposed-only cohort study.
- `exposure_location` (`Tuple[float, float, float]`, optional): 3D coordinates for exposure node. Defaults to origin if `None`.
- `exposure_size` (`float`, optional): Scaling factor for exposure node size. Defaults to automatic.
- `source` (`str`): Column name for source disease (default: `'phecode_d1'`).
- `target` (`str`): Column name for target disease (default: `'phecode_d2'`).
- `phewas_phecode` (`str`): Column for phecode in PHEWAS results (default: `'phecode'`).
- `phewas_number` (`str`): Column for case counts (default: `'N_cases_exposed'`).
- `system_col` (`str`): Column for disease system (default: `'system'`).
- `col_disease_pair` (`str`): Column for pair identifier (default `'name_disease_pair'`).
- `filter_phewas_col` (`str`): Column for PHEWAS significance filter.
- `filter_comorbidity_col` (`str`): Column for comorbidity significance filter.
- `filter_trajectory_col` (`str`): Column for trajectory significance filter.
- `**kwargs`  
  - `SYSTEM` (`List[str]`, optional): List of systems to visualize; defaults to all from PHEWAS if `None`.
  - `COLOR` (`List[str]`, optional): Colors corresponding to systems; default palette used if `None`.

---

#### Methods

##### `three_dimension_plot`

```python
three_dimension_plot(
    self,
    path: str,
    max_radius: float = 180.0,
    min_radius: float = 35.0,
    line_color: str = 'black',
    line_width: float = 1.0,
    size_reduction: float = 0.5,
    cluster_reduction_ratio: float = 0.4,
    cluster_weight: str = 'comorbidity_beta',
    layer_distance: float = 40.0,
    layout_width: float = 900.0,
    layout_height: float = 900.0,
    font_style: str = 'Times New Roman',
    font_size: float = 15.0
) -> None
```

Generate and save a 3D interactive HTML visualization.

**Parameters:**
- `path`: File path to save the HTML visualization
- `max_radius`: Maximum radial distance for node placement (default: `180.0`)
- `min_radius`: Minimum radial distance for node placement (default: `35.0`)
- `line_color`: Color for trajectory lines (default: `"black"`)
- `line_width`: Width for trajectory lines (default: `1.0`)
- `size_reduction`: Scaling factor for node sizes (default: `0.5`)
- `cluster_reduction_ratio`: Cluster compression factor for layout (default: `0.4`)
- `cluster_weight`: Edge weight metric used for clustering (default: `"comorbidity_beta"`)
- `layer_distance`: Vertical distance between layers (default: `40.0`)
- `layout_width`: Figure width in pixels (default: `900.0`)
- `layout_height`: Figure height in pixels (default: `900.0`)
- `font_style`: Font family for text elements (default: `'Times New Roman'`)
- `font_size`: Base font size in points (default: `15.0`)

---

##### `comorbidity_network_plot`

```python
comorbidity_network_plot(
    self,
    path: str,
    max_radius: float = 180.0,
    min_radius: float = 35.0,
    size_reduction: float = 0.5,
    cluster_reduction_ratio: float = 0.4,
    cluster_weight: str = 'comorbidity_beta',
    line_width: float = 1.0,
    line_color: str = 'black',
    layer_distance: float = 40.0,
    font_style: str = 'Times New Roman'
) -> None
```

Generate and save a 2D HTML visualization of the comorbidity network.

**Parameters:**
- `path`: Output file path for saving HTML visualization
- `max_radius`: Maximum radial position for nodes (default: `90.0`)
- `min_radius`: Minimum radial position for nodes (default: `35.0`)
- `size_reduction`: Scaling factor for node sizes (default: `0.5`)
- `cluster_reduction_ratio`: Compression factor for cluster layout (default: `0.4`)
- `cluster_weight`: Edge weight metric for clustering (default: `"comorbidity_beta"`)
- `line_width`: Width of comorbidity lines (default: `1.0`)
- `line_color`: Color of comorbidity lines (default: `"black"`)
- `layer_distance`: Distance between concentric circles (default: `40.0`)
- `font_style`: Font family for text elements (default: `"Times New Roman"`)

---

##### `trajectory_plot`

```python
trajectory_plot(
    self,
    path: str,
    cluster_weight: str = 'comorbidity_beta',
    source: str='phecode_d1',
    target: str='phecode_d2',
    dpi: float=500
) -> None
```

Generate and save trajectory plots per cluster as (.png files).

**Parameters:**
- `path`: Directory path to save output images
- `cluster_weight`: Edge weight metric used for clustering (default: `"comorbidity_beta"`)
- `source`: Column name representing source nodes (disease onset points) in trajectory data (default: `'phecode_d1'`)
- `target`: Column name representing target nodes (subsequent disease points) in trajectory data (default: `'phecode_d2'`)
- `dpi`: Image resolution in dots per inch for output files (default: `500`)
---

##### `phewas_plot`

```python
phewas_plot(
    self,
    path: str,
    system_font_size: float=17,
    disese_font_size: float=10,
    col_coef: str = 'phewas_coef',
    col_system: str = 'system',
    col_se: str = 'phewas_se',
    col_disease: str = 'disease',
    is_exposure_only: bool = False,
    col_exposure: str = 'N_cases_exposed',
    dpi: float=200
) -> None
```

Creates a polar bar plot visualizing disease associations across different disease categories (systems)

**Parameters:**
- `path`: Output file path for saving the plot
- `system_font_size`: Font size for disease system/category labels (default: `17`)
- `disease_font_size`: Font size for disease labels (default: `10`)
- `col_coef`: Column name for effect size coefficients (default: `"phewas_coef"`)
- `col_system`: Column name for disease system/category (default: `"system"`)
- `col_se`: Column name for standard errors (default: `"phewas_se"`)
- `col_disease`: Column name for disease names (default: `"disease"`)
- `is_exposure_only`: Identifier of exposure (default: `False`)
- `col_exposure`: Column name for exposure number (default: `"N_cases_exposed"`)
- `dpi`: Image resolution in dots per inch for output files (default: `200`)