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

## 3. Three-dimensional disease network analysis

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

The `disease_network_pipeline` will return a tuple of DataFrames: in order, these correspond to **phewas_result**, **com_strength_result** (comorbidity strength), **com_network_result** (comorbidity network), **binomial_result**, and **trajectory_result**. In addition to returning these, the function writes out certain results and logs to files in `output_dir` for your records.

**Result of PheWAS analysis**

| Variable Name           | Type      | Description |
|------------------------|-----------|-------------|
| `phecode`              | String    | Disease code (Phecode) used in PheWAS analysis |
| `disease`              | String    | Disease name corresponding to the Phecode |
| `system`               | String    | Phecode disease system corresponding to the Phecode (e.g., infectious diseases) |
| `sex`                  | String    | Sex-specificity of the disease (e.g., Both, Male, Female) |
| `N_cases_exposed`      | Integer   | Number of individuals diagnosed with the disease in the exposed group |
| `describe`             | String    | Descriptions of the model fitting state and removed covariates with reasons |
| `exposed_group`        | String    | Incidence rate (unit: per 1,000 person-years) in the exposed group |
| `unexposed_group`      | String    | Incidence rate (unit: per 1,000 person-years) in the unexposed group |
| `phewas_coef`          | Float     | Estimated coefficient from the model |
| `phewas_se`            | Float     | Standard error of the estimated coefficient |
| `phewas_p`             | Float     | P-value indicating statistical significance of the coefficient |
| `phewas_p_significance`| Boolean   | Indicates whether the result is statistically significant based on adjusted p-value (True/False) |
| `phewas_p_adjusted`    | Float     | Adjusted p-value accounting for multiple comparisons |

**Result of comorbidity strength estimation**

| Variable Name            | Type    | Description |
|-------------------------|---------|-------------|
| `phecode_d1`            | Integer | Phecode for disease 1 in the disease pair |
| `phecode_d2`            | Integer | Phecode for disease 2 in the disease pair |
| `name_disease_pair`     | String  | Name of the disease pair (format: "D1-D2") |
| `N_exposed`             | Integer | Total number of individuals in exposed group |
| `n_total`               | Integer | Number of exposed individuals included in the sub-cohort that meet the sex-specificity eligibility criteria for both diseases and after excluding those with history of either disease 1, disease 2, or related diseases |
| `n_d1d2_diagnosis`      | Integer | Number of individuals diagnosed with both diseases |
| `n_d1_diagnosis`        | Integer | Number of individuals diagnosed with disease 1 |
| `n_d2_diagnosis`        | Integer | Number of individuals diagnosed with disease 2 |
| `n_d1d2_nontemporal`    | Integer | Number of individuals diagnosed with both disease 1 and disease 2 but without defined temporal order (i.e., the time interval between the two diagnosis is smaller than or equal to `min_interval_days` or larger `max_interval_days`) |
| `n_d1d2_temporal`       | Integer | Number of individuals diagnosed with disease 1 followed by disease 2 in a defined temporal order (i.e., the time interval between the two diagnosis is larger than `min_interval_days` and smaller than or equal to`max_interval_days`) |
| `n_d2d1_temporal`       | Integer | Number of individuals diagnosed with disease 2 followed by disease 1 in a defined temporal order (i.e., the time interval between the two diagnosis is larger than `min_interval_days` and smaller than or equal to`max_interval_days`) |
| `phi_coef`              | Float   | Phi coefficient (φ), Pearson’s correlations for two binary variables |
| `phi_p`                 | Float   | P-value for Phi coefficient significance |
| `RR`                    | Float   | Relative risk of observing both conditions in the same individual relative to expectation |
| `RR_p`                  | Float   | P-value for relative risk |
| `phi_p_adjusted`        | Float   | Adjusted P-value for Phi coefficient (multiple comparisons) |
| `RR_p_adjusted`         | Float   | Adjusted P-value for relative risk (multiple comparisons) |
| `phi_p_significance`    | Boolean | Whether the Phi is statistically significant based on adjusted p-value |
| `RR_p_significance`     | Boolean | Whether the RR is statistically significant based on adjusted p-value |
| `disease_d1`            | String  | Name of disease 1 |
| `system_d1`             | String  | Phecode disease system related to disease 1 |
| `sex_d1`                | String  | Sex-specificity of disease 1 |
| `disease_d2`            | String  | Name of disease 2 |
| `system_d2`             | String  | Phecode disease system related to disease 2 |
| `sex_d2`                | String  | Sex-specificity of disease 2 |

**Result of binomial test**

| Variable Name               | Type    | Description |
|----------------------------|---------|-------------|
| `phecode_d1`               | Float   | Phecode for disease 1 in the temporal disease pair |
| `phecode_d2`               | Float   | Phecode for disease 2 in the temporal disease pair |
| `name_disease_pair`        | String  | Name of the temporal disease pair (e.g., D1->D2) |
| `n_d1d2_nontemporal`       | Float   | Number of individuals diagnosed with both disease 1 and disease 2 but without defined temporal order (i.e., the time interval between the two diagnosis is smaller than or equal to `min_interval_days` or larger `max_interval_days`) |
| `n_d1d2_temporal`          | Float   | Number of individuals diagnosed with disease 1 followed by disease 2 in a defined temporal order (i.e., the time interval between the two diagnosis is larger than `min_interval_days` and smaller than or equal to`max_interval_days`) |
| `n_d2d1_temporal`          | Float   | Number of individuals diagnosed with disease 2 followed by disease 1 in a defined temporal order (i.e., the time interval between the two diagnosis is larger than `min_interval_days` and smaller than or equal to`max_interval_days`) |
| `binomial_p`               | Float   | P-value from the binomial test for directionality |
| `binomial_proportion`      | Float   | Proportion of successful outcomes in the binomial test |
| `binomial_proportion_ci`   | String  | Confidence interval for the binomial proportion |
| `disease_d1`               | String  | Name of disease 1 |
| `system_d1`                | String  | Phecode disease system for disease 1 |
| `sex_d1`                   | String  | Sex-specificity of disease 1 |
| `disease_d2`               | String  | Name of disease 2 |
| `system_d2`                | String  | Phecode disease system for disease 2 |
| `sex_d2`                   | String  | Sex-specificity of disease 2 |
| `binomial_p_significance` | Boolean | Indicates whether the result is statistically significant based on adjusted p-value |
| `binomial_p_adjusted`     | Float   | Adjusted p-value for multiple comparisons |

**Result of comorbidity network analysis**

| Variable Name                | Type    | Description |
|-----------------------------|---------|-------------|
| `phecode_d1`                | Float   | Phecode for disease 1 in the non-temporal disease pair |
| `phecode_d2`                | Float   | Phecode for disease 2 in the non-temporal disease pair |
| `name_disease_pair`         | String  | Name of the non-temporal disease pair (e.g., "D1-D2") |
| `N_exposed`                 | Integer | Total number of individuals in exposed group |
| `n_total`                   | Integer | Number of exposed individuals included in the sub-cohort that meet the sex-specificity eligibility criteria for both diseases and after excluding those with history of either disease 1, disease 2, or related diseases |
| `n_exposed/n_cases`         | String  | Number of exposed individuals (individuals with diagnosis of  D1) among cases (individuals with diagnosis of D2) |
| `n_exposed/n_controls`      | String  | Number of exposed individuals (individuals with diagnosis of D1) among controls (individuals without diagnosis of D2) |
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
| `sex_d1`                    | String  | Sex-specificity of the disease 1 |
| `disease_d2`                | String  | Name of the disease 2 |
| `system_d2`                 | String  | Phecode disease system for the disease 2 |
| `sex_d2`                    | String  | Sex-specificity of the disease 2 |
| `comorbidity_p_significance`| Boolean | Whether the result is statistically significant based on adjusted p-value |
| `comorbidity_p_adjusted`    | Float   | Adjusted p-value accounting for multiple comparisons |
| **Columns for RPCN method** |||
| `alpha`     | Float   | Hyperparameter used for l1-norm (Weight multiplying the l1 penalty term) |
| **Columns for PCN_PCA method** |||
| `pc_sum_variance_explained` | Float   | The cumulative proportion of variance that is accounted for by a selected number of principal components in a Principal Component Analysis (sum of explained variance for principal components) |

**Result of disease trajectory analysis**

| Variable Name               | Type    | Description |
|-----------------------------|---------|-------------|
| `phecode_d1`                | Float   | Phecode for disease 1 in the temporal disease pair |
| `phecode_d2`                | Float   | Phecode for disease 2 in the temporal disease pair |
| `name_disease_pair`         | String  | Name of the temporal disease pair (e.g., "D1 → D2") |
| `N_exposed`                 | Integer | Total number of individuals in exposed group |
| `n_total`                   | Integer | Number of exposed individuals included in the nested case-control dataset, where eligible cases (with diagnosis of D2) are all selected and matched with specified number of controls using incidence density sampling |
| `n_exposed/n_cases`         | String  | Number of exposed individuals (individuals with diagnosis of  D1) among cases (individuals with diagnosis of D2) |
| `n_exposed/n_controls`      | String  | Number of exposed individuals (individuals with diagnosis of  D1) among cases (individuals with diagnosis of D2) |
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
| `sex_d1`                    | String  | Sex-specificity of the disease 1 |
| `disease_d2`                | String  | Name of the disease 2 |
| `system_d2`                 | String  | Phecode disease system for the disease 2 |
| `sex_d2`                    | String  | Sex-specificity of the disease 2 |
| `trajectory_p_significance` | Boolean | Whether the result is statistically significant based on adjusted p-value |
| `trajectory_p_adjusted`     | Float   | Adjusted p-value accounting for multiple comparisons |
| **Columns for RPCN method** |||
| `alpha`     | Float   | Hyperparameter used for l1-norm (weight multiplying the l1 penalty term) |
| **Additions for PCN_PCA method** |||
| `pc_sum_variance_explained` | Float   | The cumulative proportion of variance in a dataset that is accounted for by a selected number of principal components in a Principal Component Analysis (sum of explained variance for principal components) |

### 3.2 Step-by-step analysis

In many cases, you might want to run each part of the analysis separately — for example, to inspect intermediate outputs, adjust parameters for individual steps, or run alternative filtering between steps. ***DiNetxify*** allows this by exposing individual functions for each analysis stage. Below we illustrate a step-by-step approach performing the same overall analysis as the one-step pipeline above. We will reuse the `data` object already loaded.

#### 3.2.1 PheWAS Analysis

To perform a PheWAS (Phenome-Wide Association Study) using ***DiNetxify***, use the `dnt.phewas()` function. This function will identify which diseases (phecodes) are significantly associated with the exposure of interest. In a cohort or matched cohort design, it runs Cox proportional hazards models (stratified by matching in the latter) and reports hazard ratios (HR) and p-values for each disease comparing exposed vs unexposed. In an exposed-only design, since there is no comparison group, this function simply flags diseases that exceed a certain incidence threshold.

The example below demonstrates a PheWAS in a matched cohort study, adjusting for covariates and using parallel processing:

```python
# Reminder: if using n_process > 1, wrap calls in if __name__ == "__main__":

phewas_result = dnt.phewas(  
    data=data,                                             # our DiseaseNetworkData object  
    covariates=['BMI', 'age'],                             # adjust for BMI and age (and sex by default)  
    proportion_threshold=0.01,                             # require at least 1% of exposed have the disease  
    n_process=2,                                           # use 2 processes for parallel model fitting  
    correction='bonferroni',                               # multiple testing correction method  
    cutoff=0.05,                                           # significance threshold  
    system_inc=None,                                       # (optional) include only certain systems  
    system_exl=None,                                       # (optional) exclude certain systems  
    phecode_inc=None,                                      # (optional) include only certain phecodes  
    phecode_exl=None,                                      # (optional) exclude certain phecodes  
    log_file=None,                                         # (optional) log file prefix  
    lifelines_disable=False                                # whether to disable lifelines (Cox model library)  
)  

```

In this call:

- We specified `proportion_threshold=0.01` which means a phecode must have at least 1% of exposed individuals as cases to be considered. (This is an alternative to using `n_threshold`; you generally use one or the other, not both.)
- We left inclusion/exclusion lists as None, meaning we are analyzing all diseases.
- We enabled bonferroni correction on p-values and set a cutoff of 0.05 for significance.
- `lifelines_disable=False` means we use the lifelines Cox model (which is slower but handles certain edge cases better); if set to True, a custom faster method might be used but potentially less robust.

Running `dnt.phewas()` returns a DataFrame (`phewas_result`). In a cohort design, this DataFrame will include columns such as: phecode, disease name, number of cases in exposed and unexposed, hazard ratio (HR) and its confidence interval, p-value, and adjusted p-value (q-value). For a matched cohort, similar output with stratified Cox results. For an exposed-only design, the output will likely include counts and perhaps a placeholder measure, primarily to identify which diseases meet the incidence threshold.

By examining `phewas_result`, you can see which diseases are significantly associated with your exposure of interest. Typically, for downstream analysis, you might focus on diseases with HR > 1 and q-value < cutoff (if you want only positive associations), which is exactly what the `keep_positive_associations=True` setting in the pipeline would enforce.

After obtaining `phewas_result`, if you want to apply a different multiple testing correction or filter in a custom way, ***DiNetxify*** provides convenience functions. For example, you can use `dnt.phewas_multipletests()` to adjust p-values in `phewas_result` using a specified method (if you didn’t do it inside `phewas()` or want to try a different method):

```python
# Example: applying a different correction (e.g., FDR) to the PheWAS results  
phewas_result = dnt.phewas_multipletests(  
    df=phewas_result,  
    correction='fdr_bh',   # Benjamini-Hochberg FDR  
    cutoff=0.05  
) 
```

This will add/update columns in `phewas_result` for adjusted p-values and significance flags according to the chosen method and cutoff. (The `df` parameter is just the DataFrame from `phewas()`, so you can call this on any similar results DataFrame)

#### 3.2.2 Disease pair generation

After identifying which diseases are associated with the exposure (through PheWAS), the next step is to generate all possible disease pairs from that set of diseases for further analysis. In ***DiNetxify***, this is done using the `DiseaseNetworkData.disease_pair()` method, which operates on the `data` object. It will produce all combinations of the significant diseases (or all diseases if you choose) for each individual, distinguishing between *temporal* pairs (D1 before D2) and *non-temporal* co-occurrence. You can then filter these pairs by applying criteria like minimum co-occurrence counts or correlation as done in the next step.

The `disease_pair()` method requires the PheWAS result DataFrame to know which diseases to consider. It also allows specifying time interval constraints. For example:

```python
# Generate disease pairs for further analysis  
pair_data = data.disease_pair(  
    phewas_result=phewas_result,     # DataFrame from PheWAS (to get the list of diseases)  
    min_interval_days=0,             # minimum days between diagnoses for D1->D2 (0 = no minimum)  
    max_interval_days=float('inf'),  # maximum days between diagnoses to consider (inf = no limit)  
    force=False                      # if data already has pair info, this prevents overwrite unless True  
) 
```

- **phewas_result** – the DataFrame of PheWAS results. Typically you would filter this to only include diseases you want to carry forward (e.g., all with p < 0.05 and HR > 1 if focusing on positive associations). If you pass the full `phewas_result`, the method might by default filter internally to significant ones depending on implementation. Check the content of `phewas_result` and perhaps subset it if needed.
- **min_interval_days**, **max_interval_days** – these define the time window for considering temporal relationships. As defined earlier, here we’ve left them at 0 and infinity which means we include all occurrences and do not impose a maximum gap. If you wanted to only consider, say, disease pairs where events occur within 5 years of each other, you could set `max_interval_days=1825` (approximately 5*365).
- **force** – similar to earlier methods, if False it will not recompute pairs if they were already computed before for this data object (to avoid unnecessary re-processing). Use True to force regeneration.

The result `pair_data` (if returned, or it might store internally and return None – consult the actual function’s behavior) would contain information about each disease pair found in each person’s data. Typically, however, you won’t directly use this raw list; instead, the subsequent function `comorbidity_strength()` will use the `data` object which now has the disease pair info to calculate metrics.

> Note: In the one-step pipeline, disease_pair generation is handled internally; here we are making it explicit

#### 3.2.3 Comorbidity strength estimation

The next step is to assess the **comorbidity strength** of each disease pair – essentially measuring how strongly the two diseases are associated with each other in a cross-sectional sense (regardless of time order). ***DiNetxify***’s `comorbidity_strength()` function calculates statistics like **relative risk (RR)** and **Phi coefficient (Φ)** for each pair, and can perform filtering and significance testing.

Using the `data` object (which now has disease pairs from the previous step), we can run:

```python
# Reminder: if using n_process > 1, wrap calls in if __name__ == "__main__":
com_strength_result = dnt.comorbidity_strength(  
    data=data,  
    proportion_threshold=None,    # Alternatively, could require a certain prevalence  
    n_threshold=100,             # Only consider pairs that co-occur in at least 100 exposed individuals (same as we used above)  
    n_process=1,                 # Single process (this function is usually fast; can set >1 if needed)  
    log_file=None,               # Log file (optional)  
    correction_phi='bonferroni', # Correction method for Phi coefficient p-values  
    cutoff_phi=0.05,             # Significance cutoff for Phi  
    correction_RR='bonferroni',  # Correction method for RR p-values  
    cutoff_RR=0.05               # Significance cutoff for RR  
)  

```

Important points about `comorbidity_strength()`:

- It uses the disease pairs present in `data` (so ensure you called `disease_pair()` first).
- It focuses typically on the *exposed group* for calculations (since the design is cohort-based; if an exposed-only design, then it’s just that group).
- **n_threshold** here serves to filter out infrequent pairs (at least 100 co-occurrences as we set). This is analogous to `n_threshold_comorbidity` in the pipeline.
- It will compute metrics like:
  - **RR (Relative Risk)** of the two diseases co-occurring in the same person, compared to what would be expected if independent (often simplified as co-occurrence probability vs product of marginal probabilities).
  - **Φ (Phi correlation)** which is basically the Pearson correlation for two binary variables (disease present/absent).
  - P-values for these metrics (likely via chi-square tests or similar) and then it applies corrections (bonferroni in our case) to give adjusted p-values.
- The output `com_strength_result` will be a DataFrame where each row is a disease pair, with columns such as: disease1, disease2, count of individuals with both (possibly separate counts in exposed/unexposed), RR, RR_p-value, RR_q-value, Phi, Phi_p-value, Phi_q-value, etc., plus maybe some identifiers.

After obtaining `com_strength_result`, you might want to filter it to retain only those pairs that meet certain criteria. The typical default (reflected in pipeline when `keep_positive_associations=True`) is to keep only pairs with RR > 1, Phi > 0, and significance. However, you can examine this DataFrame and decide. A convenience function `comorbidity_strength_multipletests()` exists if you need to re-adjust p-values with a different method.

#### 3.2.4 Binomial test

For each disease pair, if we want to determine if there is a **preferred temporal order** (i.e., does D1 tend to occur before D2 more often than vice versa), ***DiNetxify*** uses a binomial test. Essentially, for every individual who has both diseases, we check how many had D1 first vs D2 first, across all individuals. Under the null hypothesis of no preferred order, these are like coin flips; the binomial test checks if one order is significantly more common than the other.

The function `dnt.binomial_test()` performs this analysis. It uses the pair information in `data` (so again ensure `disease_pair()` was run). We can call it as:

```python
binomial_result = dnt.binomial_test(  
    data=data,  
    enforce_temporal_order=False,  # If True, exclude individuals with ties/non-temporal occurrences  
    min_interval_days=0,          # (same as before)  
    max_interval_days=float('inf'),  
    correction='bonferroni',      # p-value correction method  
    cutoff=0.05                   # significance threshold  
)  
```

- **enforce_temporal_order** – If True, the function will exclude any “non-temporal” pairs (i.e., individuals where the two diseases occurred on the same day or essentially simultaneously, as defined by min_interval_days) from the counts. It also will ensure that if an individual has both orders of occurrence (e.g., D1 then D2 and later D2 then D1 due to multiple episodes), those might be handled to not bias the test. Since we set it False here, we’re being more permissive (which is fine given our min_interval_days=0, meaning simultaneous is allowed and counted in whichever order they happened first).
- **min_interval_days, max_interval_days** – These should match what was used in `disease_pair()` to ensure consistency in what “temporal” means.
- **correction, cutoff** – Correction for multiple tests (there’s one binomial test per disease pair) and significance threshold for the adjusted p-value.

The `binomial_result` DataFrame will list disease pairs (likely identified by some code or name) with the number of individuals where D1 happened first vs D2 first, the binomial p-value and adjusted p-value, and possibly an indication of which order is predominant. Typically, one might filter this to pairs where p < 0.05, indicating a significant temporal bias, and also note the direction (e.g., D1 -> D2 is the more frequent sequence). This information will be used in building trajectories and in the next network analysis.

#### 3.2.5 Comorbidity network analysis

Now we move to constructing the **comorbidity network** — a graph where nodes are diseases and edges represent associations between diseases, adjusting for covariates and other diseases per the method selected. This is effectively a deeper analysis on the set of disease pairs that survived earlier filters, focusing on direct relationships.

The function `dnt.comorbidity_network()` carries out this analysis. Depending on the `method` (RPCN, PCN_PCA, or CN), it will fit either a series of regularized regressions or simpler correlations. Here’s how we might call it (mirroring our pipeline parameters):

```python
# Reminder: if using n_process > 1, wrap calls in if __name__ == "__main__":

com_network_result = dnt.comorbidity_network(  
    data=data,  
    method="RPCN",           # or "PCN_PCA" or "CN"  
    covariates=['BMI', 'age'],        # adjust for these covariates (plus sex) in each model  
    correction='bonferroni',          # correction for multiple testing of edges  
    cutoff=0.05,                      # significance threshold  
    **{  
        "alpha": None,        # RPCN-specific kwargs: if we wanted to manually set  
        "auto_penalty": True, # whether to auto-select alpha  
        "alpha_range": (1, 15),  
        "scaling_factor": 1,  
        # If using PCN_PCA: we could provide "n_PC" or "explained_variance" as kwargs  
    }  
)  

```

- **method** – The network inference method: `'RPCN'` (default, uses L1-penalized partial correlation), `'PCN_PCA'` (partial correlation with dimensionality reduction), or `'CN'` (simple correlation network).
- **covariates** – Covariates to adjust for in the disease-disease association models. This should include 'sex' if you want to adjust for sex (note: in the code, sex is likely automatically included; but to be safe you could include it, though in some contexts they say "sex is always required, use 'sex'"). We include BMI and age as in PheWAS.
- **correction, cutoff** – multiple testing correction for the edge significance and the significance threshold. Each edge (disease pair) will get a p-value (from logistic regression coefficients, etc.), and these will be corrected. We use Bonferroni and 0.05.
- **kwargs** – RPCN and PCN_PCA have additional hyperparameters:
  - For RPCN: `alpha` (L1 penalty weight; if None or auto_penalty True, the code will try to auto-tune it), `auto_penalty` (if True, find optimal alpha via AIC), `alpha_range` (range to search for alpha), `scaling_factor` (scales alpha search). We left `alpha=None` and `auto_penalty=True` which means the software will choose the best penalty.
  - For PCN_PCA: `n_PC` (number of principal components to use; default 5) and/or `explained_variance` (if set, override n_PC to take enough PCs to explain this fraction of variance). Not used in RPCN.
- The output `com_network_result` will be a DataFrame capturing which disease pairs (edges) are considered significant in this network analysis. It likely includes columns like: disease1, disease2, coefficient (beta or correlation), p-value, adjusted p-value, etc., and an indicator of significance. Additionally, it might list communities if computed (though community detection is usually a separate visualization step).

After this, you could filter `com_network_result` to only significant edges (adjusted p < 0.05) if you want to focus on the final network structure. Non-significant edges might be dropped for interpretation.

#### 3.2.6 Disease trajectory analysis

Finally, we conduct the **disease trajectory analysis**. This aims to identify sequences of diseases, i.e., if having disease D1 increases the risk of subsequently developing disease D2 (beyond just co-occurrence). This is typically done with a nested case-control approach or similar incidence density sampling that was set up by matching_var_dict and matching_n earlier. Essentially, for each candidate pair from the comorbidity analysis, we treat the first-occurring disease as an exposure and the second as an outcome in a survival analysis framework to test the temporal association (while accounting for matching and covariates).

The function `dnt.disease_trajectory()` performs this analysis. We will call it on the filtered set of disease pairs. In pipeline mode 'v2', that set is restricted to those significant in the network; in 'v1', it’s basically all pairs from the comorbidity strength step (perhaps already filtered by RR/Phi positivity). You can decide which input to give it (e.g., you might filter `com_strength_result` or use `com_network_result`). Here, assuming we want to analyze trajectories for all pairs identified as having some comorbidity strength (since we didn’t filter by network first in our pipeline example):

```python
# Reminder: if using n_process > 1, wrap calls in if __name__ == "__main__":

trajectory_result = dnt.disease_trajectory(  
    data=data,  
    covariates=['BMI', 'age'],         # adjust for these covariates (and sex implicitly)  
    matching_var_dict={'sex': 'exact'}, # matching variables and criteria (as used earlier)  
    matching_n=2,                      # number of matched controls per case  
    enforce_time_interval=False,       # whether to enforce min/max intervals in analysis  
    correction='bonferroni',           # p-value correction method  
    cutoff=0.05                        # significance threshold  
)  

```

- **covariates** – Covariates to include in the models (again, include things like BMI, age; sex is always considered). These models are typically conditional logistic regressions or Cox models within matched sets, depending on how implemented.
- **matching_var_dict** and **matching_n** – These should match how we set them before generating the controls (the pipeline or earlier code ensures consistent matching was applied to identify controls; by providing the same parameters, we ensure it uses those controls). In our example, we matched on sex exactly and allowed up to 2 controls per case.
- **enforce_time_interval** – If True, the function would ensure that only disease pairs within the specified min/max intervals are considered in the modeling. Since we set it False, it will include all pairs (the min/max were already applied when generating pairs if at all).
- **correction, cutoff** – multiple testing correction for the trajectory analysis p-values and the significance threshold.

The output `trajectory_result` will be a DataFrame of disease pairs with metrics indicating temporal association significance. Likely, for each pair (D1 -> D2), it provides an estimated effect size (perhaps an odds ratio or hazard ratio for D1 leading to D2), confidence interval, p-value, and adjusted p-value, plus an indicator if it’s significant. Only pairs that passed previous filtering are analyzed, and among those, some will show a significant temporal relationship. In essence, this final result pinpoints which disease pairs have a directionality (D1 significantly predisposes to D2).

After completing these steps, you now have:

- `phewas_result`: diseases associated with exposure
- `com_strength_result`: disease pairs with comorbidity strength metrics
- `binomial_result`: which pairs have a dominant order (not yet incorporating covariates, just raw order bias)
- `com_network_result`: network edges that are significant controlling for others (direct relationships)
- `trajectory_result`: disease pairs with significant temporal relationships (risk from one to the other)

These correspond to the outputs of the one-step `disease_network_pipeline`. You can proceed to interpret them, and also to visualize them using the Plot class as described next.

## 4. Visualization

After performing the analyses, ***DiNetxify*** offers visualization tools to help interpret the results. The visualizations tie together the findings from PheWAS, comorbidity network, and disease trajectory analyses in a coherent way. The primary class for visualization is `dnt.visualization.Plot`, which takes the result DataFrames and generates interactive plots.

### 4.1 Initializing the Plot object

To create visualizations, first initialize a `Plot` object with your analysis results. You will pass in the PheWAS, comorbidity network, and trajectory results DataFrames, and optionally specify how the “exposure” (the primary risk factor or cohort definition) should appear in the plots. For example:

```python
from DiNetxify.visualization import Plot  

# Suppose phewas_result, com_network_result, trajectory_result are obtained from above  
result_plot = Plot(  
    phewas_result=phewas_result,  
    comorbidity_result=com_network_result,  
    trajectory_result=trajectory_result,  
    exposure_name='Short LTL',    # Name of the exposure (for labeling the exposure node). Use None for exposed-only cohorts.  
    exposure_size=15,            # Relative size scaling for the exposure node (to make it prominent). None for exposed-only cohorts.  
    exposure_location=(0, 0, 0)  # 3D coordinates for the exposure node. If None, it defaults to (0,0,0).  
)  

# If this were an exposed-only cohort (no explicit exposure variable), you would set:  
result_plot = Plot(  
    phewas_result=phewas_result,  
    comorbidity_result=com_network_result,  
    trajectory_result=trajectory_result,  
    exposure_name=None,     # No separate exposure node  
    exposure_size=None,  
    exposure_location=None  
)  
```

- **phewas_result** – DataFrame from the PheWAS analysis (must include columns for phecode identifier, disease name/system, case counts, significance etc.).
- **comorbidity_result** – DataFrame from the comorbidity network analysis. It should include columns for disease pairs (D1, D2 identifiers), some unique pair name or ID, the association metrics (like beta or correlation), and a significance indicator (True/False for whether that pair is significant in the network).
- **trajectory_result** – DataFrame from the disease trajectory analysis. It should have columns for the disease pairs (D1, D2, typically oriented as source -> target), the effect sizes (like OR or HR), and a significance indicator for temporal association (True/False for adjusted p-value significance).
- **exposure_name** – A label for the exposure of interest (the factor that defines exposed vs unexposed in the cohort). In our examples, the exposure was “Short LTL” (short leukocyte telomere length) in one of the case studies, hence the example. If you are analyzing something like “smoking” or “diabetes” as the exposure, you’d put that. For an exposed-only study, use None (because there isn’t a separate exposure node to highlight).
- **exposure_location** – The (x, y, z) coordinates where the exposure node should be placed in the 3D plot. By default, if None, the exposure node will be placed at the origin (0,0,0). This is relevant only for 3D plotting; if exposure is None, this is ignored.
- **exposure_size** – A scaling factor for the exposure node’s size in the network visualization. Increase this to make the exposure node larger relative to disease nodes (to emphasize it). If None, in an exposed-only design, the exposure node is not present.

The `Plot` class will internally verify that the required columns exist in the input DataFrames (for example, it expects certain default column names like `'phecode_d1'`, `'phecode_d2'` for pair identifiers, `'phecode'` for disease codes in PheWAS, `'system'` for disease system category, `'name_disease_pair'` for a unique pair name, `'..._significance'` for significance flags, etc. If you did not change column names, it uses these defaults). You can override these defaults by passing optional parameters if your DataFrame uses different column names (we’ll mention those as needed in each plot function below).

Now that `result_plot` is created, we can generate specific plots. The `Plot` class provides multiple methods for different visualizations. Each produces either an interactive HTML file or a static image, as noted.

### 4.2 PheWAS plot

The `result_plot.phewas_plot()` function creates a summary plot of the PheWAS results – essentially a visual depiction of which diseases were associated with the exposure. In a standard or matched cohort study, this is typically a Manhattan-style plot or bar plot showing hazard ratios (HRs) for each significant disease. In an exposed-only cohort, since we don’t have HRs, the plot might show the number of cases of each disease (to highlight which are most frequent). The plot also differentiates diseases by their category/system (often using color coding).

For example:

```python
# Generate a PheWAS plot  
result_plot.phewas_plot(  
    path="/your/project/path/phewas_plot.png",  # output file path (supports .png, .svg, .jpg)  
    is_exposure_only=False                     # False for cohort/matched designs; True if this is an exposed-only cohort  
)  
```
This will save an image to the specified path. You can use `.png` for a static image (good for publications), or `.svg` for a vector graphic, etc. If you want to just display it in a Jupyter notebook, you could omit the path and it might display inline, but typically you provide a path to save.

Parameters for `phewas_plot()` include:

- **path** – File path including filename and extension where the plot will be saved. Ensure the extension is one of the supported image types.
- **is_exposure_only** – Boolean flag; set to `True` if your analysis is an exposed-only cohort (so the plot knows it shouldn’t expect HRs and will plot counts instead). For our example (standard/matched cohort), it’s `False`.

**Optional parameters:** (if your DataFrame columns differ from defaults or if you want to adjust aesthetics)

- **col_coef** – Name of the column in `phewas_result` that contains the effect size (e.g., HR or OR). *(Default: 'phewas_coef')*.
- **col_se** – Name of the column for standard error of the effect size (used to plot error bars). *(Default: 'phewas_se')*.
- **col_system** – Column name for the disease system/category. *(Default: 'system')*.
- **col_disease** – Column name for the disease description/name. *(Default: 'disease')*.
- **col_exposure** – Column name for number of cases in exposed group (used in exposed-only plot). *(Default: 'N_cases_exposed')*.
- **disease_font_size** (parameter is actually `disese_font_size` due to a minor naming typo in code) – Font size for disease labels in the plot (if labels are shown). *(Default: 10)*.
- **system_font_size** – Font size for the disease system labels on the plot. *(Default: 17)*.
- **dpi** – Resolution of the output image (dots per inch). Higher DPI gives a higher resolution image. *(Default: 200)*.

Using these, the PheWAS plot will highlight which diseases came out as significantly associated. Typically, you’ll see something like a scatter of points or bars, colored by system, maybe labeled for the top hits. In our dummy data, since everything is random, the plot isn’t meaningful medically, but with real data this can quickly show you the pattern of associations.

### 4.3 Comorbidity network plot

The `result_plot.comorbidity_network_plot()` function creates an interactive network visualization of the comorbidity relationships. It clusters diseases into communities based on their network connections (using the Louvain community detection algorithm) and plots them, often in a circular layout where each community occupies a sector. The nodes (diseases) are colored by their disease system, and edges represent significant associations from the comorbidity network analysis. This plot is typically output as an HTML file because it’s interactive (you can hover to see details, zoom, etc.).

For example:

```python
# Generate an interactive comorbidity network plot  
result_plot.comorbidity_network_plot(  
    path="/your/project/path/comorbidity_network.html"  
)  
```
This will save an HTML file which you can open in a web browser to explore. Each node (disease) might be labeled or have tooltip info (like disease name, maybe prevalence), and edges might have tooltips for correlation values. Nodes within the same community are grouped together. This visualization helps identify clusters of diseases that frequently co-occur beyond what would be expected.

Optional parameters for `comorbidity_network_plot()` include layout and styling options:

- **max_radius** – Maximum radial distance (in pixels) from the center that nodes can be placed. This essentially controls the size of the outer circle of the plot. *(Default: 90.0)*.
- **min_radius** – Minimum radial distance (pixels) from center for node placement (the inner boundary of the network). *(Default: 35.0)*.
- **layer_distance** – Radial distance between concentric layers of nodes (if nodes are layered by some criteria, e.g., significance or degree). *(Default: 40.0)*.
- **size_reduction** – A scaling factor (0 to 1) for node sizes to ensure they fit nicely (smaller values make nodes proportionally smaller). *(Default: 0.5)*.
- **line_width** – Width (pixels) of the lines (edges) connecting nodes. *(Default: 1.0)*.
- **line_color** – Color of the edges. Can specify a named color (e.g., `'steelblue'`), hex code (`'#4682B4'`), or RGB tuple. *(Default: 'black')*.

These parameters allow fine-tuning the appearance if needed (for example, if node labels overlap, you might reduce node sizes or adjust radii). Usually, the defaults produce a clear separation of communities.

### 4.4 Disease trajectory plot

The `result_plot.trajectory_plot()` function generates visualizations for disease trajectories. Typically, it will produce one plot per community of diseases (as identified in the network analysis) to show how diseases progress within that community. For each significant disease pair (temporal relationship), an arrow or directed edge is drawn from the antecedent disease to the consequent disease. The output is often a set of static image files (e.g., one PNG per community) because it might be easier to print or inspect individually.

For example:

```python
# Generate disease trajectory plots (one per community)  
result_plot.trajectory_plot(  
    path="/your/project/path/trajectory_plots/"  
)  
```

> **Important:** Here, `path` is a directory (you should include the trailing slash). The function will then save multiple files in this directory, named perhaps by community or numbered sequentially. Ensure the directory exists or the function might attempt to create it.

Optional parameters for `trajectory_plot()` include:

- **source** – Column name in the DataFrames for the source disease (D1). *(Default: 'phecode_d1')*.
- **target** – Column name for the target disease (D2). *(Default: 'phecode_d2')*.
- **dpi** – Image resolution (like in PheWAS plot). *(Default: 500 for these plots, to ensure high clarity since many arrows/labels might be present.)*
- **cluster_weight** – This parameter specifies which edge weight from the network to use when arranging the layout. *(Default: 'comorbidity_beta')*, meaning it might use the beta coefficient from the comorbidity network as the weight for clustering layout. Usually, you won’t need to change this.

The trajectory plots will illustrate sequences. Within each community, you might see something like a directed acyclic graph (hopefully acyclic if the data suggests a direction). The arrow from disease A to B indicates A tends to precede B. By generating one plot per community, it keeps the graphs manageable and interpretable. If a community has no significant trajectories, it might be empty or skipped.

After running this, check the output directory for files. For instance, you might find files like `community_1.png`, `community_2.png`, etc., each showing the trajectory network for that community of diseases.

> *(Note: In interactive use, the function might print or log which communities are being plotted and the number of images saved.)*

### 4.5 Three-dimensional plot

The `result_plot.three_dimension_plot()` function is a special visualization that combines the comorbidity network and trajectory information into a single 3D interactive figure. It essentially places the comorbidity network on one plane (say, the horizontal plane) and the trajectory connections on a vertical plane, giving a three-dimensional view where you can see both types of relationships simultaneously. The exposure node (if any) is usually at the center, and disease nodes are arranged around. This plot is typically an interactive HTML as well, since you’d want to rotate and zoom in 3D.

Example usage:

```python
# Generate a combined 3D network plot  
result_plot.three_dimension_plot(  
    path="/your/project/path/combined_network.html"  
)  
```

This saves an HTML file with the interactive 3D visualization. When you open it, you can drag to rotate the network in 3D space. One plane (e.g., viewed from top-down) will show the clusters and connections akin to the comorbidity network, and the other plane (viewed from the side) will show the directed trajectories. Where they intersect, you get a full picture of how diseases group and progress.

Optional parameters for `three_dimension_plot()` include various layout settings:

- **max_radius** – Maximum radius for node placement from the center (similar to the 2D network plot). *(Default: 180.0)* (since in 3D we might allow a larger radius).
- **min_radius** – Minimum radius from center. *(Default: 35.0)*.
- **layer_distance** – Distance between layers in the radial direction (similar concept to above). *(Default: 40.0)*.
- **layout_width** – Width of the overall figure in pixels. *(Default: 900.0)*.
- **layout_height** – Height of the figure in pixels. *(Default: 900.0)*.
- **line_color** – Color for trajectory lines (since in 3D plot, maybe comorbidity edges are one style and trajectory edges another). *(Default: 'black')*.
- **line_width** – Width of trajectory lines. *(Default: 1.0)*.
- **size_reduction** – Node size scaling factor (0.1–1.0). *(Default: 0.5)*.
- **cluster_reduction_ratio** – A factor (0.1–1.0) to compress or spread out clusters in the 3D space. Lower means clusters are more tightly grouped. *(Default: 0.4)*.
- **cluster_weight** – Which weight to use for determining cluster layout (as in trajectory_plot, usually 'comorbidity_beta'). *(Default: 'comorbidity_beta')*.
- **font_style** – Font family for text elements (node labels, etc.). *(Default: 'Times New Roman')*.
- **font_size** – Base font size for text. *(Default: 15.0)*.

The 3D plot is a bit advanced and may require a good computer/browser to manipulate if there are many nodes, but it provides a unique integrated view.

With all these visualization tools, even a beginner user can not only run the analysis with ***DiNetxify*** but also see the results in intuitive forms, which can greatly aid interpretation and presentation. Next, we provide a quick reference to the API of the main classes and functions for convenience.

## API Reference

Below is a concise reference for ***DiNetxify***’s classes and functions, summarizing their signatures and parameters. This is useful when writing your own scripts or if you need to quickly recall how to call a function.

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

A class for handling disease network data creation and operations, for use in ***DiNetxify*** module.

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
- `binomial_test_result` (*pd.DataFrame, default=None*): DataFrame containing binomial test analysis results produced by the `DiNetxify.binomial_test` function.
- `method` (*str, default='RPCN'*): Specifies the comorbidity network analysis method to use. Choices are: - 'RPCN: Regularized Partial Correlation Network. - 'PCN_PCA: Partial Correlation Network with PCA. - 'CN': Correlation Network. **Additional Options for RPCN:** - 'alpha' : non-negative scalar The weight multiplying the l1 penalty term for other diseases covariates. Ignored if 'auto_penalty' is enabled. - 'auto_penalty' : bool, default=True If 'True', automatically determine the optimal 'alpha' based on model AIC value. - 'alpha_range' : tuple, default=(1,15) When 'auto_penalty' is True, search the optimal 'alpha' in this range. - 'scaling_factor' : positive scalar, default=1 The scaling factor for the alpha when 'auto_penalty' is True. **Additional Options for PCN_PCA:** - 'n_PC' : int, default=5 Fixed number of principal components to include in each model. - 'explained_variance' : float Determines the number of principal components based on the cumulative explained variance. Overrides 'n_PC' if specified.
- `covariates` (*list, default=None*): List of phenotypic covariates to include in the model. By default, includes ['sex'] and all covariates specified in the `DiNetxify.DiseaseNetworkData.phenotype_data()` function. To include the required variable sex as a covariate, always use 'sex' instead of its original column name. For other covariates specified in the `DiNetxify.DiseaseNetworkData.phenotype_data()` function, use their original column names.
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
- `comorbidity_strength_result` (*pd.DataFrame*): DataFrame containing comorbidity strength analysis results produced by the `DiNetxify.comorbidity_strength()` function.
- `binomial_test_result` (*pd.DataFrame*): DataFrame containing binomial test analysis results produced by the `DiNetxify.binomial_test()` function.
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