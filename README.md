# DiseaseNetPy

DiseaseNetPy is a Python package designed for comprehensive disease network analysis and visualization. This novel disease network analysis approach (three-dimensional disease network analysis) integrates and refines existing disease trajectory and comorbidity network analysis methods. This new approach enhances disease association verification by incorporating regularized partial correlations. It also facilitates robust identification and visualization of disease clusters (i.e., groups of depression-associated diseases with high within-group connectivity) through both non-temporal (illustrated by the x-axis and y-axis) and temporal (z-axis) dimensions.

## Table of Contents

- [Installation](#installation)
  - [1. Input data preparation](#1-input-data-preparation)
    - [1.1 Requirements for input data](#11-requirements-for-input-data)
    - [1.2 Description of the dummy dataset provided](#12-description-of-the-dummy-dataset-provided)
  - [2. Data Harmonization](#2-data-harmonization)
    - [2.1 Initializing the data object](#21-initializing-the-data-object)
    - [2.2 Load phenotype data](#22-load-phenotype-data)
    - [2.3 Load medical records data](#23-load-medical-records-data)
  - [3. Data Analysis](#3-data-analysis)
    - [3.1 Quick Pipeline](#31-quick-pipeline)
    - [3.2 PheWAS Analysis](#32-phewas-analysis)
    - [3.3 Disease Pair Construction](#33-disease-pair-construction)
    - [3.4 Comorbidity Strength Estimation](#34-comorbidity-strength-estimation)
    - [3.5 Binomial Test](#35-binomial-test)
    - [3.6 Comorbidity Network Analysis](#36-comorbidity-network-analysis)
    - [3.7 Trajectory Analysis](#37-trajectory-analysis)
  - [4. Visualization](#4-visualization)
    - [4.1 Initializing the plot object](#41-initializing-the-plot-object)
    - [4.2 PheWAS Plot](#42-phewas-plot)
    - [4.3 Comorbidity Network Plot](#43-comorbidity-network-plot)
    - [4.4 Disease Trajectory Plot](#44-disease-trajectory-plot)
    - [4.5 Three Dimension Plot](#45-three-dimension-plot)
- API Reference
  - Classes
    - [DiseaseNetworkData](#diseasenetworkdata)
    - [Plot](#plot)
  - Functions
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
    - [phewas_plot](#phewas_plot)
    - [comorbidity_network_plot](#comorbidity_network_plot)
    - [trajectory_plot](#trajectory_plot)
    - [three_dimension_plot](#three_dimension_plot)
- [Issues reporting and recommendations](#trajectory_multipletests)
- [License](#license)

## Installation

You can install DiseaseNetPy via `pip`. Ensure you have Python 3.13 .

```bash
pip install diseasenetpy

#optional packages
pip install lifelines #for phewas analysis with lifelines_disable set to False
```
## 1. Input data preparation

### 1.1 Requirements for input data

DiseaseNetPy performs 3D disease network analysis on cohort data derived from electronic health records (EHR). We currently support three study designs: standard cohort, matched cohort, and exposed-only cohort. Standard and matched cohort studies are suitable for investigating the disease network of individuals with a specific disease (e.g., depression) or exposure (e.g., smoking). The exposed-only cohort design is suitable for investigating the disease network in the whole population or a subset (e.g., older individuals) without a comparison group.

To begin using DiseaseNetPy, two datasets are required: a **phenotype data** file recording each individual's basic and additional characteristics, and a **medical records data** file extracted from an EHR database that stores diagnosis codes and dates for all cohort individuals over the study period. Specific requirements for these datasets are as follows:

- **Phenotype data**: An on-disk CSV (or TSV) file with headers, listing each participant with the following required and optional columns:
  
  - **Participant ID** – unique identifier for each individual.
  - **Index date** – start date of follow-up (e.g., date of exposure or baseline).
  - **End date** – end of follow-up (e.g., last visit, death, or study end).
  - **Exposure** – exposure status (1 for exposed, 0 for unexposed) for cohort and matched cohort designs (omit this column for the exposed-only design).
  - **Match ID** – group identifier for matching sets (only for matched cohort design).
  - **Sex** – sex of the individual (coded as 1 for female, 0 for male).
  - **Additional covariates (optional)** – any number of covariates useful for adjustment or matching (e.g., age, BMI, education).
  
  For required columns, missing data are not allowed, and dates must be specified using format codes according to the 1989 C standard. The **Sex** and **Exposure** columns must be coded as specified. There is no limit on the number of additional covariates. Covariate types (binary, categorical, or continuous) are automatically determined based on value type and distribution and transformed accordingly (e.g., via one-hot encoding). Missing values are allowed: for categorical variables, they form a separate NA category; for continuous variables, individuals with missing values are dropped.
  
- **Medical records data**: One or more CSV (or TSV) files with headers, containing diagnosis records for participants. Each file should include the following columns:
  
  - **Participant ID** – must match the IDs in the phenotype file.
  - **Diagnosis code** – a diagnosis identifier (ICD-10, ICD-9, or another coding system).
  - **Date of diagnosis** – the date when the diagnosis was recorded.
  
  Unlike the phenotype file, the medical records files should be in a record-per-row format, allowing multiple rows per participant, consistent with typical EHR structures. Records matching participant IDs in the phenotype file and occurring within the specified follow-up period (i.e., before the end date) are loaded; all others are discarded. Therefore, there is **no need** to pre-filter medical records - you should provide complete diagnosis information for all cohort individuals. Specifically, it is **not recommended** to filter data based on first occurrence of each code or diagnosis date.
  
  Each medical records file must use a single diagnosis code version; currently supported versions are WHO or CM versions of ICD-9 and ICD-10. Other code systems must be converted to a supported format. ICD-10 codes may be formatted with a decimal point (e.g., `F32.1`) or without one (e.g., `F321`); ICD-9 codes may use a decimal format (e.g., `9.9`) or a non-decimal “short” format (e.g., `0099`). The **Date of diagnosis** column must use a consistent date format specified by the 1989 C standard.

### 1.2 Description of the dummy dataset provided

The package includes example datasets demonstrating the required data format for disease network analysis:

- **dummy_phenotype.csv**: 60,000 records (3.72 MB)

  - **ID** – unique identifier for each individual (e.g., 1001).
  - **date_start** – start date of follow-up (e.g., 2016/10/13).
  - **date_end** – end of follow-up (e.g., 2022/8/2).
  - **exposure** – exposure status (e.g., 1 for exposed, 0 for unexposed).
  - **group_id** – group identifier for matching sets (e.g., "group_0").
  - **sex** – sex of the individual (e.g., 0 or 1).
  - **age** – age for year (e.g., 65.2).
  - **BMI** – the BMI level for each individual (e.g., "c1").

The dummy phenotype dataset (**dummy_phenotype.csv**) contains 60,000 synthetic records (3.72 MB) simulating real-world cohort data for disease network analysis. Each record includes essential participant information: a unique identifier (**ID**), follow-up period dates (**date_start** and **date_end**), exposure status (**exposure** coded as 0/1), and matching group identifiers (**group_id**) for matched cohort studies. The dataset also captures demographic characteristics including biological sex (**sex** coded as 0/1), precise age in years (**age** as a discrete variable), and BMI categories (**BMI** as a categorical variable). The covariates of dummy phenotype dataset are **sex**, **age**, and **BMI**.

- **dummy_EHR_ICD9.csv**: 10,188 records (227 kB)

  - **ID** – unique identifier for each individual (e.g., 1001).
  - **dia_date** – date of diagnosis (e.g., 2016/10/13).
  - **diag_icd9** – ICD9 of diagnosis  (e.g., "E950").

The **dummy_EHR_ICD9.csv** file contains 10,188 simulated electronic health records (227 kB) with ICD-9 coded diagnoses, structured with three key fields: patient identifier (**ID**), diagnosis date (**dia_date**), and corresponding ICD-9 code (**diag_icd9**). This standardized format demonstrates the required medical records structure for disease network analysis, where each row represents a discrete diagnosis event (e.g., patient 1001's `E950` coded diagnosis on 2016/10/13).

- **dummy_EHR_ICD10.csv**: 1,048,576 records (36.2 MB)

  - **ID** – unique identifier for each individual (e.g., 1001).
  - **dia_date** – date of diagnosis (e.g., 2014/9/23).
  - **diag_icd10** – ICD10 of diagnosis  (e.g., "L905").

The **dummy_EHR_ICD10.csv** file contains 1,048,576 simulated electronic health records (36.2 MB) with ICD-10 coded diagnoses, structured with three key fields: patient identifier (**ID**), diagnosis date (**dia_date**), and corresponding ICD-10 code (**diag_icd10**). This standardized format demonstrates the required medical records structure for disease network analysis, where each row represents a discrete diagnosis event (e.g., patient 1001's `L905` coded diagnosis on 2014/9/23).

## 2. Data Harmonization

Data harmonization loads and merges phenotype and medical records data into a single `DiseaseNetworkData` object for subsequent analysis, ensuring consistent coding (e.g., mapping diagnosis codes to phecodes) and standardized formatting.

### 2.1 Initializing the data object

First, import DiseaseNetPy and initialize an empty `DiseaseNetworkData` object, specifying the study design, phecode level, and any optional parameters if needed:

```python
import diseasenetpy as dnt

# For a standard cohort study
data = dnt.DiseaseNetworkData(
    study_design='cohort',           # 'cohort', 'matched cohort', or 'exposed-only cohort'
    phecode_level=1,                 # use phecode level 1 (3-digit ICD grouping)
)

# For a matched cohort study
data = dnt.DiseaseNetworkData(
    study_design='matched cohort',   # 'cohort', 'matched cohort', or 'exposed-only cohort'
    phecode_level=1,
)

# For an exposed-only cohort study
data = dnt.DiseaseNetworkData(
    study_design='exposed-only cohort',  # 'cohort', 'matched cohort', or 'exposed-only cohort'
    phecode_level=1,
)
```

- **study_design** – the cohort design: `'cohort'`, `'matched cohort'`, or `'exposed-only cohort'`.
- **phecode_level** – the level of phecode used for analysis; level 1 provides broader categories (~585 conditions), while level 2 offers more detail (~1257 conditions). Level 1 is recommended for smaller datasets to maintain statistical power.

#### Optional parameters:

- **min_required_icd_codes** – the minimum number of ICD codes mapping to a phecode required for it to be considered valid; default is 1. For example, setting this to 2 requires at least two records mapping to phecode 250.2 (Type 2 diabetes) for a participant to be considered diagnosed. Ensure your medical records include complete data (not limited to first occurrences) when using this parameter.
- **date_fmt** – format of date fields (Index date and End date) in the phenotype data. Default is `'%Y-%m-%d'` (year-month-day, e.g., 2005-12-01).
- **phecode_version** – currently only `'1.2'` is supported.

### 2.2 Load phenotype data

After initializing the data object, use the `phenotype_data()` method to load your cohort phenotype file by providing the file path, a dictionary mapping required columns, and a list of additional covariate names.

The following example codes show how to load the dummy phenotype dataset under different study designs. Although the file is originally formatted for a matched cohort study, you can adapt it for other designs: omitting the **Match ID** column loads it as a standard cohort study (ignoring the matching), while omitting both **Match ID** and **Exposure** columns loads it as an exposed‑only cohort study - treating all participants as exposed (i.e., representing the entire population).

```python
#--------- Load phenotype data for a matched cohort study
col_dict = {
    'Participant ID': 'ID',        # maps to ID column in file
    'Exposure': 'exposure',        # 0/1 exposure status
    'Sex': 'sex',                  # 0 for male, 1 for female
    'Index date': 'date_start',    # start of follow-up
    'End date': 'date_end',        # end of follow-up
    'Match ID': 'group_id'         # matching group identifier
}
vars_lst = ['age', 'BMI']          # covariates to be used in the analysis
data.phenotype_data(
    phenotype_data_path=r"/test/data/dummy_phenotype.csv",
    column_names=col_dict,
    covariates=vars_lst
)

#--------- Load phenotype data for a traditional cohort study
col_dict = {
    'Participant ID': 'ID',
    'Exposure': 'exposure',
    'Sex': 'sex',
    'Index date': 'date_start',
    'End date': 'date_end'
}
vars_lst = ['age', 'BMI']
data.phenotype_data(
    phenotype_data_path=r"/test/data/dummy_phenotype.csv",
    column_names=col_dict,
    covariates=vars_lst
)

#--------- Load phenotype data for an exposed-only cohort study
col_dict = {
    'Participant ID': 'ID',
    'Sex': 'sex',
    'Index date': 'date_start',
    'End date': 'date_end'
}
vars_lst = ['age', 'BMI']
data.phenotype_data(
    phenotype_data_path=r"/test/data/dummy_phenotype.csv",
    column_names=col_dict,
    covariates=vars_lst
)
```

- **phenotype_data_path** – path to your phenotype data file (CSV or TSV).
- **column_names** – dictionary mapping the required variable names (e.g., `Participant ID`, `Index date`, `End date`, `Sex`) to the corresponding headers in your file. Include `Exposure` for cohort and matched cohort designs, and `Match ID` for matched cohort designs.
- **covariates** – list of additional covariate names. Provide an empty list if none. The method automatically detects and converts variable types. Records with missing values in continuous variables are removed, while missing values in categorical variables form an NA category.

#### Optional parameters:

- **is_single_sex** – set to `True` if the dataset contains only one sex. Default is `False`.
- **force** – set to `True` to overwrite existing data in the object. Default is `False`, which raises an error if data already exist.

#### After data loading:

After loading data, you can inspect basic information (e.g., number of individuals, average follow-up time) by printing the object:

```python
print(data)
# This will output something like (e.g., for a matched cohort study):
"""
DiseaseNetPy.DiseaseNetworkData

Study design: matched cohort

Phentype data
Total number of individuals: 60,000 (10,000 exposed and 50,000 unexposed)
The average group size is: 6.00
Average follow-up years: 10.44 (exposed) and 10.46 (unexposed)

Warning: 102 exposed individuals and 440 unexposed individuals have negative or zero follow-up time.
Consider removing them before merge.
"""
```

Additionally, you can generate a basic descriptive Table 1 for all variables in your phenotype data using the `Table1()` method and save it to an Excel file:

```python
table_1 = data.Table1()
print(table_1)
# Example output (e.g., for a matched cohort study):
"""
                   Variable exposure=1 (n=10,000) exposure=0 (n=50,000)                       Test and p-value
0        _age (median, IQR)   57.08 (48.91-65.32)   57.05 (48.87-65.35)  Mann-Whitney U test p-value=9.824e-01      
1   follow_up (median, IQR)     9.18 (5.77-13.70)     9.22 (5.80-13.75)  Mann-Whitney U test p-value=6.806e-01      
2                sex (n, %)
3                sex=Female        5,045 (50.45%)       25,225 (50.45%)
4                  sex=Male        4,955 (49.55%)       24,775 (49.55%)     Chi-squared test p-value=1.000e+00      
5                BMI (n, %)
6                    BMI=c2        1,945 (19.45%)       10,170 (20.34%)
7                    BMI=c4        2,022 (20.22%)       10,022 (20.04%)
8                    BMI=c5        2,002 (20.02%)       10,031 (20.06%)
9                    BMI=c1        1,999 (19.99%)        9,952 (19.90%)
10                   BMI=c3        2,032 (20.32%)        9,825 (19.65%)     Chi-squared test p-value=2.552e-01 
"""

table_1.to_excel(r"/test/data/Table1.xlsx")  # Save Table 1 to an Excel file
```

### 2.3 Load medical records data

After loading the phenotypic data, use the `merge_medical_records()` method to load your one or more medical records files by providing the file path, a format of ICD code, and mapping required columns, and a dictionary mapping required columns.

The following example codes show how to load the dummy EHR ICD9/ICD10 dataset.

```python
# Merge with the first medical records file (dummy_EHR_ICD10.csv)
data.merge_medical_records(
    medical_records_data_path=r"/test/data/dummy_EHR_ICD10.csv",  # Path to first medical records file
    diagnosis_code='ICD-10-WHO',                                  # Diagnosis code type
    column_names={
        'Participant ID': 'ID',                                   # Participant ID column in medical records
        'Diagnosis code': 'diag_icd10',                           # Diagnosis code column
        'Date of diagnosis': 'dia_date'                           # Diagnosis date column
    }
)

# Merge with the second medical records file (dummy_EHR_ICD9.csv)
data.merge_medical_records(
    medical_records_data_path=r"/test/data/dummy_EHR_ICD9.csv",  
    diagnosis_code="ICD-9-WHO",                                  
    column_names={
        'Participant ID':'ID',                                   
        'Diagnosis code':'diag_icd9',                            
        'Date of diagnosis':'dia_date'                           
    }
)
```

- **medical_records_data_path** – path to your medical records data files (CSV or TSV).
- **diagnosis_code** – Diagnosis ICD code type used in the medical records data (e.g., `ICD-9-CM`, `ICD-9-WHO`, `ICD-10-CM`, `ICD-10-WHO`).
- **column_names** – dictionary mapping the required variable names (e.g., `Participant ID`, `Diagnosis code`, `Date of diagnosis`) to the corresponding headers in your file.

#### Optional parameters:

- **date_fmt** – The format of the date fields in your medical records data. Defaults to the same format as phenotype data if not specified.
- **chunksize** – Number of rows per chunk to read, useful for large datasets. Default=1_000_000.

#### After data loading:

Within loading data, you can inspect basic information (e.g., number of records, number of phecode mapping) by the printing information:

```python
"""
1,000,000 records read (1,000,000 included after filltering on participant ID), 0 records with missing values excluded.
1,668,795 records read (1,668,795 included after filltering on participant ID), 0 records with missing values excluded.
Total: 1,668,795 diagnosis records processed, 0 records with missing values were excluded.
1,286,386 diagnosis records mapped to phecode without truncating.
0 diagnosis records mapped to phecode after truncating to 4 digits.
72,073 diagnosis records mapped to phecode after truncating to 3 digits.
302,908 diagnosis records not mapped to any phecode.
Phecode diagnosis records successfully merged (18,486 invalid records were not merged, typically with diagnosis date later than date of follow-up end)

1 medical records data already merged, merging with a new one.
10,188 records read (10,188 included after filltering on participant ID), 0 records with missing values excluded.
Total: 10,188 diagnosis records processed, 0 records with missing values were excluded.
9,711 diagnosis records mapped to phecode without truncating.
0 diagnosis records mapped to phecode after truncating to 4 digits.
266 diagnosis records mapped to phecode after truncating to 3 digits.
211 diagnosis records not mapped to any phecode.
Phecode diagnosis records successfully merged (0 invalid records were not merged, typically with diagnosis date later than date of follow-up end)
"""
```

## 3. Data Analysis

Data analysis is based on a `DiseaseNetworkData` object to subsequent analysis, including PheWAS analysis, disease pair generation, comorbidity strength estimation, binomial testing, comorbidity network analysis, and trajectory analysis. The result format of each analysis is `pd.DataFrame`.

### 3.1 Quick Pipeline

This guide walks you through a typical workflow using DiseaseNetPy for a matched cohort study design. The process involves data preparation, PheWAS analysis, comorbidity strength estimation, binomial testing, comorbidity network analysis, and trajectory analysis.

```python
# Reminder:
# When using multiprocessing, ensure that the code is enclosed within the following block.
# This prevents entering a never ending loop of new process creation.

if __name__ == "__main__":
    dnt.disease_network_pipeline(
      data=data,                               # DiseaseNetworkData object
      n_process=2,                             # Number of parallel processes
      n_threshold_phewas=100,                  # Minimum number of cases to include in phewas analysis
      n_threshold_comorbidity=100,             # Minimum number of cases to include in comorbidity strength analysis
      output_dir="/your/project/path"          # Path to save results
      project_prefix="disease network"         # Prefix for naming output files and intermediate data
      keep_positive_associations=False         # control hazard ratio (HR) > 1 in the Phewas analysis and positive correlations in the comorbidity strength estimation
      save_intermediate_data=False             # control to save intermediate DiseaseNetworkData objects created by the `DiseaseNetPy.DiseaseNetworkData.disease_pair` function
      system_exl=[
        'symptoms', 
        'others', 
        'injuries & poisonings'
      ],                                       # List of phecode systems to exclude from the analysis
      pipeline_mode="v1",                      # Specifies the analysis order
      method="RPCN",                           # The method to use for the comorbidity network and disease trajectory analysis
      covariates=[
        'BMI', 
        'age'
      ],                                      # List of covariates to adjust for in the PheWAS, comorbidity network and disease trajectory analysis
      matching_var_dict={
        'sex':'exact'
      },                                      # Specifies the matching variables and the criteria used for incidence density sampling
      matching_n=2,                           # Specifies the maximum number of matched controls for each case
      min_interval_days=0,                    # Minimum required time interval (in days) between diagnosis dates when constructing temporal D1 → D2 disease pair for each individual
      max_interval_days=np.inf,               # Maximum allowed time interval (in days) between diagnosis dates when constructing temporal and non-temporal D1-D2 disease pair for each individual
      enforce_temporal_order=False,           # control to exclude individuals with non-temporal D1-D2 pair when performing the binomial test
      correction='bonferroni',                # Method for p-value correction from the statsmodels.stats.multitest.multipletests
      cutoff=0.05                             # The significance threshold for adjusted p-values
    )

```

- **data** – The DiseaseNetworkData object.
- **n_process** – Specifies the number of parallel processes to use for the disease network analysis. Multiprocessing is enabled when `n_process` is set to a value greater than one.
- **n_threshold_phewas** – The minimum number of cases within the exposed group required for a phecode to be included in the PheWAS analysis. This parameter will be passed to the phewas function. See the phewas function for more information.
- **n_threshold_comorbidity** – The minimum number of individuals in the exposed group in which a disease pair must co-occur (temporal or non-temporal) to be included in the comorbidity strength estimation. This parameter will be passed to the comorbidity_strength function. See the comorbidity_strength function for more information.
- **project_prefix** – Prefix for naming output files and intermediate data.

#### Optional parameters:

- **keep_positive_associations** – set to `True` if retains only diseases with hazard ratio (HR) > 1 from the PheWAS analysis. Default is `False`.
- **save_intermediate_data** – set to `True` to intermediate DiseaseNetworkData objects created by the `DiseaseNetPy.DiseaseNetworkData.disease_pair` function are saved to disk. Default is `False`.
- **system_exl** – List of phecode systems to exclude from the analysis. Default is `None`.
- **pipeline_mode** – Specifies the analysis order. Available modes: - 'v1': PheWAS → comorbidity strength → binomial test → (comorbidity network analysis/disease trajectory analysis) - 'v2': PheWAS → comorbidity strength → comorbidity network analysis → binomial test → disease trajectory analysis. In 'v1', the binomial test does not depend on results from the comorbidity network analysis; thus, disease trajectory and comorbidity network analyses can be conducted independently. In 'v2', the binomial test is performed only on disease pairs identified as significant by the comorbidity network analysis, making the disease trajectory analysis dependent on these results.
- **method** – The method to use for the comorbidity network and disease trajectory analysis. Available methods are: - 'RPCN: Regularized Partial Correlation Network. - 'PCN_PCA: Partial Correlation Network with PCA. - 'CN': Correlation Network. This parameter will be passed to the comorbidity_network and disease_trajectory function. See the these two functions for more information.
- **covariates** – List of covariates to adjust for in the PheWAS, comorbidity network and disease trajectory analysis. Default is `None`.
- **matching_var_dict** – Specifies the matching variables and the criteria used for incidence density sampling. Default is `{'sex':'exact'}`.
- **matching_n** – Specifies the maximum number of matched controls for each case. This parameter will be passed to the disease_trajectory function. Default is `2`.
- **min_interval_days** – Minimum required time interval (in days) between diagnosis dates when constructing temporal D1 → D2 disease pair for each individual. This parameter will be passed to the DiseaseNetPy.DiseaseNetworkData.disease_pair function. See the disease_pair function for more information. Default is `0`.
- **max_interval_days** – Maximum allowed time interval (in days) between diagnosis dates when constructing temporal and non-temporal D1-D2 disease pair for each individual. This parameter will be passed to the DiseaseNetPy.DiseaseNetworkData.disease_pair function. See the disease_pair function for more information. Default is `np.inf`.
- **enforce_temporal_order** – set to `True` to exclude individuals with non-temporal D1-D2 pair when performing the binomial test; also applies the specified minimum and maximum time intervals when performing disease trajectory analysis. See 'enforce_temporal_order' parameter in binomial_test function and 'enforce_time_interval' parameter in disease_trajectory function. Default is `False`.
- **correction** – Method for p-value correction from the statsmodels.stats.multitest.multipletests.          
    Available methods are: 
        none : no correction
        bonferroni : one-step correction
        sidak : one-step correction
        holm-sidak : step down method using Sidak adjustments
        holm : step-down method using Bonferroni adjustments
        simes-hochberg : step-up method (independent)
        hommel : closed method based on Simes tests (non-negative)
        fdr_bh : Benjamini/Hochberg (non-negative)
        fdr_by : Benjamini/Yekutieli (negative)
        fdr_tsbh : two stage fdr correction (non-negative)
        fdr_tsbky : two stage fdr correction (non-negative)
        See https://www.statsmodels.org/dev/generated/statsmodels.stats.multitest.multipletests.html for more details.
    Default is `bonferroni`.
- **cutoff** – The significance threshold for adjusted p-values. Default is `0.05`.

#### Kwargs Parameters:

##### RPCN Method Parameters:

- **alpha** – The weight multiplying the l1 penalty term for other diseases covariates. Ignored if 'auto_penalty' is enabled.
- **auto_penalty** – If 'True', automatically determines the best 'alpha' based on model AIC value. Default is `True`.
- **alpha_range** – When 'auto_penalty' is True, search the optimal 'alpha' in this range. Default is `(1,15)`.
- **scaling_factor** – The scaling factor for the alpha when 'auto_penalty' is True. Default is `1`.

##### PCN_PCA Method Parameters:

- **n_PC** – Fixed number of principal components to include in each model. Default is `5`.
- **explained_variance** – Cumulative explained variance threshold to determine the number of principal components. Overrides 'n_PC' if specified.

### 3.2 PheWAS Analysis
The first step of data analysis is PheWAS analysis, aiming to filter significantly occurring diseases. Conducts Phenome-wide association studies (PheWAS) using the specified DiseaseNetworkData object.

```python
# Reminder:
# When using multiprocessing, ensure that the code is enclosed within the following block.
# This prevents entering a never ending loop of new process creation.
if __name__ == "__main__":
    phewas_result = dnt.phewas(
        data=data,                                            # DiseaseNetworkData object
        proportion_threshold=0.01,                            # Minimum proportion of cases to include
        n_process=2,                                          # Number of parallel processes
        system_exl=[                                          # Phecode systems to exclude
            'symptoms', 
            'others', 
            'injuries & poisonings', 
            'pregnancy complications'
        ],                                                    # exclude phecode diseases system
        covariates=[
          'age', 
          'social', 
          'BMI', 
          'smoking', 
          'drinking'
        ],                                                    # Covariates to adjust for
        correction='bonferroni'                               # Method to Multiple test
        lifelines_disable=True,                               # Disable lifelines for faster computation
        log_file='/your/project/path/dep.log'                 # Path to log file
    )

    # Save the entire PheWAS results to a CSV file
    phewas_result.to_csv('/your/project/path/dep_phewas.csv')  # Path to save PheWAS results

    # Filter results for diseases with Hazard Ratio (HR) > 1 if necessary
    phewas_result = phewas_result[phewas_result['phewas_coef'] > 0]

    # Redo p-value adjustment using False Discovery Rate (FDR) Benjamini-Hochberg method
    phewas_result = dnt.phewas_multipletests(
        df=phewas_result,                                     # DataFrame with PheWAS results
        correction='fdr_bh',                                  # P-value correction method
        cutoff=0.05                                           # Significance threshold
    )
```

- **data** – The DiseaseNetworkData object.

#### Optional Parameters:

- **covariates** – List of phenotypic covariates to include in the model. Default is `None`.
- **proportion_threshold** – The minimum proportion of cases within the exposed group required for a phecode to be included in the PheWAS analysis. If the proportion of cases is below this threshold, the phecode is excluded from the analysis. `proportion_threshold` and `n_threshold` are mutually exclusive. Default is `None`.
- **n_threshold** - The minimum number of cases within the exposed group required for a phecode to be included in the PheWAS analysis. If the number of cases is below this threshold, the phecode is excluded.  
*Note:* `n_threshold` and `proportion_threshold` are mutually exclusive. Default is `None`.
- **n_process** - Number of parallel processes to use for analysis. Multiprocessing is enabled when set to greater than one. Default is `1`.
- **correction** - Method for p-value correction (from `statsmodels.stats.multitest.multipletests`).  
Available methods:  
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
*Reference:* [statsmodels documentation](https://www.statsmodels.org/dev/generated/statsmodels.stats.multitest.multipletests.html)  
Default is `bonferroni`.
- **cutoff** - Significance threshold for adjusted PheWAS p-values. Default is `0.05`.
- **system_inc** - Phecode systems to include in analysis. *Note:* Mutually exclusive with `system_exl`.  
Available systems:  
circulatory, congenital anomalies, dermatologic, digestive, endocrine/metabolic, genitourinary, hematopoietic, infectious diseases, injuries & poisonings, mental disorders, musculoskeletal, neoplasms, neurological, pregnancy complications, respiratory, sense organs, symptoms, others  
Default is `None`.
- **system_exl** - Phecode systems to exclude from analysis. *Note:* Mutually exclusive with `system_inc`.  
Available systems: Same as `system_inc` Default is `None`.
- **phecode_inc** - Specific phecodes to include in analysis. *Note:* Mutually exclusive with phecode_exl. Default is `None`.
- **phecode_exl** - Specific phecodes to exclude from analysis. *Note:* Mutually exclusive with phecode_inc. Default is `None`.
- **log_file** - Path and prefix for log file. If None, logs are written to temporary directory with prefix DiseaseNet_. Default is `None`.
- **lifelines_disable** - Whether to disable lifelines. Lifelines provide more robust fitting but require longer computation time. Default is `False`.

#### After PheWAS Analysis:

After PheWAS analysis, we get a `phewas_result` of `pd.DataFrame` format. And then use `phewas_result` to disease pair construction

### 3.3 Disease Pair Construction

The second step of data analysis is Disease pair construction, aiming to filter significantly occurring diseases. Conducts Phenome-wide association studies (PheWAS) using the specified DiseaseNetworkData object.

```python
data.disease_pair(
    phewas_result=phewas_result,               # Filtered PheWAS results
    min_interval_days=30,                      # Minimum interval between diagnoses (30 days here)
    max_interval_days=365.25*5,                # Maximum interval between diagnoses (5 years here)
    force=True                                 # Overwrite existing data if necessary (be cautious)
)

# Save the updated data object with disease pairs
data.save('/your/project/path/dep_withtra')    # Path to save the updated data object
```

- **phewas_result** - `pd.DataFrame` containing PheWAS analysis results produced by the `DiseaseNetPy.phewas` function.

#### Optional Parameters:

- **min_interval_days** - Minimum required time interval (in days) between diagnosis dates when constructing temporal D1→D2 disease pairs. Individuals with D1 and D2 diagnoses interval ≤ this value are considered to have non-temporal pairs. Default is `0`.
- **max_interval_days** - Maximum allowed time interval (in days) between diagnosis dates when constructing disease pairs. Individuals with interval > this value are excluded from temporal analysis. Default is `np.inf`.
- **force** - If `True`, overwrites existing data attributes. If `False`, raises error when data exists. Default is `False`.
- **n_process** - Number of processes for parallel processing. Values >1 enable multiprocessing. Default is `1`.

#### Kwargs Parameters:

- **phecode_col** - Column name for phecode identifiers in `phewas_result`. Default is `'phecode'`.  
- **significance_col** - Column name for PheWAS significance values. Default is `'phewas_p_significance'`.

#### After Disease Pair Construction:

After disease pair construction, we get a new `DiseaseNetworkData` object. And then use `DiseaseNetworkData` to comorbidity strength estimation.

### 3.4 Comorbidity Strength Estimation

The third step of data analysis is comorbidity strength estimation, aiming to estimate the comorbidity strength of disease pairs.

```python
# Reminder:
# When using multiprocessing, ensure that the code is enclosed within the following block.
# This prevents entering a never ending loop of new process creation.

if __name__ == "__main__":
    com_strength_result = dnt.comorbidity_strength(
        data=data,                                     # DiseaseNetworkData object
        proportion_threshold=0.001,                    # Minimum proportion for comorbidity
        n_process=2,                                   # Number of parallel processes
        log_file='/your/project/path/dep.log'          # Path to log file
    )

    # Save the comorbidity strength results to a CSV file
    com_strength_result.to_csv('/your/project/path/dep_com_strength.csv')  # Path to save comorbidity strength results

    # Further filter based on Relative Risk (RR) > 1 and phi-correlation > 0 if necessary
    com_strength_result = com_strength_result[
        (com_strength_result['phi'] > 0) & (com_strength_result['RR'] > 1)
    ]

    # Adjust p-values for comorbidity strength using FDR Benjamini-Hochberg method
    com_strength_result = dnt.comorbidity_strength_multipletests(
        df=com_strength_result,                         # DataFrame with comorbidity strength results
        correction_phi='fdr_bh',                        # P-value correction for phi-correlation
        correction_RR='fdr_bh',                         # P-value correction for Relative Risk
        cutoff_phi=0.05,                                # Significance threshold for phi-correlation
        cutoff_RR=0.05                                  # Significance threshold for Relative Risk
    )
```

- **data** - DiseaseNetworkData object containing the processed disease network data.

#### Optional Parameters:

- **proportion_threshold** - Minimum proportion of exposed individuals required for disease pair co-occurrence to be included in analysis. Disease pairs below this threshold are excluded. *Note:* Mutually exclusive with `n_threshold`. Default is `None`.
- **n_threshold** - Minimum number of exposed individuals required for disease pair co-occurrence to be included in analysis. Disease pairs below this threshold are excluded. *Note:* Mutually exclusive with proportion_threshold. Default is `None`.
- **n_process** - Number of parallel processes for analysis. Values >1 enable multiprocessing. Default is `1`.
- **correction_phi** - P-value correction method for phi-correlation (from statsmodels.stats.multitest). Available methods: `none`, `bonferroni`, `sidak`, `holm-sidak`, `holm`, `simes-hochberg`, `hommel`, `fdr_bh`, `fdr_by`, `fdr_tsbh`, `fdr_tsbky`. Default is `bonferroni`.
- **cutoff_phi** - Significance threshold for adjusted phi-correlation p-values. Default is `0.05`.
- **correction_RR** - P-value correction method for relative risk (same methods as correction_phi). Default is `bonferroni`.
- **cutoff_RR** - Significance threshold for adjusted RR p-values. Default is `0.05`.
- **log_file** - Path/prefix for log file. If `None`, uses temporary directory with prefix `DiseaseNet_`. Default is `None`.

#### After Comorbidity Strength Estimation:

After comorbidity strength estimation, we get `com_strength_result` of `pd.DataFrame` format. And then use `com_strength_result` to binomial test.

### 3.5 Binomial Test

The 4th step of data analysis is binomial test, aiming to filter non-temperal/temperal disease pairs.

```python
binomial_result = dnt.binomial_test(
    data=data,                                        # DiseaseNetworkData object
    comorbidity_strength_result=com_strength_result,  # Comorbidity strength results
    n_process=1,                                      # Number of CPU cores (1 to disable multiprocessing)
    enforce_temporal_order=True,                      # Enforce temporal order in testing
    log_file='/your/project/path/dep.log'             # Path to log file
)

# Adjust p-values for binomial test using FDR Benjamini-Hochberg method
binomial_result = dnt.binomial_multipletests(
    df=binomial_result,                              # DataFrame with binomial test results
    correction='fdr_bh',                             # P-value correction method
    cutoff=0.05                                      # Significance threshold
)

# Save the binomial test results to a CSV file
binomial_result.to_csv('/your/project/path/dep_binomial.csv')  # Path to save binomial test results
```

- **data** - DiseaseNetworkData object containing processed disease network data.

#### Optional Parameters:

- **comorbidity_strength_result** - DataFrame containing comorbidity strength analysis results from `DiseaseNetPy.comorbidity_strength`.
- **comorbidity_network_result** - DataFrame containing comorbidity network analysis results from `DiseaseNetPy.comorbidity_network`. When provided, limits binomial test to significant disease pairs. Default is `None`.
- **n_process** - Number of parallel processes. Note: Multiprocessing is disabled for this analysis. Default is `1`.
- **correction** - P-value correction method for binomial tests. Available methods: `none`, `bonferroni`, `sidak`, `holm-sidak`, `holm`, `simes-hochberg`, `hommel`, `fdr_bh`, `fdr_by`, `fdr_tsbh`, `fdr_tsbky`. Default is `bonferroni`.
- **cutoff** - Significance threshold for adjusted binomial p-values. Default is `0.05`.
- **log_file** - Path/prefix for log file. If `None`, uses temporary directory with prefix `DiseaseNet_`. Default is `None`.
- **enforce_temporal_order** - If `True`, excludes individuals with non-temporal D1-D2 pairs. If `False`, includes all individuals. Default is `False`.

#### Kwargs Parameters:

- **phecode_d1_col** - Column for disease 1 phecode. Default is `phecode_d1`.  
- **phecode_d2_col** - Column for disease 2 phecode. Default is `phecode_d2`.
- **n_nontemporal_col** - Column for non-temporal pair counts. Default is `n_d1d2_nontemporal`.  
- **n_temporal_d1d2_col** - Column for D1→D2 temporal counts. Default is `n_d1d2_temporal`.  
- **n_temporal_d2d1_col** - Column for D2→D1 temporal counts. Default is `n_d2d1_temporal`.  
- **significance_phi_col** - Column for phi-correlation significance. Default is `phi_p_significance`.  
- **significance_RR_col** - Column for RR significance. Default is `RR_p_significance`.
- **significance_coef_col** - Column for comorbidity significance. Default is `comorbidity_p_significance`.

#### After Binomial Test:

After binomial test, we get `binomial_result` of `pd.DataFrame` format. And then use `binomial_result` to comorbidity network analysis, and trajectory analysis.

### 3.6 Comorbidity Network Analysis

The 5th step of data analysis is comorbidity network analysis, aiming to construct the comorbidity network.

```python
# Reminder:
# When using multiprocessing, ensure that the code is enclosed within the following block.
# This prevents entering a never ending loop of new process creation.

if __name__ == "__main__":
    comorbidity_result = dnt.comorbidity_network(
        data=data,                                       # DiseaseNetworkData object
        comorbidity_strength_result=com_strength_result, # Comorbidity strength results
        binomial_test_result=binomial_result,            # Binomial test results
        n_process=2,                                     # Number of parallel processes
        covariates=['age', 'BMI'],                       # Covariates to adjust for
        method='CN',                                     # Analysis method ('CN', 'PCN_PCA', 'RPCN')
        log_file='/your/project/path/dep.log'            # Path to log file
    )

    # Adjust p-values for comorbidity network analysis using FDR Benjamini-Hochberg method
    comorbidity_result = dnt.comorbidity_multipletests(
        df=comorbidity_result,                            # DataFrame with comorbidity network results
        correction='fdr_bh',                              # P-value correction method
        cutoff=0.05                                       # Significance threshold
    )

    # Save the comorbidity network analysis results to a CSV file
    comorbidity_result.to_csv('/your/project/path/dep_comorbidity.csv')  # Path to save comorbidity network results
```

- **data** - DiseaseNetworkData object containing processed disease network data.
- **comorbidity_strength_result** - DataFrame containing comorbidity strength analysis results from `DiseaseNetPy.comorbidity_strength`.

#### Optional Parameters:

- **binomial_test_result** - DataFrame containing binomial test analysis results from `DiseaseNetPy.binomial_test`. Default is `None`.
- **method** - Comorbidity network analysis method to use. Options: `RPCN`: Regularized Partial Correlation Network. `PCN_PCA`: Partial Correlation Network with PCA. `CN`: Correlation Network. Default is `RPCN`.
- **covariates** - List of phenotypic covariates to include. Always use `sex` for sex covariate. Default includes `sex` and all covariates from `phenotype_data()`. Default is `None`.
- **n_process** - Number of parallel processes. Values >1 enable multiprocessing. Default is `1`.
- **correction** - P-value correction method. Available: `none`, `bonferroni`, `sidak`, `holm-sidak`, `holm`, `simes-hochberg`, `hommel`, `fdr_bh`, `fdr_by`, `fdr_tsbh`, `fdr_tsbky`. Default is `bonferroni`.
- **cutoff** - Significance threshold for adjusted p-values. Default is `0.05`.
- **log_file** - Path/prefix for log file. If `None`, uses temporary directory with prefix `DiseaseNet_`. Default is `None`.

#### Kwargs Parameters:

  - **phecode_d1_col** - Column for disease 1 phecode. Default is `phecode_d1`.  
  - **phecode_d2_col** - Column for disease 2 phecode. Default is `phecode_d2`.  
  - **significance_phi_col** - Column for phi-correlation. Default is `phi_p_significance`.  
  - **significance_RR_col** - Column for RR. Default is `RR_p_significance`.  
  - **significance_binomial_col** - Column for binomial test. Default is `binomial_p_significance`.

  #### RPCN Method:  

  - **alpha** - L1 penalty weight. Default is `None`  
  - **auto_penalty** - Auto-determine alpha. Default is `True`  
  - **alpha_range** - Alpha search range. Default is `(1,15)`  
  - **scaling_factor** - Alpha scaling factor. Default is `1`

  #### PCN_PCA Method:  

  - **n_PC** - Number of principal components. Default is `5`  
  - **explained_variance** - Variance threshold. Default is `None`

### 3.7 Trajectory Analysis

The 6th step of data analysis is trajectory analysis, aiming to construct disease trajectory.

```python
# Reminder:
# When using multiprocessing, ensure that the code is enclosed within the following block.
# This prevents entering a never ending loop of new process creation.

if __name__ == "__main__":
    trajectory_result = dnt.disease_trajectory(
        data=data,                                       # DiseaseNetworkData object
        comorbidity_strength_result=com_strength_result, # Comorbidity strength results
        binomial_test_result=binomial_result,            # Binomial test results
        method='RPCN',                                   # Trajectory analysis method ('CN', 'PCN_PCA', 'RPCN')
        n_process=2,                                     # Number of parallel processes
        matching_var_dict={'age': 2, 'sex': 'exact'},    # Matching variables and criteria
        matching_n=5,                                    # Number of matched controls per case
        enforce_time_interval=False,                     # Enforce time interval in trajectory analysis
        covariates=['age', 'BMI'],                       # Covariates to adjust for
        log_file='/your/project/path/dep.log'            # Path to log file
    )

    # Adjust p-values for trajectory analysis using FDR Benjamini-Hochberg method
    trajectory_result = dnt.trajectory_multipletests(
        df=trajectory_result,                             # DataFrame with trajectory analysis results
        correction='fdr_bh',                              # P-value correction method
        cutoff=0.05                                       # Significance threshold
    )

    # Save the trajectory analysis results to a CSV file
    trajectory_result.to_csv('/your/project/path/dep_trajectory.csv')  # Path to save trajectory analysis results
```

- **data** - DiseaseNetworkData object containing processed disease network data.
- **comorbidity_strength_result** - DataFrame containing comorbidity strength analysis results from `DiseaseNetPy.comorbidity_strength`.
- **binomial_test_result** - DataFrame containing binomial test analysis results from `DiseaseNetPy.binomial_test`.

#### Optional Parameters:

- **method** - Comorbidity network analysis method: `RPCN`: Regularized Partial Correlation Network. `PCN_PCA`: Partial Correlation Network with PCA. `CN`: Correlation Network. Default is `RPCN`.
- **matching_var_dict** - Dictionary specifying matching variables and criteria: Categorical/binary: `{'var':'exact'}`. Continuous: `{'var':max_diff}` (scalar > 0). Always use `'sex'` for sex matching. Default is `{'sex':'exact'}`.
- **matching_n** - Maximum matched controls per case. Default is `2`.
- **max_n_cases** - Maximum D2 cases to include (random sampling if exceeded). Default is `np.inf`.
- **global_sampling** - `True` for single sampling across all pairs, `False` for per-pair sampling. Default is `False`.
- **covariates** - List of phenotypic covariates (use `'sex'` for sex). Exclude matching variables. Default is `None`.
- **n_process** - Number of parallel processes (>1 enables multiprocessing). Default is `1`.
- **correction** - P-value correction method: `none`, `bonferroni`, `sidak`, `holm-sidak`, `holm`, `simes-hochberg`, `hommel`, `fdr_bh`, `fdr_by`, `fdr_tsbh`, `fdr_tsbky`. Default is `bonferroni`.
- **cutoff** - Significance threshold for adjusted p-values. Default is `0.05`.
- **log_file** - Log file path/prefix. If `None`, uses temp dir with `DiseaseNet_` prefix. Default is `None`.

#### Kwargs Parameters:

##### RPCN Method:

- **alpha**: L1 penalty weight (ignored if auto_penalty). Default is `None`
- **auto_penalty**: Auto-determine optimal alpha. Default is `True`
- **alpha_range**: Alpha search range. Default is `(1,15)`
- **scaling_factor**: Alpha scaling factor. Default is `1`

##### PCN_PCA Method:

- **n_PC**: Principal components count. Default is `5`
- **explained_variance**: Variance threshold (overrides n_PC). Default is `None`

##### Analysis Option:

- **enforce_time_interval** - Apply min/max time intervals for D2 outcome determination. Default is `True`.

##### Column Mappings:

  - **phecode_d1_col**: Disease 1 phecode column. Default is `phecode_d1`
  - **phecode_d2_col**: Disease 2 phecode column. Default is `phecode_d2`
  - **significance_phi_col**: Phi-correlation column. Default is `phi_p_significance`
  - **significance_RR_col**: RR column. Default is `RR_p_significance`
  - **significance_binomial_col**: Binomial test column. Default is `binomial_p_significance`

## 4. Visualization

### 4.1 Initializing the plot object

First, from DiseaseNetPy.visualization import plot and initialize an empty `Plot` object, specifying the result of PheWAS analysis, comorbidity network analysis, disease trajectory analysis, and any optional parameters if needed:

```python
# Create Plot object
from diseasenetpy.visualization import Plot

# cohort/matched cohort
result_plot = Plot(
    phewas_result=phewas_result,                         # DataFrame with PheWAS results
    comorbidity_network_result=comorbidity_result,       # DataFrame with comorbidity network results
    disease_trajectory_result=trajectory_result,         # DataFrame with trajectory analysis results
    exposure=495.2,                                      # Phecode of exposure. Default to None, means that it's a exposed-only study
    exposure_size=15,                                    # Relative size scaling factor for the exposure node in visualizations
    exposure_location=(0,0,0),                           # Custom 3D coordinates (x,y,z) for positioning the exposure node
    source='phecode_d1',                                 # Column name indicating source/antecedent diseases in disease pair
    target='phecode_d2',                                 # Column name indicating target/consequent diseases in disease pair
    phewas_phecode='phecode',                            # Column name of pd.DataFrame for phecode disease name in the PHEWAS result
    phewas_number='N_cases_exposed',                     # Column name of pd.DataFrame for case counts in the PHEWAS result
    system_col='system',                                 # Column name of pd.DataFrame for disease system classifications
    col_disease_pair='name_disease_pair',                # Column name of pd.DataFrame for disease pair identifiers
    filter_phewas_col='phewas_p_significance',           # Column name of pd.DataFrame for PHEWAS significance filter
    filter_comorbidity_col='comorbidity_p_significance', # Column name of pd.DataFrame for comorbidity significance filter
    filter_trajectory_col='trajectory_p_significance',   # Column name of pd.DataFrame for trajectory significance filter
)

# or exposed-only cohort
result_plot = Plot(
    phewas_result=phewas_result,                         # DataFrame with PheWAS results
    comorbidity_network_result=comorbidity_result,       # DataFrame with comorbidity network results
    disease_trajectory_result=trajectory_result,         # DataFrame with trajectory analysis results
    exposure=None,                                       # Phecode of exposure. Default to None, means that it's a exposed-only study
    exposure_size=None,                                  # Relative size scaling factor for the exposure node in visualizations
    exposure_location=None,                              # Custom 3D coordinates (x,y,z) for positioning the exposure node
    source='phecode_d1',                                 # Column name indicating source/antecedent diseases in disease pair
    target='phecode_d2',                                 # Column name indicating target/consequent diseases in disease pair
    phewas_phecode='phecode',                            # Column name of pd.DataFrame for phecode disease name in the PHEWAS result
    phewas_number='N_cases_exposed',                     # Column name of pd.DataFrame for case counts in the PHEWAS result
    system_col='system',                                 # Column name of pd.DataFrame for disease system classifications
    col_disease_pair='name_disease_pair',                # Column name of pd.DataFrame for disease pair identifiers
    filter_phewas_col='phewas_p_significance',           # Column name of pd.DataFrame for PHEWAS significance filter
    filter_comorbidity_col='comorbidity_p_significance', # Column name of pd.DataFrame for comorbidity significance filter
    filter_trajectory_col='trajectory_p_significance',   # Column name of pd.DataFrame for trajectory significance filter
)
```

- **comorbidity_result** - DataFrame containing comorbidity network analysis results with: Non-temporal disease pairs (D1-D2). Association metrics (beta coefficients, p-values). Significance indicators (True/False).
- **trajectory_result** - DataFrame containing temporal disease trajectory analysis with: Temporal disease pairs (source→target). Temporal association metrics. Significance indicators (True/False).
- **phewas_result** - DataFrame containing PheWAS analysis results with: Phecode diseases. Effect sizes (hazard ratios). Case counts. Disease system classifications

#### Optional Parameters:

- **exposure** - Phecode identifier for primary exposure variable. Highlights exposure-disease relationships. Default is `None` (exposed-only cohort).
- **exposure_location** - Custom 3D coordinates (x,y,z) for exposure node positioning. Default is `None` (auto-positioned at (0,0,0)).
- **exposure_size** - Relative size scaling factor for exposure node. Default is `None` (exposed-only cohort).
- **source**- Source disease column. Default is `phecode_d1`.
- **target**: Target disease column. Default is `phecode_d2`.
- **phewas_phecode**: Phecode column. Default `phecode`.
- **phewas_number**: Case count column. Default `N_cases_exposed`.
- **system_col**: Disease system column. Default `system`.
- **col_disease_pair**: Disease pair identifier column. Default `name_disease_pair`.
- **filter_phewas_col**: PheWAS significance column. Default `phewas_p_significance`.
- **filter_comorbidity_col**: Comorbidity significance column. Default `comorbidity_p_significance`.
- **filter_trajectory_col**: Trajectory significance column. Default `trajectory_p_significance`.

#### Kwargs Parameters:

- **SYSTEM** - List of phecode systems to visualize. Available systems (17 total):
  ['neoplasms', 'genitourinary', 'digestive', 'respiratory',
   'infectious diseases', 'mental disorders', 'musculoskeletal',
   'hematopoietic', 'dermatologic', 'circulatory system',
   'neurological', 'endocrine/metabolic', 'sense organs',
   'injuries & poisonings', 'congenital anomalies', 'symptoms',
   'others']
- **COLOR** - Custom colors for systems (hex, RGB, or named colors). Must match SYSTEM order. Default is 
    ['#F46D5A', '#5DA5DA', '#5EBCD1', '#C1D37F',
    '#CE5A57', '#A5C5D9', '#F5B36D', '#7FCDBB', '#ED9A8D',
    '#94B447', '#8C564B', '#E7CB94', '#8C9EB2', '#E0E0E0', "#F1C40F",
    '#9B59B6', '#4ECDC4', '#6A5ACD']

#### After initializing the plot object:

After initializing the plot object, use `result_plot` of `Plot` object to visualization.

### 4.2 PheWAS Plot

Generates a circular PheWAS (Phenome-Wide Association Study) plot. Creates a polar bar plot visualizing disease associations across different disease categories (systems). For a cohort/matched cohort study, the figure shows hazard ratios between exposure and outcome diseases. While the figure shows exposed number of diseases.

```python
# phewas plot
result_plot.phewas_plot(
    path="/your/project/path/",                          # Output file path for saving the plot
    col_coef="phewas_coef",                              # Column name for effect size coefficients (default: "phewas_coef")
    col_system="system",                                 # Column name for disease system/category (default: "system")
    col_se="phewas_se",                                  # Column name for standard errors (default: "phewas_se")
    col_disease="disease",                               # Column name for disease names (default: "disease")
    is_exposure_only=False,                              # Identifier of exposure (default: False)
    col_exposure='N_cases_exposed'                       # Column name for exposure number (default: "N_cases_exposed")
)
```
- **path** - Output file path for saving the plot (including filename and extension)

#### Optional Parameters:

- **col_coef** - Column name containing effect size coefficients (e.g., hazard ratios or odds ratios). Default is `phewas_coef`
- **col_system** - Column name containing disease system/category classifications. Default is `"system"`
- **col_se** - Column name containing standard errors for effect sizes. Default is `phewas_se`
- **col_disease** - Column name containing disease names/descriptions. Default is `"disease"`
- **is_exposure_only** - Boolean flag indicating whether the plot is for an exposure-only cohort study. Default is `False`
- **col_exposure** - Column name containing case counts for exposed individuals. Default is `"N_cases_exposed"`

#### After PheWAS plot:

After generating the PheWAS plot, the visualization will be exported as an image file in your preferred format (.jpg, .svg, or .png).

### 4.3 Comorbidity Network Plot

Generate a 2D visualization of the comorbidity network.

```python
# comorbidity network visualization
result_network.comorbidity_network_plot(
    path="/your/project/path/",                          # Output file path for saving HTML visualization
    max_radius=180.0,                                    # Maximum radial position for nodes (default: 180.0)
    min_radius=35.0,                                     # Minimum radial position for nodes (default: 35.0)
    size_reduction=0.5,                                  # Scaling factor for node sizes (default: 0.5)
    cluster_reduction_ratio=0.4,                         # Compression factor for cluster layout (default: 0.4)
    cluster_weight="comorbidity_beta",                   # Edge weight metric for clustering (default: "comorbidity_beta")
    line_width=1.0,                                      # Width of comorbidity lines (default: 1.0)
    line_color="black",                                  # Color of comorbidity lines (default: "black")
    layer_distance=40.0,                                 # Distance between concentric circles (default: 40.0)
    font_style="Times New Roman"                         # Font family for text elements (default: "Times New Roman")
)
```
- **path** - Output file path for saving the interactive HTML visualization.  

#### Optional Parameters:

- **max_radius** - Maximum radial position for nodes (in pixels). Controls outer boundary of the network. Default is `90.0`.
- **min_radius** - Minimum radial position for nodes (in pixels). Controls inner boundary. Default is `35.0`
- **layer_distance** - Spacing between concentric circles (in pixels). Affects radial grouping. Default is `40.0`.
- **size_reduction** - Scaling factor for node diameters (0-1 range). Adjusts visual prominence. Default is `0.5`
- **line_width** - Stroke width for comorbidity connections (in pixels). Default is `1.0`
- **line_color** - Color specification for comorbidity lines. Accepts: Named colors (e.g., `steelblue`). HEX codes (e.g., `#4682B4`). RGB tuples (e.g., `(70, 130, 180)`). Default is `black`
- **cluster_reduction_ratio** - Compression factor (0-1) for cluster tightness. Lower values create denser groupings. Default is `0.4`.
- **cluster_weight** - Edge attribute used for clustering calculations. Typically the association strength metric. Default is `comorbidity_beta`
- **font_style** - Font family for all text elements. Use web-safe fonts or loaded font families. Default is `Times New Roman`.

#### After Comorbidity Network Plot:

After generating the comorbidity network plot, the visualization will be exported as an image file in your preferred format (.html).

### 4.4 Disease Trajectory Plot

Creates 2D network plots showing disease trajectories within each cluster, with nodes positioned hierarchically based on trajectory relationships. Each cluster is saved as a separate image file.

```python
# Disease trajectory visualization
result_network.trajectory_plot(
    path="/your/project/path/",                          # Directory path to save output images
    cluster_weight="comorbidity_beta",                   # Edge weight metric used for clustering (default: "comorbidity_beta")
)
```

- **path** - Directory path where output visualization images will be saved. *Note:* Include trailing slash for proper path resolution (e.g., `/output/plots/`)

#### Optional Parameters:

- **cluster_weight** Specifies the edge weight metric used for network clustering calculations. Default is `comorbidity_beta`.

#### After Disease Trajectory Plot:

After generating the disease trajectory plot, the visualization will be exported as an image file in your preferred format (.jpg, .svg, or .png).

### 4.5 Three Dimension Plot

Generates and saves a 3D visualization of comorbidity and disease trajectory networks.

```python
# three-dimension visualization
result_network.plot_3d(
    path="/your/project/path/",                          # File path to save the HTML visualization
    max_radius=180.0,                                    # Maximum radial distance for node placement (default: 180.0)
    min_radius=35.0,                                     # Minimum radial distance for node placement (default: 35.0)
    line_color="black",                                  # Color for trajectory lines (default: "black")
    line_width=1.0,                                      # Width for trajectory lines (default: 1.0)
    size_reduction=0.5,                                  # Scaling factor for node sizes (default: 0.5)
    cluster_reduction_ratio=0.4,                         # Cluster compression factor for layout (default: 0.4)
    cluster_weight="comorbidity_beta",                   # Edge weight metric used for clustering (default: "comorbidity_beta")
    layer_distance=40.0,                                 # Vertical distance between layers (default: 40.0)
    layout_width=900.0,                                  # Figure width in pixels (default: 900.0)
    layout_height=900.0,                                 # Figure height in pixels (default: 900.0)
    font_style='Times New Roman',                        # Font family for text elements (default: 'Times New Roman')
    font_size=15.0,                                      # Base font size in points (default: 15.0)
)
```

- **path** - Absolute or relative file path to save the interactive HTML visualization.

#### Optional Parameters:
- **max_radius** - Maximum radial distance (in pixels) from center for node placement. Default is `180.0`.
- **min_radius** - Minimum radial distance (in pixels) from center for node placement. Default is `35.0`.
- **layer_distance** - Vertical spacing (in pixels) between concentric layers. Default is `40.0`.
- **layout_width** - Total figure width in pixels. Default is `900.0`.
- **layout_height** - Total figure height in pixels. Default is `900.0`.
- **line_color** - Color specification for trajectory pathways. Accepts: CSS color names (e.g., `steelblue`). HEX codes (e.g., `#4682B4`). RGB tuples (e.g., `(70, 130, 180)`). Default is `black`.
- **line_width** - Stroke width (in pixels) for trajectory lines. Default is `1.0`.
- **size_reduction** - Multiplicative scaling factor for node diameters (0.1-1.0). Default is `0.5`.
- **cluster_reduction_ratio** - Compression factor (0.1-1.0) for cluster density. Default is `0.4`.
- **cluster_weight** - Edge metric used for clustering calculations. Default is `comorbidity_beta`.
- **font_style** - Font family for all text elements. Use web-safe fonts. Default is `Times New Roman`
- **font_size** - Base font size in points for all text elements. Default is `15.0`.

#### After Disease Trajectory Plot:

After generating the disease trajectory plot, the visualization will be exported as an image file in your preferred format (.html).

## API Reference

### Classes

#### `DiseaseNetworkData`

```python
class DiseaseNetworkData()
```

A class for handling disease network data creation and operations, for use in the DiseaseNetPy module.

**Parameters:**

- `study_design : str`
  - Specify the type of study design, either `"cohort"`, `"matched cohort"`, or `"registry"`.
- `phecode_level : int`
  - The level of phecode to use for analysis:
    - Level 1: Corresponds to 3-digit ICD-10 codes with a total of 585 medical conditions.
    - Level 2: Corresponds to 4-digit ICD-10 codes with a total of 1257 medical conditions.
  - Recommendation:
    - Level 2 offers a more granular analysis suitable for larger studies.
    - Level 1 is recommended for smaller studies to maintain statistical power.
- `min_required_icd_codes : int`
  - The level of phecode to use for analysis, where level 1 (with a total of 585 medical conditions) corresponds to 3-digit ICD-10 codes and level 2 (a total of 1257 medical conditions) to 4-digit ICD-10 codes.
  - Recommendation:
    - For larger studies, level 2 phecodes may enhance result interpretation.
    - For smaller studies, level 1 is recommended to maintain statistical power.
- `date_fmt : str, default='%Y-%m-%d'`
  - The format of the date fields in your phenotype and medical records data.
- `phecode_version : str, default='1.2'`
  - The version of the phecode system used for converting diagnosis codes. Currently, only version 1.2 is supported.

------

#### `Plot`

```python
class Plot()
```

This class integrates and visualizes disease relationships from three complementary analyses:
1. Phenome-Wide Association Study (PHEWAS) results
2. Comorbidity network analysis
3. Disease trajectory analysis
**Parameters:**

- `comorbidity_result : pd.DataFrame`
  - Result dataframe from comorbidity network analysis containing:
    - Non-temporal disease pairs (D1-D2)
    - Association metrics (e.g., beta coefficients, p-values)
    - Significance identifier (True or False)
  
- `trajectory_result : pd.DataFrame`
  - Result dataframe from temporal disease trajectory analysis containing:
    - Temporal disease pairs (source→target)
    - Temporal association metrics (e.g., beta coefficients, p-values)
    - Significance identifier (True or False)
  
- `phewas_result : pd.DataFrame`
  - Result dataframe from PHEWAS analysis containing:
    - Phecode disease identifiers
    - Effect sizes (e.g., hazard ratios)
    - Case counts
    - Disease system classifications

- `exposure : float, optional`
  - Phecode identifier for the primary exposure variable of interest.
  - Used to highlight exposure-disease relationships in visualizations.
  - Default: None (exposed-only cohort)
  
- `exposure_location : Tuple[float], optional`
  - Custom 3D coordinates (x,y,z) for positioning the exposure node.
  - Default: None (automatically positioned at origin (0,0,0))
  
- `exposure_size : float, optional`
  - Relative size scaling factor for the exposure node.
  - Default: None (uses automatic sizing)

- `source : str, default='phecode_d1'`
  - Column name for source/antecedent diseases in disease pairs
  
- `target : str, default='phecode_d2'`
  - Column name for target/consequent diseases in disease pairs

- `phewas_phecode : str, default='phecode'`
  - Column name for phecode disease identifiers in PHEWAS results
  
- `phewas_number : str`
  - Column name for case counts in PHEWAS results
  
- `system_col : str`
  - Column name for disease system classifications
  
- `col_disease_pair : str`
  - Column name for disease pair identifiers
  
- `filter_phewas_col : str`
  - Column name for PHEWAS significance filter
  
- `filter_comorbidity_col : str`
  - Column name for comorbidity significance filter
  
- `filter_trajectory_col : str`
  - Column name for trajectory significance filter

- `SYSTEM : List[str], optional`
  - Subset of phecode disease systems to visualize. 
  - Available systems (17 total):
    ```python
    ['neoplasms',
     'genitourinary',
     'digestive',
     'respiratory',
     'infectious diseases',
     'mental disorders',
     'musculoskeletal',
     'hematopoietic',
     'dermatologic',
     'circulatory system',
     'neurological',
     'endocrine/metabolic',
     'sense organs',
     'injuries & poisonings',
     'congenital anomalies',
     'symptoms',
     'others']
    ```
  - Default: All systems present in PHEWAS results

- `COLOR : List[str], optional`
  - Custom colors for disease systems (hex/RGB/tuple format).
  - Must match length of SYSTEM parameter.
  - Default color palette:
    ```python
    ['#F46D5A', '#5DA5DA', '#5EBCD1', '#C1D37F',
     '#CE5A57', '#A5C5D9', '#F5B36D', '#7FCDBB',
     '#ED9A8D', '#94B447', '#8C564B', '#E7CB94',
     '#8C9EB2', '#E0E0E0', '#F1C40F', '#9B59B6',
     '#4ECDC4', '#6A5ACD']
    ```

**Notes:**
- All input dataframes must use consistent phecode identifiers
- Significant results are filtered using specified significance columns
- Node sizes are proportional to case counts by default
- Color schemes are automatically assigned by disease system
- Column mappings retain defaults if analysis outputs use standard names

##### Methods

###### `phenotype_data`

```python
phenotype_data(
  self, 
  phenotype_data_path:str, 
  column_names:dict, 
  is_single_sex: bool=False, 
  covariates:list, 
  force:bool=False
)
```

Merges phenotype and medical records data into the main data attribute.

**Parameters:**

- `phenotype_data_path : str`

  - File path containing phenotype data in CSV or TSV format with a header row.

- `column_names : dict`

  - Maps required variable names to their corresponding identifiers in the dataset.

  - Expected Keys:

    - `'Participant ID'`, `'Index date'`, `'End date'`, `'Exposure'` (for cohort and matched cohort study), `'Sex'`, and `'Match ID'` (for matched cohort study).

  - Example:

    ```python
    column_names = {
        'Participant ID': 'eid',
        'Exposure': 'status',
        'Sex': 'sex',
        'Index date': 'index_date',
        'End date': 'final_date',
        'Match ID': 'group'
    }
    ```

  - Coding Requirements:

    - `'Exposure'`: 0 (unexposed), 1 (exposed)
    - `'Sex'`: 1 (female), 0 (male)

  - Date Format:

    - Must follow `%Y-%m-%d` unless specified otherwise.

- `covariates : list`

  - List of additional covariate variable names, e.g., `['age', 'BMI']`.
  - If no additional covariates are included, provide an empty list.
  - Handles missing values by removing individuals with missing continuous variables and categorizing those with missing categorical variables separately.

- `is_single_sex : bool, default=False`

  - Boolean flag indicating if the dataset contains only one sex (male or female)

- `force : bool, default=False`

  - If `True`, overwrites existing data attributes even if they contain data.

**Returns:**

- `None`
  - Modifies the object's main data attribute in-place.

------

###### `merge_medical_records`

```python
merge_medical_records(
    self, 
    medical_records_data_path:str, 
    diagnosis_code:str,
    column_names:dict, 
    date_fmt:str=None, 
    chunksize:int=1000000
)
```

Merges the loaded phenotype data with one or more medical records data.

**Parameters:**

- `medical_records_data_path : str`

  - File path containing medical records data in CSV or TSV format.
  - Required Columns:
    - Participant ID, Diagnosis code (specified type), Date of diagnosis.

- `diagnosis_code : str`

  - Diagnosis code type used in the medical records data.
  - Valid Options: `'ICD-9-CM'`, `'ICD-9-WHO'`, `'ICD-10-CM'`, `'ICD-10-WHO'`.
  - Note: Mixing different ICD code types within the same file is not allowed.

- `column_names : dict`

  - Maps required variable names to their corresponding identifiers in the medical records dataset.

  - Required Keys: `'Participant ID'`, `'Diagnosis code'`, `'Date of diagnosis'`.

  - Example:

    ```python
    column_names = {
        'Participant ID': 'eid',
        'Diagnosis code': 'ICD',
        'Date of diagnosis': 'date'
    }
    ```

  - Date Format:

    - Must follow `%Y-%m-%d` unless specified otherwise.

- `date_fmt : str, optional`

  - The format of the date fields in your medical records data.
  - Defaults to the same format as phenotype data if not specified.

- `chunksize : int, default=1_000_000`

  - Number of rows per chunk to read, useful for large datasets.

**Returns:**

- `None`
  - Modifies the object's main data attribute in-place.

------

###### `modify_phecode_level(phecode_level: int)`

Modifies the phecode level setting.

**Parameters:**

- `phecode_level : int`
  - The level of phecode to use for analysis:
    - Level 1: 3-digit ICD-10 codes (585 conditions).
    - Level 2: 4-digit ICD-10 codes (1257 conditions).

**Returns:**

- `None`
  - Modifies the object's main data attribute in-place.

------

###### `Table1`

Generates a descriptive Table 1 summary of phenotype data.

```python
Table1(
    continuous_stat_mode: str='auto'
)
```

**Parameters:**

- `continuous_stat_mode : str` (default: `'auto'`)  
  Specifies the statistical display method for continuous variables:  
  - `'auto'`: Automatically selects summary statistics based on normality test (Shapiro-Wilk or similar)  
  - `'normal'`: Forces normal distribution display (mean ± standard deviation)  
  - `'nonnormal'`: Forces non-parametric display (median [interquartile range])  

**Returns:**

- `pd.DataFrame`  
  A formatted summary table containing:  
  - All variables in the phenotype dataset  
  - Appropriate descriptive statistics per variable  
  - Group comparisons when applicable (p-values from t-test/ANOVA or Mann-Whitney/Kruskal-Wallis)  

------

###### `disease_pair`

```python
disease_pair(
    self, 
    phewas_result:pd.DataFrame, 
    min_interval_days:int=0, 
    max_interval_days:int=np.inf, 
    force:bool=False,
    n_process:int=1,
    **kwargs
)
```

Constructs disease pairs based on PheWAS results.

**Parameters:**

- `phewas_result : pd.DataFrame`
  - DataFrame containing PheWAS analysis results.
- `min_interval_days : int/float, default=0`
  - Minimum required time interval (in days) between diagnosis dates for temporal pairs.
- `max_interval_days : int/float, default=np.inf`
  - Maximum allowed time interval (in days) between diagnosis dates for temporal pairs.
- `force : bool, default=False`
  - If `True`, overwrites existing data attributes.
- `n_process : int, default=1`
  - Number of processes to use for parallel processing.
  - Multiprocessing is enabled when `n_process` is set to a value greater than one.
- `**kwargs`
  - Additional parameters:
    - `phecode_col : str, default='phecode'`
    - `significance_col : str, default='phewas_p_significance'`

**Returns:**

- `None`
  - Modifies the object's main data attribute in-place.

------

###### `load`

```python
load(self, file:str, force:bool=False)
```

Load data from a `.pkl.gz` (gzip-compressed pickle) file and restore the attributes to this DiseaseNetworkData object.

**Parameters:**

- `file : str`
  - The filename with or without `.pkl.gz` extension.
- `force : bool, default=False`
  - If `True`, overwrites existing data attributes.

**Returns:**

- `None`

------

###### `save`

```python
save(self, file:str)
```

Saves the DiseaseNetworkData object's attributes to a `.pkl.gz` (gzip-compressed pickle) file

**Parameters:**

- `file : str`
  - The filename with or without `.pkl.gz` extension.

**Returns:**

- `None`

------

###### `load_npz`

```python
save_npz(self, file:str, force: bool=False)
```

Load data from a `.npz` file and restore the attributes to this DiseaseNetworkData object.

**Parameters:**

- `file : str`
  - The filename with or without `.npz` extension.
- `force : bool, default=False`
  - If `True`, overwrites existing data attributes.

**Returns:**

- `None`  
  The method modifies the filesystem but returns nothing  

------

###### `save_npz`

```python
save_npz(self, file:str)
```

Saves the `DiseaseNetworkData` object's attributes to a compressed NumPy `.npz` file.

**Parameters:**

- `file : str`  
  Target filename/path for saving the data:  
  - If extension is omitted, `.npz` will be automatically appended  
  - Supports both relative and absolute paths  

**Returns:**

- `None`  
  The method modifies the filesystem but returns nothing  


### Functions

#### `phewas`

```python
dnt.phewas(
    data:DiseaseNetworkData, 
    covariates:list=None,
    proportion_threshold:float=None, 
    n_threshold:int=None,
    n_process:int=1, 
    correction:str='bonferroni', 
    cutoff:float=0.05,
    system_inc:list=None, 
    system_exl:list=None,
    phecode_inc:list=None, 
    phecode_exl:list=None, 
    log_file:str=None,
    lifelines_disable:bool=False
) -> pd.DataFrame
```

Conducts Phenome-wide Association Studies (PheWAS) using the specified DiseaseNetworkData object.

**Parameters:**

- `data : DiseaseNetworkData`
  - The DiseaseNetworkData object containing your study data.
- `covariates : list, default=None`
  - List of phenotypic covariates to include in the model.
  - Defaults to include `'sex'` and all covariates specified in `phenotype_data()`.
  - Use `'sex'` to include sex as a covariate.
- `proportion_threshold : float`
  - Minimum proportion of cases within the exposed group required for a phecode to be included.
  - Mutually exclusive with `n_threshold`.
- `n_threshold : int`
  - Minimum number of cases within the exposed group required for a phecode to be included.
  - Mutually exclusive with `proportion_threshold`.
- `n_process : int, default=1`
  - Number of parallel processes to use. Set greater than one to enable multiprocessing.
- `correction : str, default='bonferroni'`
  - Method for p-value correction.
  - Options: `none`, `bonferroni`, `sidak`, `holm-sidak`, `holm`, `simes-hochberg`, `hommel`, `fdr_bh`, `fdr_by`, `fdr_tsbh`, `fdr_tsbky`.
- `cutoff : float, default=0.05`
  - Significance threshold for adjusted PheWAS p-values.
- `system_inc : list, default=None`
  - List of phecode systems to include.
  - Mutually exclusive with `system_exl`.
  - Eligible Systems:
    - `circulatory system`, `congenital anomalies`, `dermatologic`, `digestive`,
    - `endocrine/metabolic`, `genitourinary`, `hematopoietic`, `infectious diseases`,
    - `injuries & poisonings`, `mental disorders`, `musculoskeletal`, `neoplasms`,
    - `neurological`, `pregnancy complications`, `respiratory`, `sense organs`, `symptoms`, `others`.
- `system_exl : list, default=None`
  - List of phecode systems to exclude.
  - Mutually exclusive with `system_inc`.
- `phecode_inc : list, default=None`
  - Specific phecodes to include.
  - Mutually exclusive with `phecode_exl`.
- `phecode_exl : list, default=None`
  - Specific phecodes to exclude.
  - Mutually exclusive with `phecode_inc`.
- `log_file : str, default=None`
  - Path and prefix for the log file. Defaults to temporary files directory with prefix `DiseaseNet_`.
- `lifelines_disable : bool, default=False`
  - If `True`, disables the use of lifelines. Lifelines are more robust but require longer fitting time.

**Returns:**

- `pd.DataFrame`
  - DataFrame containing the results of the PheWAS analysis.

------

#### `phewas_multipletests`

```python
dnt.phewas_multipletests(
    df:pd.DataFrame, 
    correction:str='bonferroni', 
    cutoff:float=0.05
) -> pd.DataFrame
```

Adjusts PheWAS p-values for multiple comparisons using specified correction methods.

**Parameters:**

- `df : pd.DataFrame`
  - DataFrame containing PheWAS results.
- `correction : str, default='bonferroni'`
  - Method for p-value correction.
  - Options: `none`, `bonferroni`, `sidak`, `holm-sidak`, `holm`, `simes-hochberg`, `hommel`, `fdr_bh`, `fdr_by`, `fdr_tsbh`, `fdr_tsbky`.
- `cutoff : float, default=0.05`
  - Significance threshold for adjusted PheWAS p-values.

**Returns:**

- `pd.DataFrame`
  - DataFrame with adjusted p-values.

------

#### `comorbidity_strength`

```python
dnt.comorbidity_strength(
    data:DiseaseNetworkData, 
    proportion_threshold:float=None, 
    n_threshold:int=None, 
    n_process:int=1, 
    log_file:str=None, 
    correction_phi:str='bonferroni', 
    cutoff_phi:float=0.05, 
    correction_RR:str='bonferroni', 
    cutoff_RR:float=0.05
) -> pd.DataFrame
```

Conducts comorbidity strength estimation among exposed individuals for all possible disease pairs.

**Parameters:**

- `data : DiseaseNetworkData`
  - The DiseaseNetworkData object containing your study data.
- `proportion_threshold : float`
  - Minimum proportion of individuals in the exposed group where a disease pair must co-occur to be included.
  - Mutually exclusive with `n_threshold`.
- `n_threshold : int`
  - Minimum number of individuals in the exposed group where a disease pair must co-occur to be included.
  - Mutually exclusive with `proportion_threshold`.
- `n_process : int, default=1`
  - Number of parallel processes to use. Set greater than one to enable multiprocessing.
- `log_file : str, default=None`
  - Path and prefix for the log file.
- `correction_phi : str, default='bonferroni'`
  - P-value correction method for phi-correlation.
  - Options: `none`, `bonferroni`, `sidak`, `holm-sidak`, `holm`, `simes-hochberg`, `hommel`, `fdr_bh`, `fdr_by`, `fdr_tsbh`, `fdr_tsbky`.
- `cutoff_phi : float, default=0.05`
  - Significance threshold for adjusted phi-correlation p-values.
- `correction_RR : str, default='bonferroni'`
  - P-value correction method for Relative Risk (RR).
  - Options: Same as `correction_phi`.
- `cutoff_RR : float, default=0.05`
  - Significance threshold for adjusted RR p-values.

**Returns:**

- `pd.DataFrame`
  - DataFrame containing comorbidity strength estimation results.

------

#### `comorbidity_strength_multipletests`

```python
dnt.comorbidity_strength_multipletests(
    df:pd.DataFrame, 
    correction_phi:str='bonferroni', 
    cutoff_phi:float=0.05, 
    correction_RR:str='bonferroni', 
    cutoff_RR:float=0.05
) -> pd.DataFrame
```

Adjusts comorbidity strength p-values (phi-correlation and RR) for multiple comparisons.

**Parameters:**

- `df : pd.DataFrame`
  - DataFrame containing comorbidity strength results.
- `correction_phi : str, default='bonferroni'`
  - P-value correction method for phi-correlation.
- `cutoff_phi : float, default=0.05`
  - Significance threshold for adjusted phi-correlation p-values.
- `correction_RR : str, default='bonferroni'`
  - P-value correction method for Relative Risk (RR).
- `cutoff_RR : float, default=0.05`
  - Significance threshold for adjusted RR p-values.

**Returns:**

- `pd.DataFrame`
  - DataFrame with adjusted comorbidity strength results.

------

#### `binomial_test`

```python
dnt.binomial_test(
    data:DiseaseNetworkData, 
    comorbidity_strength_result:pd.DataFrame, 
    comorbidity_network_result: pd.DataFrame=None,
    n_process:int=1, 
    log_file:str=None, 
    correction:str='bonferroni', 
    cutoff:float=0.05, 
    enforce_temporal_order:bool=False, 
    **kwargs
) -> pd.DataFrame
```

Conducts binomial tests for disease pairs with significant comorbidity strength to identify significant temporal orders.

**Parameters:**

- `data : DiseaseNetworkData`
  - The DiseaseNetworkData object containing your study data.
- `comorbidity_strength_result : pd.DataFrame`
  - DataFrame containing comorbidity strength analysis results.
- `n_process : int, default=1`
  - Number of parallel processes to use. Set to `1` to disable multiprocessing.
- `log_file : str, default=None`
  - Path and prefix for the log file.
- `correction : str, default='bonferroni'`
  - P-value correction method.
  - Options: `none`, `bonferroni`, `sidak`, `holm-sidak`, `holm`, `simes-hochberg`, `hommel`, `fdr_bh`, `fdr_by`, `fdr_tsbh`, `fdr_tsbky`.
- `cutoff : float, default=0.05`
  - Significance threshold for adjusted binomial p-values.
- `enforce_temporal_order : bool, default=False`
  - If `True`, excludes individuals with non-temporal D1-D2 pairs.
- `**kwargs`
  - Additional parameters:
    - `phecode_d1_col : str, default='phecode_d1'`
    - `phecode_d2_col : str, default='phecode_d2'`
    - `n_nontemporal_col : str, default='n_d1d2_nontemporal'`
    - `n_temporal_d1d2_col : str, default='n_d1d2_temporal'`
    - `n_temporal_d2d1_col : str, default='n_d2d1_temporal'`
    - `significance_phi_col : str, default='phi_p_significance'`
    - `significance_RR_col : str, default='RR_p_significance'`

**Returns:**

- `pd.DataFrame`
  - DataFrame containing binomial test results.

------

#### `binomial_multipletests`

```python
dnt.binomial_multipletests(
    df:pd.DataFrame,
    correction:str='bonferroni', 
    cutoff:float=0.05
) -> pd.DataFrame
```

Adjusts binomial test p-values for multiple comparisons.

**Parameters:**

- `df : pd.DataFrame`
  - DataFrame containing binomial test results.
- `correction : str, default='bonferroni'`
  - P-value correction method.
  - Options: `none`, `bonferroni`, `sidak`, `holm-sidak`, `holm`, `simes-hochberg`, `hommel`, `fdr_bh`, `fdr_by`, `fdr_tsbh`, `fdr_tsbky`.
- `cutoff : float, default=0.05`
  - Significance threshold for adjusted binomial p-values.

**Returns:**

- `pd.DataFrame`
  - DataFrame with adjusted binomial test results.

------

#### `comorbidity_network`

```python
dnt.comorbidity_network(
    data:DiseaseNetworkData,
    comorbidity_strength_result:pd.DataFrame, 
    binomial_test_result:pd.DataFrame=None, 
    method:str='RPCN', 
    covariates:list=None, 
    n_process:int=1, 
    log_file:str=None, 
    correction:str='bonferroni', 
    cutoff:float=0.05, 
    **kwargs
) -> pd.DataFrame
```

Performs comorbidity network analysis on disease pairs with significant comorbidity strength.

**Parameters:**

- `data : DiseaseNetworkData`
  - The DiseaseNetworkData object containing your study data.
- `comorbidity_strength_result : pd.DataFrame`
  - DataFrame containing comorbidity strength analysis results.
- `binomial_test_result : pd.DataFrame`
  - DataFrame containing binomial test analysis results.
- `method : str, default='RPCN'`
  - Specifies the comorbidity network analysis method.
  - Options:
    - `'RPCN'`: Regularized Partial Correlation Network.
    - `'PCN_PCA'`: Partial Correlation Network with PCA.
    - `'CN'`: Correlation Network.
  - Method Details:
    - RPCN:
      - Utilizes L1-regularized logistic regression to estimate partial correlations.
      - Includes phenotypic variables and other diseases as covariates.
      - Supports automatic penalty determination.
    - PCN_PCA:
      - Applies logistic regression with principal components of other diseases.
    - CN:
      - Uses standard logistic regression adjusting only for phenotypic variables.
- `covariates : list, default=None`
  - List of phenotypic covariates to include.
  - Defaults to `['sex']` and all covariates specified in `phenotype_data()`.
- `n_process : int, default=1`
  - Number of parallel processes to use.
- `log_file : str, default=None`
  - Path and prefix for the log file.
- `correction : str, default='bonferroni'`
  - P-value correction method.
  - Options: `none`, `bonferroni`, `sidak`, `holm-sidak`, `holm`, `simes-hochberg`, `hommel`, `fdr_bh`, `fdr_by`, `fdr_tsbh`, `fdr_tsbky`.
- `cutoff : float, default=0.05`
  - Significance threshold for adjusted comorbidity network analysis p-values.
- `**kwargs`
  - Additional parameters based on method:
    - RPCN:
      - `alpha : float`
      - `auto_penalty : bool, default=True`
      - `alpha_range : tuple, default=(1,15)`
      - `scaling_factor : float, default=1`
    - PCN_PCA:
      - `n_PC : int, default=5`
      - `explained_variance : float`
  - Additional keyword arguments to define required columns:
    - `phecode_d1_col : str, default='phecode_d1'`
    - `phecode_d2_col : str, default='phecode_d2'`
    - `significance_phi_col : str, default='phi_p_significance'`
    - `significance_RR_col : str, default='RR_p_significance'`
    - `significance_binomial_col : str, default='binomial_p_significance'`

**Returns:**

- `pd.DataFrame`
  - DataFrame containing comorbidity network analysis results.

------

#### `comorbidity_multipletests`

```python
dnt.comorbidity_multipletests(
    df:pd.DataFrame, 
    correction:str='bonferroni', 
    cutoff:float=0.05
) -> pd.DataFrame:
```

Adjusts comorbidity network analysis p-values for multiple comparisons.

**Parameters:**

- `df : pd.DataFrame`
  - DataFrame containing comorbidity network analysis results.
- `correction : str, default='bonferroni'`
  - P-value correction method.
  - Options: `none`, `bonferroni`, `sidak`, `holm-sidak`, `holm`, `simes-hochberg`, `hommel`, `fdr_bh`, `fdr_by`, `fdr_tsbh`, `fdr_tsbky`.
- `cutoff : float, default=0.05`
  - Significance threshold for adjusted p-values.

**Returns:**

- `pd.DataFrame`
  - DataFrame with adjusted comorbidity network analysis results.

------

#### `disease_trajectory`

```python
dnt.disease_trajectory(
    data:DiseaseNetworkData, 
    comorbidity_strength_result:pd.DataFrame, 
    binomial_test_result:pd.DataFrame, 
    method:str='RPCN', 
    matching_var_dict:dict={'sex':'exact'}, 
    matching_n:int=2, 
    max_n_cases:int=np.inf,
    global_sampling: bool=False, 
    covariates:list=None, 
    n_process:int=1, 
    log_file:str=None, 
    correction:str='bonferroni', 
    cutoff:float=0.05, 
    **kwargs
) -> pd.DataFrame
```

Performs temporal comorbidity network (disease trajectory) analysis to identify pairs with confirmed temporal comorbidity associations.

**Parameters:**

- `data : DiseaseNetworkData`

  - The DiseaseNetworkData object containing your study data.

- `comorbidity_strength_result : pd.DataFrame`

  - DataFrame containing comorbidity strength analysis results.

- `binomial_test_result : pd.DataFrame`

  - DataFrame containing binomial test analysis results.

- `method : str, default='RPCN'`

  - Specifies the trajectory analysis method.
  - Options:
    - `'RPCN'`: Regularized Partial Correlation Network.
    - `'PCN_PCA'`: Partial Correlation Network with PCA.
    - `'CN'`: Correlation Network.

- `matching_var_dict : dict, default={'sex': 'exact'}`

  - Specifies matching variables and criteria for incidence density sampling.

  - Example:

    ```python
    matching_var_dict = {'age': 2, 'sex': 'exact'}
    ```

  - Criteria:

    - Categorical/Binary: `'exact'`
    - Continuous: Scalar indicating maximum allowed difference.

- `matching_n : int, default=2`

  - Maximum number of matched controls per case.

- `max_n_cases : int, default=np.inf`
  
  - Specifies the maximum number of D2 cases to include in the analysis.
  - If the number of D2 cases exceeds this value, a random sample of cases will be selected.

- `global_sampling : bool, default=False`

  - Indicates whether to perform independent incidence density sampling for each D1→D2 pair (if False),
    or to perform a single incidence density sampling for all Dx→D2 pairs with separate regression models for each D1→D2 pair (if True).
  - Global sampling is recommended when processing large datasets, though it might reduce result heterogeneity.

- `covariates : list, default=None`

  - List of phenotypic covariates to include.
  - Defaults to all covariates specified in `phenotype_data()`.

- `n_process : int, default=1`

  - Number of parallel processes to use.

- `log_file : str, default=None`

  - Path and prefix for the log file.

- `correction : str, default='bonferroni'`

  - P-value correction method.
  - Options: `none`, `bonferroni`, `sidak`, `holm-sidak`, `holm`, `simes-hochberg`, `hommel`, `fdr_bh`, `fdr_by`, `fdr_tsbh`, `fdr_tsbky`.

- `cutoff : float, default=0.05`

  - Significance threshold for adjusted trajectory analysis p-values.

- `**kwargs`

  - Additional parameters:
    - `enforce_time_interval : bool, default=True`
    - Column names:
      - `phecode_d1_col : str, default='phecode_d1'`
      - `phecode_d2_col : str, default='phecode_d2'`
      - `significance_phi_col : str, default='phi_p_significance'`
      - `significance_RR_col : str, default='RR_p_significance'`
      - `significance_binomial_col : str, default='binomial_p_significance'`
    - Method-Specific Parameters:
      - RPCN:
        - `alpha : float`
        - `auto_penalty : bool, default=True`
        - `alpha_range : tuple, default=(1,15)`
        - `scaling_factor : float, default=1`
      - PCN_PCA:
        - `n_PC : int, default=5`
        - `explained_variance : float`

**Returns:**

- `pd.DataFrame`
  - DataFrame containing trajectory analysis results.

------

#### `trajectory_multipletests`

```python
dnt.trajectory_multipletests(
    df:pd.DataFrame, 
    correction:str='bonferroni', 
    cutoff:float=0.05
) -> pd.DataFrame
```

Adjusts trajectory analysis p-values for multiple comparisons.

**Parameters:**

- `df : pd.DataFrame`
  - DataFrame containing trajectory analysis results.
- `correction : str, default='bonferroni'`
  - P-value correction method.
  - Options: `none`, `bonferroni`, `sidak`, `holm-sidak`, `holm`, `simes-hochberg`, `hommel`, `fdr_bh`, `fdr_by`, `fdr_tsbh`, `fdr_tsbky`.
- `cutoff : float, default=0.05`
  - Significance threshold for adjusted p-values.

**Returns:**

- `pd.DataFrame`
  - DataFrame with adjusted trajectory analysis results.

------

#### `phewas_plot`

```python
def phewas_plot(
    self,
    path: str,
    col_coef: Optional[str]="phewas_coef",
    col_system: Optional[str]="system",
    col_se: Optional[str]="phewas_se",
    col_disease: Optional[str]="disease",
    is_exposure_only: Optional[bool]=False,
    col_exposure: Optional[str]='N_cases_exposed'
) -> None:
```

Generates a circular PheWAS (Phenome-Wide Association Study) plot. Creates a polar bar plot visualizing disease associations across different disease categories (systems) with effect size visualization.

**Parameters:**

  - `path : str`
    - Output file path for saving the plot (PNG format)

  - `col_coef : str, default="phewas_coef"`
    - Column name containing effect size coefficients (log hazard ratios)

  - `col_system : str, default="system"`
    - Column name containing disease system/category classifications

  - `col_se : str, default="phewas_se"`
    - Column name containing standard errors of effect estimates

  - `col_disease : str, default="disease"`
    - Column name containing disease names/labels

  - `is_exposure_only : bool, default=False`
    - Whether to filter for exposure-only cohort results

  - `col_exposure : str, default="N_cases_exposed"`
    - Column name containing case counts for exposed individuals

**Plot Features:**
- Polar coordinate visualization with:
  - Outer ring showing individual disease associations
  - Inner segments grouping by disease system
  - Color gradient indicating effect size (red=positive, green=negative)
  - Automatic text rotation for optimal label readability
- System-level estimates calculated using random effects models
- High-resolution output (1200 DPI)

**Example:**
```python
network.phewas_plot(
    "results/phewas_plot.png",
    col_coef="beta",
    col_system="category",
    col_disease="condition"
)
```

**Notes:**
- Positive associations shown in red (hazard ratio > 1)
- Negative associations shown in green (hazard ratio < 1)
- Plot includes:
  - Color bar legend for effect sizes
  - System category labels
  - Individual disease labels
- Recommended minimum 8"×8" size for readability

**Returns:**

- `None`  
  The method generates the plot but returns nothing  

------

#### `comorbidity_network_plot`

```python
def comorbidity_network_plot(
    self, 
    path :str,
    max_radius: Optional[float]=180.0,
    min_radius: Optional[float]=35.0,
    size_reduction: Optional[float]=0.5,
    cluster_reduction_ratio: Optional[float]=0.4,
    cluster_weight: Optional[str]="comorbidity_beta",
    line_width: Optional[float]=1.0,
    line_color: Optional[str]="black",
    layer_distance: Optional[float]=40.0,
    font_style: Optional[str]="Times New Roman"
) -> None:
```

Generates a 2D visualization of the comorbidity network. Creates an interactive plot showing disease comorbidities with system-based clustering.

**Parameters:**

  - `path : str`
    - Output file path for saving HTML visualization (must end with .html)

  - `max_radius : float, default=90.0`
    - Maximum radial position for nodes (arbitrary units)
    
  - `min_radius : float, default=35.0`
    - Minimum radial position for nodes (arbitrary units)

  - `size_reduction : float, default=0.5`
    - Scaling factor for node sizes (0-1 range recommended)

  - `cluster_reduction_ratio : float, default=0.4`
    - Compression factor for cluster layout (higher values = more compact)

  - `cluster_weight : str, default="comorbidity_beta"`
    - Edge weight metric used for clustering. Options:
      - `"comorbidity_beta"`: Beta coefficients from comorbidity analysis
      - `"phi_coefficient"`: Phi correlation coefficients
      - `"RR"`: Risk ratios

  - `line_width : float, default=1.0`
    - Width of comorbidity connection lines (pixels)

  - `line_color : str, default="black"`
    - Color of comorbidity lines (supports CSS color names/hex codes)

  - `layer_distance : float, default=40.0`
    - Distance between concentric circles for systems (arbitrary units)

  - `font_style : str, default="Times New Roman"`
    - Font family for all text elements

**Example:**
```python
network.comorbidity_network_plot(
    "results/comorbidity_network.html",
    max_radius=120,
    line_color="#555555",
    cluster_weight="RR",
    font_style="Arial"
)
```

**Returns:**

  - `None`  
    The method generates the plot but returns nothing  

------

#### `trajectory_plot`

```python
def trajectory_plot(
    self, 
    path: str,
    cluster_weight: Optional[str]="comorbidity_beta",
) -> None:
```

Generates hierarchical trajectory visualizations for disease clusters. Creates 2D network plots showing temporal disease relationships within each identified cluster, with nodes positioned based on trajectory patterns.

**Parameters:**

- `path : str`
  - Directory path where cluster visualizations will be saved
  - Will create directory if it doesn't exist
  - Output files will be named `cluster_<number>.png`

- `cluster_weight : str, default="comorbidity_beta"`
  - Metric used for determining cluster relationships. Options:
    - `"comorbidity_beta"`: Beta coefficients from comorbidity analysis
    - `"phi_coefficient"`: Phi correlation values  
    - `"RR"`: Risk ratios
    - `"OR"`: Odds ratios

**Visualization Features:**
- Hierarchical node layout showing temporal progression
- Nodes:
  - Colored by disease system category
  - Sized by disease significance (p-value)
  - Exposure disease (if specified) highlighted in gray
- Edges:
  - Width proportional to trajectory strength
  - Direction indicates temporal sequence (source→target)
- Per-cluster images with consistent scaling

**Workflow:**
1. Cluster Analysis:
    - Performs community detection if not already completed
    - Uses specified `cluster_weight` metric
2. Trajectory Processing:
    - Filters for statistically significant trajectories
    - Calculates edge weights
3. Visualization Generation:
    - Creates separate plot for each cluster
    - Computes hierarchical node positions
    - Renders nodes and edges with matplotlib
4. Output:
    - Saves PNG images (600 DPI) to specified directory
    - Names files sequentially (cluster_1.png, cluster_2.png, etc.)

**Example:**
```python
# Basic usage
network.trajectory_plot("results/trajectory_plots")

# Custom edge weighting
network.trajectory_plot(
    "output/figures",
    cluster_weight="RR"
)
```

**Notes:**
- Requires pre-computed trajectory analysis results
- Recommended minimum cluster size: 5 diseases
- For publication-quality figures:
  - Use vector formats (PDF/SVG) by changing file extension
  - Adjust figure size via matplotlib rcParams
- Color scheme matches system categories from comorbidity plot
- Exposure node (if present) uses special marker style

**Returns:**

  - `None`  
    The method generates the plot but returns nothing  

------

#### `three_dimension_plot`

```python
def three_dimension_plot(
    self, 
    path: str,
    max_radius: Optional[float]=180.0, 
    min_radius: Optional[float]=35.0,
    line_color: Optional[str]="black", 
    line_width: Optional[float]=1.0,
    size_reduction: Optional[float]=0.5,
    cluster_reduction_ratio: Optional[float]=0.4,
    cluster_weight: str="comorbidity_beta",
    layer_distance: Optional[float]=40.0,
    layout_width: Optional[float]=900.0,
    layout_height: Optional[float]=900.0,
    font_style: Optional[str]='Times New Roman',
    font_size: Optional[float]=15.0,
) -> None:
```
Generates an interactive 3D visualization of disease trajectories and comorbidities.

  Creates a spherical plot showing temporal disease relationships with:
  - Hierarchical clustering in 3D space
  - Dynamic node positioning based on trajectory patterns
  - System-based color coding

**Parameters:**

  - `path : str`
    - Output file path for HTML visualization (must end with .html)
    - Example: "results/3d_network.html"

  - `max_radius : float, default=180.0`
    - Maximum distance from center for node placement (arbitrary units)

  - `min_radius : float, default=35.0`
    - Minimum distance from center for node placement (arbitrary units)

  - `line_color : str, default="black"`
    - Color for trajectory paths (supports hex/RGB/color names)

  - `line_width : float, default=1.0`
    - Width of trajectory lines (pixels)

  - `size_reduction : float, default=0.5`
    - Node size scaling factor (0.1-1.0 recommended)

  - `cluster_reduction_ratio : float, default=0.4`
    - Cluster compactness factor (higher = more dense)

  - `cluster_weight : str, default="comorbidity_beta"`
    - Edge weighting metric for clustering. Options:
      - `"comorbidity_beta"`: Regression coefficients
      - `"phi_coefficient"`: Phi correlations
      - `"RR"`: Risk ratios
      - `"OR"`: Odds ratios

  - `layer_distance : float, default=40.0`
    - Vertical spacing between disease systems (arbitrary units)

  - `layout_width : float, default=900.0`
    - Figure width in pixels (minimum 600 recommended)

  - `layout_height : float, default=900.0`
    - Figure height in pixels (minimum 600 recommended)

  - `font_style : str, default='Times New Roman'`
    - Font family for all text elements

  - `font_size : float, default=15.0`
    - Base font size in points (10-20 recommended)

**Visualization Features:**
  - 3D spherical coordinate system:
    - Nodes positioned by disease system clusters
    - Vertical dimension represents temporal progression
  - Interactive elements:
    - Zoom/rotate/pan functionality
    - Hover tooltips showing disease details
    - Toggleable legend
  - Visual encoding:
    - Node color = disease system (17 categories)
    - Node size = disease significance
    - Edge width = trajectory strength
    - Special marker for exposure disease (if specified)

**Workflow:**
  1. Data Preparation:
      - Performs cluster analysis if missing
      - Calculates 3D node positions
  2. Visualization:
      - Creates Plotly 3D scatter plot
      - Adds trajectory paths as 3D lines
      - Configures camera and lighting
  3. Output:
      - Generates standalone HTML file
      - Embeds all required JavaScript

**Example:**
  ```python
  # Basic usage
  network.three_dimension_plot("results/3d_network.html")

  # Customized visualization
  network.three_dimension_plot(
      path="output/3d_plot.html",
      max_radius=150,
      line_color="#3498db",
      cluster_weight="RR",
      size_reduction=0.7,
      font_style="Arial"
  )
  ```

**Notes:**
- Requires Plotly (v5.0+) for visualization
- Optimal browser: Chrome/Firefox (WebGL acceleration)
- For publication:
  - Use higher resolution (increase layout dimensions)
  - Consider PNG screenshot for static versions
- Performance scales with node count (~100 nodes recommended)
- Color scheme consistent with 2D visualizations

**Returns:**

  - `None`  
    The method generates the plot but returns nothing  

------

## Issues reporting and recommendations

Please contact:
Can Hou: houcan@wchscu.cn

Haowen Liu: haowenliu81@gmail.com

## License

DiseaseNetPy is released under [The GNU General Public License v3.0](https://www.gnu.org/licenses/gpl-3.0.html)