# Three-dimensional disease network analysis using DiseaseNetPy
## Table of contents

- [1. Input data preparation](#1-input-data-preparation)  
  - [1.1 Requirements for input data](#11-requirements-for-input-data)  
  - [1.2 Dummy dataset overview](#12-dummy-dataset-overview)  
- [2. Data Harmonization](#2-data-harmonization)  
  - [2.1 Initializing the data object](#21-initializing-the-data-object)  
  - [2.2 Load phenotype data](#22-load-phenotype-data)  
  - [2.3 Load medical records data](#23-load-medical-records-data)  
- [3. Data Analysis](#3-data-analysis)  
  - [3.1 Pipeline](#31-pipeline)  
  - [3.2 PheWAS Analysis](#32-phewas-analysis)  
  - [3.3 Disease pair construction](#33-disease-pair-construction)  
  - [3.4 Comorbidity strength estimation](#34-comorbidity-strength-estimation)  
  - [3.5 Binomial test](#35-binomial-test)  
  - [3.6 Comorbidity network analysis](#36-comorbidity-network-analysis)  
  - [3.7 Trajectory analysis](#37-trajectory-analysis)  
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
    - [concat](#concat)  
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

DiseaseNetPy enables 3D disease network analysis on cohort data from electronic health records (EHR) and offers three study designs: the standard cohort, which compares individuals with a specific disease (e.g., depression) or exposure (e.g., smoking) against the general population; the matched cohort, which pairs subjects on key characteristics to reduce bias; and the exposed-only cohort, which examines disease networks within a defined subgroup (e.g., older adults) without a comparison group.

To begin using DiseaseNetPy, two datasets are required: a **phenotype data** file recording each individual's basic and additional characteristics, and a **medical records data** file extracted from an EHR database that stores diagnosis codes and dates for all cohort individuals over the study period. Specific requirements for these datasets are as follows:

- **Phenotype data**: An on-disk CSV (or TSV) file with headers, listing each participant with the following required and optional columns:

  - **Participant ID**: Unique identifier for each individual.  
  - **Index date**: Start of follow-up (e.g., date of exposure or baseline).  
  - **End date**: End of follow-up (e.g., last visit, death, or study completion).  
  - **Exposure**: Binary indicator (1 = exposed, 0 = unexposed) for standard and matched cohort designs (omit for exposed-only cohorts).  
  - **Match ID**: Identifier for matched sets (only for matched cohort designs).  
  - **Sex**: Biological sex (1 = female, 0 = male).  
  - **Additional covariates (optional)**: Any number of extra variables for adjustment or matching (e.g., age, BMI, education).  

  For all required columns, missing values are not permitted and dates must follow the datetime formats (e.g., “%Y-%m-%d”); the **Sex** and **Exposure** fields must use the specified 1/0 coding, and you may include an unlimited number of additional covariates. Column types—binary, categorical, or continuous—are detected automatically and transformed as needed (e.g., one-hot encoding), with missing categorical values assigned to a separate “NA” category.

- **Medical records data**: One or more CSV or TSV files (with a header row), each listing diagnosis events for study participants. Every record must include:  

  - **Participant ID**: The unique identifier matching the phenotype dataset.  
  - **Diagnosis code**: A standardized code (e.g., ICD-10, ICD-9).  
  - **Date of diagnosis**: Date of diagnosis or event in datetime formats (e.g., “%Y-%m-%d”).
  
  Unlike the **Phenotype data**, the **medical records data** should be in a record-per-row format, allowing multiple rows per participant. Records matching **participant ID** in the **phenotype data** and occurring within the specified follow-up period (i.e., before the end date) are loaded. Therefore, there is **no need** to pre-filter medical records - you should provide complete diagnosis information for all cohort individuals. Specifically, it is **not recommended** to filter data based on first occurrence of each code or diagnosis date.

  Each **medical records data** must use a single diagnosis code version; currently supported versions are WHO or CM versions of ICD-9 and ICD-10. Other code systems must be converted to a supported format.

### 1.2 Dummy dataset overview

A dummy dataset is provided to help you become familiar with the required input format and to run through the full analysis workflow before applying it to your own cohort. It simulates a matched‐cohort study of 10 000 exposed individuals and 50 000 matched unexposed individuals, together with their entire follow-up EHR data.

**Caution:** All participant characteristics and diagnosis records in this dataset are randomly generated. Although the ICD-9 and ICD-10 codes correspond to real‐world classifications, and the analysis may produce apparently significant associations based on that, these results do **not** reflect any true medical findings.

- The dataset consists of three CSV files, located in the `tests/data` directory:
  - **`dummy_phenotype.csv`**
     Baseline characteristics for all 60 000 individuals, containing:
    - **ID**: unique participant identifier
    - **date_start**, **date_end**: follow-up start and end dates
    - **exposure**: exposure status (0 = unexposed, 1 = exposed)
    - **group_id**: matching group identifier
    - **sex**: biological sex (1 for female and 0 for male)
    - **age**: baseline age (years)
    - **BMI**: body‐mass index category
  - **`dummy_EHR_ICD9.csv`**
     Simulated EHR diagnoses coded using ICD-9 (n = 10,188 records), with columns:
    - **ID**: participant identifier
    - **dia_date**: diagnosis date
    - **diag_icd9**: ICD-9 diagnosis code
  - **`dummy_EHR_ICD10.csv`**
     Simulated EHR diagnoses coded using ICD-10 (n = 1,048,576 records), with columns:
    - **ID**: participant identifier
    - **dia_date**: diagnosis date
    - **diag_icd10**: ICD-10 diagnosis code

## 2. Data Harmonization

Data harmonization loads and merges phenotype and medical records data into a single `DiseaseNetworkData` object for subsequent analysis, ensuring consistent coding (e.g., mapping diagnosis codes to phecodes) and standardized formatting.

### 2.1 Initializing the data object

First, import `DiseaseNetworkData` from DiseaseNetPy and instantiate it with your chosen study design, phecode level, and any optional parameters to create an data object.

```python
import diseasenetpy as dnt

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

- **study_design** – the cohort design: `'cohort'`, `'matched cohort'`, or `'exposed-only cohort'`. Default is `'cohort'`.
- **phecode_level** – the level of phecode used for analysis; level 1 provides broader categories (~585 conditions), while level 2 offers more detail (~1257 conditions). Level 1 is recommended for smaller datasets to maintain statistical power. The level of phecode: `1`, or `2`. Default is `1`.

#### Optional parameters:

- **min_required_icd_codes** – the minimum number of ICD codes mapping to a phecode required for it to be considered valid. Default is `1`. For example, setting this to 2 requires at least two records mapping to phecode 250.2 (Type 2 diabetes) for a participant to be considered diagnosed. Ensure your medical records include complete data (not limited to first occurrences) when using this parameter.
- **date_fmt** – format of date fields (Index date and End date) in the phenotype data. Default is `'%Y-%m-%d'` (year-month-day, e.g., 2005-12-01).
- **phecode_version** – currently only `1.2` is supported. Default is `1.2`.

### 2.2 Load phenotype data

After initializing the data object, use the `phenotype_data()` method to load your one cohort phenotype file by providing the file path, a dictionary mapping required columns, and a list of additional covariate names.

The following example codes show how to load the dummy phenotype dataset under different study designs. Although the file (**dummy_phenotype.csv**) is originally formatted for a matched cohort study, you can adapt it for other designs: omitting the **Match ID** column loads it as a standard cohort study (ignoring the matching), while omitting both **Match ID** and **Exposure** columns loads it as an exposed‑only cohort study - treating all participants as exposed (i.e., representing the entire population).

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
vars_lst = ['age', 'BMI']  
data.phenotype_data(
    phenotype_data_path=r"/test/data/dummy_phenotype.csv",
    column_names=col_dict,
    covariates=vars_lst
)

# Load phenotype data for a traditional cohort study
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

# Load phenotype data for an exposed-only cohort study
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

- **phenotype_data_path** – path to your one **Phenotype data** file (CSV or TSV).
- **column_names** – dictionary mapping the required variable names (e.g., `Participant ID`, `Index date`, `End date`, `Sex`) to the corresponding headers in your file. Include `Exposure` for cohort and matched cohort designs, and `Match ID` for matched cohort designs.
- **covariates** – list of additional covariate names. Provide an empty list if none. The method automatically detects and converts variable types. Records with missing values in continuous variables are removed, while missing values in categorical variables form an NA category.

#### Optional parameters:

- **is_single_sex** – set to `True` if the dataset contains only one sex. Default is `False`.
- **force** – set to `True` to overwrite existing data in the object. Default is `False`, which raises an error if data already exist.

#### After data loading:

After loading data, you can inspect basic information (e.g., number of individuals, average follow-up time) by printing the data object:

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

Additionally, you can generate a basic descriptive table 1 (`pd.DataFrame`) for all variables in your phenotype data using the `Table1()` method and save it to multiple formats (e.g., .csv/.tsv/.xlsx):

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
# For example: save Table 1 to an Excel file
table_1.to_excel(r"/test/data/Table1.xlsx")  
```

### 2.3 Load medical records data

After loading the phenotypic data, use the `merge_medical_records()` method to load your one or more medical records files by providing the file path, a format of ICD code, and mapping required columns, and a dictionary mapping required columns. The following example codes show how to load the dummy EHR ICD9/ICD10 dataset.

```python
# Merge with the first medical records file (dummy_EHR_ICD10.csv)
data.merge_medical_records(
    medical_records_data_path=r"/test/data/dummy_EHR_ICD10.csv", 
    diagnosis_code='ICD-10-WHO',                                  
    column_names={
        'Participant ID': 'ID',                                   
        'Diagnosis code': 'diag_icd10',                           
        'Date of diagnosis': 'dia_date'        
    }
)

# Merge with the second medical records file (dummy_EHR_ICD9.csv)
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

- **medical_records_data_path** – path to your one **medical records data** file (CSV or TSV).
- **diagnosis_code** – Diagnosis ICD code type used in the medical records data (e.g., `ICD-9-CM`, `ICD-9-WHO`, `ICD-10-CM`, `ICD-10-WHO`).
- **column_names** – dictionary mapping the required variable names (e.g., `Participant ID`, `Diagnosis code`, `Date of diagnosis`) to the corresponding headers in your file.

#### Optional parameters:

- **date_fmt** – The format of the date fields in your medical records data. Defaults to the same format as phenotype data if not specified.
- **chunksize** – Number of rows per chunk to read, useful for large datasets. Default=1_000_000.

#### During data loading:

During loading data, you can inspect basic information (e.g., number of records, number of phecode mapping) by the printing information:

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

#### After data loading:

After data loading, you can inspect basic information of **medical records data** (e.g., number of ICD code mapped to phecodes, average number of disease diagnosis) by printing the object:

```python
print(data)
# This will output something like (e.g., for a matched cohort study):
"""
Merged Medical records
2 medical records data with 1,678,983 diagnosis records were merged (0 with missing values).
Average number of disease diagnosis during follow-up: 18.99 (exposed) and 7.31 (unexposed)
Average number of disease diagnosis before follow-up: 8.40 (exposed) and 3.46 (unexposed)

Warning: 102 exposed individuals and 440 unexposed individuals have negative or zero follow-up time.
Consider removing them before merge.
Warning: 18.15% of ICD-10-WHO codes were not mapped to phecodes for file d:\GitHubWarehouse\test\src/data/dummy_EHR_ICD10.csv.
Warning: 2.07% of ICD-9-WHO codes were not mapped to phecodes for file d:\GitHubWarehouse\test\src/data/dummy_EHR_ICD9.csv.
"""
```

## 3. Data Analysis

Data analysis is based on a `DiseaseNetworkData` object to subsequent analysis, including PheWAS analysis, disease pair generation, comorbidity strength estimation, binomial testing, comorbidity network analysis, and disease trajectory analysis. The result format of each analysis is `pd.DataFrame`.

### 3.1 Pipeline

In this section, we can use the `disease_network_pipeline` function to complete the entire disease network analysis workflow. Compared to a step-by-step analysis, the `disease_network_pipeline` function requires less code and avoids repetitive operations. However, it does not allow for precise control over the parameter details of each analysis step. The following code applies to all research designs.

```python
# Reminder:
# When using multiprocessing, ensure that the code is enclosed within the following block.
# This prevents entering a never ending loop of new process creation.
from diseasenetpy import disease_network_pipeline

if __name__ == "__main__":
    zipped_result = disease_network_pipeline(
      data=data,                               
      n_process=2,                            
      n_threshold_phewas=100,                  
      n_threshold_comorbidity=100,             
      output_dir="/your/project/path",          
      project_prefix="disease network",         
      keep_positive_associations=False,
      save_intermediate_data=False,
      system_exl=[
        'symptoms', 
        'others', 
        'injuries & poisonings'
      ],                                      
      pipeline_mode="v1",                     
      method="RPCN",                        
      covariates=[
        'BMI', 
        'age'
      ],                                     
      matching_var_dict={
        'sex':'exact'
      },                                  
      matching_n=2,                          
      min_interval_days=0,                    
      max_interval_days=np.inf,               
      enforce_temporal_order=False,           
      correction='bonferroni',                
      cutoff=0.05                      
    )
```

- **data** – The `DiseaseNetworkData` object.
- **n_process** – Specifies the number of parallel processes to use for the disease network analysis. Multiprocessing is enabled when `n_process` is set to a value greater than one.
- **n_threshold_phewas** – The minimum number of cases within the exposed group required for a phecode to be included in the PheWAS analysis. See the phewas function for more information.
- **n_threshold_comorbidity** – The minimum number of individuals in the exposed group in which a disease pair must co-occur (temporal or non-temporal) to be included in the comorbidity strength estimation. See the comorbidity_strength function for more information.
- **output_dir** - Directory path to store output files generated by the pipeline.
- **project_prefix** – Prefix for naming output files and intermediate data.
- **keep_positive_associations** – set to `True` if retains only diseases with hazard ratio (HR) > 1 from the PheWAS analysis. Default is `False`.
- **save_intermediate_data** – set to `True` to intermediate DiseaseNetworkData objects created by the `DiseaseNetPy.DiseaseNetworkData.disease_pair` function are saved to disk. Default is `False`.
- **system_exl** – List of phecode systems to exclude from the analysis. Default is `None`.
- **pipeline_mode** – Specifies the analysis order. 
  - Available modes: 
    - **v1**: PheWAS → comorbidity strength → binomial test → (comorbidity network analysis/disease trajectory analysis) 
    - **v2**: PheWAS → comorbidity strength → comorbidity network analysis → binomial test → disease trajectory analysis. 
  - In **v1**, the binomial test does not depend on results from the comorbidity network analysis; thus, disease trajectory and comorbidity network analyses can be conducted independently. In **v2**, the binomial test is performed only on disease pairs identified as significant by the comorbidity network analysis, making the disease trajectory analysis dependent on these results.
- **method** – The method to use for the comorbidity network and disease trajectory analysis. 
  - Available methods are: 
    - **RPCN**: Regularized Partial Correlation Network. There are four `**Kwarg` parameters.
      - **alpha** – The weight multiplying the l1 penalty term for other diseases covariates. Ignored if 'auto_penalty' is enabled.
      - **auto_penalty** – If 'True', automatically determines the best 'alpha' based on model AIC value. Default is `True`.
      - **alpha_range** – When 'auto_penalty' is True, search the optimal 'alpha' in this range. Default is `(1,15)`.
      - **scaling_factor** – The scaling factor for the alpha when 'auto_penalty' is True. Default is `1`.
    - **PCN_PCA**: Partial Correlation Network with PCA. There are two `**Kwarg` parameters.
      - **n_PC** – Fixed number of principal components to include in each model. Default is `5`.
      - **explained_variance** – Cumulative explained variance threshold to determine the number of principal components. Overrides 'n_PC' if specified.
    - **CN**: Correlation Network. This parameter will be passed to the `comorbidity_network` and `disease_trajectory` function. See the these two functions for more information.
  - The correlation network (CN) represents the most straightforward approach, utilizing a logistic regression model with user-defined covariates to assess pairwise disease associations. Building upon this foundation, the regularized partial correlation network (RPCN) (Epskamp et al., 2016) employs a regularized logistic regression framework that integrates both user-specified covariates and historical disease data. Notably, we further develop the RPCNPCA variant by incorporating principal component analysis (PCA) for covariate dimensionality reduction prior to network construction.
- **covariates** – List of covariates to adjust for in the PheWAS, comorbidity network and disease trajectory analysis. Default is `None`.
- **matching_var_dict** – Specifies the matching variables and the criteria used for incidence density sampling. Default is `{'sex':'exact'}`.
- **matching_n** – Specifies the maximum number of matched controls for each case. This parameter will be passed to the disease_trajectory function. Default is `2`.
- **min_interval_days** – Minimum required time interval (in days) between diagnosis dates when constructing temporal D1 → D2 disease pair for each individual. This parameter will be passed to the DiseaseNetPy.DiseaseNetworkData.disease_pair function. See the disease_pair function for more information. Default is `0`.
- **max_interval_days** – Maximum allowed time interval (in days) between diagnosis dates when constructing temporal and non-temporal D1-D2 disease pair for each individual. This parameter will be passed to the DiseaseNetPy.DiseaseNetworkData.disease_pair function. See the disease_pair function for more information. Default is `np.inf`.
- **enforce_temporal_order** – set to `True` to exclude individuals with non-temporal D1-D2 pair when performing the binomial test; also applies the specified minimum and maximum time intervals when performing disease trajectory analysis. See 'enforce_temporal_order' parameter in binomial_test function and 'enforce_time_interval' parameter in disease_trajectory function. Default is `False`.
- **correction** – Method for p-value correction from the `statsmodels.stats.multitest.multipletests`.          
    Available methods are: 
      - none : no correction
      - bonferroni : one-step correction
      - sidak : one-step correction
      - holm-sidak : step down method using Sidak adjustments
      - holm : step-down method using Bonferroni adjustments
      - simes-hochberg : step-up method (independent)
      - hommel : closed method based on Simes tests (non-negative)
      - fdr_bh : Benjamini/Hochberg (non-negative)
      - fdr_by : Benjamini/Yekutieli (negative)
      - fdr_tsbh : two stage fdr correction (non-negative)
      - fdr_tsbky : two stage fdr correction (non-negative)
    See https://www.statsmodels.org/dev/generated/statsmodels.stats.multitest.multipletests.html for more details. Default is `bonferroni`.
- **cutoff** – The significance threshold for adjusted p-values. Default is `0.05`.

#### After pipeline:

Following the pipeline execution, we obtain a zipped data structure containing results from five analysis steps: PheWAS analysis, comorbidity strength estimation, binomial test, comorbidity network analysis, and disease trajectory analysis. These results can be unpacked into individual components and displayed. There are some descriptions of each variable of each result.

```python
# Firstly, unpack the results
phewas_result = zipped_result[0]
com_strength_result = zipped_result[1]
com_network_result = zipped_result[2]
binomial_result = zipped_result[3]
trajectory_result = zipped_result[4]
```

```python
# And then, print the result
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

**Result of PheWAS analysis**

| Variable Name           | Type      | Description |
|------------------------|-----------|-------------|
| `phecode`              | String    | Disease code (Phecode) used in PheWAS |
| `disease`              | String    | Disease name corresponding to the Phecode |
| `system`               | String    | Phecode disease system corresponding to the Phecode (e.g., infectious diseases) |
| `sex`                  | String    | Sex-specific for the disease (e.g., Both, Male, Female) |
| `N_cases_exposed`      | Integer   | Number of cases among the exposed group |
| `describe`             | String    | Description of the analysis or model (e.g., covariates used) |
| `exposed_group`        | String    | Time at risk in the exposed group |
| `unexposed_group`      | String    | Time at risk in the unexposed group |
| `phewas_coef`          | Float     | Estimated coefficient from the model |
| `phewas_se`            | Float     | Standard error of the estimated coefficient |
| `phewas_p`             | Float     | P-value indicating statistical significance of the coefficient |
| `phewas_p_significance`| Boolean   | Indicates whether the result is statistically significant (True/False) |
| `phewas_p_adjusted`    | Float     | Adjusted p-value accounting for multiple comparisons |

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

**Result of comorbidity strength estimation**

| Variable Name            | Type    | Description |
|-------------------------|---------|-------------|
| `phecode_d1`            | Integer | Phecode for disease 1 in the disease pair |
| `phecode_d2`            | Integer | Phecode for disease 2 in the disease pair |
| `name_disease_pair`     | String  | Concatenated identifier of the disease pair (format: d1-d2) |
| `N_exposed`             | Integer | Number of individuals in exposed group (All individuals is esposed group in the exposed-only cohort) |
| `n_total`               | Integer | Total number of individuals in the sub-cohort |
| `n_d1d2_diagnosis`      | Integer | Number of individuals diagnosed with both diseases |
| `n_d1_diagnosis`        | Integer | Number of individuals diagnosed with disease 1 |
| `n_d2_diagnosis`        | Integer | Number of individuals diagnosed with disease 2 |
| `n_d1d2_nontemporal`    | Integer | Co-diagnoses without considering temporal order |
| `phi_coef`              | Float   | Phi coefficient measuring co-occurrence strength |
| `phi_p`                 | Float   | P-value for Phi coefficient significance |
| `RR`                    | Float   | Relative risk estimate for the disease pair |
| `RR_p`                  | Float   | P-value for relative risk |
| `phi_p_adjusted`        | Float   | Adjusted P-value for Phi coefficient (multiple comparisons) |
| `RR_p_adjusted`         | Float   | Adjusted P-value for relative risk (multiple comparisons) |
| `phi_p_significance`    | Boolean | Whether the Phi is statistically significant |
| `RR_p_significance`     | Boolean | Whether the RR is statistically significant |
| `disease_d1`            | String  | Name of disease 1 |
| `system_d1`             | String  | Phecode disease system related to disease 1 |
| `sex_d1`                | String  | Sex-specific for disease 1 |
| `disease_d2`            | String  | Name of disease 2 |
| `system_d2`             | String  | Phecode disease system related to disease 2 |
| `sex_d2`                | String  | Sex-specific for disease 2 |

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

**Result of binomial test**

| Variable Name               | Type    | Description |
|----------------------------|---------|-------------|
| `phecode_d1`               | Float   | Phecode for the first disease in the temporal disease pair |
| `phecode_d2`               | Float   | Phecode for the second disease in the temporal disease pair |
| `name_disease_pair`        | String  | Identifier for the disease pair, showing direction (e.g., d1->d2) |
| `n_d1d2_nontemporal`       | Float   | Count of individuals diagnosed with both diseases (non-temporal) |
| `n_d1d2_temporal`          | Float   | Count of individuals where disease 1 was diagnosed before disease 2 |
| `n_d2d1_temporal`          | Float   | Count of individuals where disease 2 was diagnosed before disease 1 |
| `binomial_p`               | Float   | P-value from the binomial test for directionality (d1 before d2) |
| `binomial_proportion`      | Float   | Proportion of d1 before d2 among all temporal diagnoses |
| `binomial_proportion_ci`   | String  | Confidence interval for the binomial proportion (format: lower-upper) |
| `disease_d1`               | String  | Name of disease 1 |
| `system_d1`                | String  | Phecode disease system for disease 1 |
| `sex_d1`                   | String  | Sex-specific for disease 1 |
| `disease_d2`               | String  | Name of disease 2 |
| `system_d2`                | String  | Phecode disease system for disease 2 |
| `sex_d2`                   | String  | Sex-specific for disease 2 |
| `binomial_p_significance` | Boolean | Indicates whether the result is statistically significant |
| `binomial_p_adjusted`     | Float   | Adjusted p-value for multiple comparisons |

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

**Result of comorbidity network analysis**

| Variable Name                | Type    | Description |
|-----------------------------|---------|-------------|
| `phecode_d1`                | Float   | Phecode for the first disease in the disease pair |
| `phecode_d2`                | Float   | Phecode for the second disease in the disease pair |
| `name_disease_pair`         | String  | Identifier for the disease pair (format: d1-d2) |
| `N_exposed`                 | Integer | Number of individuals in the exposed group |
| `n_total`                   | Integer | Total number of individuals in the sub-cohort |
| `n_exposed/n_cases`         | String  | Number of exposed among cases (format: exposed/cases) |
| `n_exposed/n_controls`      | String  | Number of exposed among controls (format: exposed/controls) |
| `comorbidity_network_method`| String  | Method used for comorbidity network analysis |
| `describe`                  | String  | Description of the model fitting (e.g., covariates removed) |
| `co_vars_list`              | String  | List of covariates used in the model |
| `co_vars_zvalues`           | String  | Z-values for each covariate in the model |
| `comorbidity_beta`          | Float   | Estimated coefficient from the comorbidity model |
| `comorbidity_se`            | Float   | Standard error of the estimated coefficient |
| `comorbidity_p`             | Float   | P-value for the comorbidity coefficient |
| `comorbidity_aic`           | Float   | Akaike Information Criterion for the model |
| `disease_d1`                | String  | Name of the first disease |
| `system_d1`                 | String  | Phecode disease system for the first disease |
| `sex_d1`                    | String  | Sex-specific for the first disease |
| `disease_d2`                | String  | Name of the second disease |
| `system_d2`                 | String  | Phecode disease system for the second disease |
| `sex_d2`                    | String  | Sex-specific for the second disease |
| `comorbidity_p_significance`| Boolean | Whether the result is statistically significant |
| `comorbidity_p_adjusted`    | Float   | Adjusted p-value accounting for multiple comparisons |

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

**Result of disease trajectory analysis**

| Variable Name               | Type    | Description |
|-----------------------------|---------|-------------|
| `phecode_d1`                | Float   | Phecode for the first disease in the disease pair |
| `phecode_d2`                | Float   | Phecode for the second disease in the disease pair |
| `name_disease_pair`         | String  | Identifier for the disease pair (format: d1->d2) |
| `N_exposed`                 | Integer | Number of individuals in the exposed group |
| `n_total`                   | Integer | Total number of individuals in the sub-cohort |
| `n_exposed/n_cases`         | String  | Number of exposed among cases (format: exposed/cases) |
| `n_exposed/n_controls`      | String  | Number of exposed among controls (format: exposed/controls) |
| `trajectory_method`         | String  | Method used for disease trajectory analysis |
| `describe`                  | String  | Description of model fitting and covariates used/removed |
| `co_vars_list`              | String  | List of covariates included in the model |
| `co_vars_zvalues`           | String  | Z-values for each covariate in the model |
| `trajectory_beta`           | Float   | Estimated coefficient from the model |
| `trajectory_se`             | Float   | Standard error of the estimated coefficient |
| `trajectory_p`              | Float   | P-value for the coefficient |
| `trajectory_aic`            | Float   | Akaike Information Criterion for the model |
| `disease_d1`                | String  | Name of the first disease |
| `system_d1`                 | String  | Phecode disease system for the first disease |
| `sex_d1`                    | String  | Sex-specific for the first disease |
| `disease_d2`                | String  | Name of the second disease |
| `system_d2`                 | String  | Phecode disease system for the second disease |
| `sex_d2`                    | String  | Sex-specific for the second disease |
| `trajectory_p_significance` | Boolean | Whether the result is statistically significant |
| `trajectory_p_adjusted`     | Float   | Adjusted p-value accounting for multiple comparisons |

#### Save the results:
All analysis outputs are standardized `pd.DataFrame` objects, enabling flexible data serialization through pandas' native I/O methods. The results can be exported to multiple formats including but not limited to:
- Tabular formats (.csv/.tsv)
- Spreadsheet files (.xlsx/.ods)
- Binary formats (.feather/.parquet)
- Database interfaces (SQL)

```python
# For example: save the results to .csv file
phewas_result.to_csv("/your/project/path/phewas_result.csv")
com_strength_result.to_csv("/your/project/path/com_strength_result.csv")
com_network_result.to_csv("/your/project/path/com_network_result.csv")
binomial_result.to_csv("/your/project/path/binomial_result.csv")
trajectory_result.to_csv("/your/project/path/trajectory_result.csv")
```

### 3.2 PheWAS Analysis

The analysis begins with Phenome-Wide Association Studies (PheWAS) using the `DiseaseNetworkData` object. For cohort/matched-cohort designs, this step assesses correlations between exposure and outcome diseases, filtering those with strong associations. In exposure-only cohorts, it calculates disease incidence rates and applies minimum incidence thresholds.

```python
# Reminder:
# When using multiprocessing, ensure that the code is enclosed within the following block.
# This prevents entering a never ending loop of new process creation.
if __name__ == "__main__":
    phewas_result = dnt.phewas(
        data=data,                                            
        proportion_threshold=0.01,                            
        n_process=2,                                          
        system_exl=[                                          
            'symptoms', 
            'others', 
            'injuries & poisonings', 
            'pregnancy complications'
        ],                                                   
        covariates=[
          'age', 
          'sex',
          'BMI',
        ],                                                    
        correction='bonferroni',                               
        lifelines_disable=True,                               
        log_file='/your/project/path/dep.log'                 
    )
```

- **data** – The `DiseaseNetworkData` object.
- **proportion_threshold** – The minimum proportion of cases within the exposed group required for a phecode to be included in the PheWAS analysis. If the proportion of cases is below this threshold, the phecode is excluded from the analysis. `proportion_threshold` and `n_threshold` are mutually exclusive. Default is `None`.
- **n_process** - Number of parallel processes to use for analysis. Multiprocessing is enabled when set to greater than one. Default is `1`.
- **system_exl** - Phecode systems to exclude from analysis. *Note:* Mutually exclusive with `system_inc`. Available systems: Same as `system_inc`. Default is `None`.
- **covariates** – List of phenotypic covariates to include in the model. Default is `None`.
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
- **lifelines_disable** - Whether to disable lifelines. Lifelines provide more robust fitting but require longer computation time. Default is `False`.
- **log_file** - Path and prefix for log file. If None, logs are written to temporary directory with prefix `DiseaseNet_phewas_`. Default is `None`.

#### Optional parameters:

- **cutoff** - Significance threshold for adjusted PheWAS p-values. Default is `0.05`.
- **system_inc** - Phecode systems to include in analysis. *Note:* Mutually exclusive with `system_exl`. Available systems: `circulatory`, `congenital anomalies`, `dermatologic`, `digestive`, `endocrine/metabolic`, `genitourinary`, `hematopoietic`, `infectious diseases`, `injuries & poisonings`, `mental disorders`, `musculoskeletal`, `neoplasms`, `neurological`, `pregnancy complications`, `respiratory`, `sense organs`, `symptoms`, `others`. Default is `None`.
- **phecode_inc** - Specific phecodes to include in analysis. *Note:* Mutually exclusive with phecode_exl. Default is `None`.
- **phecode_exl** - Specific phecodes to exclude from analysis. *Note:* Mutually exclusive with phecode_inc. Default is `None`.
- **n_threshold** - The minimum number of cases within the exposed group required for a phecode to be included in the PheWAS analysis. If the number of cases is below this threshold, the phecode is excluded. *Note:* `n_threshold` and `proportion_threshold` are mutually exclusive. Default is `None`.

#### After PheWAS analysis:

In cohort/matched-cohort studies, we typically focus on outcome diseases where the exposure disease results in increased incidence (i.e., Hazard ratio (HR) > 1). These selected outcome diseases with HR > 1 are retained for downstream analysis. Additionally, alternative p-value adjustment methods can be implemented to enhance statistical rigor.

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

#### Save the result of PheWAS analysis:

Same to the pipeline, we can save the result to multiple formats including but not limited to:
- Tabular formats (.csv/.tsv)
- Spreadsheet files (.xlsx/.ods)
- Binary formats (.feather/.parquet)
- Database interfaces (SQL)

```python
# For example: save the entire PheWAS results to a CSV file
phewas_result.to_csv('/your/project/path/dep_phewas.csv')  
```

### 3.3 Disease pair construction

Following diseases filtering, we construct pairwise disease combinations through exhaustive permutation of the remaining diseases. These generated disease pairs are subsequently archived within the `DiseaseNetworkData` object using its native `disease_pair() `method, with configurable temporal constraints (e.g., minimum/maximum allowable time intervals between disease onset events) to enforce clinically meaningful temporal relationships.

```python
# Disease pair construction of the `DiseaseNetworkData` object
data.disease_pair(
    phewas_result=phewas_result,               
    min_interval_days=30,                      
    max_interval_days=365.25*5,                
    force=True                              
)
```

- **phewas_result** - `pd.DataFrame` containing PheWAS analysis results produced by the `DiseaseNetPy.phewas` function.
- **min_interval_days** - Minimum required time interval (in days) between diagnosis dates when constructing temporal D1→D2 disease pairs. Individuals with D1 and D2 diagnoses interval ≤ this value are considered to have non-temporal pairs. Default is `0`.
- **max_interval_days** - Maximum allowed time interval (in days) between diagnosis dates when constructing disease pairs. Individuals with interval > this value are excluded from temporal analysis. Default is `np.inf`.
- **force** - If `True`, overwrites existing data attributes. If `False`, raises error when data exists. Default is `False`.

#### Optional Parameters:

- **n_process** - Number of processes for parallel processing. Values >1 enable multiprocessing. Default is `1`.
- **phecode_col** - Column name for phecode identifiers in `phewas_result`. Default is `'phecode'`.  
- **significance_col** - Column name for PheWAS significance values. Default is `'phewas_p_significance'`.

#### After disease pair construction

Following disease pair construction, the `DiseaseNetworkData` object can be serialized into compressed file formats (e.g., .npz for NumPy-based storage or .pkl.gz for gzipped Python object persistence) to enable reproducible experimentation.

```python
# For example: save the updated data object with disease pairs
data.save('/your/project/path/dep_withtra')
```

### 3.4 Comorbidity strength estimation

Following disease pair generation via the `disease_pair()` method, we quantify comorbidity strength metrics (Relative risk, and phi-correlation) for each disease pairs (D₁–D₂) using the `DiseaseNetworkData` object, which serves as input for downstream analysis.

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

- **data** - DiseaseNetworkData object containing the processed disease network data.
- **proportion_threshold** - Minimum proportion of exposed individuals required for disease pair co-occurrence to be included in analysis. Disease pairs below this threshold are excluded. *Note:* Mutually exclusive with `n_threshold`. Default is `None`.
- **n_process** - Number of parallel processes for analysis. Values >1 enable multiprocessing. Default is `1`.
- **log_file** - Path/prefix for log file. If `None`, uses temporary directory with prefix `DiseaseNet_com_strength_`. Default is `None`.

#### Optional Parameters:

- **n_threshold** - Minimum number of exposed individuals required for disease pair co-occurrence to be included in analysis. Disease pairs below this threshold are excluded. *Note:* Mutually exclusive with proportion_threshold. Default is `None`.
- **correction_phi** - P-value correction method for phi-correlation (from statsmodels.stats.multitest). Available methods:  
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
- **cutoff_phi** - Significance threshold for adjusted phi-correlation p-values. Default is `0.05`.
- **correction_RR** - P-value correction method for relative risk (same methods as `correction_phi`). Default is `bonferroni`.
- **cutoff_RR** - Significance threshold for adjusted RR p-values. Default is `0.05`.

#### After Comorbidity Strength Estimation:

Following comorbidity strength estimation, we apply filtering criteria based on relative risk (RR) > 1 and phi-correlation > 0 to ensure statistically significant correlations and identify associations where the exposure may increase disease risk. The filtered disease pairs will be analysed in the next step. The details of result is shows in the pipeline. Additionally, alternative p-value adjustment methods can be implemented to enhance statistical rigor.

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

#### Save the result of comorbidity strength estimation:

Same to the pipeline, we can save the result to multiple formats including but not limited to:
- Tabular formats (.csv/.tsv)
- Spreadsheet files (.xlsx/.ods)
- Binary formats (.feather/.parquet)
- Database interfaces (SQL)

```python
# For example: save the comorbidity strength estimation results to a CSV file
com_strength_result.to_csv('/your/project/path/dep_com_strength.csv')  
```

### 3.5 Binomial test

This phase performs a binomial test to identify non-temporal disease pairs (D₁-D₂ co-occurrence) and temporal disease pairs (D₁→D₂) based on the result of comorbidity strength estimation. Additionally, alternative p-value adjustment methods can be implemented to enhance statistical rigor.

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
- **comorbidity_strength_result** - DataFrame containing comorbidity strength analysis results from `DiseaseNetPy.comorbidity_strength()`.
- **n_process** - Number of parallel processes. Note: Multiprocessing is disabled for this analysis. Default is `1`.
- **enforce_temporal_order** - If `True`, excludes individuals with non-temporal D1-D2 pairs. If `False`, includes all individuals. Default is `False`.
- **log_file** - Path/prefix for log file. If `None`, uses temporary directory with prefix `DiseaseNet_binomial_test_`. Default is `None`.

#### Optional Parameters:

- **comorbidity_network_result** - DataFrame containing comorbidity network analysis results from `DiseaseNetPy.comorbidity_network()`. When provided, limits binomial test to significant disease pairs. Default is `None`.
- **correction** - P-value correction method for binomial tests. 
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
- **cutoff** - Significance threshold for adjusted binomial p-values. Default is `0.05`.
- **phecode_d1_col** - Column for disease 1 phecode. Default is `'phecode_d1'`.  
- **phecode_d2_col** - Column for disease 2 phecode. Default is `'phecode_d2'`.
- **n_nontemporal_col** - Column for non-temporal pair counts. Default is `'n_d1d2_nontemporal'`.  
- **n_temporal_d1d2_col** - Column for D1→D2 temporal counts. Default is `'n_d1d2_temporal'`.  
- **n_temporal_d2d1_col** - Column for D2→D1 temporal counts. Default is `'n_d2d1_temporal'`.  
- **significance_phi_col** - Column for phi-correlation significance. Default is `'phi_p_significance'`.  
- **significance_RR_col** - Column for RR significance. Default is `'RR_p_significance'`.
- **significance_coef_col** - Column for comorbidity significance. Default is `'comorbidity_p_significance'`.

#### After binomial test:

After the binomial test, we get a `binomial_result` (`pd.DataFrame`) including the non-temporal disease pairs (D₁-D₂ co-occurrence) and temporal disease pairs (D₁→D₂). The details of result is shows in the pipeline. Additionally, alternative p-value adjustment methods can be implemented to enhance statistical rigor.

```python
# Redo p-value adjustment using False Discovery Rate (FDR) with other methods
# For example: Benjamini-Hochberg method
binomial_result = dnt.binomial_multipletests(
    df=binomial_result,                              
    correction='fdr_bh',                             
    cutoff=0.05                                      
)
```

#### Save the result of binomial test:

Same to the pipeline, we can save the result to multiple formats including but not limited to:
- Tabular formats (.csv/.tsv)
- Spreadsheet files (.xlsx/.ods)
- Binary formats (.feather/.parquet)
- Database interfaces (SQL)

```python
# For example: save the binomial test results to a CSV file
binomial_result.to_csv('/your/project/path/dep_binomial.csv')
```

### 3.6 Comorbidity network analysis

Building upon the non-temporal disease pairs derived from the `com_strength_result()` and `binomial_result()`, we perform unconditional regression analyses on each pair to identify statistically significant associations for constructing the comorbidity network. This stage incorporates three distinct analysis methods (`RPCN`, `PCN_PCA`, and `CN`) with methodological specifics outlined in the pipeline's protocol documentation.

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

- **data** - DiseaseNetworkData object containing processed disease network data.
- **comorbidity_strength_result** - DataFrame containing comorbidity strength analysis results from `DiseaseNetPy.comorbidity_strength()`.
- **binomial_test_result** - DataFrame containing binomial test analysis results from `DiseaseNetPy.binomial_test()`. Default is `None`.
- **n_process** - Number of parallel processes. Values > 1 enable multiprocessing. Default is `1`.
- **covariates** - List of phenotypic covariates to include. If `None`, `sex` includes in **covariates**. Default is `None`.
- **method** - Comorbidity network analysis method to use. 
  - Options: 
    - `RPCN`: Regularized Partial Correlation Network. 
        - **alpha** - L1 penalty weight. Default is `None`  
        - **auto_penalty** - Auto-determine alpha. Default is `True`  
        - **alpha_range** - Alpha search range. Default is `(1,15)`  
        - **scaling_factor** - Alpha scaling factor. Default is `1`
    - `PCN_PCA`: Partial Correlation Network with PCA. 
        - **n_PC** - Number of principal components. Default is `5`  
        - **explained_variance** - Variance threshold. Default is `None`
    - `CN`: Correlation Network. Default is `RPCN`.
- **log_file** - Path/prefix for log file. If `None`, uses temporary directory with prefix `DiseaseNet_comorbidity_network_`. Default is `None`.

#### Optional parameters:

- **correction** - P-value correction method. 
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
- **cutoff** - Significance threshold for adjusted p-values. Default is `0.05`.
- **phecode_d1_col** - Column for disease 1 phecode. Default is `phecode_d1`.  
- **phecode_d2_col** - Column for disease 2 phecode. Default is `phecode_d2`.  
- **significance_phi_col** - Column for phi-correlation. Default is `phi_p_significance`.  
- **significance_RR_col** - Column for RR. Default is `RR_p_significance`.  
- **significance_binomial_col** - Column for binomial test. Default is `binomial_p_significance`.

#### After comorbidity network analysis:

After comorbidity network analysis, we get a `comorbidity_result` (`pd.DataFrame`) including non-temporal disease pairs. The details of result is shows in the pipeline. Additionally, alternative p-value adjustment methods can be implemented to enhance statistical rigor.

```python
# Redo p-value adjustment using False Discovery Rate (FDR) with other methods
# For example: Benjamini-Hochberg method
comorbidity_result = dnt.comorbidity_multipletests(
    df=comorbidity_result,                            
    correction='fdr_bh',                              
    cutoff=0.05                                      
)
```

#### Save the result of comorbidity network analysis:

Same to the pipeline, we can save the result to multiple formats including but not limited to:
- Tabular formats (.csv/.tsv)
- Spreadsheet files (.xlsx/.ods)
- Binary formats (.feather/.parquet)
- Database interfaces (SQL)

```python
# For example: save the comorbidity network analysis results to a CSV file
comorbidity_result.to_csv('/your/project/path/dep_comorbidity.csv')
```

#### 3.7 Trajectory analysis

Building upon the temporal disease pairs derived from `com_strength_result` (comorbidity strength estimation) and `binomial_result` (binomial test), we perform conditional regression analyses for each temporal disease pair (D₁→D₂) to identify statistically validated temporal disease pair, which are subsequently integrated into a directed disease trajectory network. This stage incorporates three analysis methods (`RPCN`, `PCN_PCA`, and `CN`) with methodological specifics outlined in the pipeline's protocol documentation.

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
- **comorbidity_strength_result** - DataFrame containing comorbidity strength analysis results from `DiseaseNetPy.comorbidity_strength()`.
- **binomial_test_result** - DataFrame containing binomial test analysis results from `DiseaseNetPy.binomial_test()`.
- **method** - Comorbidity network analysis method: 
 - `RPCN`: Regularized Partial Correlation Network. 
    - **alpha**: L1 penalty weight (ignored if auto_penalty). Default is `None`
    - **auto_penalty**: Auto-determine optimal alpha. Default is `True`
    - **alpha_range**: Alpha search range. Default is `(1,15)`
    - **scaling_factor**: Alpha scaling factor. Default is `1`
 - `PCN_PCA`: Partial Correlation Network with PCA. 
    - **n_PC**: Principal components count. Default is `5`
    - **explained_variance**: Variance threshold (overrides n_PC). Default is `None`
 - `CN`: Correlation Network. Default is `RPCN`.
- **matching_var_dict** - Dictionary specifying matching variables and criteria: Categorical/binary: `{'var':'exact'}`. Continuous: `{'var':max_diff}` (scalar > 0). Always use `'sex'` for sex matching. Default is `{'sex':'exact'}`.
- **matching_n** - Maximum matched controls per case. Default is `2`.
- **covariates** - List of phenotypic covariates (use `sex` for sex). Exclude matching variables. Default is `None`.
- **n_process** - Number of parallel processes (>1 enables multiprocessing). Default is `1`.
- **log_file** - Log file path/prefix. If `None`, uses temp dir with `DiseaseNet_trajectory_` prefix. Default is `None`.
- **enforce_time_interval** - Apply min/max time intervals for D2 outcome determination. Default is `True`.

#### Optional Parameters:

- **max_n_cases** - Maximum D2 cases to include (random sampling if exceeded). Default is `np.inf`.
- **global_sampling** - `True` for single sampling across all pairs, `False` for per-pair sampling. Default is `False`.
- **correction** - P-value correction method. 
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
- **cutoff** - Significance threshold for adjusted p-values. Default is `0.05`.
- **phecode_d1_col**: Disease 1 phecode column. Default is `'phecode_d1'`
- **phecode_d2_col**: Disease 2 phecode column. Default is `'phecode_d2'`
- **significance_phi_col**: Phi-correlation column. Default is `'phi_p_significance'`
- **significance_RR_col**: RR column. Default is `'RR_p_significance'`
- **significance_binomial_col**: Binomial test column. Default is `'binomial_p_significance'`

#### After trajectory analysis:

After trajectory analysis, we get a `trajectory_result` (`pd.DataFrame`) including temporal disease pairs. The details of result is shows in the pipeline. Additionally, alternative p-value adjustment methods can be implemented to enhance statistical rigor.

```python
# Redo p-value adjustment using False Discovery Rate (FDR) with other methods
# For example: Benjamini-Hochberg method
trajectory_result = dnt.trajectory_multipletests(
    df=trajectory_result,                             
    correction='fdr_bh',                            
    cutoff=0.05                                       
)
```

#### Save the result of trajectory analysis:

Same to the pipeline, we can save the result to multiple formats including but not limited to:
- Tabular formats (.csv/.tsv)
- Spreadsheet files (.xlsx/.ods)
- Binary formats (.feather/.parquet)
- Database interfaces (SQL)

```python
# For example: save the trajectory analysis analysis results to a CSV file
trajectory_result.to_csv('/your/project/path/dep_trajectory.csv')
```

## 4. Visualization

### 4.1 Initializing the plot object

After `import the Plot class from DiseaseNetPy.visualization` and initializing a `Plot` object with results from PheWAS analysis, comorbidity network analysis, and disease trajectory analysis (along with optional parameters). The Plot object provides four visualization method: PheWAS Plot, Comorbidity Network plot, Disease Trajectory plot, and Three Dimension Plot.

```python
# Create Plot object
from diseasenetpy.visualization import Plot

# cohort/matched cohort
result_plot = Plot(
    phewas_result=phewas_result,                         
    comorbidity_network_result=comorbidity_result,       
    disease_trajectory_result=trajectory_result,         
    exposure=495.2,                                      
    exposure_size=15,                                    
    exposure_location=(0,0,0),                           
)

# or exposed-only cohort
result_plot = Plot(
    phewas_result=phewas_result,                         
    comorbidity_network_result=comorbidity_result,       
    disease_trajectory_result=trajectory_result,         
    exposure=None,                                       
    exposure_size=None,                                 
    exposure_location=None,                              
)
```

- **comorbidity_result** - DataFrame containing comorbidity network analysis results with: Non-temporal disease pairs (D1-D2). Association metrics (beta coefficients, p-values). Significance indicators (True/False).
- **trajectory_result** - DataFrame containing temporal disease trajectory analysis with: Temporal disease pairs (source→target). Temporal association metrics. Significance indicators (True/False).
- **phewas_result** - DataFrame containing PheWAS analysis results with: Phecode diseases. Effect sizes (hazard ratios). Case counts. Disease system classifications
- **exposure** - Phecode identifier for primary exposure variable. Highlights exposure-disease relationships. Default is `None` (exposed-only cohort).
- **exposure_location** - Custom 3D coordinates (x,y,z) for exposure node positioning. Default is `None` (auto-positioned at (0,0,0)).
- **exposure_size** - Relative size scaling factor for exposure node. Default is `None` (exposed-only cohort).

#### Optional Parameters:

- **source**- Source disease column. Default is `'phecode_d1'`.
- **target**: Target disease column. Default is `'phecode_d2'`.
- **phewas_phecode**: Phecode column. Default `'phecode'`.
- **phewas_number**: Case count column. Default `'N_cases_exposed'`.
- **system_col**: Disease system column. Default `'system'`.
- **col_disease_pair**: Disease pair identifier column. Default `'name_disease_pair'`.
- **filter_phewas_col**: PheWAS significance column. Default `'phewas_p_significance'`.
- **filter_comorbidity_col**: Comorbidity significance column. Default `'comorbidity_p_significance'`.
- **filter_trajectory_col**: Trajectory significance column. Default `'trajectory_p_significance'`.
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

#### After initializing the Plot object:

After initializing the plot object, use `result_plot` object to visualization.

### 4.2 PheWAS plot

Generates a circular PheWAS (Phenome-Wide Association Study) plot. Creates a polar bar plot visualizing disease associations across different disease categories (systems). For a cohort/matched cohort study, the figure shows hazard ratios between exposure and outcome diseases. While the figure shows exposed number of diseases in the exposed-only cohort.

```python
# phewas plot
result_plot.phewas_plot(
    path="/your/project/path/",
    is_exposure_only=False,                         
)
```
- **path** - Output file path for saving the plot (including filename and extension)
- **is_exposure_only** - Boolean flag indicating whether the plot is for an exposure-only cohort study. Default is `False`

#### Optional parameters:

- **col_coef** - Column name containing effect size coefficients (e.g., hazard ratios or odds ratios). Default is `'phewas_coef'`
- **col_system** - Column name containing disease system/category classifications. Default is `"system"`
- **col_se** - Column name containing standard errors for effect sizes. Default is `'phewas_se'`
- **col_disease** - Column name containing disease names/descriptions. Default is `"disease"`
- **col_exposure** - Column name containing case counts for exposed individuals. Default is `"N_cases_exposed"`

#### After PheWAS plot:

After generating the PheWAS plot, the visualization will be exported as an image file in your preferred format (.jpg, .svg, and .png).

### 4.3 Comorbidity network plot

Generate a 2D visualization of the comorbidity network. Illustrates phecode disease co-occurrence patterns through network visualization, with each community using the Louvain clustering algorithm.

```python
# comorbidity network visualization
result_network.comorbidity_network_plot(
    path="/your/project/path/"
)
```
- **path** - Output file path for saving the interactive HTML visualization.  

#### Optional parameters:

- **max_radius** - Maximum radial position for nodes (in pixels). Controls outer boundary of the network. Default is `90.0`.
- **min_radius** - Minimum radial position for nodes (in pixels). Controls inner boundary. Default is `35.0`
- **layer_distance** - Spacing between concentric circles (in pixels). Affects radial grouping. Default is `40.0`.
- **size_reduction** - Scaling factor for node diameters (0-1 range). Adjusts visual prominence. Default is `0.5`
- **line_width** - Stroke width for comorbidity connections (in pixels). Default is `1.0`
- **line_color** - Color specification for comorbidity lines. Accepts: Named colors (e.g., `steelblue`). HEX codes (e.g., `#4682B4`). RGB tuples (e.g., `(70, 130, 180)`). Default is `black`
- **cluster_reduction_ratio** - Compression factor (0-1) for cluster tightness. Lower values create denser groupings. Default is `0.4`.
- **cluster_weight** - Edge attribute used for clustering calculations. Typically the association strength metric. Default is `comorbidity_beta`
- **font_style** - Font family for all text elements. Use web-safe fonts or loaded font families. Default is `Times New Roman`.

#### After comorbidity network plot:

After generating the comorbidity network plot, the visualization will be exported as an image file in your preferred format (.html).

### 4.4 Disease trajectory plot

Creates 2D network plots showing disease trajectories within each cluster, with nodes positioned hierarchically based on trajectory relationships. Each cluster is saved as a separate image file.

```python
# Disease trajectory visualization
result_network.trajectory_plot(
    path="/your/project/path/"
)
```

- **path** - Directory path where output visualization images will be saved. *Note:* Include trailing slash for proper path resolution (e.g., `/output/plots/`)

#### Optional parameters:

- **cluster_weight** Specifies the edge weight metric used for network clustering calculations. Default is `comorbidity_beta`.

#### After disease trajectory plot:

After generating the disease trajectory plot, the visualization will be exported as image files in your preferred format (.jpg, .svg, and .png).

### 4.5 Three dimension plot

Generates and saves a 3D visualization of comorbidity and disease trajectory networks. Combines trajectory and network visualizations by displaying comorbidity networks in the x-y plane while simultaneously showing disease trajectories in the x-z plane. Phecode diseases of the same cluster will gather together, and lines of plot are expressed as the trajectory.

```python
# three dimension visualization
result_network.three_dimension_plot(
    path="/your/project/path/"
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

#### After disease trajectory plot:

After generating the disease trajectory plot, the visualization will be exported as an image file in your preferred format (.html).

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

A class for handling disease network data creation and operations, for use in DiseaseNetPy module.

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

##### `concat` (classmethod)
```python
@classmethod
concat(
    cls,
    first_data: DiseaseNetworkData,
    second_data: DiseaseNetworkData,
    duplicates: str = 'raise'
) -> DiseaseNetworkData
```
Concatenate two `DiseaseNetworkData` objects.

**Parameters:**
- `first_data` (`DiseaseNetworkData`): First object to concatenate.
- `second_data` (`DiseaseNetworkData`): Second object to concatenate.
- `duplicates` (`str`): Duplicate handling. Choices:
  - `raise`: Error on duplicates.
  - `first`: Keep first object’s records.
  - `second`: Keep second object’s records.
  Defaults to `'raise'`.

**Returns:**
- `DiseaseNetworkData`

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
- `comorbidity_strength_result` (*pd.DataFrame*): DataFrame containing comorbidity strength analysis results produced by the 'DiseaseNetPy.comorbidity_strength' function.
- `comorbidity_network_result` (*pd.DataFrame, default=None*): DataFrame containing comorbidity network analysis results produced by the 'DiseaseNetPy.comorbidity_network' function. When provided, the binomial test is limited to disease pairs deemed significant in the comorbidity network analysis.
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
- `comorbidity_strength_result` (*pd.DataFrame*): DataFrame containing comorbidity strength analysis results produced by the 'DiseaseNetPy.comorbidity_strength' function.
- `binomial_test_result` (*pd.DataFrame, default=None*): DataFrame containing binomial test analysis results produced by the 'DiseaseNetPy.binomial_test' function.
- `method` (*str, default='RPCN'*): Specifies the comorbidity network analysis method to use. Choices are: - 'RPCN: Regularized Partial Correlation Network. - 'PCN_PCA: Partial Correlation Network with PCA. - 'CN': Correlation Network. **Additional Options for RPCN:** - 'alpha' : non-negative scalar The weight multiplying the l1 penalty term for other diseases covariates. Ignored if 'auto_penalty' is enabled. - 'auto_penalty' : bool, default=True If 'True', automatically determine the optimal 'alpha' based on model AIC value. - 'alpha_range' : tuple, default=(1,15) When 'auto_penalty' is True, search the optimal 'alpha' in this range. - 'scaling_factor' : positive scalar, default=1 The scaling factor for the alpha when 'auto_penalty' is True. **Additional Options for PCN_PCA:** - 'n_PC' : int, default=5 Fixed number of principal components to include in each model. - 'explained_variance' : float Determines the number of principal components based on the cumulative explained variance. Overrides 'n_PC' if specified.
- `covariates` (*list, default=None*): List of phenotypic covariates to include in the model. By default, includes ['sex'] and all covariates specified in the 'DiseaseNetPy.DiseaseNetworkData.phenotype_data()' function. To include the required variable sex as a covariate, always use 'sex' instead of its original column name. For other covariates specified in the 'DiseaseNetPy.DiseaseNetworkData.phenotype_data()' function, use their original column names.
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
- `comorbidity_strength_result` (*pd.DataFrame*): DataFrame containing comorbidity strength analysis results produced by the 'DiseaseNetPy.comorbidity_strength' function.
- `binomial_test_result` (*pd.DataFrame*): DataFrame containing binomial test analysis results produced by the 'DiseaseNetPy.binomial_test' function.
- `method` (*str, default='RPCN'*): Specifies the comorbidity network analysis method to use. Choices are: 
  - `'RPCN'`: Regularized Partial Correlation Network.
  - `'PCN_PCA'`: Partial Correlation Network with PCA. 
  - `'CN'`: Correlation Network. 

- `matching_var_dict` (*dict, default={'sex':'exact'}*): Specifies the matching variables and the criteria used for incidence density sampling. For categorical and binary variables, the matching criteria should always be `'exact'`. For continuous variables, provide a scalar greater than 0 as the matching criterion, indicating the maximum allowed difference when matching. To include the required variable sex as a matching variable, always use `'sex'` instead of its original column name. For other covariates specified in the `DiseaseNetPy.DiseaseNetworkData.phenotype_data()` function, use their original column names.
- `matching_n` (*int, default=2*): Specifies the maximum number of matched controls for each case.
- `max_n_cases` (*int, default=np.inf*): Specifies the maximum number of D2 cases to include in the analysis. If the number of D2 cases exceeds this value, a random sample of cases will be selected.
- `global_sampling` (*bool, default=False*): Indicates whether to perform independent incidence density sampling for each D1→D2 pair (if False), or to perform a single incidence density sampling for all Dx→D2 pairs with separate regression models for each D1→D2 pair (if True). Global sampling is recommended when processing large datasets, though it might reduce result heterogeneity.
- `covariates` (*list, default=None*): List of phenotypic covariates to include in the model. By default, includes all covariates specified in the `DiseaseNetPy.DiseaseNetworkData.phenotype_data()` function. Categorical and binary variables used for matching should not be included as covariates. Continuous variables used for matching can be included as covariates, but caution is advised. To include the required variable sex as a covariate, always use `sex` instead of its original column name. For other covariates specified in the `DiseaseNetPy.DiseaseNetworkData.phenotype_data()` function, use their original column names.
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
  - `enforce_time_interval` : bool, default=True If set to True, applies the specified minimum and maximum time intervals when determining the D2 outcome among individuals diagnosed with D1. These time interval requirements should be defined using the `DiseaseNetPy.DiseaseNetworkData.disease_pair()` function.
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
    exposure: float = None,
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
    SYSTEM: List[str] = None,
    COLOR: List[str] = None
)
```

A class for integrating and visualizing disease relationships from PHEWAS, comorbidity network, and trajectory analyses.

**Constructor Parameters:**
- `comorbidity_result` (`pd.DataFrame`): Non-temporal disease pairs with association metrics and significance flag.
- `trajectory_result` (`pd.DataFrame`): Temporal disease pairs (source→target) with metrics and significance flag.
- `phewas_result` (`pd.DataFrame`): PheWAS results including phecode, effect sizes, case counts, and system classifications.
- `exposure` (`float`, optional): Phecode for primary exposure; highlights exposure node. Defaults to `None`.
- `exposure_location` (`Tuple[float, float, float]`, optional): 3D coordinates for exposure node. Defaults to origin if `None`.
- `exposure_size` (`float`, optional): Scaling factor for exposure node size. Defaults to automatic.
- `source` (`str`): Column name for source disease (default `'phecode_d1'`).
- `target` (`str`): Column name for target disease (default `'phecode_d2'`).
- `phewas_phecode` (`str`): Column for phecode in PHEWAS results (default `'phecode'`).
- `phewas_number` (`str`): Column for case counts (default `'N_cases_exposed'`).
- `system_col` (`str`): Column for disease system (default `'system'`).
- `col_disease_pair` (`str`): Column for pair identifier (default `'name_disease_pair'`).
- `filter_phewas_col` (`str`): Column for PHEWAS significance filter.
- `filter_comorbidity_col` (`str`): Column for comorbidity significance filter.
- `filter_trajectory_col` (`str`): Column for trajectory significance filter.
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
- `max_radius`: Maximum radial distance for node placement (default: 180.0)
- `min_radius`: Minimum radial distance for node placement (default: 35.0)
- `line_color`: Color for trajectory lines (default: "black")
- `line_width`: Width for trajectory lines (default: 1.0)
- `size_reduction`: Scaling factor for node sizes (default: 0.5)
- `cluster_reduction_ratio`: Cluster compression factor for layout (default: 0.4)
- `cluster_weight`: Edge weight metric used for clustering (default: "comorbidity_beta")
- `layer_distance`: Vertical distance between layers (default: 40.0)
- `layout_width`: Figure width in pixels (default: 900.0)
- `layout_height`: Figure height in pixels (default: 900.0)
- `font_style`: Font family for text elements (default: 'Times New Roman')
- `font_size`: Base font size in points (default: 15.0)

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
- `max_radius`: Maximum radial position for nodes (default: 90.0)
- `min_radius`: Minimum radial position for nodes (default: 35.0)
- `size_reduction`: Scaling factor for node sizes (default: 0.5)
- `cluster_reduction_ratio`: Compression factor for cluster layout (default: 0.4)
- `cluster_weight`: Edge weight metric for clustering (default: "comorbidity_beta")
- `line_width`: Width of comorbidity lines (default: 1.0)
- `line_color`: Color of comorbidity lines (default: "black")
- `layer_distance`: Distance between concentric circles (default: 40.0)
- `font_style`: Font family for text elements (default: "Times New Roman")

---

##### `trajectory_plot`

```python
trajectory_plot(
    self,
    path: str,
    cluster_weight: str = 'comorbidity_beta'
) -> None
```

Generate and save trajectory plots per cluster as (.png/.svg/.jpg files).

**Parameters:**
- `path`: Directory path to save output images
- `cluster_weight`: Edge weight metric used for clustering (default: "comorbidity_beta")

---

##### `phewas_plot`

```python
phewas_plot(
    self,
    path: str,
    col_coef: str = 'phewas_coef',
    col_system: str = 'system',
    col_se: str = 'phewas_se',
    col_disease: str = 'disease',
    is_exposure_only: bool = False,
    col_exposure: str = 'N_cases_exposed'
) -> None
```

Creates a polar bar plot visualizing disease associations across different disease categories (systems)

**Parameters:**
- `path`: Output file path for saving the plot
- `col_coef`: Column name for effect size coefficients (default: "phewas_coef")
- `col_system`: Column name for disease system/category (default: "system")
- `col_se`: Column name for standard errors (default: "phewas_se")
- `col_disease`: Column name for disease names (default: "disease")
- `is_exposure_only`: Identifier of exposure (default: False)
- `col_exposure`: Column name for exposure number (default: "N_cases_exposed")