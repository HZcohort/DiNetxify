# DiseaseNetPy

DiseaseNetPy is a Python package designed for comprehensive disease network analysis and visualization. This novel disease network analysis approach (three-dimensional disease network analysis) integrates and refines existing disease trajectory and comorbidity network analysis methods. This new approach enhances disease association verification by incorporating regularized partial correlations. It also facilitates robust identification and visualization of disease clusters (i.e., groups of depression-associated diseases with high within-group connectivity) through both non-temporal (illustrated by the x-axis and y-axis) and temporal (z-axis) dimensions.

## Table of Contents

- [Installation](#installation)
- [Study design](#study-design)
- [Data of test](#data-of-test)
  - [Dataset Overview](#dataset-overview)
  - [Data Description](#data-description)
  - [Usage Guide](#usage-guide)
  - [Important Notes](#important-notes)
- Quick Start
  - [Test](#test)
  - [Workflows and Example: Matched Cohort Study Design](#workflows-and-example-matched-cohort-study-design)
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

### Study design
This package suports three different study, such as matched cohort study, cohort study, and exposed-only cohort

1. Matched Cohort Study
A matched cohort study is an observational research design where exposed and unexposed groups are matched based on specific  covariates (e.g., age, sex, comorbidities) to reduce bias and improve comparability. 

2. Cohort Study
A cohort study is an observational research design that follows groups of individuals (cohorts) over time to assess the association between exposures (e.g., risk factors, treatments) and outcomes (e.g., disease incidence, mortality). 

3. Exposed-only Cohort
The exposed-only cohort design is a variation of cohort studies where only individuals exposed to a risk factor are followed over time, and their outcomes are compared to expected population rates (external controls) rather than an internal unexposed group.

## Data of test
### Dataset Overview
There are three dummy data (one phnotypic data, and two medical records data) which are suitable for chort, matched cohort, exposed-only cohort study.
- **Format**: Both are CSV files
- **Size**: 
  - phnotypic data: 60,000 records | 3.72 MB
  - ICD9 medical records data: 10,188 records  | 227 kB
  - ICD10 medical records data: 1,048,576 records  | 36.2 MB

### Data Description
- **Fields**:
  phnotypic data

  | Column | Type | Description | Example |
  |--------|------|-------------|---------|
  | ID     | int  | Unique identifier of one individual | 1001 |
  | group_id | str | Identifier of group (matched cohort) | "group_0" |
  | exposure | int  | Exposure label (cohort/matched cohort) | 0 or 1 |
  | date_start | datetime | Start time of follow-up | 2016/10/13 |
  | date_end | datetime  | End time of follow-up | 2022/8/2 |
  | age | float  | Age for year | 71.23251 |
  | sex | int  | Identifier of sex | 0 or 1 |
  | BMI | str  | The level of BMI | c1 |

  medical records data

  | Column | Type | Description | Example |
  |--------|------|-------------|---------|
  | ID     | int  | Unique identifier of one individual | 1001 |
  | dia_date | datetime | Time of diagnosis | 2018/10/01 |
  | diag_icd10/diag_icd9 | str  | ICD9/ICD10 medical records | "L905"/"E950" |

### Usage Guide
1. **Download**
```bash
git clone {https://github.com/HZcohort/DiseaseNetPy.git}
# or
wget {https://github.com/HZcohort/DiseaseNetPy.git}
```

2. **Loading Example**
  In the tests folder, there are three python scripts (cohort.py, mathced cohort.py, and exposed-only cohort.py) which are suitable for three study designs (cohort/mathced cohort/exposed-only cohort) respectively.
```bash
cd tests
# test of cohort study
python cohort.py

# or test of matched cohort study
python matched cohort.py

# or test of exposed-only cohort study
python exposed-only cohort.py
```

3. **Recommended Uses**:
   - python 3.13.0
   - A multi-core processor with more than 8 cores

### Important Notes
⚠️ **Usage Restrictions**:
- [✔] Commercial use allowed
- [x] Attribution required
- [✔] Redistribution prohibited

📌 **Data Characteristics**:
- There are no missing values
- Covariates (sex, bmi, age)

## Quick Start
This guide walks you through a typical workflow using DiseaseNetPy for a matched cohort study design. The process involves data preparation, PheWAS analysis, comorbidity strength estimation, binomial testing, comorbidity network analysis, and trajectory analysis. The difference in each study design lies in Step 1 (creating the DiseaseNetworkData object). The data are based on test dataset.

Step 1: Create DiseaseNetworkData object
Define required columns and other covariates columns

Example for matched cohort

```python
import diseasenetpy as dnt
col_dict = {
    'Participant ID': 'ID',                # Maps the participant identifier
    'Exposure': 'exposure',                # Defines exposure status (0 or 1)
    'Sex': 'sex',                          # Indicates sex (1 for female, 0 for male)
    'Index date': 'date_start',            # Start date of the study
    'End date': 'date_end',                # End date of the study
    'Match ID': 'group_id'                 # Identifier for matching group
}
vars_lst = ['age', 'BMI']  # List of covariates to be used

# Initialize the data object with study design and phecode level
data = dnt.DiseaseNetworkData(
    study_design='matched cohort',          # Type of study design
    phecode_level=1,                        # Level of phecode (1 or 2)
    date_fmt='%Y-%m-%d'                     # Date format in data files
)
```

Example for cohort

```python
import diseasenetpy as dnt
col_dict = {
    'Participant ID': 'ID',                # Maps the participant identifier
    'Exposure': 'exposure',                # Defines exposure status (0 or 1)
    'Sex': 'sex',                          # Indicates sex (1 for female, 0 for male)
    'Index date': 'date_start',            # Start date of the study
    'End date': 'date_end',                # End date of the study
}
vars_lst = ['age', 'BMI']  # List of covariates to be used

# Initialize the data object with study design and phecode level
data = dnt.DiseaseNetworkData(
    study_design='cohort',          # Type of study design
    phecode_level=1,                # Level of phecode (1 or 2)
    date_fmt='%Y-%m-%d'             # Date format in data files
)
```

Example for exposed-only cohort

```python
import diseasenetpy as dnt
col_dict = {
    'Participant ID': 'ID',                # Maps the participant identifier
    'Sex': 'sex',                          # Indicates sex (1 for female, 0 for male)
    'Index date': 'date_start',            # Start date of the study
    'End date': 'date_end',                # End date of the study
}
vars_lst = ['age', 'BMI']  # List of covariates to be used

# Initialize the data object with study design and phecode level
data = dnt.DiseaseNetworkData(
    study_design='exposed-only cohort',    # Type of study design
    phecode_level=1,                       # Level of phecode (1 or 2)
    date_fmt='%Y-%m-%d'                    # Date format in data files
)
```

Step 2: Data harmonization

```python
# Load the phenotype CSV file into the data object
data.phenotype_data(
    phenotype_data_path='/test/data/dummy_cohort.csv',  # Path to phenotype data
    column_names=col_dict,                              # Column mappings
    covariates=vars_lst                                 # Covariates to include
)

# Merge with the first medical records file (CSV)
data.merge_medical_records(
    medical_records_data_path='/test/data/dummy_EHR_ICD10.csv',  # Path to first medical records file
    diagnosis_code='ICD-10-WHO',                                 # Diagnosis code type
    column_names={
        'Participant ID': 'ID',                                  # Participant ID column in medical records
        'Diagnosis code': 'diag_icd10',                          # Diagnosis code column
        'Date of diagnosis': 'dia_date'                          # Diagnosis date column
    }
)

# Merge with the second medical records file (CSV)
data.merge_medical_records(
    medical_records_data_path="/test/data/dummy_EHR_ICD9.csv",  # Path to first medical records file
    diagnosis_code="ICD-9-WHO",                                 # Diagnosis code type
    column_names={
        'Participant ID':'ID',                                  # Participant ID column in medical records
        'Diagnosis code':'diag_icd9',                           # Diagnosis code column
        'Date of diagnosis':'dia_date'                          # Diagnosis date column
    }
)
```

Step 3: Use the pipeline disease network

Reminder:
When using multiprocessing, ensure that the code is enclosed within the following block.
This prevents entering a never ending loop of new process creation.

```python
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

### Basic workflows and examples
The difference in each study design lies in Step 1 (creating the DiseaseNetworkData object).

Step 1: Create DiseaseNetworkData object

Example for matched cohort
```python
# Define required columns and other covariates columns
import diseasenetpy as dnt
col_dict = {
    'Participant ID': 'ID',                # Maps the participant identifier
    'Exposure': 'exposure',                # Defines exposure status (0 or 1)
    'Sex': 'sex',                          # Indicates sex (1 for female, 0 for male)
    'Index date': 'date_start',            # Start date of the study
    'End date': 'date_end',                # End date of the study
    'Match ID': 'group_id'                 # Identifier for matching group
}
vars_lst = ['age', 'BMI']  # List of covariates to be used

# Initialize the data object with study design and phecode level
data = dnt.DiseaseNetworkData(
    study_design='matched cohort',          # Type of study design
    phecode_level=1,                        # Level of phecode (1 or 2)
    date_fmt='%Y-%m-%d'                     # Date format in data files
)
```

Example for cohort

```python
import diseasenetpy as dnt
col_dict = {
    'Participant ID': 'ID',                # Maps the participant identifier
    'Exposure': 'exposure',                # Defines exposure status (0 or 1)
    'Sex': 'sex',                          # Indicates sex (1 for female, 0 for male)
    'Index date': 'date_start',            # Start date of the study
    'End date': 'date_end',                # End date of the study
}
vars_lst = ['age', 'BMI']  # List of covariates to be used

# Initialize the data object with study design and phecode level
data = dnt.DiseaseNetworkData(
    study_design='cohort',          # Type of study design
    phecode_level=1,                # Level of phecode (1 or 2)
    date_fmt='%Y-%m-%d'             # Date format in data files
)
```

Example for exposed-only cohort

```python
import diseasenetpy as dnt
col_dict = {
    'Participant ID': 'ID',                # Maps the participant identifier
    'Sex': 'sex',                          # Indicates sex (1 for female, 0 for male)
    'Index date': 'date_start',            # Start date of the study
    'End date': 'date_end',                # End date of the study
}
vars_lst = ['age', 'BMI']  # List of covariates to be used

# Initialize the data object with study design and phecode level
data = dnt.DiseaseNetworkData(
    study_design='exposed-only cohort',    # Type of study design
    phecode_level=1,                       # Level of phecode (1 or 2)
    date_fmt='%Y-%m-%d'                    # Date format in data files
)
```

Step 2: Data harmonization

```python
# Load the phenotype CSV file into the data object
data.phenotype_data(
    phenotype_data_path='/test/data/dummy_cohort.csv',  # Path to phenotype data
    column_names=col_dict,                              # Column mappings
    covariates=vars_lst                                 # Covariates to include
)

# Merge with the first medical records file (CSV)
data.merge_medical_records(
    medical_records_data_path='/test/data/dummy_EHR_ICD10.csv',  # Path to first medical records file
    diagnosis_code='ICD-10-WHO',                                 # Diagnosis code type
    column_names={
        'Participant ID': 'ID',                                  # Participant ID column in medical records
        'Diagnosis code': 'diag_icd10',                          # Diagnosis code column
        'Date of diagnosis': 'dia_date'                          # Diagnosis date column
    }
)

# Merge with the second medical records file (CSV)
data.merge_medical_records(
    medical_records_data_path="/test/data/dummy_EHR_ICD9.csv",  # Path to first medical records file
    diagnosis_code="ICD-9-WHO",                                 # Diagnosis code type
    column_names={
        'Participant ID':'ID',                                  # Participant ID column in medical records
        'Diagnosis code':'diag_icd9',                           # Diagnosis code column
        'Date of diagnosis':'dia_date'                          # Diagnosis date column
    }
)

# Describe the basic information fo phenotype data and medical data
data.Table1(
    continuous_stat_mode="auto"                # method to pre-processing continuous variable
)

# Save the data object for later use
data.save('/your/project/path/dep')  # Path to save the data object
```

Step 3: PheWAS Analysis
Reminder:
When using multiprocessing, ensure that the code is enclosed within the following block.
This prevents entering a never ending loop of new process creation.

```python
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
    df=phewas_result,                         # DataFrame with PheWAS results
    correction='fdr_bh',                      # P-value correction method
    cutoff=0.05                               # Significance threshold
)
```

Step 4: Generate Disease Pair for Each Individual and Update the Data Object

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

Step 5: Comorbidity Strength Estimation
Reminder:
When using multiprocessing, ensure that the code is enclosed within the following block.
This prevents entering a never ending loop of new process creation.

```python
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

Step 6: Binomial Test

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

Step 7: Comorbidity Network Analysis
Reminder:
When using multiprocessing, ensure that the code is enclosed within the following block.
This prevents entering a never ending loop of new process creation.

```python
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

Step 8: Trajectory Analysis
Reminder:
When using multiprocessing, ensure that the code is enclosed within the following block.
This prevents entering a never ending loop of new process creation.

```python
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

Step 9: Result visualization (three-dimension visualization, comorbidity network visualization, trajectory visualization, PheWAS visualization)

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

# trajectory visualization
result_network.trajectory_plot(
    path="/your/project/path/",                          # Directory path to save output images
    cluster_weight="comorbidity_beta",                   # Edge weight metric used for clustering (default: "comorbidity_beta")
)

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