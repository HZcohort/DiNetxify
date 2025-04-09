# DiseaseNetPy

DiseaseNetPy is a Python package designed for comprehensive disease network analysis and visualization. This novel disease network analysis approach (three-dimensional disease network analysis) integrates and refines existing disease trajectory and comorbidity network analysis methods. This new approach enhances disease association verification by incorporating regularized partial correlations. It also facilitates robust identification and visualization of disease clusters (i.e., groups of depression-associated diseases with high within-group connectivity) through both non-temporal (illustrated by the x-axis and y-axis) and temporal (z-axis) dimensions.

## Table of Contents

- [Installation](#installation)
- Quick Start
  - [Example: Matched Cohort Study Design](#example-matched-cohort-study-design)
- API Reference
  - Classes
    - [DiseaseNetworkData](#diseasenetworkdata)
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
- [Issues reporting and recommendations](#trajectory_multipletests)
- [License](#license)

## Installation

You can install DiseaseNetPy via `pip`. Ensure you have Python 3.7 or higher.

```bash
pip install diseasenetpy

#optional packages
pip install lifelines #for phewas analysis with lifelines_disable set to False

```

## Quick Start

This guide walks you through a typical workflow using DiseaseNetPy for a matched cohort study design. The process involves data preparation, PheWAS analysis, comorbidity strength estimation, binomial testing, comorbidity network analysis, and trajectory analysis.

```python
import diseasenetpy as dnt
# Step 1: Create DiseaseNetworkData object
# Define required columns and other covariates columns
col_dict = {
    'Participant ID': 'new_index',          # Maps the participant identifier
    'Exposure': 'outcome',                  # Defines exposure status (0 or 1)
    'Sex': 'sex',                           # Indicates sex (1 for female, 0 for male)
    'Index date': 'date_start',             # Start date of the study
    'End date': 'time_end',                 # End date of the study
    'Match ID': 'match_2'                   # Identifier for matching group
}
vars_lst = ['age', 'social', 'BMI', 'smoking', 'drinking']  # List of covariates to be used

# Initialize the data object with study design and phecode level
data = dnt.DiseaseNetworkData(
    study_design='matched cohort',          # Type of study design
    phecode_level=1,                        # Level of phecode (1 or 2)
    date_fmt='%Y-%m-%d'                     # Date format in data files
)

# Load the phenotype CSV file into the data object
data.phenotype_data(
    phenotype_data_path='/your/project/path/phenotype.csv',  # Path to phenotype data
    column_names=col_dict,                                   # Column mappings
    covariates=vars_lst                                      # Covariates to include
)

# Merge with the first medical records file (CSV)
data.merge_medical_records(
    medical_records_data_path='/your/project/path/inp_1.csv',  # Path to first medical records file
    diagnosis_code='ICD-10-WHO',                               # Diagnosis code type
    column_names={
        'Participant ID': 'eid',                               # Participant ID column in medical records
        'Diagnosis code': 'diag_icd10',                        # Diagnosis code column
        'Date of diagnosis': 'date'                            # Diagnosis date column
    }
)

# Step 2: Use the pipeline disease network
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
        'smoking', 
        'drinking'
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
### Test
In the test folder, we provide three comprehensive examples demonstrating different study designs to help users understand and implement their analyses effectively. Each example includes relevant code, output, and visualizations.

1. Example 1: Matched Cohort Study Design
A matched cohort study is an observational research design where exposed and unexposed groups are matched based on specific  covariates (e.g., age, sex, comorbidities) to reduce bias and improve comparability. 

2. Example 2: Cohort Study Design
A cohort study is an observational research design that follows groups of individuals (cohorts) over time to assess the association between exposures (e.g., risk factors, treatments) and outcomes (e.g., disease incidence, mortality). 

3. Example 3: Exposed-only Cohort Design
The exposed-only cohort design is a variation of cohort studies where only individuals exposed to a risk factor are followed over time, and their outcomes are compared to expected population rates (external controls) rather than an internal unexposed group.

### Workflows and Example: Matched Cohort Study Design

```python
import diseasenetpy as dnt

# Step 1: Create DiseaseNetworkData object
# Define required columns and other covariates columns
col_dict = {
    'Participant ID': 'new_index',          # Maps the participant identifier
    'Exposure': 'outcome',                  # Defines exposure status (0 or 1)
    'Sex': 'sex',                           # Indicates sex (1 for female, 0 for male)
    'Index date': 'date_start',             # Start date of the study
    'End date': 'time_end',                 # End date of the study
    'Match ID': 'match_2'                   # Identifier for matching group
}
vars_lst = ['age', 'social', 'BMI', 'smoking', 'drinking']  # List of covariates to be used

# Initialize the data object with study design and phecode level
data = dnt.DiseaseNetworkData(
    study_design='matched cohort',          # Type of study design
    phecode_level=1,                        # Level of phecode (1 or 2)
    date_fmt='%Y-%m-%d'                     # Date format in data files
)

# Load the phenotype CSV file into the data object
data.phenotype_data(
    phenotype_data_path='/your/project/path/phenotype.csv',  # Path to phenotype data
    column_names=col_dict,                                   # Column mappings
    covariates=vars_lst                                      # Covariates to include
)

# Merge with the first medical records file (CSV)
data.merge_medical_records(
    medical_records_data_path='/your/project/path/inp_1.csv',  # Path to first medical records file
    diagnosis_code='ICD-10-WHO',                                # Diagnosis code type
    column_names={
        'Participant ID': 'eid',                                # Participant ID column in medical records
        'Diagnosis code': 'diag_icd10',                        # Diagnosis code column
        'Date of diagnosis': 'date'                            # Diagnosis date column
    }
)

# Merge with the second medical records file (CSV)
data.merge_medical_records(
    medical_records_data_path='/your/project/path/inp_2.csv',  # Path to second medical records file
    diagnosis_code='ICD-10-WHO',                                # Diagnosis code type
    column_names={
        'Participant ID': 'eid',                                # Participant ID column in medical records
        'Diagnosis code': 'diag_icd10',                        # Diagnosis code column
        'Date of diagnosis': 'date'                            # Diagnosis date column
    }
)

# Describe the basic information fo phenotype data and medical data
data.Table1(
    continuous_stat_mode="auto"                # method to pre-processing continuous variable
)

# Save the data object for later use
data.save('/your/project/path/dep')  # Path to save the data object

# Step 2: PheWAS Analysis
# Reminder:
# When using multiprocessing, ensure that the code is enclosed within the following block.
# This prevents entering a never ending loop of new process creation.
if __name__ == "__main__":
    phewas_result = dnt.phewas(
        data=data,                                             # DiseaseNetworkData object
        proportion_threshold=0.01,                            # Minimum proportion of cases to include
        n_process=2,                                          # Number of parallel processes
        system_exl=[                                           # Phecode systems to exclude
            'symptoms', 'others', 'injuries & poisonings', 'pregnancy complications'
        ],
        covariates=['age', 'social', 'BMI', 'smoking', 'drinking'],  # Covariates to adjust for
        lifelines_disable=True,                                # Disable lifelines for faster computation
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

# Step 3: Generate Disease Pair for Each Individual and Update the Data Object
data.disease_pair(
    phewas_result=phewas_result,               # Filtered PheWAS results
    min_interval_days=30,                      # Minimum interval between diagnoses (30 days here)
    max_interval_days=365.25*5,                # Maximum interval between diagnoses (5 years here)
    force=True                                 # Overwrite existing data if necessary (be cautious)
)

# Save the updated data object with disease pairs
data.save('/your/project/path/dep_withtra')    # Path to save the updated data object

# Step 4: Comorbidity Strength Estimation
# Reminder:
# When using multiprocessing, ensure that the code is enclosed within the following block.
# This prevents entering a never ending loop of new process creation.
if __name__ == "__main__":
    com_strength_result = dnt.comorbidity_strength(
        data=data,                                     # DiseaseNetworkData object
        proportion_threshold=0.001,                   # Minimum proportion for comorbidity
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
    df=com_strength_result,                        # DataFrame with comorbidity strength results
    correction_phi='fdr_bh',                       # P-value correction for phi-correlation
    correction_RR='fdr_bh',                         # P-value correction for Relative Risk
    cutoff_phi=0.05,                                # Significance threshold for phi-correlation
    cutoff_RR=0.05                                  # Significance threshold for Relative Risk
)

# Step 5: Binomial Test
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

# Step 6: Comorbidity Network Analysis
# Reminder:
# When using multiprocessing, ensure that the code is enclosed within the following block.
# This prevents entering a never ending loop of new process creation.
if __name__ == "__main__":
    comorbidity_result = dnt.comorbidity_network(
        data=data,                                       # DiseaseNetworkData object
        comorbidity_strength_result=com_strength_result, # Comorbidity strength results
        binomial_test_result=binomial_result,           # Binomial test results
        n_process=2,                                   # Number of parallel processes
        covariates=['social', 'BMI', 'smoking', 'drinking', 'sex'],  # Covariates to adjust for
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

# Step 7: Trajectory Analysis
# Reminder:
# When using multiprocessing, ensure that the code is enclosed within the following block.
# This prevents entering a never ending loop of new process creation.
if __name__ == "__main__":
    trajectory_result = dnt.disease_trajectory(
        data=data,                                       # DiseaseNetworkData object
        comorbidity_strength_result=com_strength_result, # Comorbidity strength results
        binomial_test_result=binomial_result,           # Binomial test results
        method='RPCN',                                   # Trajectory analysis method ('CN', 'PCN_PCA', 'RPCN')
        n_process=2,                                     # Number of parallel processes
        matching_var_dict={'age': 2, 'sex': 'exact'},    # Matching variables and criteria
        matching_n=5,                                    # Number of matched controls per case
        enforce_time_interval=False,                     # Enforce time interval in trajectory analysis
        covariates=['social', 'BMI', 'smoking', 'drinking'],  # Covariates to adjust for
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

# Step 8: Result visualization (three-dimension visualization, comorbidity network visualization, significant trajectory visualization)
# Create ThreeDimensionalDiseaseNetwork object
result_network = dnt.visualization.ThreeDimensionalDiseaseNetwork(
  comorbidity_network_result=comorbidity_result,       # DataFrame with comorbidity network results
  disease_trajectory_result=trajectory_result,         # DataFrame with trajectory analysis results
  phewas_result=phewas_result,                         # DataFrame with PheWAS results
  exposure_disease=594.0,                              # Phecode of exposure. Default to 9999, means that it's a registry study
  exposure_disease_location=(0,0,0),                   # Three-dimensional coordinate. Default to (0,0,0), it's a initial point to plot
  exposure_disease_size=1,                             # Size of initial point (phecode of exposure) to plot. Default to 1
  source="phecode_d1",                                 # Column name of D1 (D1->D2). Defaults to 'phecode_d1'
  target="phecode_d2",                                 # Column name of D2 (D1->D2). Defaults to 'phecode_d2'
)

# three-dimension visualization
result_network.plot_3d(
  max_radius=45,                                       # Maximum of radius in one sector(cluster)
  min_radius=15,                                       # Minimum of radius in one sector(cluster)
  plot_method="full",                                  # Method of plot ("full", "half", "compact")
  line_color="black",                                  # Color of lines connected with each node(phecode) in the "full" plot method 
  line_width=2,                                        # Size of lines connected with each node(phecode)
  layer_distance=20,                                   # Distance of two adjoining layers
  file_name="/your/project/path"                       # Path to save three-dimension visualization
  layout_width=900,                                    # Width of layout in the figure. Defaults to 900
  layout_height=900,                                   # Height of layout in the figure. Defaults to 900
  font_style="Times New Roman",                        # Font style of layout in the figure. Defaults to "Times New Roman"
  font_size=15,                                        # Font size of layout in the figure. Defaults to 15
  location_method="random"                             # Method to calculate the three-dimension location of nodes(phecodes). Defaults to 'random'
)

# comorbidity network visualization
result_network.comorbidity_network_plot(
  max_radius=45,                                       # Maximum of radius in one sector(cluster)
  min_radius=15,                                       # Minimum of radius in one sector(cluster)
  line_width=1,                                        # Size of lines connected with each node(phecode). Defaults to 1
  source="phecode_d1",                                 # Column name of D1 (D1->D2). Defaults to 'phecode_d1'
  target="phecode_d2",                                 # Column name of D2 (D1->D2). Defaults to 'phecode_d2'
  location_method="random",                            # Method to calculate the three-dimension location of nodes(phecodes). Defaults to 'random'
  line_color="black"                                   # Color of lines connected with each node(phecode)
)

# significant trajectory visualization
result_network.incluster_trajectory_plot(
  distance=5,                                          # Distance of each nodes(phecodes) in x-y plane
  layer_distance=5,                                    # Distance of two adjoining layers
  line_width=1,                                        # Size of lines connected with each node(phecode)
  line_color="black",                                  # Color of lines connected with each node(phecode)
  max_radius=45,                                       # Maximum of radius in one sector(cluster)
  min_radius=15,                                       # Minimum of radius in one sector(cluster)
  location_method="random",                            # Method to calculate the three-dimension location of nodes(phecodes). Defaults to 'random'
  source="phecode_d1",                                 # Column name of D1 (D1->D2). Defaults to 'phecode_d1'
  target="phecode_d2"                                  # Column name of D2 (D1->D2). Defaults to 'phecode_d2'
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
- `date_fmt : str, default='%Y-%m-%d'`
  - The format of the date fields in your phenotype and medical records data.
- `phecode_version : str, default='1.2'`
  - The version of the phecode system used for converting diagnosis codes. Currently, only version 1.2 is supported.

------

##### Methods

###### `phenotype_data`

```python
phenotype_data(self, phenotype_data_path:str, column_names:dict, 
               covariates:list, force:bool=False)
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

- `force : bool, default=False`

  - If `True`, overwrites existing data attributes even if they contain data.

**Returns:**

- `None`
  - Modifies the object's main data attribute in-place.

------

###### `merge_medical_records`

```python
merge_medical_records(self, medical_records_data_path:str, diagnosis_code:str,
                      column_names:dict, date_fmt:str=None, chunksize:int=1000000)
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

###### `disease_pair`

```python
disease_pair(self, phewas_result:pd.DataFrame, min_interval_days:int=0,
             max_interval_days:int=np.inf, force:bool=False, **kwargs)
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

### Functions

#### `phewas`

```python
dnt.phewas(data:DiseaseNetworkData, covariates:list=None,
           proportion_threshold:float=None, n_threshold:int=None,
           n_process:int=1, correction:str='bonferroni', cutoff:float=0.05,
           system_inc:list=None, system_exl:list=None,
           phecode_inc:list=None, phecode_exl:list=None, log_file:str=None,
           lifelines_disable:bool=False) -> pd.DataFrame
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
dnt.phewas_multipletests(df:pd.DataFrame, correction:str='bonferroni', 
                         cutoff:float=0.05) -> pd.DataFrame
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
dnt.comorbidity_strength(data:DiseaseNetworkData, proportion_threshold:float=None, 
                         n_threshold:int=None, n_process:int=1, log_file:str=None, 
                         correction_phi:str='bonferroni', cutoff_phi:float=0.05, 
                         correction_RR:str='bonferroni', 
                         cutoff_RR:float=0.05) -> pd.DataFrame
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
dnt.comorbidity_strength_multipletests(df:pd.DataFrame, correction_phi:str='bonferroni', 
                                       cutoff_phi:float=0.05, 
                                       correction_RR:str='bonferroni', 
                                       cutoff_RR:float=0.05) -> pd.DataFrame
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
dnt.binomial_test(data:DiseaseNetworkData, comorbidity_strength_result:pd.DataFrame, 
                  n_process:int=1, log_file:str=None, 
                  correction:str='bonferroni', cutoff:float=0.05, 
                  enforce_temporal_order:bool=False, **kwargs) -> pd.DataFrame
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
dnt.binomial_multipletests(df:pd.DataFrame, correction:str='bonferroni', 
                           cutoff:float=0.05) -> pd.DataFrame
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
dnt.comorbidity_network(data:DiseaseNetworkData,comorbidity_strength_result:pd.DataFrame, 
                        binomial_test_result:pd.DataFrame, method:str='RPCN', 
                        covariates:list=None, n_process:int=1, log_file:str=None, 
                        correction:str='bonferroni', 
                        cutoff:float=0.05, **kwargs) -> pd.DataFrame
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
dnt.comorbidity_multipletests(df:pd.DataFrame, correction:str='bonferroni', 
                              cutoff:float=0.05) -> pd.DataFrame:
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
dnt.disease_trajectory(data:DiseaseNetworkData, comorbidity_strength_result:pd.DataFrame, 
                       binomial_test_result:pd.DataFrame, 
                       method:str='RPCN', matching_var_dict:dict={'sex':'exact'}, 
                       matching_n:int=2, covariates:list=None,
                       n_process:int=1, log_file:str=None, correction:str='bonferroni', 
                       cutoff:float=0.05, **kwargs) -> pd.DataFrame
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
dnt.trajectory_multipletests(df:pd.DataFrame, correction:str='bonferroni', 
                             cutoff:float=0.05) -> pd.DataFrame
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

## Issues reporting and recommendations

Please contact:
Can Hou: houcan@wchscu.cn

Haowen Liu: haowenliu81@gmail.com

## License

DiseaseNetPy is released under [The GNU General Public License v3.0](https://www.gnu.org/licenses/gpl-3.0.html)