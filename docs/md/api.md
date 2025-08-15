# API Reference

Below is a concise reference for ***DiNetxify***’s classes and functions, summarizing their signatures and parameters. This is useful when writing your own scripts or if you need to quickly recall how to call a function.

## Class `DiseaseNetworkData`

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
- `min_required_icd_codes` (`int`): The minimum number of ICD codes mapping to a specific phecode required for the phecode to be considered valid. For example, if set to 2, a single diagnosis record will not be sufficient to count as an occurrence. Ensure that your medical record are complete (i.e., not limited to only the first occurrence for each code) when using this parameter. Defaults to `1`.
- `date_fmt` (`str`): The format of the date fields in your phenotype and medical record data. Defaults to `'%Y-%m-%d'`.
- `phecode_version` (`str`): The version of the phecode system used for converting diagnosis codes. Version 1.2 is the official version of the phecode system, with mapping files available for ICD-9-CM, ICD-9-WHO, ICD-10-CM, and ICD-10-WHO codes. While option 1.3a is provided, it’s an unofficial version and not recommended for general use. Defaults to `'1.2'`.

---

### Instance Methods

#### `phenotype_data`

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

#### `Table1`

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

#### `merge_medical_records`

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

Load one or more medical record datasets.

**Parameters:**

- `medical_records_data_path` (`str`): Path to CSV/TSV medical record file.
- `diagnosis_code` (`str`): Code type: `'ICD-9-CM'`, `'ICD-9-WHO'`, `'ICD-10-CM'`, or `'ICD-10-WHO'`.
- `column_names` (`dict`): Mapping for dataset columns. Required keys: `'Participant ID'`, `'Diagnosis code'`, `'Date of diagnosis'`.
- `date_fmt` (`str`): Date format (defaults to phenotype data format). Defaults to `None`.
- `chunksize` (`int`): Rows per chunk for large files. Defaults to `1000000`.

**Returns:**

- `None`

---

#### `get_attribute`

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

#### `medical_records_to_dataframe`

```python
concat(
    self, 
    phecode_list: list,
    medical_history: bool=False
) -> DiseaseNetworkData
```

Convert stored medical record into a tidy pandas DataFrame.

**Parameters:**

- `phecode_list` (`list`): List of phecodes to extract from the medical record. Only phecodes valid for the current phecode_level are accepted.
- `medical_history` (`bool`): Include a binary history column for each phecode if set to True. Default to `False`
      

**Returns:**

- `pd.DataFrame`

---

#### `modify_phecode_level`

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

#### `disease_pair`

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

#### `save`

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

#### `load`

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

#### `save_npz`

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

#### `load_npz`

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

## Analysis Functions

### Function: `disease_network_pipeline`

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

### Function: `phewas`

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

### Function: `phewas_multipletests`

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

### Function: `comorbidity_strength`

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


### Function: `comorbidity_strength_multipletests`

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

### Function: `binomial_test`

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


### Function: `binomial_multipletests`

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


### Function: `comorbidity_network`

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


### Function: `comorbidity_multipletests`

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


### Function: `disease_trajectory`

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

### Function: `trajectory_multipletests`

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

## Class `Plot`

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

### Instance Methods

#### `three_dimension_plot`

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

#### `comorbidity_network_plot`

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

#### `trajectory_plot`

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

#### `phewas_plot`

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