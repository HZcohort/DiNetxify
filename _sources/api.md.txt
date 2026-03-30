# API Reference

Below is a concise reference for ***DiNetxify*** classes and functions, based on the current package implementation.

## Class `DiseaseNetworkData`

```python
class DiseaseNetworkData(
    study_design: str = "cohort",
    phecode_level: int = 1,
    min_required_icd_codes: int = 1,
    date_fmt: str = "%Y-%m-%d",
    phecode_version: str = "1.2",
)
```

Create a data container for phenotype data, medical records, phecode mappings, and downstream pair-construction results.

**Parameters:**

- `study_design` (`str`): One of `"cohort"`, `"matched cohort"`, or `"exposed-only cohort"`.
- `phecode_level` (`int`): Phecode granularity, either `1` or `2`.
- `min_required_icd_codes` (`int`): Minimum mapped ICD count required before a phecode is counted as present.
- `date_fmt` (`str`): Default date format used for phenotype data and, unless overridden, medical-record data.
- `phecode_version` (`str`): One of `"1.2"` or `"1.3a"`. Version `"1.2"` is the recommended general-purpose option.

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
    force: bool = False,
) -> None
```

Load phenotype data into the object.

**Parameters:**

- `phenotype_data_path` (`str`): Path to a CSV or TSV phenotype file.
- `column_names` (`dict`): Mapping from DiNetxify-required field names to dataset columns. Required keys depend on `study_design`.
- `covariates` (`list`): Additional phenotype variables to load.
- `is_single_sex` (`bool`): Set to `True` if the cohort contains only one sex.
- `force` (`bool`): Overwrite existing phenotype and medical-record data if `True`.

**Returns:**

- `None`

---

#### `Table1`

```python
Table1(
    self,
    continuous_stat_mode: str = "auto",
) -> pd.DataFrame
```

Generate a Table 1 style summary of the phenotype data.

**Parameters:**

- `continuous_stat_mode` (`str`): One of `"auto"`, `"normal"`, or `"nonnormal"`.

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
    chunksize: int = 1000000,
    diagnosis_code_exclusion: list = [],
) -> None
```

Load and merge one medical-record file into the object.

**Parameters:**

- `medical_records_data_path` (`str`): Path to a CSV or TSV diagnosis file.
- `diagnosis_code` (`str`): One of `"ICD-9-CM"`, `"ICD-9-WHO"`, `"ICD-10-CM"`, or `"ICD-10-WHO"`.
- `column_names` (`dict`): Mapping for `"Participant ID"`, `"Diagnosis code"`, and `"Date of diagnosis"`.
- `date_fmt` (`str | None`): Date format for this file. If `None`, use the object's `date_fmt`.
- `chunksize` (`int`): Number of rows processed per chunk.
- `diagnosis_code_exclusion` (`list`): Diagnosis codes to exclude before phecode mapping.

**Returns:**

- `None`

---

#### `get_attribute`

```python
get_attribute(
    self,
    attr_name: str,
) -> Any
```

Access selected internal metadata from the object.

**Supported attribute names:**

- `warning_phenotype`
- `phenotype_statistics`
- `phenotype_info`
- `warning_medical_records`
- `medical_records_statistics`
- `medical_records_info`
- `module_dir`
- `significant_phecodes`

**Returns:**

- Requested value (`Any`)

---

#### `concat`

```python
concat(
    cls,
    first_data: "DiseaseNetworkData",
    second_data: "DiseaseNetworkData",
    duplicates: str = "raise",
) -> "DiseaseNetworkData"
```

Class method reserved for concatenating two `DiseaseNetworkData` objects.

**Current status:**

- Present in the API, but currently raises `NotImplementedError`.

---

#### `modify_phecode_level`

```python
modify_phecode_level(
    self,
    phecode_level: int,
) -> None
```

Switch between phecode level `1` and `2`. This cannot be done after trajectory data have already been generated.

**Returns:**

- `None`

---

#### `disease_pair`

```python
disease_pair(
    self,
    phewas_result: pd.DataFrame,
    min_interval_days: int = 0,
    max_interval_days: float = np.inf,
    force: bool = False,
    n_process: int = 1,
    **kwargs,
) -> None
```

Construct temporal and non-temporal disease pairs from significant PheWAS phecodes among exposed individuals.

**Parameters:**

- `phewas_result` (`pd.DataFrame`): Result table from `phewas()`.
- `min_interval_days` (`int | float`): Minimum gap required for a temporal D1 -> D2 relationship.
- `max_interval_days` (`int | float`): Maximum gap allowed before a pair is treated as non-temporal.
- `force` (`bool`): Overwrite existing trajectory data if `True`.
- `n_process` (`int`): Number of processes used for pair construction.
- `**kwargs`: Optional column-name overrides:
  - `phecode_col` (default `'phecode'`)
  - `significance_col` (default `'phewas_p_significance'`)

**Returns:**

- `None`

---

#### `medical_records_to_dataframe`

```python
medical_records_to_dataframe(
    self,
    phecode_list: list,
    medical_history: bool = False,
) -> pd.DataFrame
```

Export selected phecodes from the stored medical records into a participant-level DataFrame.

**Parameters:**

- `phecode_list` (`list`): Phecodes to extract.
- `medical_history` (`bool`): If `True`, add `<phecode>_history` indicators.

**Returns:**

- `pd.DataFrame`

---

#### `save`

```python
save(
    self,
    file: str,
) -> None
```

Save the object as a gzip-compressed pickle file. The `.pkl.gz` suffix is appended automatically if needed.

---

#### `load`

```python
load(
    self,
    file: str,
    force: bool = False,
) -> None
```

Load the object from a gzip-compressed pickle file.

---

#### `save_npz`

```python
save_npz(
    self,
    file: str,
) -> None
```

Save the object as a compressed NumPy archive. The output file is written as `.npz`.

---

#### `load_npz`

```python
load_npz(
    self,
    file: str,
    force: bool = False,
) -> None
```

Load the object from a compressed NumPy archive.

## Analysis Functions

For p-value correction arguments throughout the analysis API, DiNetxify accepts `'none'` or any method supported by `statsmodels.stats.multitest.multipletests`.

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
    pipeline_mode: str = "v1",
    method: str = "RPCN",
    covariates: list = None,
    matching_var_dict: dict = {"sex": "exact"},
    matching_n: int = 2,
    min_interval_days: int = 0,
    max_interval_days: float = np.inf,
    enforce_temporal_order: bool = False,
    correction: str = "bonferroni",
    cutoff=0.05,
    **kwargs,
) -> tuple
```

Run the main workflow:

`PheWAS -> disease_pair -> comorbidity_strength -> binomial/comorbidity_network -> disease_trajectory`

**Returns:**

- A 5-tuple in this order:
  - `phewas_result`
  - `com_strength_result`
  - `com_network_result`
  - `binomial_result`
  - `trajectory_result`

**Method-specific kwargs:**

- For `method="RPCN"`:
  - `auto_penalty` (`bool`, default `True`)
  - `alpha` (`float`, required if `auto_penalty=False`)
  - `alpha_range` (`tuple`, default `(1, 15)`)
  - `scaling_factor` (`float`, default `1`)
- For `method="PCN_PCA"`:
  - `n_PC` (`int`, default `5`)
  - `explained_variance` (`float`)

---

### Function: `phewas`

```python
phewas(
    data: DiseaseNetworkData,
    covariates: list = None,
    proportion_threshold: float = None,
    n_threshold: int = None,
    n_process: int = 1,
    correction: str = "bonferroni",
    cutoff: float = 0.05,
    system_inc: list = None,
    system_exl: list = None,
    phecode_inc: list = None,
    phecode_exl: list = None,
    log_file: str = None,
    lifelines_disable: bool = False,
) -> pd.DataFrame
```

Run a phecode-wide association scan.

**Notes:**

- `n_threshold` and `proportion_threshold` are mutually exclusive.
- For `cohort` and `matched cohort`, PheWAS fits Cox models.
- For `exposed-only cohort`, significance is based on the case-count threshold rather than a model-based p-value.

---

### Function: `phewas_multipletests`

```python
phewas_multipletests(
    df: pd.DataFrame,
    correction: str = "bonferroni",
    cutoff: float = 0.05,
) -> pd.DataFrame
```

Apply multiple-testing correction to the PheWAS result table.

---

### Function: `comorbidity_strength`

```python
comorbidity_strength(
    data: DiseaseNetworkData,
    proportion_threshold: float = None,
    n_threshold: int = None,
    n_process: int = 1,
    log_file: str = None,
    correction_phi: str = "bonferroni",
    cutoff_phi: float = 0.05,
    correction_RR: str = "bonferroni",
    cutoff_RR: float = 0.05,
) -> pd.DataFrame
```

Estimate disease-pair strength among exposed individuals using phi correlation and relative risk.

**Notes:**

- Requires disease pairs to have already been built with `DiseaseNetworkData.disease_pair()`.
- `n_threshold` and `proportion_threshold` are mutually exclusive.

---

### Function: `comorbidity_strength_multipletests`

```python
comorbidity_strength_multipletests(
    df: pd.DataFrame,
    correction_phi: str = "bonferroni",
    cutoff_phi: float = 0.05,
    correction_RR: str = "bonferroni",
    cutoff_RR: float = 0.05,
) -> pd.DataFrame
```

Apply multiple-testing correction to `phi_p` and `RR_p`.

---

### Function: `binomial_test`

```python
binomial_test(
    data: DiseaseNetworkData,
    comorbidity_strength_result: pd.DataFrame,
    comorbidity_network_result: pd.DataFrame = None,
    n_process: int = 1,
    log_file: str = None,
    correction: str = "bonferroni",
    cutoff: float = 0.05,
    enforce_temporal_order: bool = False,
    **kwargs,
) -> pd.DataFrame
```

Test whether one temporal direction is more common than the reverse direction for disease pairs that are significant in comorbidity strength.

**Parameters of interest:**

- `comorbidity_network_result`: Optional filter table. If supplied, only disease pairs retained by the network result are tested.
- `enforce_temporal_order`: If `True`, exclude non-temporal D1-D2 pairs when forming the binomial test counts.

**Notes:**

- Multiprocessing is currently disabled for this function.
- `**kwargs` can be used to override relevant input-column names.

---

### Function: `binomial_multipletests`

```python
binomial_multipletests(
    df: pd.DataFrame,
    correction: str = "bonferroni",
    cutoff: float = 0.05,
) -> pd.DataFrame
```

Apply multiple-testing correction to `binomial_p`.

---

### Function: `comorbidity_network`

```python
comorbidity_network(
    data: DiseaseNetworkData,
    comorbidity_strength_result: pd.DataFrame,
    binomial_test_result: pd.DataFrame = None,
    method: str = "RPCN",
    covariates: list = None,
    n_process: int = 1,
    log_file: str = None,
    correction: str = "bonferroni",
    cutoff: float = 0.05,
    **kwargs,
) -> pd.DataFrame
```

Fit pairwise non-temporal comorbidity models.

**Supported methods:**

- `'CN'`: correlation network
- `'RPCN'`: regularized partial correlation network
- `'PCN_PCA'`: partial correlation network with principal components

**Supported kwargs:**

- Column-name overrides for the input result tables:
  - `phecode_d1_col`
  - `phecode_d2_col`
  - `significance_phi_col`
  - `significance_RR_col`
  - `significance_binomial_col`
- Method-specific options:
  - `alpha`, `auto_penalty`, `alpha_range`, `scaling_factor`
  - `n_PC`, `explained_variance`
- `enforce_time_interval` (`bool`, default `True`)

---

### Function: `comorbidity_multipletests`

```python
comorbidity_multipletests(
    df: pd.DataFrame,
    correction: str = "bonferroni",
    cutoff: float = 0.05,
) -> pd.DataFrame
```

Apply multiple-testing correction to `comorbidity_p`.

---

### Function: `disease_trajectory`

```python
disease_trajectory(
    data: DiseaseNetworkData,
    comorbidity_strength_result: pd.DataFrame,
    binomial_test_result: pd.DataFrame,
    method: str = "RPCN",
    matching_var_dict: dict = {"sex": "exact"},
    matching_n: int = 2,
    max_n_cases: int = np.inf,
    global_sampling: bool = False,
    covariates: list = None,
    n_process: int = 1,
    log_file: str = None,
    correction: str = "bonferroni",
    cutoff: float = 0.05,
    **kwargs,
) -> pd.DataFrame
```

Fit temporal disease-trajectory models using nested case-control sampling.

**Parameters of interest:**

- `matching_var_dict` (`dict`): Matching criteria for trajectory sampling. Use `'exact'` for categorical variables; use a positive numeric tolerance for continuous variables.
- `matching_n` (`int`): Maximum number of matched controls per case.
- `max_n_cases` (`int | np.inf`): Optional cap on the number of D2 cases.
- `global_sampling` (`bool`): If `True`, sample once per unique D2 and fit separate D1 models within that sampled set.

**Supported kwargs:**

- Column-name overrides for the input result tables:
  - `phecode_d1_col`
  - `phecode_d2_col`
  - `significance_phi_col`
  - `significance_RR_col`
  - `significance_binomial_col`
- Method-specific options:
  - `alpha`, `auto_penalty`, `alpha_range`, `scaling_factor`
  - `n_PC`, `explained_variance`
- `enforce_time_interval` (`bool`, default `True`)

---

### Function: `trajectory_multipletests`

```python
trajectory_multipletests(
    df: pd.DataFrame,
    correction: str = "bonferroni",
    cutoff: float = 0.05,
) -> pd.DataFrame
```

Apply multiple-testing correction to `trajectory_p`.

## Class `Plot`

Import with:

```python
from DiNetxify.visualization import Plot
```

```python
class Plot(
    phewas_result: pd.DataFrame,
    comorbidity_result: pd.DataFrame | None = None,
    trajectory_result: pd.DataFrame | None = None,
    exposure_name: str | None = None,
    exposure_location: Tuple[float] | None = None,
    exposure_size: float | None = None,
    phecode_col: str = "phecode",
    disease_col: str = "disease",
    system_col: str = "system",
    phewas_number_col: str = "N_cases_exposed",
    phewas_coef_col: str = "phewas_coef",
    phewas_se_col: str = "phewas_se",
    source_col: str = "phecode_d1",
    target_col: str = "phecode_d2",
    disease_pair_col: str = "name_disease_pair",
    comorbidity_beta_col: str = "comorbidity_beta",
    trajectory_beta_col: str = "trajectory_beta",
    phewas_significance_col: str = "phewas_p_significance",
    comorbidity_significance_col: str = "comorbidity_p_significance",
    trajectory_significance_col: str = "trajectory_p_significance",
    **kwargs,
)
```

Create a visualization object from PheWAS results, with optional comorbidity-network and trajectory result tables.

**Required inputs:**

- `phewas_result` (`pd.DataFrame`)

**Optional network inputs:**

- `comorbidity_result` (`pd.DataFrame | None`): Required for `comorbidity_network_plot()`, `three_dimension_plot()`, and `trajectory_plot()`.
- `trajectory_result` (`pd.DataFrame | None`): Required for `three_dimension_plot()` and `trajectory_plot()`.

**Optional display inputs:**

- `exposure_name`: Name of the exposure node. Set to `None` for exposed-only analyses.
- `exposure_location`: 3D location of the exposure node.
- `exposure_size`: Marker size for the exposure node.

**Optional kwargs:**

- `SYSTEM`: Ordered list of disease systems to use in the legend and color mapping.
- `COLOR`: Colors corresponding to `SYSTEM`.

---

### Instance Methods

#### `three_dimension_plot`

```python
three_dimension_plot(
    self,
    path: str,
    max_radius: float = 180.0,
    min_radius: float = 35.0,
    line_color: str = "black",
    line_width: float = 1.0,
    size_reduction: float = 0.5,
    cluster_reduction_ratio: float = 1,
    layer_distance: float = 40.0,
    layout_width: float = 900.0,
    layout_height: float = 900.0,
    font_style: str = "Times New Roman",
    font_size: float = 15.0,
) -> None
```

Generate an interactive 3D HTML disease-network plot. Requires both `comorbidity_result` and `trajectory_result`.

---

#### `comorbidity_network_plot`

```python
comorbidity_network_plot(
    self,
    path: str,
    max_radius: float = 180.0,
    min_radius: float = 35.0,
    size_reduction: float = 0.5,
    cluster_reduction_ratio: float = 1,
    line_width: float = 1.0,
    line_color: str = "black",
    layer_distance: float = 40.0,
    font_style: str = "Times New Roman",
) -> None
```

Generate an interactive 2D HTML comorbidity-network plot. Requires `comorbidity_result`.

---

#### `trajectory_plot`

```python
trajectory_plot(
    self,
    path: str,
    dpi: float = 500,
) -> None
```

Generate trajectory plots as PNG files, one cluster per image. Requires both `comorbidity_result` and `trajectory_result`.

---

#### `phewas_plot`

```python
phewas_plot(
    self,
    path: str,
    system_font_size: float = 17,
    disease_font_size: float = 10,
    HR_max: float = 2,
    incident_number_max: int = None,
    exposed_only_cohort: bool = False,
    dpi: float = 200,
) -> None
```

Generate a circular PheWAS plot as a static image.
