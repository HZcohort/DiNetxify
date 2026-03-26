# Three-dimensional disease network analysis

Once your data have been loaded into a `DiseaseNetworkData` object, DiNetxify supports two ways to run the three-dimensional disease network analysis:

1. A one-step pipeline with `disease_network_pipeline()`
2. A step-by-step workflow using the individual analysis functions

Both workflows return `pandas.DataFrame` objects that you can inspect, export, or pass to the visualization module.

## One-step analysis

`disease_network_pipeline()` runs the main analysis workflow in order:

`PheWAS -> disease pair construction -> comorbidity strength -> binomial test / comorbidity network -> disease trajectory`

Two pipeline modes are available:

- `pipeline_mode='v1'`: run the binomial test before the comorbidity network analysis
- `pipeline_mode='v2'`: run the comorbidity network analysis before the binomial test, and restrict the binomial test to network-significant pairs

Example:

```python
from DiNetxify import disease_network_pipeline

# When using multiprocessing, wrap the call in the main guard.
if __name__ == "__main__":
    (
        phewas_result,
        com_strength_result,
        com_network_result,
        binomial_result,
        trajectory_result,
    ) = disease_network_pipeline(
        data=data,
        n_process=4,
        n_threshold_phewas=100,
        n_threshold_comorbidity=100,
        output_dir="./results",
        project_prefix="disease_network",
        keep_positive_associations=True,
        save_intermediate_data=False,
        system_exl=["symptoms", "others", "injuries & poisonings"],
        pipeline_mode="v1",
        method="RPCN",
        covariates=["age", "BMI"],
        matching_var_dict={"sex": "exact"},
        matching_n=2,
        min_interval_days=0,
        max_interval_days=float("inf"),
        enforce_temporal_order=False,
        correction="bonferroni",
        cutoff=0.05,
    )
```

Important notes:

- `output_dir` must already exist.
- `n_threshold_comorbidity` must be less than or equal to `n_threshold_phewas`.
- `keep_positive_associations=True` keeps PheWAS results with `phewas_coef >= 0` for non-exposed-only designs, and keeps comorbidity pairs with `phi > 0` and `RR > 1`.
- The pipeline writes result tables and log files to `output_dir`, in addition to returning the result DataFrames.
- If `save_intermediate_data=True`, the intermediate `DiseaseNetworkData` object created after disease-pair generation is saved with the prefix `project_prefix`.

### Core pipeline arguments

- `data`: A `DiseaseNetworkData` object with phenotype data and medical records already loaded.
- `n_process`: Number of processes for the parallelized steps.
- `n_threshold_phewas`: Minimum number of exposed cases required for a phecode to enter PheWAS.
- `n_threshold_comorbidity`: Minimum number of exposed individuals with a disease pair for comorbidity strength estimation.
- `system_exl`: Optional list of phecode systems to exclude.
- `method`: One of `'RPCN'`, `'PCN_PCA'`, or `'CN'`.
- `covariates`: Covariates used in PheWAS, comorbidity network analysis, and trajectory analysis.
- `matching_var_dict`: Matching rules for trajectory analysis. Use `'sex'` for sex.
- `matching_n`: Maximum number of matched controls per case in trajectory analysis.
- `min_interval_days`, `max_interval_days`: Time-window settings used when building disease pairs.
- `enforce_temporal_order`: Passed to `binomial_test()`.
- `correction`, `cutoff`: Multiple-testing correction method and significance threshold.

### Method-specific arguments

For `method="RPCN"`:

- `auto_penalty=True` by default
- `alpha_range=(1, 15)` by default when `auto_penalty=True`
- `scaling_factor=1` by default
- If `auto_penalty=False`, you must provide `alpha`

For `method="PCN_PCA"`:

- `n_PC=5` by default
- or use `explained_variance=<float between 0 and 1>` instead of `n_PC`

## Step-by-step analysis

Use the step-by-step workflow if you want to inspect intermediate outputs or customize each stage separately.

```python
import DiNetxify as dnt
```

### 1. PheWAS

`dnt.phewas()` identifies phecodes associated with the exposure.

- In `cohort` and `matched cohort` designs, it fits Cox models.
- In `exposed-only cohort`, it marks a phecode as significant when `N_cases_exposed` reaches the threshold.

Example:

```python
phewas_result = dnt.phewas(
    data=data,
    covariates=["age", "BMI"],
    n_threshold=100,
    n_process=4,
    correction="bonferroni",
    cutoff=0.05,
    system_exl=None,
    phecode_inc=None,
    phecode_exl=None,
    log_file=None,
    lifelines_disable=False,
)
```

Notes:

- `n_threshold` and `proportion_threshold` are mutually exclusive.
- If `covariates=None`, DiNetxify uses `'sex'` plus the covariates provided in `phenotype_data()`.
- In matched cohorts, including matching variables as covariates may cause singular-model issues.

If you want to reapply multiple-testing correction later:

```python
phewas_result = dnt.phewas_multipletests(
    df=phewas_result,
    correction="fdr_bh",
    cutoff=0.05,
)
```

### 2. Disease pair construction

`DiseaseNetworkData.disease_pair()` uses the significant phecodes from the PheWAS result to build temporal and non-temporal disease pairs among exposed individuals.

Example:

```python
data.disease_pair(
    phewas_result=phewas_result,
    min_interval_days=0,
    max_interval_days=float("inf"),
    force=True,
    n_process=4,
)
```

Notes:

- `force=True` is needed if `data.trajectory` already exists and you want to rebuild the pairs.
- By default, the method expects `phewas_result` to contain `phecode` and `phewas_p_significance`.
- `min_interval_days` and `max_interval_days` are stored on the data object and reused later by trajectory analysis.

### 3. Comorbidity strength

`dnt.comorbidity_strength()` evaluates disease-pair co-occurrence strength among exposed individuals using phi-correlation and relative risk (RR).

Example:

```python
com_strength_result = dnt.comorbidity_strength(
    data=data,
    n_threshold=100,
    n_process=4,
    correction_phi="bonferroni",
    cutoff_phi=0.05,
    correction_RR="bonferroni",
    cutoff_RR=0.05,
    log_file=None,
)
```

Notes:

- `n_threshold` and `proportion_threshold` are mutually exclusive.
- This function requires disease pairs to have already been generated.
- Downstream network and trajectory analyses only keep pairs significant for both phi and RR.

If needed, you can rerun the multiple-testing correction:

```python
com_strength_result = dnt.comorbidity_strength_multipletests(
    df=com_strength_result,
    correction_phi="fdr_bh",
    cutoff_phi=0.05,
    correction_RR="fdr_bh",
    cutoff_RR=0.05,
)
```

### 4. Binomial test

`dnt.binomial_test()` tests whether one temporal direction is more common than the reverse direction for each eligible disease pair.

Example:

```python
binomial_result = dnt.binomial_test(
    data=data,
    comorbidity_strength_result=com_strength_result,
    comorbidity_network_result=None,
    n_process=1,
    correction="bonferroni",
    cutoff=0.05,
    enforce_temporal_order=False,
    log_file=None,
)
```

Notes:

- Multiprocessing is disabled for this function, so `n_process` must remain `1`.
- If you pass `comorbidity_network_result`, the test is restricted to network-significant pairs. This mirrors pipeline mode `v2`.
- If `enforce_temporal_order=True`, non-temporal D1-D2 occurrences are excluded from the test.

If needed:

```python
binomial_result = dnt.binomial_multipletests(
    df=binomial_result,
    correction="fdr_bh",
    cutoff=0.05,
)
```

### 5. Comorbidity network analysis

`dnt.comorbidity_network()` estimates adjusted disease-disease associations for the comorbidity network.

Example:

```python
com_network_result = dnt.comorbidity_network(
    data=data,
    comorbidity_strength_result=com_strength_result,
    binomial_test_result=None,
    method="RPCN",
    covariates=["age", "BMI"],
    n_process=4,
    correction="bonferroni",
    cutoff=0.05,
    auto_penalty=True,
    alpha_range=(1, 15),
    scaling_factor=1,
    log_file=None,
)
```

Notes:

- `method` can be `'RPCN'`, `'PCN_PCA'`, or `'CN'`.
- The current implementation is driven by disease pairs that are significant for both phi and RR in `comorbidity_strength_result`.
- `binomial_test_result` is accepted by the function signature for compatibility with the shared column-handling logic, but the latest code does not use it to filter the comorbidity network analysis.
- `RPCN` and `PCN_PCA` add network-disease covariates automatically; `CN` uses only phenotypic covariates.

If needed:

```python
com_network_result = dnt.comorbidity_multipletests(
    df=com_network_result,
    correction="fdr_bh",
    cutoff=0.05,
)
```

### 6. Disease trajectory analysis

`dnt.disease_trajectory()` performs the final trajectory analysis on disease pairs that are significant in both comorbidity strength and the binomial test.

Example:

```python
trajectory_result = dnt.disease_trajectory(
    data=data,
    comorbidity_strength_result=com_strength_result,
    binomial_test_result=binomial_result,
    method="RPCN",
    matching_var_dict={"sex": "exact"},
    matching_n=2,
    max_n_cases=float("inf"),
    global_sampling=False,
    covariates=["age", "BMI"],
    n_process=4,
    correction="bonferroni",
    cutoff=0.05,
    enforce_time_interval=True,
    auto_penalty=True,
    alpha_range=(1, 15),
    scaling_factor=1,
    log_file=None,
)
```

Notes:

- `matching_var_dict` controls incidence-density sampling for the nested case-control design.
- Categorical or binary matching variables should usually not also be included as covariates.
- `enforce_time_interval=True` is the default and applies the interval rules defined earlier in `data.disease_pair()`.
- `global_sampling=True` reuses incidence-density sampling across pairs with the same D2 and can be useful for larger analyses.
- `max_n_cases` can cap the number of D2 cases analyzed for each model.

If needed:

```python
trajectory_result = dnt.trajectory_multipletests(
    df=trajectory_result,
    correction="fdr_bh",
    cutoff=0.05,
)
```

## Outputs

The main result tables are:

- `phewas_result`: PheWAS results for individual phecodes
- `com_strength_result`: Comorbidity strength results based on phi and RR
- `com_network_result`: Adjusted comorbidity network results
- `binomial_result`: Temporal-order test results
- `trajectory_result`: Final disease trajectory results

These outputs can be passed directly to `DiNetxify.visualization.Plot` for visualization.
