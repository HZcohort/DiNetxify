# Three-dimensional disease network analysis

Once the data is prepared and stored in a `DiseaseNetworkData` object, DiNetxify offers two approaches to perform the disease network analysis:

1. **One-step analysis:** a comprehensive pipeline that automates the entire sequence of analyses (PheWAS → disease pair generation → comorbidity strength estimation → binomial test → comorbidity network analysis → disease trajectory analysis) with one function call. This is convenient and ensures all steps are performed in the correct order with default or specified parameters.
2. **Step-by-step analysis:** individual functions for each analysis component, allowing you to run and inspect each step separately. This approach offers more control and flexibility (e.g., to tweak parameters at each step or to examine intermediate results), at the expense of writing a bit more code.

Both approaches output their results as pandas DataFrames, which you can further analyze or export (to CSV/Excel, etc.) using pandas. The one-step pipeline minimizes redundant computations and code, but does not allow modifying certain internal parameters beyond what its arguments expose. The step-by-step approach is more verbose but lets you adjust and understand each phase of the analysis in detail. We’ll demonstrate both.

## One-step analysis

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

> **Note:** When using multiprocessing, multi-threading may not always close successfully, which can cause conflicts that significantly affect performance. We recommend disabling multi-threading with the following code (Linux):
>
> ```shell
> export OPENBLAS_NUM_THREADS=1
> export MKL_NUM_THREADS=1
> export BLIS_NUM_THREADS=1
> export OMP_NUM_THREADS=1
> export NUMEXPR_NUM_THREADS=1
> ```
>
> or the following code in Windows:
>
> ```powershell
> set OPENBLAS_NUM_THREADS=1
> set MKL_NUM_THREADS=1
> set BLIS_NUM_THREADS=1
> set OMP_NUM_THREADS=1
> set NUMEXPR_NUM_THREADS=1
> ```

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

## Step-by-step analysis

In many cases, you might want to run each part of the analysis separately — for example, to inspect intermediate outputs, adjust parameters for individual steps, or run alternative filtering between steps. ***DiNetxify*** allows this by exposing individual functions for each analysis stage. Below we illustrate a step-by-step approach performing the same overall analysis as the one-step pipeline above. We will reuse the `data` object already loaded.

### 3.2.1 PheWAS Analysis

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

### 3.2.2 Disease pair generation

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

### 3.2.3 Comorbidity strength estimation

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

### 3.2.4 Binomial test

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

### 3.2.5 Comorbidity network analysis

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

### 3.2.6 Disease trajectory analysis

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