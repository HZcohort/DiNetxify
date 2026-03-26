# Result tables description

## Result of PheWAS analysis

| Variable Name | Type | Description |
| --- | --- | --- |
| `phecode` | String/Float | Phecode included in the PheWAS result |
| `disease` | String | Disease name corresponding to the phecode |
| `system` | String | Phecode disease system corresponding to the phecode |
| `sex` | String | Sex specificity of the phecode (`Both`, `Male`, or `Female`) |
| `N_cases_exposed` | Integer | Number of exposed individuals diagnosed with the phecode during follow-up |
| `describe` | String | Model-fitting note or reason a full model estimate was not returned |
| `exposed_group` | String | Event summary for the exposed group, formatted as cases/person-years (incidence per 1,000 person-years) |
| `unexposed_group` | String | Event summary for the unexposed group, formatted as cases/person-years (incidence per 1,000 person-years) |
| `phewas_coef` | Float | Estimated effect size from the PheWAS model |
| `phewas_se` | Float | Standard error of the estimated effect size |
| `phewas_p` | Float | P-value from the PheWAS model |
| `phewas_p_significance` | Boolean | Whether the phecode is significant after multiple-testing correction |
| `phewas_p_adjusted` | Float | Adjusted P-value for `phewas_p` |

For `exposed-only cohort`, the returned table is simpler: significance is based on the case-count threshold, and model-based columns such as `phewas_coef`, `phewas_se`, `phewas_p`, and `phewas_p_adjusted` may be absent.

## Result of comorbidity strength estimation

| Variable Name | Type | Description |
| --- | --- | --- |
| `phecode_d1` | Float | Phecode for disease 1 in the disease pair |
| `phecode_d2` | Float | Phecode for disease 2 in the disease pair |
| `name_disease_pair` | String | Disease-pair label in the format `D1-D2` |
| `N_exposed` | Integer | Total number of exposed individuals in the analysis dataset |
| `n_total` | Integer | Number of exposed individuals eligible for the pair-specific analysis after excluding ineligible disease histories |
| `n_d1d2_diagnosis` | Integer | Number of eligible individuals diagnosed with both diseases |
| `n_d1_diagnosis` | Integer | Number of eligible individuals diagnosed with disease 1 |
| `n_d2_diagnosis` | Integer | Number of eligible individuals diagnosed with disease 2 |
| `n_d1d2_nontemporal` | Integer | Number of individuals with a non-temporal D1-D2 pair |
| `n_d1d2_temporal` | Integer | Number of individuals with temporal order `D1 -> D2` |
| `n_d2d1_temporal` | Integer | Number of individuals with temporal order `D2 -> D1` |
| `n_d1d2_pair` | Integer | Total number of individuals contributing a D1/D2 pair, counting temporal and non-temporal pair occurrences together |
| `description` | String | Note describing why statistics were not estimated, if applicable |
| `phi` | Float | Phi correlation for the disease pair |
| `phi_theta` | Float | Standard-error term used for inference on phi |
| `phi_p` | Float | P-value for phi correlation |
| `RR` | Float | Relative risk for co-occurrence of the disease pair |
| `RR_theta` | Float | Standard-error term used for inference on RR |
| `RR_p` | Float | P-value for relative risk |
| `phi_p_significance` | Boolean | Whether `phi_p` is significant after correction |
| `phi_p_adjusted` | Float | Adjusted P-value for `phi_p` |
| `RR_p_significance` | Boolean | Whether `RR_p` is significant after correction |
| `RR_p_adjusted` | Float | Adjusted P-value for `RR_p` |
| `disease_d1` | String | Disease name for `phecode_d1` |
| `system_d1` | String | Disease system for `phecode_d1` |
| `sex_d1` | String | Sex specificity for `phecode_d1` |
| `disease_d2` | String | Disease name for `phecode_d2` |
| `system_d2` | String | Disease system for `phecode_d2` |
| `sex_d2` | String | Sex specificity for `phecode_d2` |

## Result of binomial test

| Variable Name | Type | Description |
| --- | --- | --- |
| `phecode_d1` | Float | Phecode for disease 1 in the temporal disease pair |
| `phecode_d2` | Float | Phecode for disease 2 in the temporal disease pair |
| `name_disease_pair` | String | Temporal disease-pair label, typically `D1->D2` |
| `n_d1d2_nontemporal` | Integer | Number of individuals with a non-temporal D1-D2 pair |
| `n_d1d2_temporal` | Integer | Number of individuals with temporal order `D1 -> D2` |
| `n_d2d1_temporal` | Integer | Number of individuals with temporal order `D2 -> D1` |
| `binomial_p` | Float | P-value from the binomial test for directionality |
| `binomial_proportion` | Float | Estimated directionality proportion from the binomial test |
| `binomial_proportion_ci` | String | Confidence interval for the directionality proportion |
| `binomial_p_significance` | Boolean | Whether `binomial_p` is significant after correction |
| `binomial_p_adjusted` | Float | Adjusted P-value for `binomial_p` |
| `disease_d1` | String | Disease name for `phecode_d1` |
| `system_d1` | String | Disease system for `phecode_d1` |
| `sex_d1` | String | Sex specificity for `phecode_d1` |
| `disease_d2` | String | Disease name for `phecode_d2` |
| `system_d2` | String | Disease system for `phecode_d2` |
| `sex_d2` | String | Sex specificity for `phecode_d2` |

## Result of comorbidity network analysis

| Variable Name | Type | Description |
| --- | --- | --- |
| `phecode_d1` | Float | Phecode for disease 1 in the non-temporal disease pair |
| `phecode_d2` | Float | Phecode for disease 2 in the non-temporal disease pair |
| `name_disease_pair` | String | Disease-pair label in the format `D1-D2` |
| `N_exposed` | Integer | Total number of exposed individuals in the analysis dataset |
| `n_total` | Integer | Number of individuals included in the pair-specific regression analysis |
| `n_exposed/n_cases` | String | Exposure summary among cases for the pair-specific model |
| `n_exposed/n_controls` | String | Exposure summary among controls for the pair-specific model |
| `comorbidity_network_method` | String | Method used for comorbidity network analysis (`CN`, `RPCN`, or `PCN_PCA`) |
| `describe` | String | Model-fitting description, including removed covariates when relevant |
| `co_vars_list` | String | Covariates retained in the fitted model |
| `co_vars_zvalues` | String | Z-values for covariates retained in the fitted model |
| `comorbidity_beta` | Float | Estimated effect size from the comorbidity network model |
| `comorbidity_se` | Float | Standard error of the estimated effect size |
| `comorbidity_p` | Float | P-value for the comorbidity network effect |
| `comorbidity_aic` | Float | Akaike information criterion of the fitted model |
| `comorbidity_p_significance` | Boolean | Whether `comorbidity_p` is significant after correction |
| `comorbidity_p_adjusted` | Float | Adjusted P-value for `comorbidity_p` |
| `disease_d1` | String | Disease name for `phecode_d1` |
| `system_d1` | String | Disease system for `phecode_d1` |
| `sex_d1` | String | Sex specificity for `phecode_d1` |
| `disease_d2` | String | Disease name for `phecode_d2` |
| `system_d2` | String | Disease system for `phecode_d2` |
| `sex_d2` | String | Sex specificity for `phecode_d2` |

**Method-specific columns**

| Variable Name | Type | Description |
| --- | --- | --- |
| `alpha` | Float | Present for `RPCN`; the L1 penalty level used in the final model |
| `pc_sum_variance_explained` | Float | Present for `PCN_PCA`; cumulative variance explained by retained principal components |

## Result of disease trajectory analysis

| Variable Name | Type | Description |
| --- | --- | --- |
| `phecode_d1` | Float | Phecode for disease 1 in the temporal disease pair |
| `phecode_d2` | Float | Phecode for disease 2 in the temporal disease pair |
| `name_disease_pair` | String | Temporal disease-pair label, typically `D1->D2` |
| `N_exposed` | Integer | Total number of exposed individuals in the analysis dataset |
| `n_total` | Integer | Number of individuals included in the nested case-control analysis for the pair |
| `n_exposed/n_cases` | String | Exposure summary among sampled cases |
| `n_exposed/n_controls` | String | Exposure summary among sampled controls |
| `trajectory_method` | String | Method used for disease trajectory analysis (`CN`, `RPCN`, or `PCN_PCA`) |
| `describe` | String | Model-fitting description, including removed covariates when relevant |
| `co_vars_list` | String | Covariates retained in the fitted model |
| `co_vars_zvalues` | String | Z-values for covariates retained in the fitted model |
| `trajectory_beta` | Float | Estimated effect size from the trajectory model |
| `trajectory_se` | Float | Standard error of the estimated effect size |
| `trajectory_p` | Float | P-value for the trajectory effect |
| `trajectory_aic` | Float | Akaike information criterion of the fitted model |
| `trajectory_p_significance` | Boolean | Whether `trajectory_p` is significant after correction |
| `trajectory_p_adjusted` | Float | Adjusted P-value for `trajectory_p` |
| `disease_d1` | String | Disease name for `phecode_d1` |
| `system_d1` | String | Disease system for `phecode_d1` |
| `sex_d1` | String | Sex specificity for `phecode_d1` |
| `disease_d2` | String | Disease name for `phecode_d2` |
| `system_d2` | String | Disease system for `phecode_d2` |
| `sex_d2` | String | Sex specificity for `phecode_d2` |

**Method-specific columns**

| Variable Name | Type | Description |
| --- | --- | --- |
| `alpha` | Float | Present for `RPCN`; the L1 penalty level used in the final model |
| `pc_sum_variance_explained` | Float | Present for `PCN_PCA`; cumulative variance explained by retained principal components |
