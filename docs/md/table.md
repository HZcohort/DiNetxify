# Result tables description

## Result of PheWAS analysis

| Variable Name           | Type    | Description                                                  |
| ----------------------- | ------- | ------------------------------------------------------------ |
| `phecode`               | String  | Disease code (Phecode) used in PheWAS analysis               |
| `disease`               | String  | Disease name corresponding to the Phecode                    |
| `system`                | String  | Phecode disease system corresponding to the Phecode (e.g., infectious diseases) |
| `sex`                   | String  | Sex-specificity of the disease (e.g., Both, Male, Female)    |
| `N_cases_exposed`       | Integer | Number of individuals diagnosed with the disease in the exposed group |
| `describe`              | String  | Descriptions of the model fitting state and removed covariates with reasons |
| `exposed_group`         | String  | Incidence rate (unit: per 1,000 person-years) in the exposed group |
| `unexposed_group`       | String  | Incidence rate (unit: per 1,000 person-years) in the unexposed group |
| `phewas_coef`           | Float   | Estimated coefficient from the model                         |
| `phewas_se`             | Float   | Standard error of the estimated coefficient                  |
| `phewas_p`              | Float   | P-value indicating statistical significance of the coefficient |
| `phewas_p_significance` | Boolean | Indicates whether the result is statistically significant based on adjusted p-value (True/False) |
| `phewas_p_adjusted`     | Float   | Adjusted p-value accounting for multiple comparisons         |

## Result of comorbidity strength estimation

| Variable Name        | Type    | Description                                                  |
| -------------------- | ------- | ------------------------------------------------------------ |
| `phecode_d1`         | Integer | Phecode for disease 1 in the disease pair                    |
| `phecode_d2`         | Integer | Phecode for disease 2 in the disease pair                    |
| `name_disease_pair`  | String  | Name of the disease pair (format: "D1-D2")                   |
| `N_exposed`          | Integer | Total number of individuals in exposed group                 |
| `n_total`            | Integer | Number of exposed individuals included in the sub-cohort that meet the sex-specificity eligibility criteria for both diseases and after excluding those with history of either disease 1, disease 2, or related diseases |
| `n_d1d2_diagnosis`   | Integer | Number of individuals diagnosed with both diseases           |
| `n_d1_diagnosis`     | Integer | Number of individuals diagnosed with disease 1               |
| `n_d2_diagnosis`     | Integer | Number of individuals diagnosed with disease 2               |
| `n_d1d2_nontemporal` | Integer | Number of individuals diagnosed with both disease 1 and disease 2 but without defined temporal order (i.e., the time interval between the two diagnosis is smaller than or equal to `min_interval_days` or larger `max_interval_days`) |
| `n_d1d2_temporal`    | Integer | Number of individuals diagnosed with disease 1 followed by disease 2 in a defined temporal order (i.e., the time interval between the two diagnosis is larger than `min_interval_days` and smaller than or equal to`max_interval_days`) |
| `n_d2d1_temporal`    | Integer | Number of individuals diagnosed with disease 2 followed by disease 1 in a defined temporal order (i.e., the time interval between the two diagnosis is larger than `min_interval_days` and smaller than or equal to`max_interval_days`) |
| `phi_coef`           | Float   | Phi coefficient (φ), Pearson’s correlations for two binary variables |
| `phi_p`              | Float   | P-value for Phi coefficient significance                     |
| `RR`                 | Float   | Relative risk of observing both conditions in the same individual relative to expectation |
| `RR_p`               | Float   | P-value for relative risk                                    |
| `phi_p_adjusted`     | Float   | Adjusted P-value for Phi coefficient (multiple comparisons)  |
| `RR_p_adjusted`      | Float   | Adjusted P-value for relative risk (multiple comparisons)    |
| `phi_p_significance` | Boolean | Whether the Phi is statistically significant based on adjusted p-value |
| `RR_p_significance`  | Boolean | Whether the RR is statistically significant based on adjusted p-value |
| `disease_d1`         | String  | Name of disease 1                                            |
| `system_d1`          | String  | Phecode disease system related to disease 1                  |
| `sex_d1`             | String  | Sex-specificity of disease 1                                 |
| `disease_d2`         | String  | Name of disease 2                                            |
| `system_d2`          | String  | Phecode disease system related to disease 2                  |
| `sex_d2`             | String  | Sex-specificity of disease 2                                 |

## Result of binomial test

| Variable Name             | Type    | Description                                                  |
| ------------------------- | ------- | ------------------------------------------------------------ |
| `phecode_d1`              | Float   | Phecode for disease 1 in the temporal disease pair           |
| `phecode_d2`              | Float   | Phecode for disease 2 in the temporal disease pair           |
| `name_disease_pair`       | String  | Name of the temporal disease pair (e.g., D1->D2)             |
| `n_d1d2_nontemporal`      | Float   | Number of individuals diagnosed with both disease 1 and disease 2 but without defined temporal order (i.e., the time interval between the two diagnosis is smaller than or equal to `min_interval_days` or larger `max_interval_days`) |
| `n_d1d2_temporal`         | Float   | Number of individuals diagnosed with disease 1 followed by disease 2 in a defined temporal order (i.e., the time interval between the two diagnosis is larger than `min_interval_days` and smaller than or equal to`max_interval_days`) |
| `n_d2d1_temporal`         | Float   | Number of individuals diagnosed with disease 2 followed by disease 1 in a defined temporal order (i.e., the time interval between the two diagnosis is larger than `min_interval_days` and smaller than or equal to`max_interval_days`) |
| `binomial_p`              | Float   | P-value from the binomial test for directionality            |
| `binomial_proportion`     | Float   | Proportion of successful outcomes in the binomial test       |
| `binomial_proportion_ci`  | String  | Confidence interval for the binomial proportion              |
| `disease_d1`              | String  | Name of disease 1                                            |
| `system_d1`               | String  | Phecode disease system for disease 1                         |
| `sex_d1`                  | String  | Sex-specificity of disease 1                                 |
| `disease_d2`              | String  | Name of disease 2                                            |
| `system_d2`               | String  | Phecode disease system for disease 2                         |
| `sex_d2`                  | String  | Sex-specificity of disease 2                                 |
| `binomial_p_significance` | Boolean | Indicates whether the result is statistically significant based on adjusted p-value |
| `binomial_p_adjusted`     | Float   | Adjusted p-value for multiple comparisons                    |

## Result of comorbidity network analysis

| Variable Name                  | Type    | Description                                                  |
| ------------------------------ | ------- | ------------------------------------------------------------ |
| `phecode_d1`                   | Float   | Phecode for disease 1 in the non-temporal disease pair       |
| `phecode_d2`                   | Float   | Phecode for disease 2 in the non-temporal disease pair       |
| `name_disease_pair`            | String  | Name of the non-temporal disease pair (e.g., "D1-D2")        |
| `N_exposed`                    | Integer | Total number of individuals in exposed group                 |
| `n_total`                      | Integer | Number of exposed individuals included in the sub-cohort that meet the sex-specificity eligibility criteria for both diseases and after excluding those with history of either disease 1, disease 2, or related diseases |
| `n_exposed/n_cases`            | String  | Number of exposed individuals (individuals with diagnosis of  D1) among cases (individuals with diagnosis of D2) |
| `n_exposed/n_controls`         | String  | Number of exposed individuals (individuals with diagnosis of D1) among controls (individuals without diagnosis of D2) |
| `comorbidity_network_method`   | String  | Method used for comorbidity network analysis                 |
| `describe`                     | String  | Description of the model fitting, removed covariates in the model, and reasons for removal of covariates in the model |
| `co_vars_list`                 | String  | List of covariates used in the model                         |
| `co_vars_zvalues`              | String  | Z-values for each covariate in the model                     |
| `comorbidity_beta`             | Float   | Estimated coefficient from the comorbidity model             |
| `comorbidity_se`               | Float   | Standard error of the estimated coefficient                  |
| `comorbidity_p`                | Float   | P-value for the comorbidity coefficient                      |
| `comorbidity_aic`              | Float   | Akaike information criterion for the model                   |
| `disease_d1`                   | String  | Name of the disease 1                                        |
| `system_d1`                    | String  | Phecode disease system for the disease 1                     |
| `sex_d1`                       | String  | Sex-specificity of the disease 1                             |
| `disease_d2`                   | String  | Name of the disease 2                                        |
| `system_d2`                    | String  | Phecode disease system for the disease 2                     |
| `sex_d2`                       | String  | Sex-specificity of the disease 2                             |
| `comorbidity_p_significance`   | Boolean | Whether the result is statistically significant based on adjusted p-value |
| `comorbidity_p_adjusted`       | Float   | Adjusted p-value accounting for multiple comparisons         |
| **Columns for RPCN method**    |         |                                                              |
| `alpha`                        | Float   | Hyperparameter used for l1-norm (Weight multiplying the l1 penalty term) |
| **Columns for PCN_PCA method** |         |                                                              |
| `pc_sum_variance_explained`    | Float   | The cumulative proportion of variance that is accounted for by a selected number of principal components in a Principal Component Analysis (sum of explained variance for principal components) |

## Result of disease trajectory analysis

| Variable Name                  | Type    | Description                                                  |
| ------------------------------ | ------- | ------------------------------------------------------------ |
| `phecode_d1`                   | Float   | Phecode for disease 1 in the temporal disease pair           |
| `phecode_d2`                   | Float   | Phecode for disease 2 in the temporal disease pair           |
| `name_disease_pair`            | String  | Name of the temporal disease pair (e.g., "D1 → D2")          |
| `N_exposed`                    | Integer | Total number of individuals in exposed group                 |
| `n_total`                      | Integer | Number of exposed individuals included in the nested case-control dataset, where eligible cases (with diagnosis of D2) are all selected and matched with specified number of controls using incidence density sampling |
| `n_exposed/n_cases`            | String  | Number of exposed individuals (individuals with diagnosis of  D1) among cases (individuals with diagnosis of D2) |
| `n_exposed/n_controls`         | String  | Number of exposed individuals (individuals with diagnosis of  D1) among cases (individuals with diagnosis of D2) |
| `trajectory_method`            | String  | Method used for disease trajectory analysis                  |
| `describe`                     | String  | Description of the model fitting, removed covariates in the model, and reasons for removal of covariates in the model |
| `co_vars_list`                 | String  | List of covariates included in the model                     |
| `co_vars_zvalues`              | String  | Z-values for each covariate in the model                     |
| `trajectory_beta`              | Float   | Estimated coefficient from the model                         |
| `trajectory_se`                | Float   | Standard error of the estimated coefficient                  |
| `trajectory_p`                 | Float   | P-value for the coefficient                                  |
| `trajectory_aic`               | Float   | Akaike information criterion for the model                   |
| `disease_d1`                   | String  | Name of the disease 1                                        |
| `system_d1`                    | String  | Phecode disease system for the disease 1                     |
| `sex_d1`                       | String  | Sex-specificity of the disease 1                             |
| `disease_d2`                   | String  | Name of the disease 2                                        |
| `system_d2`                    | String  | Phecode disease system for the disease 2                     |
| `sex_d2`                       | String  | Sex-specificity of the disease 2                             |
| `trajectory_p_significance`    | Boolean | Whether the result is statistically significant based on adjusted p-value |
| `trajectory_p_adjusted`        | Float   | Adjusted p-value accounting for multiple comparisons         |
| **Columns for RPCN method**    |         |                                                              |
| `alpha`                        | Float   | Hyperparameter used for l1-norm (weight multiplying the l1 penalty term) |
| **Columns for PCN_PCA method** |         |                                                              |
| `pc_sum_variance_explained`    | Float   | The cumulative proportion of variance in a dataset that is accounted for by a selected number of principal components in a Principal Component Analysis (sum of explained variance for principal components) |