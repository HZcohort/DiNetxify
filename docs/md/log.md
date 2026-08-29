# Changelog

**0.2.0 - 2026-08-29**

- Added `Plot.interactive_website()` for generating self-contained offline websites from available PheWAS, comorbidity-network, disease-trajectory, and 3D network results.

**0.1.17 - 2026-08-26**

- Added `method` and `maxiter` parameters to `phewas()` for configuring the `statsmodels` Cox optimizer.
- Fixed forwarding of method-specific RPCN and PCN-PCA parameters in `disease_network_pipeline()`.
- Updated the analysis API documentation to distinguish the PheWAS optimizer method from the pipeline network-analysis method.


**0.1.15 - 2026-07-23**

- enhance convergence warning handling in Cox model functions.


**0.1.14 - 2026-05-09**

- Fixed inconsistent module assignments across comorbidity, trajectory, and 3D network plots when both comorbidity and trajectory results are provided.
- Improved multiprocessing stability for PheWAS, comorbidity network, and disease trajectory model-fitting steps on Linux servers.
- Added `multiprocessing_start_method` support for model-fitting steps and improved progress reporting for parallel runs.


**0.1.13 - 2026-03-30**

- Updated the `Plot` module so PheWAS, comorbidity, and trajectory plots can be created with only the result tables each plot type requires.


**0.1.12 - 2026-03-11**

- Fixed a bug in the visualization module.
- Added the `get_louvain_clusters` function to the visualization module.
- Improved the wording and clarity of medical record summary text.


**0.1.11 - 2026-02-11**

- Fixed a failure in `compute_vif_sm_exact` when multi-column dependency issues were encountered.


**0.1.10 - 2025-12-09** 

- Added the ability to highlight connected edges when a node is clicked in the 3D figure.


**0.1.9 - 2025-11-13**

- Fixed incorrect parameter descriptions in the visualization module.


**0.1.8 - 2025-10-23**

- Fixed a bug in the visualization module.


**0.1.7 - 2025-10-15**

- Added the correct node legend to the comorbidity network plot.
- Updated phecode system labels in the plots.


**0.1.6 - 2025-10-06**

- Fixed several bugs in the visualization module.


**0.1.5 - 2025-10-04**

- Added support for phecode version 1.3a.


**0.1.4 - 2025-10-02**

- Fixed a bug in the comorbidity network estimation step that could produce imaginary numbers in the phi correlation.


**0.1.3 - 2025-08-28**

- Fixed a bug in `visualization.py`.


**0.1.0 - 2025-08-13**

- First stable version released.
