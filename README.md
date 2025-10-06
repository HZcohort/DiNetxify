<div align="center">
  <img src="https://github.com/HZcohort/DiNetxify/blob/main/docs/img/DiNetxify-logo.png" 
       alt="DiNetxify Logo" 
       width="300">
</div>


--------------------------------------------------------------------------------
## About *DiNetxify*

***DiNetxify*** is an open-source Python package for comprehensive three-dimensional (3D) disease network analysis of large-scale electronic health record (EHR) data. It integrates data harmonization, analysis, and visualization into a user-friendly package to uncover multimorbidity patterns and disease progression pathways. ***DiNetxify*** is optimized for efficiency (capable of handling cohorts of hundreds of thousands of patients within hours on standard hardware) and supports multiple study designs with customizable parameters and parallel computing. ***DiNetxify*** is released under GPL-3.0 license. 
![analytical framework](https://github.com/HZcohort/DiNetxify/blob/main/docs/img/framework.png)

***DiNetxify*** provides an end-to-end solution for 3D disease network analysis, featuring:

- **Integrated Workflow:** From raw EHR data to results and plots. ***DiNetxify*** guides you through data preprocessing, sequential analyses, and interactive visualizations in one coherent framework.
- **Flexibility:** Supports various cohort study designs, including standard cohort, matched cohort, and exposed-only cohort, and offers numerous parameters to tailor the analysis (e.g. significance thresholds, methods for network construction, etc.).
- **User-Friendly API:** High-level functions (e.g. a one-step pipeline) reduce coding overhead, while modular components allow fine-grained control. A dedicated data class handles data loading, cleaning, and ICD code mapping (to phecodes) automatically.
- **Comprehensive Analyses:** Combines phenome-wide association studies (PheWAS), comorbidity network analysis, and disease trajectory analysis to identify meaningful disease clusters and temporal sequences concurrently.
- **Visualization:** Built-in plotting tools generate interactive 3D network visualizations and static plots for PheWAS results, comorbidity networks, and disease trajectories, facilitating intuitive exploration of findings.

![architecture](https://github.com/HZcohort/DiNetxify/blob/main/docs/img/architecture.png)


## Installation and Quick Start

### Installation
***DiNetxify*** requires **Python 3.10+**. Install the latest release from PyPI using pip:

```bash
pip install dinetxify
```

This will install ***DiNetxify*** along with its dependencies. The required dependencies include: numpy, pandas, matplotlib, plotly, python_louvain, networkx, scikit_learn, scipy, statsmodels (>=0.14.4), and lifelines (optional).
### Quick start
To begin using ***DiNetxify***:
1. **Install the package:** Use the pip command above to install ***DiNetxify*** in your environment (Linux or Windows).

2. **Initialize and load data:** Import ***DiNetxify*** and create a `DiseaseNetworkData` object with your chosen study design. Then load your cohort’s phenotype and medical records data into this object. The package will handle data validation and ICD-to-phecode mapping for you. You can download our test [dummy data](https://github.com/HZcohort/DiNetxify/tree/main/tests/data) and run the following code:

   ```python
   import DiNetxify as dnt
   
   # Define required columns and other covariates columns
   col_dict = {'Participant ID': 'ID','Exposure': 'exposure','Sex': 'sex','Index date': 'date_start','End date': 'date_end'}
   vars_lst = ['age', 'BMI']
   # Initialize the data object with study design and phecode level
   data = dnt.DiseaseNetworkData(study_design="cohort",phecode_level=1,date_fmt="%Y-%m-%d")
   # Load the phenotype CSV file into the data object
   data.phenotype_data(phenotype_data_path="dummy_phenotype.csv",column_names=col_dict,covariates=vars_lst)
   # Merge with the first medical records file (CSV)
   data.merge_medical_records(medical_records_data_path="dummy_EHR_ICD9.csv",diagnosis_code="ICD-9-WHO",column_names={'Participant ID':'ID','Diagnosis code':'diag_icd9','Date of diagnosis':'dia_date'})
   data.merge_medical_records(medical_records_data_path="dummy_EHR_ICD10.csv",diagnosis_code="ICD-10-WHO",column_names={'Participant ID':'ID','Diagnosis code':'diag_icd10','Date of diagnosis':'dia_date'})
   ```
   
   
   
3. **Run the analysis:** Utilize the high-level pipeline function to perform the entire 3D network analysis on your `DiseaseNetworkData`:

   ```python
   from DiNetxify import disease_network_pipeline
   
   # When using multiprocessing, ensure that the code is enclosed within the following block.
   # This prevents entering a never ending loop of new process creation.
   if __name__ == "__main__":
       results = disease_network_pipeline(data=data, n_process=4,
                                          n_threshold_phewas=100,
                                          n_threshold_comorbidity=100,
                                          output_dir="./results/",
                                          project_prefix="my_analysis")
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



For a detailed tutorial on using ***DiNetxify***, see our documentation at [https://hzcohort.github.io/DiNetxify/](https://hzcohort.github.io/DiNetxify/)

## Citation

If you use this software in your research, please cite the following papers:

1. [Disease clusters and their genetic determinants following a diagnosis of depression: analyses based on a novel three-dimensional disease network approach](https://www.nature.com/articles/s41380-025-03120-y) ([PMID: 40681841](https://pubmed.ncbi.nlm.nih.gov/40681841/))
1. [DiNetxify: a Python package for three-dimensional disease network analysis based on electronic health record data](https://www.medrxiv.org/content/10.1101/2025.08.19.25333629v1)

## Contact

- Can Hou: [houcan@wchscu.cn](mailto:houcan@wchscu.cn)
- Haowen Liu: [haowenliu81@gmail.com](mailto:haowenliu81@gmail.com)
