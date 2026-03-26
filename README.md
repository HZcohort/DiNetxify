<div align="center">
  <img src="https://github.com/HZcohort/DiNetxify/blob/main/docs/img/DiNetxify-logo.png"
       alt="DiNetxify Logo"
       width="300">
</div>

--------------------------------------------------------------------------------
## About *DiNetxify*

***DiNetxify*** is an open-source Python package for three-dimensional (3D) disease network analysis of large-scale electronic health record (EHR) data. It integrates data harmonization, analysis, and visualization in a single workflow for studying multimorbidity patterns and disease progression pathways. The package supports `cohort`, `matched cohort`, and `exposed-only cohort` designs, and is released under the GPL-3.0 license.

![analytical framework](https://github.com/HZcohort/DiNetxify/blob/main/docs/img/framework.png)

***DiNetxify*** provides:

- **Integrated workflow:** from cohort data and diagnosis records to analysis outputs and plots.
- **Flexible study designs:** support for standard cohorts, matched cohorts, and exposed-only cohorts.
- **Modular analysis:** use the one-step pipeline or run PheWAS, comorbidity, and trajectory analyses separately.
- **Built-in visualization:** generate static and interactive figures directly from the result tables.

![architecture](https://github.com/HZcohort/DiNetxify/blob/main/docs/img/architecture.png)

## Installation

***DiNetxify*** requires **Python 3.10+**.

```bash
pip install dinetxify
```

Core dependencies include `numpy`, `pandas`, `matplotlib`, `plotly`, `python_louvain`, `networkx`, `scikit_learn`, `scipy`, `statsmodels>=0.14.4`, `lifelines>=0.27.0`, and `tqdm`.

## Quick Start

### 1. Load phenotype data and medical records

You can try the package with the dummy data under [`tests/data`](./tests/data).

```python
import DiNetxify as dnt

col_dict = {
    "Participant ID": "ID",
    "Exposure": "exposure",
    "Sex": "sex",
    "Index date": "date_start",
    "End date": "date_end",
}

covariates = ["age", "BMI"]

data = dnt.DiseaseNetworkData(
    study_design="cohort",
    phecode_level=1,
    date_fmt="%Y-%m-%d",
)

data.phenotype_data(
    phenotype_data_path="tests/data/dummy_phenotype.csv",
    column_names=col_dict,
    covariates=covariates,
)

data.merge_medical_records(
    medical_records_data_path="tests/data/dummy_EHR_ICD9.csv",
    diagnosis_code="ICD-9-WHO",
    column_names={
        "Participant ID": "ID",
        "Diagnosis code": "diag_icd9",
        "Date of diagnosis": "dia_date",
    },
)

data.merge_medical_records(
    medical_records_data_path="tests/data/dummy_EHR_ICD10.csv",
    diagnosis_code="ICD-10-WHO",
    column_names={
        "Participant ID": "ID",
        "Diagnosis code": "diag_icd10",
        "Date of diagnosis": "dia_date",
    },
)
```

### 2. Run the one-step analysis pipeline

`disease_network_pipeline()` returns five result tables:

- `phewas_result`
- `com_strength_result`
- `com_network_result`
- `binomial_result`
- `trajectory_result`

Example:

```python
from DiNetxify import disease_network_pipeline

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
        project_prefix="my_analysis",
        keep_positive_associations=True,
        method="RPCN",
        covariates=["age", "BMI"],
        matching_var_dict={"sex": "exact"},
        matching_n=2,
        correction="bonferroni",
        cutoff=0.05,
    )
```

Notes:

- `output_dir` must already exist.
- When using multiprocessing, keep the call inside `if __name__ == "__main__":`.
- `method` can be `'RPCN'`, `'PCN_PCA'`, or `'CN'`.

### 3. Visualize the results

```python
from DiNetxify.visualization import Plot

plot = Plot(
    phewas_result=phewas_result,
    comorbidity_result=com_network_result,
    trajectory_result=trajectory_result,
    exposure_name="Exposure",
    exposure_location=(0, 0, 0),
    exposure_size=15,
)

plot.phewas_plot("results/phewas_plot.png")
plot.comorbidity_network_plot("results/comorbidity_network.html")
plot.three_dimension_plot("results/three_dimension_network.html")
```

For exposed-only cohorts, set `exposure_name=None`, `exposure_location=None`, and `exposure_size=None` when creating `Plot`.

## Documentation

Full documentation is available at:

[https://hzcohort.github.io/DiNetxify/](https://hzcohort.github.io/DiNetxify/)

It includes guides for:

- data preparation and harmonization
- one-step and step-by-step 3D analysis
- visualization
- table generation
- API reference

## Citation

If you use this software in your research, please cite:

1. [DiNetxify: a python package for three-dimensional disease network analysis based on electronic health record data](https://link.springer.com/article/10.1007/s10654-025-01360-4) ([PMID: 41579291](https://pubmed.ncbi.nlm.nih.gov/41579291/))
2. [Disease clusters and their genetic determinants following a diagnosis of depression: analyses based on a novel three-dimensional disease network approach](https://www.nature.com/articles/s41380-025-03120-y) ([PMID: 40681841](https://pubmed.ncbi.nlm.nih.gov/40681841/))

## Contact

- Can Hou: [houcan@wchscu.cn](mailto:houcan@wchscu.cn)
- Haowen Liu: [haowenliu81@gmail.com](mailto:haowenliu81@gmail.com)
