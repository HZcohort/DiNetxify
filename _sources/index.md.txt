### ***DiNetxify***

![DiNetxify logo](./img/DiNetxify-logo.png){width=400px}

***DiNetxify*** is an open-source Python package for comprehensive three-dimensional (3D) disease network analysis of large-scale electronic health record (EHR) data. It integrates data harmonization, statistical analysis, and visualization in a single workflow to uncover multimorbidity patterns and disease progression pathways. The package supports `cohort`, `matched cohort`, and `exposed-only cohort` study designs, provides both step-by-step and one-step analysis workflows, and is released under the GPL-3.0 license.

![DiNetxify framework](./img/framework.png){width=1200px}

## Installation

***DiNetxify*** requires **Python 3.10+**. Install the latest release from PyPI with:

```bash
pip install dinetxify
```

Core dependencies include `numpy`, `pandas`, `scipy`, `statsmodels`, `lifelines`, `matplotlib`, `plotly`, `networkx`, `python_louvain`, `scikit_learn`, and `tqdm`.

## Source code and issue report

Source code, release history, and issue tracking are available on GitHub: [HZcohort/DiNetxify](https://github.com/HZcohort/DiNetxify).

## Citation

If you use this software in your research, please cite:

1. [DiNetxify: a python package for three-dimensional disease network analysis based on electronic health record data](https://link.springer.com/article/10.1007/s10654-025-01360-4) ([PMID: 41579291](https://pubmed.ncbi.nlm.nih.gov/41579291/))
2. [Disease clusters and their genetic determinants following a diagnosis of depression: analyses based on a novel three-dimensional disease network approach](https://www.nature.com/articles/s41380-025-03120-y) ([PMID: 40681841](https://pubmed.ncbi.nlm.nih.gov/40681841/))

## Contact

- **Can Hou**: [houcan@wchscu.cn](mailto:houcan@wchscu.cn)
- **Haowen Liu**: [haowenliu81@gmail.com](mailto:haowenliu81@gmail.com)

```{toctree}
:maxdepth: 2
:caption: Documentation

data_prep
data_harm
3d_analysis
visual
table
api
log
```
