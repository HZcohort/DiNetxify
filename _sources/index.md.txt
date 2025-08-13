### ***DiNetxify***

![DiNetxify logo](./img/DiNetxify-logo.png){width=400px}

***DiNetxify*** is an open-source Python package for comprehensive three-dimensional (3D) disease network analysis of large-scale electronic health record (EHR) data. It integrates data harmonization, analysis, and visualization into a user-friendly package to uncover multimorbidity patterns and disease progression pathways. ***DiNetxify*** is optimized for efficiency (capable of handling cohorts of hundreds of thousands of patients within hours on standard hardware) and supports multiple study designs with customizable parameters and parallel computing. ***DiNetxify*** is released under GPL-3.0 license. 

![DiNetxify logo](./img/framework.png){width=1200px}

## Installation

***DiNetxify*** requires **Python 3.10+**. Install the latest release from PyPI using pip:

```bash
pip install dinetxify
```

## Source code and issue report

Available on Github, [HZcohort/DiNetxify](https://github.com/HZcohort/DiNetxify). Please report bugs and issues there. 

## Citation

If you use this software in your research, please cite the following papers:

1. [Disease clusters and their genetic determinants following a diagnosis of depression: analyses based on a novel three-dimensional disease network approach](https://www.nature.com/articles/s41380-025-03120-y) ([PMID: 40681841](https://pubmed.ncbi.nlm.nih.gov/40681841/))

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