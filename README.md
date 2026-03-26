# Borborema Magmatism Toolkit

[![Python](https://img.shields.io/badge/Python-3.9%2B-blue.svg)]()
[![License](https://img.shields.io/badge/License-MIT-green.svg)]()

Python toolkit for geochemical, isotopic, and geochronological analysis of post-collisional magmatic systems, developed for the study of Ediacaran–Cambrian granitoid magmatism in the Rio Grande do Norte Domain, Borborema Province (NE Brazil).

This toolkit supports integrated studies of:

- whole-rock geochemistry  
- radiogenic isotopes  
- zircon U–Pb geochronology  

---

## Reproducibility

All figures and results presented in this repository can be fully reproduced using the provided datasets and scripts.

---

# Conceptual basis

The codes apply classical igneous petrology models to quantify mantle partial melting, basement assimilation, and fractional crystallization processes.

The classification into six magmatic suites follows Nascimento et al. (2015) and represents the post-collisional magmatism of the Borborema Province during the final stages of West Gondwana assembly.

Assimilation–fractional crystallization (AFC) models follow DePaolo (1981), allowing quantification of mantle and crustal contributions.

Sm–Nd isotopic modelling uses CHUR (Jacobsen & Wasserburg, 1980) and depleted mantle evolution (DePaolo, 1984), while enriched mantle reservoirs (EM1, EM2) follow Zindler & Hart (1986).

---

# Geological context

## Borborema Province

![Figure 1](figures/BORBOREMA.png)

## Rio Grande do Norte Domain

![Figure 2](figures/DRN.png)

---

# Overview

- whole-rock PCA  
- AFC modelling  
- Sm–Nd isotopic modelling  
- zircon U–Pb geochronology  

---

# Requirements

Python 3.9+

```
numpy
pandas
matplotlib
scikit-learn
scipy
```

---

# Installation

```
git clone https://github.com/cataclase/borborema-magmatism-toolkit
cd borborema-magmatism-toolkit
pip install -r requirements.txt
pip install -e .
```

---

# Quick Start

```
python run_all_figures.py
```

---

# Repository structure

```
borborema-magmatism-toolkit/
├── README.md
├── requirements.txt
├── pyproject.toml
├── run_all_figures.py
├── examples/
│   └── Example_DRN.xlsx
├── figures/
│   ├── PCA.png
│   ├── AFC_SrSr.png
│   ├── SrNd.png
│   └── UPB_geochronology.png
└── borborema/
    ├── __init__.py
    ├── wr_pca.py
    ├── afc_model.py
    ├── sr_nd_models.py
    ├── upb_geochronology.py
    ├── datasets.py
    └── data_cleaning.py
```

---

# Author

Caio Tavares  

---

# License

MIT License
