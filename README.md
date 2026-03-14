# Borborema Magmatism Toolkit

[![Python](https://img.shields.io/badge/Python-3.9%2B-blue.svg)]()
[![License](https://img.shields.io/badge/License-MIT-green.svg)]()
[![DOI](https://zenodo.org/badge/DOI/10.xxxx/zenodo.xxxxxx.svg)]()

Python toolkit for geochemical, isotopic, and geochronological analysis of magmatic systems in the Borborema Province (NE Brazil).

This toolkit was developed to support integrated studies of whole-rock geochemistry, radiogenic isotopes, and zircon U–Pb geochronology in post-collisional granitoids.

---

## Overview

The toolkit integrates analytical workflows commonly used in igneous petrology and isotope geochemistry, including:

- whole-rock PCA
- AFC modelling
- Sm–Nd isotopic panels
- zircon U–Pb age plots

Applications include studies of the magmatic evolution of the Rio Grande do Norte Domain in the Borborema Province.

---

## Repository structure

```text
borborema-magmatism-toolkit/
├── README.md
├── requirements.txt
├── run_all_figures.py
│
├── sample_data/
│   ├── wr_pca_exemple.csv
│   ├── afc_model_exemple.csv
│   ├── sr_nd_models_exemple.csv
│   └── upb_geochronology_exemple.csv
│
├── figures/
│   ├── pca_plot.png
│   ├── afc_model.png
│   ├── sm_nd_panel.png
│   └── upb_ages.png
│
└── borborema/
    ├── __init__.py
    ├── wr_pca.py
    ├── afc_model.py
    ├── sr_nd_models.py
    ├── upb_geochronology.py
    ├── datasets.py
    └── data_cleaning.py
