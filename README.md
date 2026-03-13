# Borborema Magmatism Toolkit

[![Python](https://img.shields.io/badge/Python-3.9%2B-blue.svg)]()
[![License](https://img.shields.io/badge/License-MIT-green.svg)]()
[![DOI](https://zenodo.org/badge/DOI/10.xxxx/zenodo.xxxxxx.svg)]()

Python toolkit for geochemical, isotopic, and geochronological analysis of magmatic systems, developed for the study of Ediacaran–Cambrian granitoid magmatism in the Borborema Province (NE Brazil).

This toolkit was developed to support integrated studies of whole-rock geochemistry, radiogenic isotopes, and zircon U–Pb geochronology in post-collisional granitoids.


---

# Overview

The toolkit integrates multiple analytical workflows commonly used in igneous petrology and isotope geochemistry.

Applications include studies of the magmatic evolution of the Rio Grande do Norte Domain in the Borborema Province.

---

# Requirements

The toolkit requires:

- Python 3.9+
- numpy
- pandas
- matplotlib
- scikit-learn
- scipy

# Repository structure

borborema-magmatism-toolkit

```
├── README.md
├── requirements.txt
├── run_all_figures.py
│
├── examples
│   ├── wr_pca_example.csv
│   ├── afc_model_example.csv
│   ├── sr_nd_models_example.csv
│   └── upb_geochronology_example.csv
│
├── figures
│   ├── pca_plot.png
│   ├── afc_model.png
│   ├── sm_nd_evolution.png
│   └── upb_ages.png
│
└── borborema
    ├── __init__.py
    ├── wr_pca.py
    ├── afc_model.py
    ├── sr_nd_models.py
    ├── upb_geochronology.py
    ├── datasets.py
    └── data_cleaning.py
```

---

# Installation

Clone the repository:

```
git clone https://github.com/cataclase/borborema-magmatism-toolkit
cd borborema-magmatism-toolkit
```

Install dependencies:

```
pip install -r requirements.txt
```

---

# Whole-rock geochemistry – PCA

Principal Component Analysis allows identification of geochemical trends and clustering of magmatic suites.

![PCA](figures/pca_plot.png)

Example:

```
from borborema.wr_pca import run_pca
```

---

# Rb–Sr isotope modelling – AFC

Assimilation–Fractional Crystallization modelling following DePaolo (1981).

![AFC](figures/afc_model.png)

Example:

```
from borborema.afc_model import monte_carlo_afc
```

---

# Sm–Nd isotopic evolution

Visualization of εNd(t) versus age with comparison to the depleted mantle evolution curve.

![Sm-Nd](figures/sm_nd_evolution.png)

Example:

```
from borborema.sr_nd_models import plot_eNd_age
```

---

# U–Pb geochronology

Visualization of zircon crystallization ages and associated analytical uncertainties.

![U-Pb](figures/upb_ages.png)

Example:

```
from borborema.upb_geochronology import plot_upb_ages
```

---

# Reproducing all figures

All figures used in the examples can be generated automatically:

```
python run_all_figures.py
```

Figures will be saved in:

```
figures/
```

---

# Scientific context

The toolkit was developed to investigate:

• Post-collisional magmatism
• Crust–mantle interaction
• Magmatic differentiation processes
• Isotopic evolution of granitoid systems

The workflows implemented here were applied to magmatic suites of the Borborema Province.

---

# Code availability

All scripts used for geochemical modelling and figure generation are available in this repository.

---

# Author

Caio Tavares
Mine Geologist | Igneous Petrology | Isotope Geochemistry

---

# License

MIT License


