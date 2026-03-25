# Borborema Magmatism Toolkit

[![Python](https://img.shields.io/badge/Python-3.9%2B-blue.svg)]()
[![License](https://img.shields.io/badge/License-MIT-green.svg)]()

Python toolkit for geochemical, isotopic, and geochronological analysis of post-collisional magmatic systems, developed for the study of Ediacaran–Cambrian granitoid magmatism in the Rio Grande do Norte Domain of Borborema Province (NE Brazil).

This toolkit was developed to support integrated studies of:

- whole-rock geochemistry
- radiogenic isotopes
- zircon U–Pb geochronology

---

# Conceptual basis

The codes listed above apply classical igneous petrology models to qualify and quantify the processes of mantle partial melting, basement assimilation, and fractional crystallization. The classification into six magmatic suites adopted here is based on Nascimento et al. (2015) and represents the post-collisional magmatism of the Borborema Province, exposed in the crystalline basement of the Rio Grande do Norte Domain. These rocks record the final stages of West Gondwana assembly during the Ediacaran.

Assimilation–fractional crystallization (AFC) models after DePaolo (1981) were applied to quantify the relative contributions of mantle melting and crustal assimilation in the magmatic suites of the Rio Grande do Norte Domain, using natural compositions as proxies for primitive magma and assimilated crustal material.

εNd(t) versus age models are constructed using isotopic evolution curves for CHUR (Jacobsen and Wasserburg, 1980) and the depleted mantle (DM) following DePaolo (1984). Theoretical compositional fields for enriched mantle reservoirs (EM1 and EM2) are incorporated based on the global mantle component framework of Zindler and Hart (1986), allowing discrimination between depleted, enriched, and crustally contaminated sources.


## Regional geological context

### Borborema Province

At the regional scale, the Borborema Province comprises a complex Neoproterozoic orogenic system structured by major shear zones, basement domains, supracrustal belts, and widespread granitoid plutonism. Within this framework, the Rio Grande do Norte Domain preserves one of the most expressive records of Ediacaran post-collisional magmatism, making it an ideal natural laboratory to investigate the relationships among mantle source enrichment, crustal assimilation, magma differentiation, and lithospheric evolution.

![Regional geological framework of the Borborema Province](figures/BORBOREMA.png)

*Regional geological framework of the Borborema Province in northeastern Brazil, showing its major tectonic domains, lineaments, basement inliers, and the distribution of Brasiliano-Pan Africano plutonism. The dashed rectangle indicates the Rio Grande do Norte Domain shown in detail below.*

### Rio Grande do Norte Domain

![Post-collisional magmatic suites of the Rio Grande do Norte Domain](figures/DRN.png)

*Geological map of the Rio Grande do Norte Domain highlighting the spatial distribution of the six post-collisional magmatic suites investigated in this workflow.*


# Overview

The toolkit integrates analytical workflows commonly used in igneous petrology and isotope geochemistry.

Applications include studies of the magmatic evolution of the **Rio Grande do Norte Domain** in the Borborema Province.

Implemented workflows include:

- whole-rock PCA
- AFC modelling
- Sm–Nd isotopic panels
- zircon U–Pb age distributions

---

# Requirements

Python 3.9+

Required packages:

```text
numpy
pandas
matplotlib
scikit-learn
scipy
```

Install all dependencies with:

```bash
pip install -r requirements.txt
```

---

# Repository structure

```
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
│   ├── PCA.png
│   ├── AFC_SrSr.png
│   ├── SrNd.png
│   └── UPB_geochronology.png
│
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

# Quick Start

Run all example workflows and reproduce the figures:

---

# Whole-rock geochemistry – PCA

Principal Component Analysis allows identification of geochemical trends and clustering of magmatic suites.

![PCA](figures/PCA.png)

Example:

```python
import pandas as pd
from borborema.wr_pca import run_pca_from_dataframe

df = pd.read_csv("sample_data/wr_pca_exemple.csv", encoding="utf-8-sig")
df.columns = df.columns.str.replace("ï»¿", "", regex=False)

variables = [
    "SiO2", "MgO", "CaO", "Na2O", "K2O",
    "Rb", "Ba", "Sr", "Nb", "Zr", "Y", "Th", "La", "Ce"
]

fig, results = run_pca_from_dataframe(
    df,
    variables=variables,
    series_col="Suite"
)

print("Explained variance:", results["explained_variance"])
print("Groups:", results["groups"])
```

---

# Rb–Sr isotope modelling – AFC

Assimilation–Fractional Crystallization modelling following DePaolo (1981).

![AFC](figures/AFC_SrSr.png)

Example:

```python
import pandas as pd
from borborema.afc_model import run_afc_from_dataframe

df = pd.read_csv("sample_data/afc_model_exemple.csv", encoding="latin1")

fig, results = run_afc_from_dataframe(
    df,
    Sr_m=1356,
    R_m=0.71008,
    Sr_c=96.8,
    R_c=0.7808,
    series_col="SERIE",
    iterations=10000,
    random_state=42
)

print("Best r:", results["best_r"])
print("Best D:", results["best_D"])
```

---

# Sm–Nd isotopic evolution

Visualization of εNd(t) versus age with comparison to the depleted mantle evolution curve.

![Sm-Nd](figures/SrNd.png)

Example:

```python
import pandas as pd
from borborema.sr_nd_models import run_sr_nd_from_dataframe

df = pd.read_csv("sample_data/sr_nd_models_exemple.csv", encoding="latin1")

fig, results = run_sr_nd_from_dataframe(df)

print(results)
```

---

# U–Pb geochronology

Visualization of zircon crystallization ages and associated analytical uncertainties.

![U-Pb](figures/UPB_geochronology.png)

Example:

```python
import pandas as pd
from borborema.upb_geochronology import plot_upb_ages

df = pd.read_csv("sample_data/upb_geochronology_exemple.csv", encoding="utf-8-sig")

fig = plot_upb_ages(df)
```

---

# Reproducing all figures

All figures used in the examples can be generated automatically:

```
python -m borborema.run_all_figures
```

Figures will be saved in:

```
figures/
```

---

# Testing in Google Colab

You can test the toolkit in Google Colab by uploading the repository files and running the example workflows directly.

Example:

```python
import pandas as pd
import matplotlib.pyplot as plt
from borborema.sr_nd_models import run_sr_nd_from_dataframe

df = pd.read_csv("sample_data/sr_nd_models_exemple.csv", encoding="latin1")

fig, results = run_sr_nd_from_dataframe(df)

print(results)
plt.show()
```

To generate all example figures in Colab:

```python
%cd /content
!python -m borborema.run_all_figures
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
# How to Cite

If you use this toolkit in scientific work, please cite:

Tavares, C. (2026)  
Borborema Magmatism Toolkit: Python tools for geochemical and isotopic analysis of granitoid systems.  
GitHub repository.  
https://github.com/cataclase/borborema-magmatism-toolkit

---
# License

MIT License
