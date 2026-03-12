# Borborema Magmatism Toolkit

Python toolkit for geochemical and isotopic analysis of post-collisional magmatism in the **Rio Grande do Norte Domain (DRN), Borborema Province, NE Brazil**.

This repository contains scripts used to investigate the geochemical and isotopic evolution of Ediacaran–Cambrian granitoids.

## Modules

### Sr isotopes
Assimilation–Fractional Crystallization (AFC) modeling of Sr isotopes.

Python script for modeling **Sr isotope evolution using an Assimilation–Fractional Crystallization (AFC) model**.

The script generates a **Sr vs ⁸⁷Sr/⁸⁶Sr diagram** and estimates the best AFC parameters using a Monte Carlo approach.

This workflow was developed to investigate the evolution of **post-collisional granitoid magmatism from the Rio Grande do Norte Domain (Borborema Province, NE Brazil)**.

The resulting figures follow graphical standards commonly used in the journal *Lithos*.

---

# Features

* Monte Carlo AFC parameter estimation
* Automatic recognition of dataset column names
* Plot of:

  * Magmatic suites
  * Archean crust
  * Paleoproterozoic crust
  * Parental magma
  * Country rock
  * Melt fraction along AFC trajectory
* Publication-quality figures

Figures are exported as:

PNG
SVG
PDF

---

# Installation

Clone the repository:

```bash
git clone https://github.com/YOUR_USERNAME/afc-sr-model.git
cd afc-sr-model
```

Install dependencies:

```bash
pip install -r requirements.txt
```

---

# Usage

Run the script with a dataset:

```bash
python afc_sr_model.py example_data.csv
```

or

```bash
python afc_sr_model.py example_data.xlsx
```

The script will output the best AFC parameters:

```
Best r = ...
Best D = ...
```

and generate the figure:

```
AFC_SrSr.png
AFC_SrSr.svg
AFC_SrSr.pdf
```

---

# Input Data Format

The dataset must contain the following columns:

| Column                   | Description                    |
| ------------------------ | ------------------------------ |
| SERIE / Serie            | Magmatic series classification |
| Sr ppm / Sr              | Strontium concentration        |
| 87Sr/86Sr(T) / 87Sr/86Sr | Strontium isotope ratio        |

Example:

```
Sample,SERIE,Sr ppm,87Sr/86Sr(T)
A1,PHKCalcAlk,512,0.7121
A2,Shos,689,0.7108
A3,Archean,120,0.885
```

---

# Output Example

The script produces a Sr–Sr isotope diagram showing:

* AFC best-fit trajectory
* melt fraction evolution
* magmatic suites
* crustal endmembers

---

# Citation

If you use this script in scientific work, please cite:

Tavares, C. (2026)
*AFC Sr–Sr modeling script for granitoid magmatic suites.*
GitHub repository.

---

# License

This project is distributed under the MIT License.

