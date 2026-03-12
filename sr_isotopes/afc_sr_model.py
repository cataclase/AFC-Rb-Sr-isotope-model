
### Future modules

The toolkit will also include:

- Sm–Nd isotope analysis
- U–Pb geochronology visualization
- Whole-rock PCA analysis
- Geochemical diagrams

## Example figure

![AFC Sr-Sr Diagram](figures/AFC_SrSr.png)

## Installation

Install dependencies:

"""
AFC Sr–Sr model plotter
Author: Caio Tavares
Description:
Generates an AFC (Assimilation–Fractional Crystallization) Sr isotope diagram
with Monte Carlo best-fit modeling.

Usage:
    python afc_sr_model.py data.csv
    python afc_sr_model.py data.xlsx
"""

import sys
import numpy as np
import pandas as pd
import re
import unicodedata
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe
from matplotlib.lines import Line2D


# ====================== DATA CLEANING ======================

def strip_accents(s):
    s = unicodedata.normalize("NFKD", s)
    return "".join(ch for ch in s if not unicodedata.combining(ch))


def normalize_series(x):

    if pd.isna(x):
        return np.nan

    s = strip_accents(str(x)).lower()
    key = re.sub("[^a-z]", "", s)

    mapping = {
        "alk": "Alk",
        "phkcalcalk": "PHKCalcAlk",
        "ehkcalcalk": "EHKCalcAlk",
        "chalk": "ChAlk",
        "shos": "Shos",
        "archean": "Archean",
        "paleoproterozoic": "Paleoproterozoic",
    }

    return mapping.get(key, x)


def to_num(col):

    return pd.to_numeric(
        pd.Series(col, dtype=str)
        .str.replace(",", ".")
        .str.replace(" ", ""),
        errors="coerce",
    )


# ====================== LOAD DATA ======================

def load_data(file_path):

    if file_path.endswith(".csv"):
        df = pd.read_csv(file_path)

    elif file_path.endswith(".xlsx"):
        df = pd.read_excel(file_path)

    else:
        raise ValueError("File must be CSV or XLSX")

    return df


# ====================== AFC MODEL ======================

def afc_fit(df):

    Sr_m, R_m = 1356, 0.71008
    Sr_c, R_c = 96.8, 0.7808

    F = np.linspace(1, 0.01, 500)

    x_data = df["Sr"].values
    y_data = df["ratio"].values

    best_error = np.inf

    for _ in range(10000):

        r = np.random.uniform(0.01, 0.5)
        D = np.random.uniform(0.5, 3.0)

        Sr_model = (Sr_m * F**D + Sr_c * r * (1 - F)) / (F + r * (1 - F))
        R_model = (R_m * Sr_m * F**D + R_c * Sr_c * r * (1 - F)) / (
            Sr_model * (F + r * (1 - F))
        )

        dist = np.abs(R_model[:, None] - y_data) + np.abs(
            Sr_model[:, None] - x_data
        ) / 1000

        error = np.sum(np.min(dist, axis=0))

        if error < best_error:
            best_error = error
            best_r = r
            best_D = D
            best_curve = (Sr_model, R_model)

    return best_r, best_D, best_curve, Sr_m, R_m, Sr_c, R_c


# ====================== PLOT ======================

def make_plot(df, best_r, best_D, best_curve, Sr_m, R_m, Sr_c, R_c):

    palette = plt.cm.tab10.colors

    series_colors = {
        "PHKCalcAlk": palette[0],
        "EHKCalcAlk": palette[1],
        "ChAlk": palette[2],
        "Shos": palette[3],
        "Alk": palette[4],
    }

    fig, ax = plt.subplots(figsize=(9, 7), facecolor="white")
    ax.set_facecolor("white")

    # AFC curve
    ax.plot(best_curve[0], best_curve[1], color="black", linewidth=3.5)

    handles = []

    for serie, color in series_colors.items():

        sub = df[df["SERIE_norm"] == serie]

        if len(sub) > 0:

            sc = ax.scatter(
                sub["Sr"],
                sub["ratio"],
                s=55,
                color=color,
                edgecolor="black",
                linewidth=0.3,
                alpha=0.8,
                label=f"{serie} (n={len(sub)})",
            )

            handles.append(sc)

    # crust
    crust_symbols = {"Archean": "s", "Paleoproterozoic": "^"}

    for crust in crust_symbols:

        sub = df[df["SERIE_norm"] == crust]

        if len(sub) > 0:

            sc = ax.scatter(
                sub["Sr"],
                sub["ratio"],
                marker=crust_symbols[crust],
                facecolor="white",
                edgecolor="black",
                s=50,
                linewidth=0.9,
                label=f"{crust} (n={len(sub)})",
            )

            handles.append(sc)

    # Endmembers
    ax.scatter(Sr_m, R_m, color="black", s=120, marker="o")
    ax.scatter(Sr_c, R_c, color="red", s=120, marker="^")

    # Labels
    txt = ax.text(
        Sr_m,
        R_m + 0.0025,
        "Parental magma",
        fontsize=10,
        ha="center",
        va="bottom",
    )

    txt.set_path_effects([pe.withStroke(linewidth=3, foreground="white")])

    ax.text(Sr_c + 15, R_c, "Country rock", fontsize=10)

    # Melt fractions
    melt_labels = [90, 80, 70, 60, 50, 40, 30, 20, 10]

    for pct in melt_labels:

        f = 1 - pct / 100

        Sr_f = (Sr_m * f**best_D + Sr_c * best_r * (1 - f)) / (f + best_r * (1 - f))
        R_f = (R_m * Sr_m * f**best_D + R_c * Sr_c * best_r * (1 - f)) / (
            Sr_f * (f + best_r * (1 - f))
        )

        txt = ax.text(Sr_f + 15, R_f, f"{pct}%", fontsize=9)

        txt.set_path_effects([pe.withStroke(linewidth=3, foreground="white")])

    # legend
    best_handle = Line2D(
        [0],
        [0],
        color="black",
        linewidth=3,
        label=f"Best AFC fit (r={best_r:.2f}, D={best_D:.2f})",
    )

    handles.append(best_handle)

    ax.legend(handles=handles, title="Magmatic suites", fontsize=9, loc="upper right")

    ax.set_xlabel("Sr (ppm)")
    ax.set_ylabel(r"$^{87}Sr/^{86}Sr$")
    ax.set_xlim(0, 1500)

    plt.tight_layout()

    plt.savefig("AFC_SrSr.png", dpi=600)
    plt.savefig("AFC_SrSr.svg")
    plt.savefig("AFC_SrSr.pdf")

    plt.show()


# ====================== MAIN ======================

def main():

    if len(sys.argv) < 2:
        print("Usage: python afc_sr_model.py data.csv")
        sys.exit()

    file_path = sys.argv[1]

    df = load_data(file_path)

    sercol = "SERIE" if "SERIE" in df.columns else "Serie"
    ratio_col = "87Sr/86Sr(T)" if "87Sr/86Sr(T)" in df.columns else "87Sr/86Sr"
    xcol = "Sr ppm" if "Sr ppm" in df.columns else "Sr"

    df["SERIE_norm"] = df[sercol].apply(normalize_series)

    df["Sr"] = to_num(df[xcol])
    df["ratio"] = to_num(df[ratio_col])

    df = df.dropna(subset=["Sr", "ratio"])

    best_r, best_D, best_curve, Sr_m, R_m, Sr_c, R_c = afc_fit(df)

    print(f"Best r = {best_r:.3f}")
    print(f"Best D = {best_D:.3f}")

    make_plot(df, best_r, best_D, best_curve, Sr_m, R_m, Sr_c, R_c)


if __name__ == "__main__":
    main()

