import os
import pandas as pd
import matplotlib.pyplot as plt

from borborema.wr_pca import run_pca_from_dataframe
from borborema.afc_model import run_afc_from_dataframe
from borborema.sr_nd_models import run_sr_nd_from_dataframe
from borborema.upb_geochronology import plot_upb_ages


os.makedirs("figures", exist_ok=True)

# =====================
# PCA
# =====================

df = pd.read_csv("sample_data/wr_pca_exemple.csv", encoding="utf-8-sig")
df.columns = df.columns.str.replace("ï»¿", "", regex=False)

variables = [
"SiO2","MgO","CaO","Na2O","K2O",
"Rb","Ba","Sr","Nb","Zr","Y","Th","La","Ce"
]

fig, _ = run_pca_from_dataframe(
df,
variables=variables,
series_col="Suite"
)

fig.savefig("figures/pca_plot.png", dpi=300)


# =====================
# AFC
# =====================

df = pd.read_csv("sample_data/afc_model_exemple.csv", encoding="latin1")

fig, _ = run_afc_from_dataframe(
df,
Sr_m=1356,
R_m=0.71008,
Sr_c=96.8,
R_c=0.7808,
series_col="SERIE",
iterations=10000
)

fig.savefig("figures/afc_model.png", dpi=300)


# =====================
# Sm-Nd
# =====================

df = pd.read_csv("sample_data/sr_nd_models_exemple.csv", encoding="latin1")

fig, _ = run_sr_nd_from_dataframe(df)

fig.savefig("figures/sm_nd_panel.png", dpi=300)


# =====================
# U-Pb
# =====================

df = pd.read_csv("sample_data/upb_geochronology_exemple.csv", encoding="utf-8-sig")

fig = plot_upb_ages(df)

fig.savefig("figures/upb_ages.png", dpi=300)

print("All figures generated in 'figures/' directory")