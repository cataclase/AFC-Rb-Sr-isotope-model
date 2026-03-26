import re
import numpy as np
import matplotlib.pyplot as plt
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from matplotlib.lines import Line2D

from .data_cleaning import clean_numeric_series, strip_accents


# ============================================================
# Suite normalization (CONSISTENTE COM TODO O PROJETO)
# ============================================================

def normalize_suite(value):
    if value is None:
        return np.nan

    try:
        if np.isnan(value):
            return np.nan
    except:
        pass

    s = strip_accents(str(value)).lower().strip()
    key = re.sub(r"[^a-z]", "", s)

    mapping = {
        "shos": "Shos",
        "shoshonitic": "Shos",

        "alk": "Alk",
        "alkaline": "Alk",

        "calcalk": "CalcAlk",
        "calcalkaline": "CalcAlk",
        "calc": "CalcAlk",

        "chalk": "ChAlk",
        "alkalinecharnockitic": "ChAlk",
        "charnockitic": "ChAlk",

        "ehkcalcalk": "EHKCalcAlk",
        "equigranularhighkcalcalkaline": "EHKCalcAlk",
        "calce": "EHKCalcAlk",

        "phkcalcalk": "PHKCalcAlk",
        "porphyritichighkcalcalkaline": "PHKCalcAlk",
        "calcp": "PHKCalcAlk",
    }

    return mapping.get(key, str(value).strip())


# ============================================================
# Cores padrão (IGUAL A TODOS OS OUTROS MÓDULOS)
# ============================================================

def _default_colors():
    return {
        "Shos": "#2ca02c",
        "Alk": "#f4a6b5",
        "EHKCalcAlk": "#1f77b4",
        "PHKCalcAlk": "#d62728",
        "CalcAlk": "#e377c2",
        "ChAlk": "#17becf",
        "Data": "gray",
    }


def _default_order():
    return [
        "Shos",
        "Alk",
        "EHKCalcAlk",
        "PHKCalcAlk",
        "CalcAlk",
        "ChAlk",
        "Data",
    ]


# ============================================================
# PCA computation
# ============================================================

def run_pca(df, variables):
    X = df[variables].to_numpy(dtype=float)

    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(X)

    pca = PCA(n_components=2)
    PC = pca.fit_transform(X_scaled)
    loadings = pca.components_.T
    explained = pca.explained_variance_ratio_

    return PC, loadings, explained


# ============================================================
# Plot PCA
# ============================================================

def plot_pca(
    PC,
    loadings,
    variables,
    suites=None,
    colors=None,
    figsize=(8, 7),
    save_path=None,
    legend_title="Suite",
    explained=None,
):
    if suites is None:
        suites = np.array(["Data"] * len(PC), dtype=object)
    else:
        suites = np.asarray([normalize_suite(x) for x in suites], dtype=object)

    detected = list(dict.fromkeys(suites))

    if colors is None:
        colors = _default_colors()

    preferred = _default_order()
    order = [s for s in preferred if s in detected]

    fig, ax = plt.subplots(figsize=figsize)

    for serie in order:
        mask = suites == serie
        n = np.sum(mask)

        if n == 0:
            continue

        ax.scatter(
            PC[mask, 0],
            PC[mask, 1],
            s=70,
            color=colors.get(serie, "gray"),
            edgecolor="black",
            linewidth=0.4,
            alpha=0.85,
        )

    # Loadings
    pc_range = max(abs(PC[:, 0]).max(), abs(PC[:, 1]).max())
    scale = pc_range * 0.85 if pc_range > 0 else 1.0

    for i, element in enumerate(variables):
        x = loadings[i, 0] * scale
        y = loadings[i, 1] * scale

        ax.arrow(
            0, 0, x, y,
            color="black",
            linewidth=1.4,
            head_width=0.18,
            length_includes_head=True,
        )

        ax.text(
            x * 1.12,
            y * 1.12,
            element,
            fontsize=9,
            ha="center",
            va="center",
        )

    if explained is not None:
        ax.set_xlabel(f"PC1 ({explained[0] * 100:.1f}%)")
        ax.set_ylabel(f"PC2 ({explained[1] * 100:.1f}%)")
    else:
        ax.set_xlabel("PC1")
        ax.set_ylabel("PC2")

    ax.axhline(0, color="gray", lw=0.8)
    ax.axvline(0, color="gray", lw=0.8)

    # Legenda com n=
    legend_elements = []

    for serie in order:
        n = np.sum(suites == serie)
        if n == 0:
            continue

        legend_elements.append(
            Line2D(
                [0], [0],
                marker="o",
                linestyle="",
                label=f"{serie} (n={n})",
                markerfacecolor=colors.get(serie, "gray"),
                markeredgecolor="black",
                markersize=8,
            )
        )

    ax.legend(handles=legend_elements, title=legend_title, loc="best")

    plt.tight_layout()

    if save_path is not None:
        fig.savefig(save_path, dpi=300, bbox_inches="tight")

    return fig


# ============================================================
# High-level wrapper
# ============================================================

def run_pca_from_dataframe(
    df,
    variables,
    series_col=None,
    colors=None,
    figsize=(8, 7),
    save_path=None,
    legend_title="Suite",
    below_detection="value",
):
    missing_cols = [col for col in variables if col not in df.columns]
    if missing_cols:
        raise ValueError(f"Missing PCA variables: {missing_cols}")

    data = df.copy()

    for col in variables:
        data[col] = clean_numeric_series(data[col], below_detection=below_detection)

    data = data.dropna(subset=variables)

    if data.empty:
        raise ValueError("No valid PCA data after cleaning")

    if series_col and series_col in data.columns:
        suites = data[series_col].to_numpy()
    else:
        suites = np.array(["Data"] * len(data))

    PC, loadings, explained = run_pca(data, variables)

    fig = plot_pca(
        PC,
        loadings,
        variables,
        suites=suites,
        colors=colors,
        figsize=figsize,
        save_path=save_path,
        legend_title=legend_title,
        explained=explained,
    )

    results = {
        "PC": PC,
        "loadings": loadings,
        "explained_variance": explained,
        "groups": list(dict.fromkeys(suites)),
        "n_points": len(data),
    }

    return fig, results