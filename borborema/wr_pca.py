import re
import unicodedata
import numpy as np
import matplotlib.pyplot as plt
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from matplotlib.lines import Line2D


def _strip_accents(text):
    text = unicodedata.normalize("NFKD", str(text))
    return "".join(ch for ch in text if not unicodedata.combining(ch))


def normalize_series(value):
    """
    Normalize suite/series names to standard toolkit labels.
    """
    if value is None:
        return np.nan

    try:
        if np.isnan(value):
            return np.nan
    except TypeError:
        pass

    s = _strip_accents(value).lower()
    key = re.sub(r"[^a-z]", "", s)

    mapping = {
        "alk": "Alk",
        "alkaline": "Alk",
        "phkcalcalk": "PHKCalcAlk",
        "phkcalcalkalk": "PHKCalcAlk",
        "ehkcalcalk": "EHKCalcAlk",
        "ehkcalcalkalk": "EHKCalcAlk",
        "calcalk": "CalcAlk",
        "chalk": "ChAlk",
        "shos": "Shos",
        "shoshonitic": "Shos"
    }

    return mapping.get(key, str(value).strip())


def _make_color_map(categories, cmap_name="tab10"):
    categories = list(categories)
    cmap = plt.get_cmap(cmap_name, len(categories))
    return {cat: cmap(i) for i, cat in enumerate(categories)}


def _default_colors():
    return {
        "PHKCalcAlk": "red",
        "EHKCalcAlk": "orange",
        "ChAlk": "deeppink",
        "Shos": "green",
        "Alk": "lightpink",
        "CalcAlk": "#FFD700"
    }


def run_pca(df, variables):
    """
    Perform Principal Component Analysis on geochemical data.

    Parameters
    ----------
    df : pandas.DataFrame
        Geochemical dataset.
    variables : list
        List of geochemical variables.

    Returns
    -------
    PC : ndarray
        PCA scores.
    loadings : ndarray
        PCA loadings.
    explained_variance : ndarray
        Variance explained by each component.
    """
    X = df[variables].values
    X_scaled = StandardScaler().fit_transform(X)

    pca = PCA(n_components=2)
    PC = pca.fit_transform(X_scaled)
    loadings = pca.components_.T

    return PC, loadings, pca.explained_variance_ratio_


def plot_pca(
    PC,
    loadings,
    variables,
    suites=None,
    colors=None,
    figsize=(8, 7),
    save_path=None,
    legend_title="Suite"
):
    """
    Plot PCA biplot.

    Parameters
    ----------
    PC : ndarray
        PCA scores.
    loadings : ndarray
        PCA loadings.
    variables : list
        Geochemical variables.
    suites : array-like, optional
        Magmatic suite labels.
    colors : dict, optional
        Color mapping for suites.
    figsize : tuple, optional
        Figure size.
    save_path : str, optional
        Path to save figure.
    legend_title : str, optional
        Legend title.

    Returns
    -------
    matplotlib.figure.Figure
    """
    if suites is None:
        suites = np.array(["Data"] * len(PC), dtype=object)
    else:
        suites = np.asarray(suites, dtype=object)
        suites = np.array([normalize_series(x) for x in suites], dtype=object)

    unique_suites = list(dict.fromkeys(suites))

    if colors is None:
        default_colors = _default_colors()
        missing = [s for s in unique_suites if s not in default_colors]
        auto_colors = _make_color_map(missing, cmap_name="tab10") if missing else {}
        colors = {**auto_colors, **default_colors}
    else:
        missing = [s for s in unique_suites if s not in colors]
        if missing:
            auto_colors = _make_color_map(missing, cmap_name="tab10")
            colors = {**colors, **auto_colors}

    fig, ax = plt.subplots(figsize=figsize)

    for serie in unique_suites:
        mask = suites == serie
        ax.scatter(
            PC[mask, 0],
            PC[mask, 1],
            s=70,
            color=colors.get(serie, "gray"),
            edgecolor="black",
            linewidth=0.4,
            alpha=0.85,
            label=serie
        )

    pc_range = max(abs(PC[:, 0]).max(), abs(PC[:, 1]).max())
    scale = pc_range * 0.85

    for i, element in enumerate(variables):
        x = loadings[i, 0] * scale
        y = loadings[i, 1] * scale

        ax.arrow(
            0, 0, x, y,
            color="black",
            linewidth=1.4,
            head_width=0.18,
            length_includes_head=True
        )

        ax.text(
            x * 1.12,
            y * 1.12,
            element,
            fontsize=9,
            ha="center",
            va="center"
        )

    ax.set_xlabel("PC1")
    ax.set_ylabel("PC2")

    ax.axhline(0, color="gray", lw=0.8)
    ax.axvline(0, color="gray", lw=0.8)

    legend_elements = [
        Line2D(
            [0], [0],
            marker="o",
            linestyle="",
            label=serie,
            markerfacecolor=colors.get(serie, "gray"),
            markeredgecolor="black",
            markersize=8
        )
        for serie in unique_suites
    ]

    ax.legend(handles=legend_elements, title=legend_title, loc="best")
    plt.tight_layout()

    if save_path is not None:
        fig.savefig(save_path, dpi=300, bbox_inches="tight")

    return fig


def run_pca_from_dataframe(
    df,
    variables,
    series_col=None,
    colors=None,
    figsize=(8, 7),
    save_path=None,
    legend_title="Suite"
):
    """
    High-level PCA workflow from a pandas DataFrame.

    Parameters
    ----------
    df : pandas.DataFrame
        Input dataframe containing geochemical variables.
    variables : list
        List of variables to include in PCA.
    series_col : str, optional
        Name of suite/series column.
    colors : dict, optional
        Color mapping for suites.
    figsize : tuple, optional
        Figure size.
    save_path : str, optional
        Path to save figure.
    legend_title : str, optional
        Legend title.

    Returns
    -------
    fig : matplotlib.figure.Figure
    results : dict
    """
    import pandas as pd

    missing_cols = [col for col in variables if col not in df.columns]
    if missing_cols:
        raise ValueError(
            f"Missing required PCA variables: {missing_cols}. "
            f"Available columns are: {list(df.columns)}"
        )

    data = df.copy()

    # numeric cleaning
    for col in variables:
        data[col] = pd.to_numeric(
            data[col].astype(str)
            .str.replace(",", ".", regex=False)
            .str.replace(" ", "", regex=False),
            errors="coerce"
        )

    data = data.dropna(subset=variables).copy()

    if data.empty:
        raise ValueError("No valid rows available after filtering missing PCA variables.")

    if series_col is not None and series_col in data.columns:
        suites = data[series_col].apply(normalize_series).to_numpy(dtype=object)
    else:
        suites = np.array(["Data"] * len(data), dtype=object)

    PC, loadings, explained = run_pca(data, variables)

    fig = plot_pca(
        PC=PC,
        loadings=loadings,
        variables=variables,
        suites=suites,
        colors=colors,
        figsize=figsize,
        save_path=save_path,
        legend_title=legend_title
    )

    results = {
        "PC": PC,
        "loadings": loadings,
        "explained_variance": explained,
        "groups": list(dict.fromkeys(suites)),
        "n_points": len(data)
    }

    return fig, results