import re
import unicodedata
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D


def afc_curve(Sr_m, R_m, Sr_c, R_c, r, D):
    """
    Compute AFC evolution curve.
    """
    F = np.linspace(1, 0.01, 500)

    Sr_model = (Sr_m * F**D + Sr_c * r * (1 - F)) / (F + r * (1 - F))
    R_model = (R_m * Sr_m * F**D + R_c * Sr_c * r * (1 - F)) / (
        Sr_model * (F + r * (1 - F))
    )

    return Sr_model, R_model


def monte_carlo_afc(
    x_data,
    y_data,
    Sr_m,
    R_m,
    Sr_c,
    R_c,
    iterations=10000,
    random_state=None
):
    """
    Estimate AFC parameters via Monte Carlo simulation.
    """
    x_data = np.asarray(x_data, dtype=float)
    y_data = np.asarray(y_data, dtype=float)

    rng = np.random.default_rng(random_state)

    best_error = np.inf
    best_r = None
    best_D = None
    best_curve = None

    for _ in range(iterations):
        r = rng.uniform(0.01, 0.5)
        D = rng.uniform(0.5, 3.0)

        Sr_model, R_model = afc_curve(Sr_m, R_m, Sr_c, R_c, r, D)

        dist = np.abs(R_model[:, None] - y_data) + np.abs(Sr_model[:, None] - x_data) / 1000
        error = np.sum(np.min(dist, axis=0))

        if error < best_error:
            best_error = error
            best_r = r
            best_D = D
            best_curve = (Sr_model, R_model)

    return best_r, best_D, best_curve


def plot_afc_model(
    x_data,
    y_data,
    Sr_m,
    R_m,
    Sr_c,
    R_c,
    iterations=10000,
    figsize=(8, 6),
    data_label="Observed data",
    curve_label="Best-fit AFC curve",
    save_path=None,
    random_state=None
):
    """
    Run AFC Monte Carlo fitting and plot the best-fit AFC model.
    """
    best_r, best_D, best_curve = monte_carlo_afc(
        x_data,
        y_data,
        Sr_m,
        R_m,
        Sr_c,
        R_c,
        iterations=iterations,
        random_state=random_state
    )

    Sr_model, R_model = best_curve

    fig, ax = plt.subplots(figsize=figsize)

    ax.scatter(
        x_data,
        y_data,
        edgecolor="black",
        label=data_label
    )

    ax.plot(
        Sr_model,
        R_model,
        linewidth=2,
        color="black",
        label=f"{curve_label} (r={best_r:.3f}, D={best_D:.3f})"
    )

    ax.set_xlabel("Sr (ppm)")
    ax.set_ylabel(r"$^{87}$Sr/$^{86}$Sr")
    ax.grid(True, linestyle="--", alpha=0.3)
    ax.legend()
    plt.tight_layout()

    if save_path is not None:
        fig.savefig(save_path, dpi=300, bbox_inches="tight")

    return fig, best_r, best_D, best_curve


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
        "ehkcalcalk": "EHKCalcAlk",
        "calcalk": "CalcAlk",
        "chalk": "ChAlk",
        "shos": "Shos",
        "shoshonitic": "Shos",
        "archean": "Archean",
        "arqueano": "Archean",
        "paleoproterozoic": "Paleoproterozoic",
        "paleoproterozoico": "Paleoproterozoic"
    }

    return mapping.get(key, str(value).strip())


def _to_numeric_clean(series):
    return (
        series.astype(str)
        .str.replace(",", ".", regex=False)
        .str.replace(" ", "", regex=False)
        .pipe(lambda s: np.array(s, dtype=object))
    )


def run_afc_from_dataframe(
    df,
    Sr_m,
    R_m,
    Sr_c,
    R_c,
    sr_col="Sr ppm",
    ratio_col="87Sr/86Sr(T)",
    series_col=None,
    iterations=10000,
    figsize=(9, 7),
    save_path=None,
    random_state=42,
    styled=True,
    annotate_endmembers=True,
    annotate_percentages=True
):
    """
    High-level AFC workflow from a pandas DataFrame.

    Parameters
    ----------
    df : pandas.DataFrame
        Input dataframe containing Sr and isotope ratio columns.
    Sr_m, R_m, Sr_c, R_c : float
        AFC end-member parameters.
    sr_col : str, optional
        Name of Sr concentration column.
    ratio_col : str, optional
        Name of isotopic ratio column.
    series_col : str, optional
        Name of suite/series column. If provided, points are classified by series.
    iterations : int, optional
        Number of Monte Carlo iterations.
    figsize : tuple, optional
        Figure size.
    save_path : str, optional
        Path to save figure.
    random_state : int, optional
        Seed for reproducibility.
    styled : bool, optional
        If True, generates a styled figure similar to the original AFC plot.
    annotate_endmembers : bool, optional
        If True, label parental magma and country rock.
    annotate_percentages : bool, optional
        If True, add assimilation and fractionation labels along the curve.

    Returns
    -------
    fig : matplotlib.figure.Figure
    results : dict
        Dictionary with best-fit parameters and model curve.
    """
    required_cols = [sr_col, ratio_col]
    missing_cols = [col for col in required_cols if col not in df.columns]
    if missing_cols:
        raise ValueError(
            f"Missing required columns: {missing_cols}. Expected columns are: {required_cols}"
        )

    data = df.copy()

    data[sr_col] = data[sr_col].astype(str).str.replace(",", ".", regex=False).str.replace(" ", "", regex=False)
    data[ratio_col] = data[ratio_col].astype(str).str.replace(",", ".", regex=False).str.replace(" ", "", regex=False)

    data[sr_col] = np.array(data[sr_col], dtype=object)
    data[ratio_col] = np.array(data[ratio_col], dtype=object)

    import pandas as pd
    data[sr_col] = pd.to_numeric(data[sr_col], errors="coerce")
    data[ratio_col] = pd.to_numeric(data[ratio_col], errors="coerce")

    if series_col is not None and series_col in data.columns:
        data["SERIE_norm"] = data[series_col].apply(normalize_series)
    else:
        data["SERIE_norm"] = "Observed data"

    data = data.dropna(subset=[sr_col, ratio_col])

    x_data = data[sr_col].to_numpy(dtype=float)
    y_data = data[ratio_col].to_numpy(dtype=float)

    best_r, best_D, best_curve = monte_carlo_afc(
        x_data=x_data,
        y_data=y_data,
        Sr_m=Sr_m,
        R_m=R_m,
        Sr_c=Sr_c,
        R_c=R_c,
        iterations=iterations,
        random_state=random_state
    )

    Sr_model, R_model = best_curve

    fig, ax = plt.subplots(figsize=figsize, facecolor="white")
    ax.set_facecolor("white")

    # Best-fit curve
    ax.plot(
        Sr_model,
        R_model,
        color="black",
        linewidth=3.0 if styled else 2.0,
        alpha=0.75,
        label=f"Best AFC fit (r={best_r:.2f}, D={best_D:.2f})"
    )

    if styled:
        series_colors = {
            "PHKCalcAlk": "red",
            "EHKCalcAlk": "orange",
            "ChAlk": "deeppink",
            "Shos": "green",
            "Alk": "lightpink",
            "CalcAlk": "#FFD700"
        }

        crust_symbols = {
            "Archean": "s",
            "Paleoproterozoic": "^"
        }

        handles = []

        # magmatic suites
        for serie, color in series_colors.items():
            sub = data[data["SERIE_norm"] == serie]
            if len(sub) > 0:
                sc = ax.scatter(
                    sub[sr_col],
                    sub[ratio_col],
                    s=55,
                    color=color,
                    edgecolor="black",
                    linewidth=0.3,
                    alpha=0.8,
                    label=f"{serie} (n={len(sub)})"
                )
                handles.append(sc)

        # crustal groups
        for crust, marker in crust_symbols.items():
            sub = data[data["SERIE_norm"] == crust]
            if len(sub) > 0:
                sc = ax.scatter(
                    sub[sr_col],
                    sub[ratio_col],
                    marker=marker,
                    facecolor="white",
                    edgecolor="black",
                    s=50,
                    linewidth=0.9,
                    label=f"{crust} (n={len(sub)})"
                )
                handles.append(sc)

        # fallback for unmatched groups
        unmatched = sorted(
            set(data["SERIE_norm"].dropna().unique()) -
            set(series_colors.keys()) -
            set(crust_symbols.keys())
        )
        for group in unmatched:
            sub = data[data["SERIE_norm"] == group]
            if len(sub) > 0:
                sc = ax.scatter(
                    sub[sr_col],
                    sub[ratio_col],
                    s=45,
                    color="gray",
                    edgecolor="black",
                    linewidth=0.3,
                    alpha=0.8,
                    label=f"{group} (n={len(sub)})"
                )
                handles.append(sc)

        # End-members
        ax.scatter(Sr_m, R_m, color="black", s=120, marker="o")
        ax.scatter(Sr_c, R_c, color="red", s=120, marker="^")

        if annotate_endmembers:
            ax.text(Sr_m, R_m + 0.0025, "Parental magma", fontsize=10, ha="center")
            ax.text(Sr_c + 15, R_c, "Country rock", fontsize=10)

        # Assimilation and fractionation labels
        if annotate_percentages:
            melt_labels = [90, 80, 70, 60, 50, 40, 30, 20, 10]

            for pct in melt_labels:
                f = pct / 100.0

                Sr_f = (Sr_m * f**best_D + Sr_c * best_r * (1 - f)) / (f + best_r * (1 - f))
                R_f = (R_m * Sr_m * f**best_D + R_c * Sr_c * best_r * (1 - f)) / (
                    Sr_f * (f + best_r * (1 - f))
                )

                assimilated = best_r * (1 - f) * 100
                fractionated = (1 - f) * 100

                ax.text(
                    Sr_f - 20,
                    R_f + 0.0055,
                    f"{assimilated:.0f}%",
                    fontsize=9,
                    color="blue"
                )

                ax.text(
                    Sr_f - 20,
                    R_f - 0.0080,
                    f"{fractionated:.0f}%",
                    fontsize=9,
                    color="black"
                )

        best_handle = Line2D(
            [0], [0],
            color="black",
            linewidth=3,
            label=f"Best AFC fit (r={best_r:.2f}, D={best_D:.2f})"
        )
        handles.append(best_handle)

        ax.legend(
            handles=handles,
            title="Magmatic suites",
            fontsize=9,
            loc="upper right"
        )

        ax.set_xlim(0, 1500)

    else:
        ax.scatter(
            x_data,
            y_data,
            s=35,
            edgecolor="black",
            label="Observed data"
        )
        ax.legend()

    ax.set_xlabel("Sr (ppm)", fontsize=12)
    ax.set_ylabel(r"$^{87}$Sr/$^{86}$Sr", fontsize=12)
    ax.grid(False if styled else True, linestyle="--", alpha=0.3)
    plt.tight_layout()

    if save_path is not None:
        fig.savefig(save_path, dpi=300, bbox_inches="tight")

    results = {
        "best_r": best_r,
        "best_D": best_D,
        "best_curve": best_curve
    }

    return fig, results