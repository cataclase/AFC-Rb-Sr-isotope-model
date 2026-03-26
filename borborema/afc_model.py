import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

from .data_cleaning import clean_numeric_series, strip_accents


def standardize_suite_name(value):
    """
    Standardize suite names to internal labels used across the toolkit.
    """
    if value is None:
        return np.nan

    try:
        if np.isnan(value):
            return np.nan
    except Exception:
        pass

    value = str(value).strip()
    key = strip_accents(value).lower().strip()

    suite_map = {
        "shoshonitic": "Shos",
        "shosho": "Shos",
        "shos": "Shos",

        "alkaline": "Alk",
        "alk": "Alk",

        "calc-alkaline": "CalcAlk",
        "calcalk": "CalcAlk",
        "calc": "CalcAlk",

        "porphyritic high-k calc-alkaline": "PHKCalcAlk",
        "phkcalcalk": "PHKCalcAlk",
        "calcp": "PHKCalcAlk",

        "equigranular high-k calc-alkaline": "EHKCalcAlk",
        "ehkcalcalk": "EHKCalcAlk",
        "calce": "EHKCalcAlk",

        "alkaline-charnockitic": "ChAlk",
        "chalk": "ChAlk",
        "charnockitic": "ChAlk",

        "archean": "Archean",
        "arqueano": "Archean",

        "paleoproterozoic": "Paleoproterozoic",
        "paleoproterozoico": "Paleoproterozoic",
    }

    return suite_map.get(key, value)


def _default_colors():
    """
    Default color palette used across the toolkit.
    """
    return {
        "Shos": "#2ca02c",
        "Alk": "#f4a6b5",
        "EHKCalcAlk": "#1f77b4",
        "PHKCalcAlk": "#d62728",
        "CalcAlk": "#e377c2",
        "ChAlk": "#17becf",
        "Archean": "#ffffff",
        "Paleoproterozoic": "#d9d9d9",
        "Observed data": "gray",
    }


def _default_labels():
    """
    Human-readable labels for legends.
    """
    return {
        "Shos": "Shos",
        "Alk": "Alk",
        "EHKCalcAlk": "Equigranular High-K Calc-Alkaline",
        "PHKCalcAlk": "Porphyritic High-K Calc-Alkaline",
        "CalcAlk": "Calc-Alkaline",
        "ChAlk": "Alkaline-Charnockitic",
        "Archean": "Archean",
        "Paleoproterozoic": "Paleoproterozoic",
        "Observed data": "Observed data",
    }


def afc_curve(Sr_m, R_m, Sr_c, R_c, r, D):
    """
    Compute AFC evolution curve.

    Parameters
    ----------
    Sr_m : float
        Initial Sr concentration in parental magma.
    R_m : float
        Initial 87Sr/86Sr ratio in parental magma.
    Sr_c : float
        Sr concentration in contaminant/country rock.
    R_c : float
        87Sr/86Sr ratio in contaminant/country rock.
    r : float
        Assimilation/crystallization ratio.
    D : float
        Bulk distribution coefficient.

    Returns
    -------
    Sr_model : ndarray
    R_model : ndarray
    F : ndarray
    """
    F = np.linspace(1.0, 0.01, 500)

    Sr_model = (Sr_m * F**D + Sr_c * r * (1 - F)) / (F + r * (1 - F))
    R_model = (R_m * Sr_m * F**D + R_c * Sr_c * r * (1 - F)) / (
        Sr_model * (F + r * (1 - F))
    )

    return Sr_model, R_model, F


def monte_carlo_afc(
    x_data,
    y_data,
    Sr_m,
    R_m,
    Sr_c,
    R_c,
    iterations=10000,
    random_state=None,
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

        Sr_model, R_model, F = afc_curve(Sr_m, R_m, Sr_c, R_c, r, D)

        dist = (
            np.abs(R_model[:, None] - y_data)
            + np.abs(Sr_model[:, None] - x_data) / 1000.0
        )
        error = np.sum(np.min(dist, axis=0))

        if error < best_error:
            best_error = error
            best_r = r
            best_D = D
            best_curve = (Sr_model, R_model, F)

    return best_r, best_D, best_curve, best_error


def _annotate_percentages(
    ax,
    Sr_model,
    R_model,
    F,
    percentages=(90, 80, 70, 60, 50, 40, 30, 20, 10),
):
    """
    Annotate crystallization percentages along the AFC curve.
    """
    crystallized = (1 - F) * 100

    for pct in percentages:
        idx = np.argmin(np.abs(crystallized - pct))
        ax.text(
            Sr_model[idx],
            R_model[idx],
            f"{pct}%",
            fontsize=8,
            ha="left",
            va="bottom",
        )


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
    random_state=None,
    annotate_endmembers=True,
    annotate_percentages=True,
):
    """
    Run AFC Monte Carlo fitting and plot the best-fit AFC model.
    """
    best_r, best_D, best_curve, best_error = monte_carlo_afc(
        x_data=x_data,
        y_data=y_data,
        Sr_m=Sr_m,
        R_m=R_m,
        Sr_c=Sr_c,
        R_c=R_c,
        iterations=iterations,
        random_state=random_state,
    )

    Sr_model, R_model, F = best_curve

    fig, ax = plt.subplots(figsize=figsize)

    ax.scatter(
        x_data,
        y_data,
        s=65,
        color="gray",
        edgecolor="black",
        linewidth=0.4,
        alpha=0.85,
        label=data_label,
        zorder=2,
    )

    ax.plot(
        Sr_model,
        R_model,
        linewidth=2.5,
        color="black",
        label=f"{curve_label} (r={best_r:.3f}, D={best_D:.3f})",
        zorder=3,
    )

    if annotate_endmembers:
        ax.scatter([Sr_m, Sr_c], [R_m, R_c], color="black", s=35, zorder=5)
        ax.text(Sr_m, R_m, " Parental magma", fontsize=9, ha="left", va="bottom")
        ax.text(Sr_c, R_c, " Country rock", fontsize=9, ha="left", va="bottom")

    if annotate_percentages:
        _annotate_percentages(ax, Sr_model, R_model, F)

    ax.set_xlabel("Sr (ppm)")
    ax.set_ylabel(r"$^{87}$Sr/$^{86}$Sr")
    ax.grid(True, linestyle="--", alpha=0.3)
    ax.legend()
    plt.tight_layout()

    if save_path is not None:
        fig.savefig(save_path, dpi=300, bbox_inches="tight")

    results = {
        "best_r": best_r,
        "best_D": best_D,
        "best_error": best_error,
        "curve": {
            "Sr_model": Sr_model,
            "R_model": R_model,
            "F": F,
        },
    }

    return fig, results


def run_afc_from_dataframe(
    df,
    Sr_m,
    R_m,
    Sr_c,
    R_c,
    sr_col="Sr",
    ratio_col="87Sr/86Sr",
    series_col="Suite",
    iterations=10000,
    figsize=(9, 7),
    save_path=None,
    random_state=42,
    styled=True,
    annotate_endmembers=True,
    annotate_percentages=True,
    below_detection="value",
):
    """
    High-level AFC workflow from a dataframe.
    """
    missing_cols = [col for col in [sr_col, ratio_col] if col not in df.columns]
    if missing_cols:
        raise ValueError(
            f"Missing required AFC columns: {missing_cols}. "
            f"Available columns are: {list(df.columns)}"
        )

    data = df.copy()
    data[sr_col] = clean_numeric_series(data[sr_col], below_detection=below_detection)
    data[ratio_col] = clean_numeric_series(
        data[ratio_col],
        below_detection=below_detection,
    )

    if series_col is not None and series_col in data.columns:
        data["SERIE_norm"] = data[series_col].apply(standardize_suite_name)
    else:
        data["SERIE_norm"] = "Observed data"

    data = data.dropna(subset=[sr_col, ratio_col]).copy()

    if data.empty:
        raise ValueError("No valid AFC rows available after cleaning numeric columns.")

    x_data = data[sr_col].to_numpy(dtype=float)
    y_data = data[ratio_col].to_numpy(dtype=float)

    fig, fit = plot_afc_model(
        x_data=x_data,
        y_data=y_data,
        Sr_m=Sr_m,
        R_m=R_m,
        Sr_c=Sr_c,
        R_c=R_c,
        iterations=iterations,
        figsize=figsize,
        save_path=None,
        random_state=random_state,
        annotate_endmembers=annotate_endmembers,
        annotate_percentages=annotate_percentages,
    )

    if styled:
        ax = fig.axes[0]
        default_colors = _default_colors()
        default_labels = _default_labels()

        if ax.get_legend() is not None:
            ax.get_legend().remove()

        plot_order = [
            "PHKCalcAlk",
            "EHKCalcAlk",
            "Shos",
            "Alk",
            "ChAlk",
            "CalcAlk",
            "Archean",
            "Paleoproterozoic",
        ]

        present_series = [s for s in plot_order if s in data["SERIE_norm"].unique()]
        other_series = [
            s for s in data["SERIE_norm"].dropna().unique()
            if s not in present_series
        ]
        final_order = present_series + list(other_series)

        # remove the gray base scatter so only styled groups remain
        if len(ax.collections) > 0:
            ax.collections[0].remove()

        for serie in final_order:
            subset = data[data["SERIE_norm"] == serie]

            facecolor = default_colors.get(serie, "gray")
            edgecolor = "black"

            ax.scatter(
                subset[sr_col],
                subset[ratio_col],
                s=65,
                facecolor=facecolor,
                edgecolor=edgecolor,
                linewidth=0.7,
                alpha=0.85,
                zorder=4,
            )

        handles = []
        for serie in final_order:
            subset = data[data["SERIE_norm"] == serie]
            if len(subset) == 0:
                continue

            markerfacecolor = default_colors.get(serie, "gray")
            markeredgecolor = "black"

            handles.append(
                Line2D(
                    [0],
                    [0],
                    marker="o",
                    linestyle="",
                    label=f"{default_labels.get(serie, serie)} (n={len(subset)})",
                    markerfacecolor=markerfacecolor,
                    markeredgecolor=markeredgecolor,
                    markersize=8,
                )
            )

        handles.append(
            Line2D(
                [0],
                [0],
                color="black",
                lw=2.5,
                label=f"Best AFC fit (r={fit['best_r']:.2f}, D={fit['best_D']:.2f})",
            )
        )

        ax.legend(handles=handles, title="Magmatic suites", loc="best")
        plt.tight_layout()

    if save_path is not None:
        fig.savefig(save_path, dpi=300, bbox_inches="tight")

    results = {
        "best_r": fit["best_r"],
        "best_D": fit["best_D"],
        "best_error": fit["best_error"],
        "n_points": len(data),
    }

    return fig, results