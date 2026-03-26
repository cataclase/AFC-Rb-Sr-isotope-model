import re
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import ticker as mticker
from matplotlib.lines import Line2D

from .data_cleaning import clean_numeric_series, strip_accents


# ============================================================
# Reference curves
# ============================================================

def dm_curve():
    """
    Generate depleted mantle εNd evolution curve.

    Returns
    -------
    tempo_ma : ndarray
        Time in Ma.
    epsilon_nd_dm : ndarray
        Depleted mantle epsilon Nd values.
    """
    tempo_ga = np.arange(0, 4.6, 0.1)

    epsilon_nd_dm = np.array([
        8.5, 8.2025, 7.91, 7.6225, 7.34, 7.0625, 6.79, 6.5225, 6.26, 6.0025,
        5.75, 5.5025, 5.26, 5.0225, 4.79, 4.5625, 4.34, 4.1225, 3.91, 3.7025,
        3.5, 3.3025, 3.11, 2.9225, 2.74, 2.5625, 2.39, 2.2225, 2.06, 1.9025,
        1.75, 1.6025, 1.46, 1.3225, 1.19, 1.0625, 0.94, 0.8225, 0.71, 0.6025,
        0.5, 0.4025, 0.31, 0.2225, 0.14, 0.0625
    ])

    tempo_ma = tempo_ga * 1000
    return tempo_ma, epsilon_nd_dm


def chur_curve():
    """
    Generate CHUR εNd curve (horizontal line at 0).
    """
    tempo_ma, _ = dm_curve()
    epsilon_chur = np.zeros_like(tempo_ma)
    return tempo_ma, epsilon_chur


def eNd_DM(x):
    """
    Interpolate εNd of depleted mantle for a given age (Ma).
    """
    tempo_ma, epsilon_dm = dm_curve()
    return np.interp(x, tempo_ma, epsilon_dm)


# ============================================================
# Utilities
# ============================================================

def _normalize_header(text):
    """
    Normalize header names for robust matching.
    """
    text = strip_accents(str(text)).lower().strip()
    text = text.replace("−", "-")
    text = re.sub(r"\s+", "", text)
    return text


def normalize_suite(value):
    """
    Normalize suite names from the example dataset and toolkit conventions.
    """
    if value is None:
        return np.nan

    try:
        if np.isnan(value):
            return np.nan
    except Exception:
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


def _display_suite_name(suite):
    """
    Human-readable labels for legends.
    """
    display = {
        "Shos": "Shos",
        "Alk": "Alk",
        "CalcAlk": "Calc-Alkaline",
        "ChAlk": "Alkaline-Charnockitic",
        "EHKCalcAlk": "Equigranular High-K Calc-Alkaline",
        "PHKCalcAlk": "Porphyritic High-K Calc-Alkaline",
        "Data": "Data",
    }
    return display.get(suite, suite)


def _default_colors():
    """
    Default suite colors used across the toolkit.
    """
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


def _default_sizes():
    return {
        "Shos": 70,
        "Alk": 95,
        "EHKCalcAlk": 90,
        "PHKCalcAlk": 90,
        "CalcAlk": 95,
        "ChAlk": 90,
        "Data": 70,
    }


def _find_column(df, candidates, required=True):
    """
    Find a column in a dataframe using multiple candidate names.
    Matching is done after header normalization.
    """
    normalized_map = {_normalize_header(col): col for col in df.columns}

    for candidate in candidates:
        key = _normalize_header(candidate)
        if key in normalized_map:
            return normalized_map[key]

    if required:
        raise ValueError(
            f"Could not find any of these columns in dataframe: {candidates}. "
            f"Available columns are: {list(df.columns)}"
        )

    return None


# ============================================================
# Main workflow
# ============================================================

def run_sr_nd_from_dataframe(
    df,
    age_col=None,
    end_col=None,
    sr_col=None,
    tdm_col=None,
    series_col=None,
    colors=None,
    order=None,
    sizes=None,
    figsize=(16, 6),
    save_path=None,
    show_projection=True,
    show_boxes=True,
    panel_labels=True,
    below_detection="value",
):
    """
    Generate a two-panel Sm-Nd isotopic figure:
    (A) εNd(t) vs age
    (B) εNd(t) vs 87Sr/86Sr(t)
    """
    # Auto-detect columns if not provided
    if age_col is None:
        age_col = _find_column(df, ["Age (MA)", "Age (Ma)", "Age", "Idade", "Estimated_Age(MA)", "UPb_Age (Ma)"])
    if end_col is None:
        end_col = _find_column(df, ["eNd_t", "eNd(t)", "εNd(t)", "epsilon Nd", "epsilonNd"])
    if sr_col is None:
        sr_col = _find_column(df, ["87Sr/86Sr(t)", "87Sr/86Sr(T)", "87Sr/86Sr", "Sr_ratio"])
    if tdm_col is None:
        tdm_col = _find_column(df, ["TDM", "T_DM", "Nd TDM", "Tdm"], required=False)
    if series_col is None:
        series_col = _find_column(df, ["Suite", "Serie", "SERIE", "Series"], required=False)

    data = df.copy()

    # Clean numeric columns
    data[age_col] = clean_numeric_series(data[age_col], below_detection=below_detection)
    data[end_col] = clean_numeric_series(data[end_col], below_detection=below_detection)
    data[sr_col] = clean_numeric_series(data[sr_col], below_detection=below_detection)

    if tdm_col is not None and tdm_col in data.columns:
        data[tdm_col] = clean_numeric_series(data[tdm_col], below_detection=below_detection)
    else:
        data["__TDM__"] = np.nan
        tdm_col = "__TDM__"

    if series_col is not None and series_col in data.columns:
        data["Suite_norm"] = data[series_col].apply(normalize_suite)
    else:
        data["Suite_norm"] = "Data"

    # Keep rows with required data for panel A
    data = data.dropna(subset=[age_col, end_col]).copy()

    if data.empty:
        raise ValueError("No valid rows available after filtering missing age and εNd(t).")

    if colors is None:
        colors = _default_colors()

    if sizes is None:
        sizes = _default_sizes()

    detected = list(data["Suite_norm"].dropna().unique())

    if order is None:
        preferred = _default_order()
        order = [s for s in preferred if s in detected] + [s for s in detected if s not in preferred]

    # Add colors for unexpected groups, if any
    missing_color_groups = [s for s in detected if s not in colors]
    if missing_color_groups:
        cmap = plt.get_cmap("tab10", len(missing_color_groups))
        for i, grp in enumerate(missing_color_groups):
            colors[grp] = cmap(i)

    # Curves
    tempo_ma_dm, epsilon_nd_dm = dm_curve()
    tempo_ma_chur, epsilon_chur = chur_curve()

    # Figure
    fig, axs = plt.subplots(1, 2, figsize=figsize)
    fig.patch.set_facecolor("white")

    # ====================================================
    # Panel A
    # ====================================================
    ax = axs[0]
    ax.set_facecolor("white")

    # CHUR
    ax.plot(tempo_ma_chur, epsilon_chur, "k--", lw=1.3)

    # DM
    ax.plot(tempo_ma_dm, epsilon_nd_dm, color="royalblue", lw=2.2)

    ax.text(350, 1, "CHUR", fontsize=12)
    ax.text(850, 7, "DM", color="royalblue", fontsize=12)

    for suite in order:
        g = data[data["Suite_norm"] == suite]
        if len(g) == 0:
            continue

        ax.scatter(
            g[age_col],
            g[end_col],
            s=sizes.get(suite, 75),
            color=colors.get(suite, "gray"),
            edgecolor="black",
            linewidth=0.4,
            alpha=0.85,
            label=f"{_display_suite_name(suite)} (n={len(g)})",
            zorder=3,
        )

    # TDM projection lines
    if show_projection and tdm_col in data.columns:
        mask = data[tdm_col].notna() & (data[tdm_col] >= data[age_col])

        for _, row in data[mask].iterrows():
            x0 = row[age_col]
            y0 = row[end_col]
            x1 = row[tdm_col]
            y1 = eNd_DM(x1)

            ax.plot(
                [x0, x1],
                [y0, y1],
                "--",
                lw=0.7,
                alpha=0.18,
                color=colors.get(row["Suite_norm"], "gray"),
            )

    # Fixed x-range and cleaner tick spacing
    ax.set_xlim(0, 3500)
    ax.set_ylim(-30, 10)

    ax.set_xlabel("U-Pb age (Ma)")
    ax.set_ylabel("εNd(t)")

    ax.xaxis.set_major_locator(mticker.MultipleLocator(500))
    ax.xaxis.set_minor_locator(mticker.MultipleLocator(100))
    ax.yaxis.set_major_locator(mticker.MultipleLocator(5))
    ax.grid(True, alpha=0.12)

    if panel_labels:
        ax.text(0.02, 0.95, "(A)", transform=ax.transAxes, fontweight="bold", fontsize=18)

    # ====================================================
    # Panel B
    # ====================================================
    ax = axs[1]
    ax.set_facecolor("white")

    data_b = data.dropna(subset=[sr_col]).copy()

    for suite in order:
        g = data_b[data_b["Suite_norm"] == suite]
        if len(g) == 0:
            continue

        ax.scatter(
            g[sr_col],
            g[end_col],
            s=sizes.get(suite, 75),
            color=colors.get(suite, "gray"),
            edgecolor="black",
            linewidth=0.4,
            alpha=0.85,
            label=f"{_display_suite_name(suite)} (n={len(g)})",
            zorder=3,
        )

    ax.axhline(0, color="black", ls="--", lw=1.2)

    if show_boxes:
        def box(ax_obj, xmin, xmax, ymin, ymax, color, label):
            rect = plt.Rectangle(
                (xmin, ymin),
                xmax - xmin,
                ymax - ymin,
                facecolor=color,
                edgecolor="gray",
                alpha=0.10,
            )
            ax_obj.add_patch(rect)
            ax_obj.text(
                (xmin + xmax) / 2,
                (ymin + ymax) / 2,
                label,
                ha="center",
                va="center",
                fontsize=15,
                fontweight="bold",
            )

        box(ax, 0.7022, 0.7036, -2, 12, "#8fb6ff", "DM")
        box(ax, 0.7045, 0.7080, -22, -4, "#ffc89f", "EM1")
        box(ax, 0.7080, 0.7150, -22, -2, "#a8ddb5", "EM2")

    ax.set_xlim(0.701, 0.716)
    ax.set_ylim(-30, 10)

    ax.set_xlabel(r"$^{87}$Sr/$^{86}$Sr(t)")
    ax.set_ylabel("εNd(t)")
    ax.grid(True, alpha=0.12)

    if panel_labels:
        ax.text(0.02, 0.95, "(B)", transform=ax.transAxes, fontweight="bold", fontsize=18)

    # Unified legend with counts
    handles = []
    labels = []

    for suite in order:
        subset = data[data["Suite_norm"] == suite]
        n = len(subset)

        if n == 0:
            continue

        handles.append(
            Line2D(
                [0], [0],
                marker="o",
                linestyle="",
                markerfacecolor=colors.get(suite, "gray"),
                markeredgecolor="black",
                markersize=9,
            )
        )
        labels.append(f"{_display_suite_name(suite)} (n={n})")

    fig.legend(
        handles,
        labels,
        loc="center left",
        bbox_to_anchor=(0.84, 0.5),
        title="Magmatic suites",
    )

    plt.tight_layout(rect=[0, 0, 0.84, 1])

    if save_path is not None:
        fig.savefig(save_path, dpi=300, bbox_inches="tight")

    results = {
        "n_points_panel_a": len(data),
        "n_points_panel_b": len(data_b),
        "groups": [g for g in order if g in data["Suite_norm"].unique()],
        "columns_used": {
            "age_col": age_col,
            "end_col": end_col,
            "sr_col": sr_col,
            "tdm_col": tdm_col if tdm_col != "__TDM__" else None,
            "series_col": series_col,
        },
    }

    return fig, results