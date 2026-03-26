import re
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

from .data_cleaning import clean_numeric_series, strip_accents


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


# ============================================================
# Main plotting functions
# ============================================================

def plot_upb_ages(
    df,
    age_col="Age (Ma)",
    error_col="Error 2S",
    pluton_col="Pluton",
    suite_col="Suite",
    sample_col=None,
    colors=None,
    order=None,
    figsize=(16, 8),
    save_path=None,
):
    """
    Plot zircon U-Pb ages with horizontal error bars.
    """
    data = df.copy()

    # Normalize text columns
    data[pluton_col] = data[pluton_col].astype(str).str.strip()
    data[suite_col] = data[suite_col].apply(normalize_suite)

    # Clean numeric columns
    data[age_col] = clean_numeric_series(data[age_col], below_detection="value")
    data[error_col] = clean_numeric_series(data[error_col], below_detection="value")

    # Remove incomplete rows
    data = data.dropna(subset=[age_col, error_col, pluton_col, suite_col]).copy()

    if data.empty:
        raise ValueError("No valid U-Pb rows available after cleaning required columns.")

    if colors is None:
        colors = _default_colors()

    detected = list(data[suite_col].dropna().unique())

    if order is None:
        preferred = _default_order()
        order = [s for s in preferred if s in detected] + [s for s in detected if s not in preferred]

    # Add colors for unexpected groups, if any
    missing_color_groups = [s for s in detected if s not in colors]
    if missing_color_groups:
        cmap = plt.get_cmap("tab10", len(missing_color_groups))
        for i, grp in enumerate(missing_color_groups):
            colors[grp] = cmap(i)

    # Sort by age descending so older ages appear at the top
    data = data.sort_values(by=age_col, ascending=False).reset_index(drop=True)
    data["index"] = range(len(data))

    fig, ax = plt.subplots(figsize=figsize)
    fig.patch.set_facecolor("white")
    ax.set_facecolor("white")

    # Plot by suite
    for suite in order:
        g = data[data[suite_col] == suite]
        if len(g) == 0:
            continue

        ax.errorbar(
            g[age_col],
            g["index"],
            xerr=g[error_col],
            fmt="o",
            color=colors.get(suite, "gray"),
            ecolor=colors.get(suite, "gray"),
            markerfacecolor=colors.get(suite, "gray"),
            markeredgecolor="black",
            markeredgewidth=0.5,
            elinewidth=1.4,
            capsize=3,
            markersize=10,
            linestyle="",
            zorder=3,
        )

    ax.set_xlabel("Zircon U-Pb age (Ma)")
    ax.set_ylabel("Pluton")

    ax.set_yticks(data["index"])
    ax.set_yticklabels(data[pluton_col])

    ax.invert_xaxis()
    ax.grid(True, linestyle="--", alpha=0.3)

    # Legend with counts
    handles = []
    for suite in order:
        g = data[data[suite_col] == suite]
        if len(g) == 0:
            continue

        handles.append(
            Line2D(
                [0], [0],
                marker="o",
                linestyle="",
                label=f"{_display_suite_name(suite)} (n={len(g)})",
                markerfacecolor=colors.get(suite, "gray"),
                markeredgecolor="black",
                markersize=10,
            )
        )

    # Legend outside the plot
    ax.legend(
        handles=handles,
        title="Magmatic suites",
        loc="center left",
        bbox_to_anchor=(1.04, 0.5),
        borderaxespad=0.0,
        frameon=True,
    )

    # Layout adjusted so legend is not cropped
    plt.tight_layout(rect=[0, 0, 0.72, 1])

    if save_path is not None:
        fig.savefig(save_path, dpi=300, bbox_inches="tight")

    results = {
        "n_points": len(data),
        "groups": [g for g in order if g in data[suite_col].unique()],
        "age_min": float(data[age_col].min()),
        "age_max": float(data[age_col].max()),
    }

    return fig, results


def run_upb_from_dataframe(
    df,
    age_col=None,
    error_col=None,
    pluton_col=None,
    suite_col=None,
    sample_col=None,
    colors=None,
    order=None,
    figsize=(16, 8),
    save_path=None,
    below_detection="value",
):
    """
    High-level U-Pb plotting workflow from a dataframe with auto-detected columns.
    """
    # Auto-detect columns if not provided
    if age_col is None:
        age_col = _find_column(df, ["Age (MA)", "Age (Ma)", "Age", "Idade", "UPb_Age (Ma)"])
    if error_col is None:
        error_col = _find_column(df, ["Error 2S", "Erro 2S", "2SE", "2σ", "2s", "UPb_Error 2S"])
    if pluton_col is None:
        pluton_col = _find_column(df, ["Pluton", "Plutons", "PLUTON", "Body", "Intrusion"])
    if suite_col is None:
        suite_col = _find_column(df, ["Suite", "Serie", "SERIE", "Series"])
    if sample_col is None:
        sample_col = _find_column(df, ["Sample_ID", "Sample ID", "Sample", "Amostra"], required=False)

    data = df.copy()

    # Clean numeric columns here for consistency
    data[age_col] = clean_numeric_series(data[age_col], below_detection=below_detection)
    data[error_col] = clean_numeric_series(data[error_col], below_detection=below_detection)

    fig, results = plot_upb_ages(
        df=data,
        age_col=age_col,
        error_col=error_col,
        pluton_col=pluton_col,
        suite_col=suite_col,
        sample_col=sample_col,
        colors=colors,
        order=order,
        figsize=figsize,
        save_path=save_path,
    )

    results["columns_used"] = {
        "age_col": age_col,
        "error_col": error_col,
        "pluton_col": pluton_col,
        "suite_col": suite_col,
        "sample_col": sample_col,
    }

    return fig, results