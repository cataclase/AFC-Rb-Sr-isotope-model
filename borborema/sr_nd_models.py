import re
import unicodedata
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import ticker as mticker


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

def _strip_accents(text):
    text = unicodedata.normalize("NFKD", str(text))
    return "".join(ch for ch in text if not unicodedata.combining(ch))


def _normalize_header(text):
    """
    Normalize header names for robust matching.
    """
    text = _strip_accents(str(text)).lower().strip()
    text = text.replace("−", "-")
    text = re.sub(r"\s+", "", text)
    return text


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


def _clean_numeric(series):
    import pandas as pd
    return pd.to_numeric(
        series.astype(str)
        .str.replace("−", "-", regex=False)
        .str.replace(",", ".", regex=False)
        .str.replace(" ", "", regex=False)
        .str.replace(r"[^0-9eE\+\-\.]", "", regex=True),
        errors="coerce"
    )


def _default_colors():
    return {
        "CalcAlk": "#FFD700",
        "Shos": "green",
        "Alk": "lightpink",
        "ChAlk": "deeppink",
        "EHKCalcAlk": "orange",
        "PHKCalcAlk": "red"
    }


def _default_order():
    return [
        "Shos",
        "PHKCalcAlk",
        "EHKCalcAlk",
        "ChAlk",
        "CalcAlk",
        "Alk"
    ]


def _default_sizes():
    return {
        "Shos": 65,
        "PHKCalcAlk": 75,
        "EHKCalcAlk": 75,
        "ChAlk": 80,
        "CalcAlk": 90,
        "Alk": 90
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
    figsize=(14, 6),
    save_path=None,
    show_projection=True,
    show_boxes=True,
    panel_labels=True
):
    """
    Generate a two-panel Sm-Nd isotopic figure:
    (A) εNd(t) vs age
    (B) εNd(t) vs 87Sr/86Sr(t)

    Parameters
    ----------
    df : pandas.DataFrame
        Input dataframe.
    age_col : str, optional
        Column with U-Pb age (Ma). If None, tries to detect automatically.
    end_col : str, optional
        Column with εNd(t). If None, tries to detect automatically.
    sr_col : str, optional
        Column with 87Sr/86Sr(t). If None, tries to detect automatically.
    tdm_col : str, optional
        Column with TDM model age. If None, tries to detect automatically.
    series_col : str, optional
        Column with magmatic suite. If None, tries to detect automatically.
    colors : dict, optional
        Mapping of suites to colors.
    order : list, optional
        Plotting order of suites.
    sizes : dict, optional
        Marker sizes by suite.
    figsize : tuple, optional
        Figure size.
    save_path : str, optional
        Base path to save the figure, without extension.
    show_projection : bool, optional
        If True, plot TDM projection lines.
    show_boxes : bool, optional
        If True, show DM/EM1/EM2 isotopic fields.
    panel_labels : bool, optional
        If True, show (A) and (B).

    Returns
    -------
    fig : matplotlib.figure.Figure
    results : dict
    """
    # Auto-detect columns if not provided
    if age_col is None:
        age_col = _find_column(df, ["Age (MA)", "Age (Ma)", "Age", "Idade"])
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
    data[age_col] = _clean_numeric(data[age_col])
    data[end_col] = _clean_numeric(data[end_col])
    data[sr_col] = _clean_numeric(data[sr_col])

    if tdm_col is not None and tdm_col in data.columns:
        data[tdm_col] = _clean_numeric(data[tdm_col])
    else:
        data["__TDM__"] = np.nan
        tdm_col = "__TDM__"

    if series_col is not None and series_col in data.columns:
        data["Suite_norm"] = data[series_col].apply(normalize_series)
    else:
        data["Suite_norm"] = "Data"

    # Keep rows with required data
    data = data[data[age_col].notna() & data[end_col].notna()].copy()

    if data.empty:
        raise ValueError("No valid rows available after filtering missing age and εNd(t).")

    if colors is None:
        colors = _default_colors()

    if order is None:
        order = _default_order()

    if sizes is None:
        sizes = _default_sizes()

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
    ax.plot(tempo_ma_dm, epsilon_nd_dm, color="royalblue", lw=2)

    ax.text(350, 1, "CHUR")
    ax.text(850, 7, "DM", color="royalblue")

    for serie in order:
        g = data[data["Suite_norm"] == serie]
        if len(g) == 0:
            continue

        z = 1 if serie == "Shos" else 3
        alpha = 0.65 if serie == "Shos" else 0.95

        ax.scatter(
            g[age_col],
            g[end_col],
            s=sizes.get(serie, 70),
            color=colors.get(serie, "gray"),
            edgecolor="black",
            linewidth=0.4,
            alpha=alpha,
            label=f"{serie} (n = {len(g)})",
            zorder=z
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
                lw=0.6,
                alpha=0.20,
                color=colors.get(row["Suite_norm"], "gray")
            )

    ax.set_xlim(0, 3200)
    ax.set_ylim(-30, 10)

    ax.set_xlabel("U–Pb age (Ma)")
    ax.set_ylabel("εNd(t)")

    ax.xaxis.set_major_locator(mticker.MultipleLocator(500))
    ax.yaxis.set_major_locator(mticker.MultipleLocator(5))
    ax.grid(True, alpha=0.12)

    if panel_labels:
        ax.text(0.02, 0.95, "(A)", transform=ax.transAxes, fontweight="bold", fontsize=18)

    # ====================================================
    # Panel B
    # ====================================================
    ax = axs[1]
    ax.set_facecolor("white")

    # Only rows with Sr for panel B
    data_b = data[data[sr_col].notna()].copy()

    for serie in order:
        g = data_b[data_b["Suite_norm"] == serie]
        if len(g) == 0:
            continue

        z = 1 if serie == "Shos" else 3
        alpha = 0.65 if serie == "Shos" else 0.95

        ax.scatter(
            g[sr_col],
            g[end_col],
            s=sizes.get(serie, 70),
            color=colors.get(serie, "gray"),
            edgecolor="black",
            linewidth=0.4,
            alpha=alpha,
            label=serie,
            zorder=z
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
                alpha=0.10
            )
            ax_obj.add_patch(rect)

            ax_obj.text(
                (xmin + xmax) / 2,
                (ymin + ymax) / 2,
                label,
                ha="center",
                va="center",
                fontsize=15,
                fontweight="bold"
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

    # Unified legend
    handles, labels = ax.get_legend_handles_labels()
    bylabel = dict(zip(labels, handles))

    fig.legend(
        bylabel.values(),
        bylabel.keys(),
        loc="center left",
        bbox_to_anchor=(0.86, 0.5),
        title="Magmatic suites"
    )

    plt.tight_layout(rect=[0, 0, 0.86, 1])

    if save_path is not None:
        fig.savefig(f"{save_path}.png", dpi=300, bbox_inches="tight")
        fig.savefig(f"{save_path}.pdf", dpi=300, bbox_inches="tight")
        fig.savefig(f"{save_path}.svg", dpi=300, bbox_inches="tight")

    results = {
        "n_points_panel_a": len(data),
        "n_points_panel_b": len(data_b),
        "groups": [g for g in order if g in data["Suite_norm"].unique()],
        "columns_used": {
            "age_col": age_col,
            "end_col": end_col,
            "sr_col": sr_col,
            "tdm_col": tdm_col if tdm_col != "__TDM__" else None,
            "series_col": series_col
        }
    }

    return fig, results