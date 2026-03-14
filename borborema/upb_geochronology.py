import matplotlib.pyplot as plt
from matplotlib.lines import Line2D


def plot_upb_ages(df, colors=None):

    fig, ax = plt.subplots(figsize=(10,8))

    df = df.copy()

    df["Suite"] = df["Suite"].astype(str).str.strip()
    df["Pluton"] = df["Pluton"].astype(str).str.strip()

    df = df.sort_values(by="Age (Ma)", ascending=False).reset_index(drop=True)

    df["index"] = range(len(df))

    # gerar cores automaticamente se não fornecidas
    suites = df["Suite"].unique()

    if colors is None:
        cmap = plt.get_cmap("tab10")
        colors = {suite: cmap(i) for i, suite in enumerate(suites)}

    # plot
    for _, row in df.iterrows():

        c = colors.get(row["Suite"], "gray")

        ax.errorbar(
            row["Age (Ma)"],
            row["index"],
            xerr=row["Error 2S"],
            fmt="o",
            color=c,
            ecolor=c,
            elinewidth=1.5,
            capsize=4
        )

    ax.set_xlabel("Zircon U–Pb age (Ma)")
    ax.set_ylabel("Pluton")

    ax.set_yticks(df["index"])
    ax.set_yticklabels(df["Pluton"])

    ax.invert_xaxis()

    ax.grid(True, linestyle="--", alpha=0.3)

    # legenda
    legend_elements = [
        Line2D(
            [0], [0],
            marker="o",
            linestyle="",
            label=suite,
            markerfacecolor=colors[suite],
            markeredgecolor="black",
            markersize=8
        )
        for suite in suites
    ]

    ax.legend(handles=legend_elements, title="Suite", loc="upper right")

    plt.tight_layout()

    return fig