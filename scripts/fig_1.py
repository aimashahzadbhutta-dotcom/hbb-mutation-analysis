from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


# --------------------------------------------------
# 1. Set project paths and plotting style
# --------------------------------------------------

BASE_DIR = Path(__file__).resolve().parent.parent
RESULTS_DIR = BASE_DIR / "results"
FIGURES_DIR = RESULTS_DIR / "figures"

FIGURES_DIR.mkdir(parents=True, exist_ok=True)

SASA_FILE = RESULTS_DIR / "beta6_comparison.csv"
NEIGHBORS_FILE = RESULTS_DIR / "beta6_neighbors_comparison.csv"

# A consistent, publication-inspired style without requiring
# a particular Matplotlib preset to be installed.
plt.rcParams.update({
    "figure.dpi": 120,
    "savefig.dpi": 300,
    "font.family": "DejaVu Sans",
    "font.size": 10,
    "axes.titlesize": 13,
    "axes.labelsize": 10,
    "xtick.labelsize": 9,
    "ytick.labelsize": 9,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "axes.grid": True,
    "grid.alpha": 0.22,
    "grid.linestyle": "--",
})


# --------------------------------------------------
# 2. Load and validate the data
# --------------------------------------------------

sasa_df = pd.read_csv(SASA_FILE)
neighbors_df = pd.read_csv(NEIGHBORS_FILE)

required_sasa = {
    "Structure", "Beta_chain", "Residue_6",
    "SASA_A2", "Neighbor_count"
}
required_neighbors = {
    "Structure", "Beta_chain", "Neighbor_chain",
    "Neighbor_residue", "Neighbor_name", "Same_chain"
}

if not required_sasa.issubset(sasa_df.columns):
    raise ValueError(
        f"Missing columns in {SASA_FILE.name}: "
        f"{required_sasa - set(sasa_df.columns)}"
    )

if not required_neighbors.issubset(neighbors_df.columns):
    raise ValueError(
        f"Missing columns in {NEIGHBORS_FILE.name}: "
        f"{required_neighbors - set(neighbors_df.columns)}"
    )

if sasa_df.empty:
    raise ValueError("The beta6 comparison CSV contains no data.")

# Use readable labels for the plotted chains.
sasa_df["Label"] = (
    sasa_df["Structure"].astype(str)
    + " | Chain "
    + sasa_df["Beta_chain"].astype(str)
)

sasa_df = sasa_df.sort_values(["Structure", "Beta_chain"])


# --------------------------------------------------
# FIGURE 1: SASA dot plot
# --------------------------------------------------

fig, ax = plt.subplots(figsize=(10, 5.5))

x = np.arange(len(sasa_df))

# Use one color per structure.
structures = sasa_df["Structure"].unique()
colors = plt.get_cmap("tab10")(
    np.linspace(0, 1, max(len(structures), 1))
)
color_map = dict(zip(structures, colors))

for structure in structures:
    mask = sasa_df["Structure"] == structure

    ax.scatter(
        x[mask.to_numpy()],
        sasa_df.loc[mask, "SASA_A2"],
        s=115,
        color=color_map[structure],
        edgecolor="black",
        linewidth=0.7,
        label=structure,
        zorder=3,
    )

# Label each measured value.
for xpos, value in zip(x, sasa_df["SASA_A2"]):
    ax.annotate(
        f"{value:.1f}",
        (xpos, value),
        xytext=(0, 9),
        textcoords="offset points",
        ha="center",
        fontsize=8,
    )

ax.set_xticks(x)
ax.set_xticklabels(sasa_df["Label"], rotation=20, ha="right")
ax.set_ylabel("Solvent-accessible surface area (Å²)")
ax.set_xlabel("Structure and beta-globin chain")
ax.set_title("Solvent Accessibility of Beta-Globin Residue 6")
ax.grid(axis="x", visible=False)
ax.legend(title="Structure", frameon=False)
ax.margins(y=0.18)

fig.tight_layout()
fig.savefig(
    FIGURES_DIR / "beta6_sasa_comparison.png",
    bbox_inches="tight",
)
plt.close(fig)


# --------------------------------------------------
# FIGURE 2: Neighbor-count lollipop plot
# --------------------------------------------------

fig, ax = plt.subplots(figsize=(10, 5.5))

y = sasa_df["Neighbor_count"].to_numpy()

for xpos, value, structure in zip(
    x, y, sasa_df["Structure"]
):
    ax.vlines(
        xpos, 0, value,
        color=color_map[structure],
        linewidth=2,
        alpha=0.8,
    )

ax.scatter(
    x,
    y,
    s=110,
    c=[color_map[s] for s in sasa_df["Structure"]],
    edgecolor="black",
    linewidth=0.7,
    zorder=3,
)

for xpos, value in zip(x, y):
    ax.annotate(
        f"{value:.0f}",
        (xpos, value),
        xytext=(0, 7),
        textcoords="offset points",
        ha="center",
        fontsize=9,
    )

ax.set_xticks(x)
ax.set_xticklabels(sasa_df["Label"], rotation=20, ha="right")
ax.set_ylabel("Number of neighboring residues")
ax.set_xlabel("Structure and beta-globin chain")
ax.set_title(
    "Local Structural Environment of Beta-Globin Residue 6"
)
ax.grid(axis="x", visible=False)
ax.set_ylim(bottom=0)
ax.margins(y=0.15)

# Structure legend
handles = [
    plt.Line2D(
        [0], [0],
        marker="o",
        linestyle="",
        markerfacecolor=color_map[s],
        markeredgecolor="black",
        label=s,
    )
    for s in structures
]
ax.legend(handles=handles, title="Structure", frameon=False)

fig.tight_layout()
fig.savefig(
    FIGURES_DIR / "beta6_neighbor_comparison.png",
    bbox_inches="tight",
)
plt.close(fig)


# --------------------------------------------------
# FIGURE 3: Neighbor-residue identity heatmap
# --------------------------------------------------

if neighbors_df.empty:
    print(
        "Warning: neighbor CSV is empty; "
        "Figure 3 cannot be generated."
    )
else:
    # Count the occurrences of each neighboring amino-acid type
    # for each structure and beta chain.
    neighbors_df["Row_label"] = (
        neighbors_df["Structure"].astype(str)
        + " | beta chain "
        + neighbors_df["Beta_chain"].astype(str)
    )

    heatmap_data = pd.crosstab(
        neighbors_df["Row_label"],
        neighbors_df["Neighbor_name"],
    )

    # Keep rows in the same order as the SASA comparison.
    row_order = sasa_df["Label"].tolist()
    row_mapping = dict(zip(sasa_df["Label"], sasa_df["Label"]))

    # Match the labels used by the two CSV files.
    heatmap_data.index = heatmap_data.index.map(
        lambda label: label.replace(" | beta chain ", " | Chain ")
    )

    existing_rows = [
        label for label in row_order if label in heatmap_data.index
    ]
    heatmap_data = heatmap_data.reindex(existing_rows)

    # Sort amino-acid columns for easier reading.
    heatmap_data = heatmap_data.reindex(
        sorted(heatmap_data.columns),
        axis=1,
    )

    fig_width = max(10, 0.55 * len(heatmap_data.columns))
    fig_height = max(4.5, 0.65 * len(heatmap_data.index))

    fig, ax = plt.subplots(figsize=(fig_width, fig_height))

    image = ax.imshow(
        heatmap_data.to_numpy(),
        aspect="auto",
        cmap="YlGnBu",
        interpolation="nearest",
    )

    ax.set_xticks(np.arange(len(heatmap_data.columns)))
    ax.set_xticklabels(
        heatmap_data.columns,
        rotation=45,
        ha="right",
    )

    ax.set_yticks(np.arange(len(heatmap_data.index)))
    ax.set_yticklabels(heatmap_data.index)

    ax.set_xlabel("Neighbor amino-acid type")
    ax.set_ylabel("Structure and beta-globin chain")
    ax.set_title("Amino-Acid Composition of the Local Environment")

    # Add count labels to cells with nonzero counts.
    values = heatmap_data.to_numpy()
    for row in range(values.shape[0]):
        for col in range(values.shape[1]):
            count = values[row, col]
            if count > 0:
                ax.text(
                    col, row, str(count),
                    ha="center",
                    va="center",
                    fontsize=8,
                    color="black",
                )

    colorbar = fig.colorbar(image, ax=ax, pad=0.02)
    colorbar.set_label("Neighbor-residue count")

    fig.tight_layout()
    fig.savefig(
        FIGURES_DIR / "beta6_neighbor_heatmap.png",
        bbox_inches="tight",
    )
    plt.close(fig)


# --------------------------------------------------
# 3. Report the generated figures
# --------------------------------------------------

print("\nFigure generation complete.")
print(f"Figures directory: {FIGURES_DIR}")

for filename in [
    "beta6_sasa_comparison.png",
    "beta6_neighbor_comparison.png",
    "beta6_neighbor_heatmap.png",
]:
    filepath = FIGURES_DIR / filename
    if filepath.exists():
        print(f"Saved: {filename}")