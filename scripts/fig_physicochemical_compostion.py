from pathlib import Path

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from Bio.SeqUtils.ProtParam import ProteinAnalysis


# ------------------------------------------
# 1. Project paths
# ------------------------------------------

BASE_DIR = Path(__file__).resolve().parent.parent
RESULTS_DIR = BASE_DIR / "results"
FIGURES_DIR = RESULTS_DIR / "figures"

FIGURES_DIR.mkdir(parents=True, exist_ok=True)


# ------------------------------------------
# 2. Reference and mutant sequences
# ------------------------------------------

# Same precursor sequence convention as Figure 2.
# Position 7: glutamic acid (E) -> valine (V).

normal_sequence = (
    "MVHLTPEEKSAVTALWGKVNVDEVGGEALGRLLVVYPWTQRFFESFGDLSTPDAVMGNPK"
    "VKAHGKKVLGAFSDGLAHLDNLKGTFATLSELHCDKLHVDPENFRLLGNVLVCVLAHHF"
    "GKEFTPPVQAAYQKVVAGVANALAHKYH"
)

mutant_sequence = normal_sequence[:6] + "V" + normal_sequence[7:]

if len(normal_sequence) != len(mutant_sequence):
    raise ValueError("Reference and mutant sequence lengths differ.")

if normal_sequence[6] != "E" or mutant_sequence[6] != "V":
    raise ValueError("Expected the E7V substitution was not found.")


# ------------------------------------------
# 3. Calculate physicochemical properties
# ------------------------------------------

def calculate_properties(sequence):
    analysis = ProteinAnalysis(sequence)

    return {
        "Molecular weight (Da)": analysis.molecular_weight(),
        "Theoretical pI": analysis.isoelectric_point(),
        "Instability index": analysis.instability_index(),
        "GRAVY hydropathy index": analysis.gravy(),
        "Aromaticity (fraction)": analysis.aromaticity(),
    }


normal_properties = calculate_properties(normal_sequence)
mutant_properties = calculate_properties(mutant_sequence)

df = pd.DataFrame({
    "Property": normal_properties.keys(),
    "Reference": normal_properties.values(),
    "Mutant": mutant_properties.values(),
})

df["Change (mutant - reference)"] = (
    df["Mutant"] - df["Reference"]
)

df.to_csv(
    RESULTS_DIR / "physicochemical_comparison.csv",
    index=False,
)


# ------------------------------------------
# 4. Create a multi-panel comparison figure
# ------------------------------------------

properties = df["Property"].tolist()

fig, axes = plt.subplots(
    nrows=len(properties),
    ncols=1,
    figsize=(10, 15),
)

colors = ["#4477AA", "#CC6677"]

for ax, (_, row) in zip(axes, df.iterrows()):
    values = [row["Reference"], row["Mutant"]]
    bars = ax.barh(
        ["Reference", "E7V mutant"],
        values,
        color=colors,
        height=0.55,
    )

    ax.set_title(row["Property"], loc="left", fontweight="bold")
    ax.set_xlabel("Calculated value")
    ax.grid(axis="x", alpha=0.25)
    ax.set_axisbelow(True)

    # Add numerical labels without forcing the two
    # different properties onto a shared scale.
    for bar, value in zip(bars, values):
        ax.annotate(
            f"{value:.4f}" if abs(value) < 1 else f"{value:.2f}",
            xy=(bar.get_width(), bar.get_y() + bar.get_height() / 2),
            xytext=(5, 0),
            textcoords="offset points",
            va="center",
            fontsize=9,
        )

    # Include zero for properties that may be negative.
    lower = min(values)
    upper = max(values)
    span = upper - lower
    padding = max(span * 0.25, abs(upper) * 0.03, 0.01)

    if lower < 0:
        ax.set_xlim(lower - padding, upper + padding)
    else:
        ax.set_xlim(0, upper + padding)

fig.suptitle(
    "Physicochemical Properties of Beta-Globin\n"
    "Reference vs. E7V Variant",
    fontsize=15,
    fontweight="bold",
    y=1.01,
)

fig.tight_layout()

fig.savefig(
    FIGURES_DIR / "physicochemical_comparison.png",
    dpi=300,
    bbox_inches="tight",
)

plt.close(fig)


# ------------------------------------------
# 5. Print the calculated results
# ------------------------------------------

print("\nPhysicochemical comparison:")
print(df.to_string(index=False))

print("\nFigure saved to:")
print(FIGURES_DIR / "physicochemical_comparison.png")

print("\nTable saved to:")
print(RESULTS_DIR / "physicochemical_comparison.csv")