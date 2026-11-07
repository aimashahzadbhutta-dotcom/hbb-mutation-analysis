from pathlib import Path
from collections import Counter

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np


# ------------------------------------------
# 1. Project paths
# ------------------------------------------

BASE_DIR = Path(__file__).resolve().parent.parent
FIGURES_DIR = BASE_DIR / "results" / "figures"
FIGURES_DIR.mkdir(parents=True, exist_ok=True)


# ------------------------------------------
# 2. Reference and mutant beta-globin sequences
# ------------------------------------------

# Full translated beta-globin sequence, including
# the initial methionine used for precursor numbering.
# The E -> V mutation is at position 7 in this sequence.

normal_sequence = (
    "MVHLTPEEKSAVTALWGKVNVDEVGGEALGRLLVVYPWTQRFFESFGDLSTPDAVMGNPK"
    "VKAHGKKVLGAFSDGLAHLDNLKGTFATLSELHCDKLHVDPENFRLLGNVLVCVLAHHF"
    "GKEFTPPVQAAYQKVVAGVANALAHKYH"
)

mutant_sequence = normal_sequence[:6] + "V" + normal_sequence[7:]


# ------------------------------------------
# 3. Calculate amino-acid composition
# ------------------------------------------

amino_acids = sorted(set(normal_sequence) | set(mutant_sequence))

normal_counts = Counter(normal_sequence)
mutant_counts = Counter(mutant_sequence)

comparison = pd.DataFrame({
    "Amino_Acid": amino_acids,
    "Reference_Count": [
        normal_counts[aa] for aa in amino_acids
    ],
    "Mutant_Count": [
        mutant_counts[aa] for aa in amino_acids
    ],
})

comparison["Reference_Percent"] = (
    comparison["Reference_Count"] / len(normal_sequence) * 100
)

comparison["Mutant_Percent"] = (
    comparison["Mutant_Count"] / len(mutant_sequence) * 100
)

comparison["Change_percentage_points"] = (
    comparison["Mutant_Percent"]
    - comparison["Reference_Percent"]
)

# Save the numerical comparison for reproducibility.
comparison.to_csv(
    BASE_DIR / "results" / "amino_acid_composition_comparison.csv",
    index=False,
)


# ------------------------------------------
# 4. Plot composition changes
# ------------------------------------------

comparison = comparison.sort_values(
    "Change_percentage_points"
)

fig, ax = plt.subplots(figsize=(10, 7))

changes = comparison["Change_percentage_points"].to_numpy()
labels = comparison["Amino_Acid"].to_numpy()

colors = [
    "#C44E52" if value < 0
    else "#4C72B0" if value > 0
    else "#B8B8B8"
    for value in changes
]

ax.barh(
    labels,
    changes,
    color=colors,
    edgecolor="white",
    height=0.7,
)

ax.axvline(0, color="black", linewidth=0.9)

ax.set_xlabel(
    "Change in amino-acid composition (percentage points)"
)
ax.set_ylabel("Amino-acid one-letter code")
ax.set_title(
    "Change in Beta-Globin Amino-Acid Composition\n"
    "Reference vs. Sickle-Cell Variant (E → V)"
)

# Keep the scale honest and make the small changes visible.
ax.set_xlim(-1.0, 1.0)
ax.grid(axis="x", linestyle="--", alpha=0.25)
ax.grid(axis="y", visible=False)

fig.tight_layout()

fig.savefig(
    FIGURES_DIR / "amino_acid_composition_comparison.png",
    dpi=300,
    bbox_inches="tight",
)

plt.close(fig)


# ------------------------------------------
# 5. Report results
# ------------------------------------------

print("Reference sequence length:", len(normal_sequence))
print("Mutant sequence length:", len(mutant_sequence))
print("\nAmino acids with changed composition:")

print(
    comparison.loc[
        comparison["Change_percentage_points"] != 0,
        [
            "Amino_Acid",
            "Reference_Count",
            "Mutant_Count",
            "Change_percentage_points",
        ],
    ].to_string(index=False)
)

print("\nFigure saved to:", FIGURES_DIR)
print("Numerical table saved to: results/amino_acid_composition_comparison.csv")

