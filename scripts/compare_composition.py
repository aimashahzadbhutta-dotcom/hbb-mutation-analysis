from Bio import SeqIO
from Bio.Data import IUPACData
from pathlib import Path
import pandas as pd
from collections import Counter

Base = Path(__file__).resolve().parent.parent
Beta_globin = SeqIO.read(Base / "data" / "Beta_globin.fasta", "fasta")
mutated_Beta_globin = SeqIO.read(Base / "data" / "mutant_hbb.fasta", "fasta")

h_seq = str(Beta_globin.seq)
m_seq = str(mutated_Beta_globin.seq)

count_h = Counter(h_seq)
count_m = Counter(m_seq)

rows = []
for aa in IUPACData.protein_letters:  # all 20, alphabetical
    h = count_h[aa]  # Counter gives 0 if the letter is missing
    m = count_m[aa]
    rows.append({
        "Amino_acid": aa,
        "Ref_count": h,
        "Mut_count": m,
        "Count_change": m - h,
        "Ref_percent": h / len(h_seq) * 100,
        "Mut_percent": m / len(m_seq) * 100,
    })

df = pd.DataFrame(rows)
df["Percent_change"] = (df["Mut_percent"] - df["Ref_percent"]).round(2)

# Checks
assert len(h_seq) == len(m_seq), "Sequence lengths differ"
assert df["Ref_count"].sum() == len(h_seq), "Reference counts don't add up"
assert df["Mut_count"].sum() == len(m_seq), "Mutant counts don't add up"
df = df.round(2)
df.to_csv(Base / "results" / "composition_comparison.csv", index=False)

print("Total length (ref, mutant):", len(h_seq), len(m_seq))
print(df[df["Count_change"] != 0])

