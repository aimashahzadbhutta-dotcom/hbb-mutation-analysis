from Bio import SeqIO
import pandas as pd
from itertools import zip_longest
from pathlib import Path
from Bio.Data import IUPACData

Base = Path(__file__).resolve().parent.parent

Beta_globin = SeqIO.read(Base / "data" / "Beta_globin.fasta", "fasta")
mutated_Beta_globin = SeqIO.read(Base / "data" / "mutant_hbb.fasta", "fasta")
h_seq = Beta_globin.seq
m_seq = mutated_Beta_globin.seq
print(f"Amino Acid at position 7 of healthy Beta Globin: {Beta_globin.seq[6]}")
print(f"Healthy Beta Globin Amino Acid sequence from positions 4 to 10: {Beta_globin.seq[3:10]}")


print("Reference Length: ", len(h_seq))
print("Mutant Length: ", len(m_seq))
Data = []
for e,(h, m) in enumerate(zip_longest(h_seq, m_seq), start=1):
    result = "Yes" if h == m else "No"
    Data.append({
    "Position" : e,
    "Reference" : h,
    "Mutant" : m,
    "Status" : result
})

df = pd.DataFrame(Data)
diff = df[df["Status"] == "No"]
diff_list = list(zip(diff["Position"], diff["Reference"], diff["Mutant"]))
if len(diff_list) == 1:
    Mutant_position, Reference_residue, Mutant_residue = diff_list[0]
else:
    Mutant_position = Reference_residue = Mutant_residue = "N/A"

valid = set(IUPACData.protein_letters)
checks = {
    "Reference length is 147": len(h_seq) == 147,
    "Mutant length is 147": len(m_seq) == 147,
    "Only standard amino acids": set(h_seq) <= valid and set(m_seq) <= valid,
    "Exactly one difference": len(diff) == 1,
    "Difference is E7V": diff_list == [(7, "E", "V")],
}
status = "PASS" if all(checks.values()) else "FAIL"

print("Number of different Amino Acids: ", len(diff))
print("Position of difference: ", Mutant_position)
print("Reference residue: ", Reference_residue)
print("Mutant residue: ", Mutant_residue)
print("Validation status:", status)

with open(Base / "results" / "mutation.txt", "w") as f:
    f.write(f"Amino Acid at position 7 of healthy Beta Globin: {Beta_globin.seq[6]}\n")
    f.write(f"Healthy Beta Globin Amino Acid sequence from positions 4 to 10: {Beta_globin.seq[3:10]}\n")
    f.write(f"Number of different Amino Acids: {len(diff)}\n")
    f.write(f"Position of difference: {Mutant_position}\n")
    f.write(f"Reference residue: {Reference_residue}\n")
    f.write(f"Mutant residue: {Mutant_residue}\n")
    f.write(f"Reference Length: {len(h_seq)}\n")
    f.write(f"Mutant Length: {len(m_seq)}")
    for name, ok in checks.items():
        f.write(f"{name}: {'PASS' if ok else 'FAIL'}\n")
        f.write(f"Overall validation status: {status}\n")