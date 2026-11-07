from pathlib import Path
from Bio import SeqIO
from collections import Counter
from Bio.SeqUtils.ProtParam import ProteinAnalysis
import pandas as pd
import Bio

Base = Path(__file__).resolve().parent.parent
Beta_globin = SeqIO.read(Base / "data" / "Beta_globin.fasta", "fasta")
mutated_Beta_globin = SeqIO.read(Base / "data" / "mutant_hbb.fasta", "fasta")

h_seq = Beta_globin.seq
m_seq = mutated_Beta_globin.seq

h_count = str(Counter(h_seq))
m_count = Counter(m_seq)

print(f"Length of Reference Sequence: {len(h_seq)} \nLength of Mutant Sequence: {len(m_seq)}")
assert len(h_seq) == len(m_seq) == 147, "Unexpected sequence length"

def properties(seq):
    p = ProteinAnalysis(seq)
    return { 
        "Length_aa" : len(seq),
        "Molecular_weight_Da" : p.molecular_weight(),
        "Isoelectric_point_pI" : p.isoelectric_point(),
        "Charge_at_pH_7.4" : p.charge_at_pH(7.4),
        "GRAVY" : p.gravy(),

    }

ref_props = properties(h_seq)
mut_props = properties(m_seq)

rows = []
for name in ref_props:
    rows.append({
        "Property" : name,
        "Reference" : ref_props[name],
        "Mutant" : mut_props[name],
        "Change" : mut_props[name] - ref_props[name]
    })

df = pd.DataFrame(rows)
df = df.round(2)
df.to_csv(Base / "results" / "sequence_properties.csv", index=False)
print(df.to_string(index=False))
print(f"\nMethod: Biopython {Bio.__version__}  ProtParam (theoretical estimates)")