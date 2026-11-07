from pathlib import Path
import pandas as pd
from Bio.PDB import PDBParser, NeighborSearch
from Bio.PDB.SASA import ShrakeRupley
from Bio.SeqUtils import seq1

Base = Path(__file__).resolve().parent.parent
struct_dir = Base / "data" / "structures" 

STRUCTURES = [
    ("2HBS_deoxy_HbS", struct_dir / "2HBS.pdb", "VAL"),
    ("2HHB_deoxy_normal", struct_dir / "2HHB.pdb", "GLU"),
]
CUTOFF = 4.0 


def chain_sequence(chain):
    return "".join(seq1(r.get_resname()) for r in chain if r.id[0] == " ")


def analyze(path, label, expected_res):
    structure = PDBParser(QUIET=True).get_structure(label, str(path))
    model = structure[0]

    # Removing water
    for chain in model:
        for res in list(chain):
            if res.id[0] == "W":
                chain.detach_child(res.id)

    # Identifying beta chains (mature beta starts VHLTP)
    beta_chains = []
    for chain in model:
        seq = chain_sequence(chain)
        print(f"[{label}] chain {chain.id}: {len(seq)} residues, starts {seq[:8]}")
        if seq.startswith("VHLTP"):
            beta_chains.append(chain)
    print(f"[{label}] beta chains: {[c.id for c in beta_chains]}")
    assert beta_chains, f"No beta chains found in {label}"

    ShrakeRupley().compute(model, level="R")

    atoms = [a for a in model.get_atoms() if a.get_parent().id[0] == " "]
    ns = NeighborSearch(atoms)

    summary, neighbors = [], []
    for chain in beta_chains:
        res6 = chain[6]  # PDB (mature) numbering: beta6
        name = res6.get_resname()
        assert name == expected_res, (
            f"{label} chain {chain.id}: residue 6 is {name}, expected {expected_res}"
        )

        near = set()
        for atom in res6:
            for other in ns.search(atom.coord, CUTOFF, level="R"):
                if other is res6:
                    continue
                same_chain = other.get_parent().id == chain.id
                if same_chain and abs(other.id[1] - 6) <= 1:
                    continue  
                near.add((other.get_parent().id, other.id[1], other.get_resname()))

        for ch_id, num, resname in sorted(near):
            neighbors.append({
                "Structure": label,
                "Beta_chain": chain.id,
                "Neighbor_chain": ch_id,
                "Neighbor_residue": num,
                "Neighbor_name": resname,
                "Same_chain": ch_id == chain.id,
            })

        summary.append({
            "Structure": label,
            "Beta_chain": chain.id,
            "Residue_6": name,
            "SASA_A2": round(res6.sasa, 1),
            "Neighbor_count": len(near),
        })

    return summary, neighbors


all_summary, all_neighbors = [], []
for label, path, expected in STRUCTURES:
    summary, neighbors = analyze(path, label, expected)
    all_summary += summary
    all_neighbors += neighbors

summary_df = pd.DataFrame(all_summary)
neighbors_df = pd.DataFrame(all_neighbors)

summary_df.to_csv(Base / "results" / "beta6_comparison.csv", index=False)
neighbors_df.to_csv(Base / "results" / "beta6_neighbors_comparison.csv", index=False)

print("\nSummary per beta chain:")
print(summary_df.to_string(index=False))
print("\nMean beta6 SASA by structure (A^2):")
print(summary_df.groupby("Structure")["SASA_A2"].mean().round(1))