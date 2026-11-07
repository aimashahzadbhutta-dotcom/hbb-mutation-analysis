from pathlib import Path
from Bio.PDB import PDBParser, NeighborSearch
from Bio.PDB.SASA import ShrakeRupley
from Bio.SeqUtils import seq1
import pandas as pd

Base = Path(__file__).resolve().parent.parent
pdb_file = Base / "data" / "structures" / "2HBS.pdb"

structure = PDBParser(QUIET=True).get_structure("2HBS", str(pdb_file))
model = structure[0]

# Remove water molecules so they do not affect contacts or surface area
for chain in model:
    for res in list(chain):
        if res.id[0] == "W":
            chain.detach_child(res.id)


def chain_sequence(chain):
    # Standard amino-acid residues only (heme and other ligands are skipped)
    return "".join(seq1(r.get_resname()) for r in chain if r.id[0] == " ")


# Mature beta-globin starts VHLTPEEK / VHLTPVEK; alpha starts VLSPADK.

beta_chains, alpha_chains = [], []
for chain in model:
    seq = chain_sequence(chain)
    print(f"Chain {chain.id}: {len(seq)} residues, starts {seq[:8]}")
    if seq.startswith("VHLTP"):
        beta_chains.append(chain)
    elif seq.startswith("VLSP"):
        alpha_chains.append(chain)

print("Alpha chains:", [c.id for c in alpha_chains])
print("Beta chains:", [c.id for c in beta_chains])
assert beta_chains, "No beta chains found"

# Residue 6 in PDB (mature) numbering = beta6 = full-length position 7
rows = []
sr = ShrakeRupley()
sr.compute(model, level="R")  # solvent-accessible surface area per residue

all_atoms = [a for a in model.get_atoms() if a.get_parent().id[0] == " "]
ns = NeighborSearch(all_atoms)

for chain in beta_chains:
    res6 = chain[6]
    name = res6.get_resname()
    print(f"Chain {chain.id}: residue 6 is {name}, SASA = {res6.sasa:.1f} A^2")
    assert name == "VAL", f"Chain {chain.id} residue 6 is {name}, expected VAL"

    near = set()
    for atom in res6:
        for other in ns.search(atom.coord, 4.0, level="R"):
            if other is res6:
                continue
            same_chain = other.get_parent().id == chain.id
            if same_chain and abs(other.id[1] - 6) <= 1:
                continue  # skip covalent neighbors in the chain
            near.add((other.get_parent().id, other.id[1], other.get_resname()))

    for ch_id, num, resname in sorted(near):
        rows.append({
            "Beta_chain": chain.id,
            "Neighbor_chain": ch_id,
            "Neighbor_residue": num,
            "Neighbor_name": resname,
            "Same_chain": ch_id == chain.id,
            "Beta6_SASA_A2": round(res6.sasa, 1),
        })

df = pd.DataFrame(rows)
df.to_csv(Base / "results" / "beta6_environment.csv", index=False)
print(df.to_string(index=False))