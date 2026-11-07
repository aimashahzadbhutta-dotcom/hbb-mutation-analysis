import pandas as pd
from pathlib import Path

Base = Path(__file__).resolve().parent.parent
my_dict = {
    "Properties" : ["Side-chain character", "Hydrophilic/ Hydrophobic", "Side-chain Structure", "Charge near physiological pH(~7.4)", "Side-chain pKa", "Residue mass", "Kyte-Doolittle hydropathy"],
    "Glutamic Acid(E)" : ["Acidic, Polar, Ionizable", "Hydrophilic", "Carboxylic acid group", "-1", "4.1 - 4.3", "129.12 Da", "-3.5"],
    "Valine(V)" : ["Non-polar, Aliphatic", "Hydrophobic", "Aliphatic isopropyl group", "0", "Not ionizable", "99.13 Da", "4.2"],
}

df = pd.DataFrame(my_dict)
df.to_csv(Base / "results" / "physicochemical_table.csv", index=False)
print(df)