from Bio import SeqIO
from Bio.Data import IUPACData
from collections import Counter

analysis = SeqIO.read("/home/aimashahzadbhutta/hbb_mutation_analysis/data/Beta_globin.fasta", "fasta")
print("Sequence ID: ", analysis.id)
print("Length of sequence: ", len(analysis.seq))
seq = analysis.seq
print("Sequence: ", seq)

unique_aa = set(seq)
valid_aa = set(IUPACData.protein_letters)
invalid = unique_aa - valid_aa

if invalid == set():
    print("Invalid amino acids: None")
else:
    print("Invalid amino acids: ", invalid)
print("Valid amino acids: ", valid_aa)
if len(invalid) == 0:
    print("Sequence is valid")
else:
    print("Sequence is invalid")

counts = Counter(seq)
print("Amino Acid count: ", counts)
with open("/home/aimashahzadbhutta/hbb_mutation_analysis/results/amino_acid_count.txt", "w") as f:
    for amino_acid, count in counts.items():
        f.write(f"{amino_acid}: {count}\n")

most_common_count = counts.most_common(1)[0][1]
all_most_common_counts = [(aa,num) for aa, num in counts.items() if num == most_common_count]

print("Most frequent Amino Acids:\n", all_most_common_counts)
with open("/home/aimashahzadbhutta/hbb_mutation_analysis/results/Most_frequent_count.txt", "w") as f:
    f.write("Most Frequent Count:\n")
    f.write(f"{all_most_common_counts}\n")

print("Percentage of each Amino Acid in Beta Globin: ")
with open("/home/aimashahzadbhutta/hbb_mutation_analysis/results/amino_acid_percentage.txt", "w") as f:
    f.write("Percentage of each Amino Acid in Beta Globin: \n")
    for aa, num in counts.items():
        print(f"{aa}:{round(((num / len(seq))*100), 2)}%")
        f.write(f"{aa}:{round(((num / len(seq))*100))}%\n")


    
