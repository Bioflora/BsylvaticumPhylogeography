import os
from Bio import SeqIO
from collections import defaultdict

input_dir = "./Epichloe/Fastas/"
output_dir = "./Epichloe/Genes/"
os.makedirs(output_dir, exist_ok=True)

genes_dict = defaultdict(list)

for fasta_file in os.listdir(input_dir):
    if fasta_file.endswith((".fasta", ".fa")):
        fasta_path = os.path.join(input_dir, fasta_file)
        individual = os.path.splitext(fasta_file)[0]
        for record in SeqIO.parse(fasta_path, "fasta"):
            genes_dict[record.id].append((individual, record.seq))

for gene, records in genes_dict.items():
    with open(os.path.join(output_dir, f"{gene}.fasta"), "w") as output_file:
        for indiv, seq in records:
            output_file.write(f">{indiv}\n{seq}\n")

print("Gene files generated in:", output_dir)
