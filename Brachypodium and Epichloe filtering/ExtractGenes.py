import os
from Bio import SeqIO

# Directory with input FASTA files
input_dir = "./Brachypodium/Fastas"
# Output directory
output_dir = "./Brachypodium/Genes"

if not os.path.exists(output_dir):
    os.makedirs(output_dir)

# Extract sequences and save them one per file
for fasta_file in os.listdir(input_dir):
    if fasta_file.endswith(".fasta"):
        input_path = os.path.join(input_dir, fasta_file)
        for record in SeqIO.parse(input_path, "fasta"):
            gene_name = record.id
            individual = fasta_file.replace(".fasta", "")
            output_path = os.path.join(output_dir, f"{individual}_{gene_name}.fasta")
            SeqIO.write(record, output_path, "fasta")

print("Genes extracted to:", output_dir)
