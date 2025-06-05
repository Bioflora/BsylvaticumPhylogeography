import os
from Bio import SeqIO

input_dir = "./Epichloe/Genes/"
output_dir = "./Epichloe/Genes_Filtered/"
os.makedirs(output_dir, exist_ok=True)

min_length = 500
max_missing = 0.5
min_seqs = 10

for fasta_file in os.listdir(input_dir):
    if fasta_file.endswith((".fasta", ".fa")):
        path = os.path.join(input_dir, fasta_file)
        records = []
        for record in SeqIO.parse(path, "fasta"):
            seq = str(record.seq).upper()
            if len(seq) > min_length and seq.count("N") / len(seq) <= max_missing:
                records.append(record)

        if len(records) >= min_seqs:
            SeqIO.write(records, os.path.join(output_dir, fasta_file), "fasta")
            print("Filtered file created:", fasta_file)
        else:
            print("Skipped (not enough sequences):", fasta_file)

print("Filtering complete. Output in:", output_dir)
