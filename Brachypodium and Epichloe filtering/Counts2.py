import os

# Count the number of genes (lines) in each FASTA file in the folder
fasta_dir = "./Brachypodium/Fastas"
output_file = "./Brachypodium/gene_counts.txt"

with open(output_file, "w") as out:
    out.write("Individual\tGeneCount\n")
    for file in os.listdir(fasta_dir):
        if file.endswith(".fasta"):
            path = os.path.join(fasta_dir, file)
            count = 0
            with open(path, "r") as f:
                for line in f:
                    if line.startswith(">"):
                        count += 1
            name = file.replace(".fasta", "")
            out.write(f"{name}\t{count}\n")

print("Gene counts written to", output_file)
