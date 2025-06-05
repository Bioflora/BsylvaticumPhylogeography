import os

stats_dir = "./Epichloe/Stats/"
output_file = "./Epichloe/Name_Genes.txt"

with open(output_file, "w") as out:
    out.write("Ind\tGenes\n")
    for stats_file in os.listdir(stats_dir):
        if stats_file.endswith("_stats.txt"):
            genes = []
            with open(os.path.join(stats_dir, stats_file), "r") as f:
                for line in f:
                    if line.startswith("Sequence Name") or "Total Statistics" in line: continue
                    parts = line.strip().split("\t")
                    if len(parts) >= 4 and float(parts[3]) < 0.1:
                        genes.append(parts[0])
            indiv = stats_file.replace("_stats.txt", "")
            out.write(f"{indiv}\t{', '.join(genes)}\n")

print("Gene lists generated at:", output_file)
