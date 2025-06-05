import os

stats_dir = "./Epichloe/Stats"
output_file = "./Epichloe/summary_stats.txt"

with open(output_file, "w") as out:
    out.write("Ind\tGenes50\tGenes25\tGenes10\n")
    for stats_file in os.listdir(stats_dir):
        if stats_file.endswith(".txt"):
            with open(os.path.join(stats_dir, stats_file), "r") as f:
                lines = f.readlines()
                if len(lines) >= 3:
                    g50 = lines[-3].split(":")[1].strip()
                    g25 = lines[-2].split(":")[1].strip()
                    g10 = lines[-1].split(":")[1].strip()
                    out.write(f"{stats_file}\t{g50}\t{g25}\t{g10}\n")

print("Summary file created at:", output_file)
