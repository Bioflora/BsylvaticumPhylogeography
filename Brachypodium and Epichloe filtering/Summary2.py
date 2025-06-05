import os

# Directory with stats files
stats_dir = "./Brachypodium/Stats2"
output_file = "./Brachypodium/summary_stats.txt"

# Output file
with open(output_file, "w") as out:
    out.write("Ind\tGenes50\tGenes25\tGenes0\n")

    for stats_file in os.listdir(stats_dir):
        if stats_file.endswith(".txt"):
            file_path = os.path.join(stats_dir, stats_file)

            with open(file_path, "r") as file:
                lines = file.readlines()
                if len(lines) >= 3:
                    less_than_50 = lines[-3].split(":")[1].strip()
                    less_than_25 = lines[-2].split(":")[1].strip()
                    less_than_0 = lines[-1].split(":")[1].strip()

                    out.write(f"{stats_file}\t{less_than_50}\t{less_than_25}\t{less_than_0}\n")

print(f"Summary generated at: {output_file}")
