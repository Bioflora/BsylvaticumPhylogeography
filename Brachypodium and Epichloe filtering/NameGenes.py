import os

# Directory with statistics
stats_dir = "./Brachypodium/Stats"
output_file = "./Brachypodium/Name_Genes.txt"

# Output file
with open(output_file, "w") as out:
    # Table header
    out.write("Ind\tGenes\n")

    # Process each stats file
    for stats_file in os.listdir(stats_dir):
        if stats_file.endswith("_stats.txt"):  # Only process stats files
            file_path = os.path.join(stats_dir, stats_file)

            genes_less_than_10 = []
            with open(file_path, "r") as file:
                for line in file:
                    if line.startswith("Sequence Name"):  # Skip header
                        continue
                    if line.strip() == "" or line.startswith("Total Statistics"):
                        break

                    parts = line.strip().split("\t")
                    if len(parts) < 4:
                        continue

                    seq_name, _, _, fraction_n = parts
                    if float(fraction_n) < 0.1:
                        genes_less_than_10.append(seq_name)

            ind_name = stats_file.replace("_stats.txt", "")
            genes_list = ", ".join(genes_less_than_10)
            out.write(f"{ind_name}\t{genes_list}\n")

print(f"Processing complete. Results saved to: {output_file}")
