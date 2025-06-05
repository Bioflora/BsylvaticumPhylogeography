import os
from Bio import SeqIO

input_dir = "./Epichloe/Fastas/"
output_dir = "./Epichloe/Stats/"

for fasta_file in os.listdir(input_dir):
    if fasta_file.endswith((".fasta", ".fa")):
        fasta_path = os.path.join(input_dir, fasta_file)
        output_path = os.path.join(output_dir, f"{os.path.splitext(fasta_file)[0]}_stats.txt")

        with open(output_path, "w") as output_file:
            output_file.write("Sequence Name\tLength\tN Count\tFraction N\n")
            total_n = total_length = count_50 = count_25 = count_10 = 0

            for record in SeqIO.parse(fasta_path, "fasta"):
                seq = str(record.seq).upper()
                n_count = seq.count("N")
                length = len(seq)
                fraction_n = n_count / length if length > 0 else 0
                total_n += n_count
                total_length += length

                if fraction_n < 0.5: count_50 += 1
                if fraction_n < 0.25: count_25 += 1
                if fraction_n < 0.1: count_10 += 1

                output_file.write(f"{record.id}\t{length}\t{n_count}\t{fraction_n:.4f}\n")

            output_file.write("\nTotal Statistics:\n")
            output_file.write(f"Total Length: {total_length}\n")
            output_file.write(f"Total N Count: {total_n}\n")
            output_file.write(f"Fraction N (Total): {total_n / total_length:.4f}\n" if total_length > 0 else "N/A\n")
            output_file.write(f"Genes with <50% N: {count_50}\n")
            output_file.write(f"Genes with <25% N: {count_25}\n")
            output_file.write(f"Genes with <10% N: {count_10}\n")

print("Processing complete. Results saved in:", output_dir)
