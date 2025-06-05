# Brachypodium and Epichloe filtering

In this section, you can find the scripts used for read assembly and split in Brachypodium / Epichloe reads, posterior assembly to each reference genome, and filtering SNPs

> This repository contains data and instructions to run scripts used in different analyses of *Brachypodium sylvaticum* data included in the paper:
>
> "A Palearctic divide, niche conservatism and host-fungal endophyte interactions shaped the phylogeography of the grass *Brachypodium sylvaticum*" by Mª Ángeles Decena & Miguel Campos.
>
> And co-authored by Diana Calderón Pardo, Valeriia Shiposha, Marina Olonova, Ernesto Pérez-Collazos and Pilar Catalán. 

## Read mapping and sequence extraction
- `01.SplitReads.sh` that contains the main script that maps against a concatenated reference using Bsylvaticum (nuclear), Bsylvaticum (plastome), Esylvatica (endophyte) reference genomes. and splits into different subfolders.
- `02.Counts.py` and `06.Counts.py` that contain scripts to summarize the amount of genes based on different thresholds and generate a table (for Brachypodium first and second for Epichloe).
- `03.Summary.py` and `07.Summary.py` that contain scripts to summarize the table generated in ` 02/06.Counts.py`.
- `04.Extract_Genes.py` and `08.Extract_Genes.py` that contain scripts to split the individual fasta sequences into gene fastas.
- `05.FilterGenes.py` and `09.FilterGenes.py` that contain scripts to filter the fastas based on missing data and depth.


## SNP calling
- `10.SNP_Calling.sh` main script that uses the bam files to generate individual vcfs and then combine into a single VCF file.
