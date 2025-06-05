# Population Genomics

In this section you can find the scripts used for Linkage Disequilibrium Filtering, ADMIXTURE analysis and Isolation by Distance/Environment

> This repository contains data and instructions to run scripts used in different analyses of *Brachypodium sylvaticum* data included in the paper:
>
> "A Palearctic divide, niche conservatism and host-fungal endophyte interactions shaped the phylogeography of the grass *Brachypodium sylvaticum*" by Mª Ángeles Decena & Miguel Campos.
>
> And co-authored by Diana Calderón Pardo, Valeriia Shiposha, Marina Olonova, Ernesto Pérez-Collazos and Pilar Catalán. 

## Linkage Disequilibrium Filtering
For LD Filtering, we need: 
- `LD_Filtering.sh` that contains the main script to filter linked variants.
- `Bsylvaticum.vcf` that contains the vcf file obtained from SNP calling.

## ADMIXTURE
For the ADMIXTURE, we need: 
- `RunAdmixture.sh` that contains the main script to run the analysis.
- `Bsylvaticum_Filtered.vcf` that contains the vcf file obtained from SNP calling and after filtering for LD.

## IBD / IBE
For the IBD/IBE analysis, we need: 
- `dbRDA.R` that contains the main script to run the analysis.
- `Climatic.tsv` that contains the climatic variables for each sample extracted from the ENM analysis.
- `Bsylvaticum_Filtered.vcf` that contains the vcf file obtained from SNP calling and after filtering for LD.
