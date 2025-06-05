# Biogeography and Dating analysis

In this section, you can find the scripts used for TreePL and BioGeoBEARS analysis

> This repository contains data and instructions to run scripts used in different analyses of *Brachypodium sylvaticum* data included in the paper:
>
> "A Palearctic divide, niche conservatism and host-fungal endophyte interactions shaped the phylogeography of the grass *Brachypodium sylvaticum*" by Mª Ángeles Decena & Miguel Campos.
>
> And co-authored by Diana Calderón Pardo, Valeriia Shiposha, Marina Olonova, Ernesto Pérez-Collazos and Pilar Catalán. 

## TreePL
For TreePL, we used the following command:

docker run --rm -v "$PWD":/TreePL naturalis/docker-treepl "/TreePL/Brachypodium_Config.txt"

And the files needed are the following: 
- `Bsylvaticum.tree` that contains the tree for the dating analysis obtained from IQTree2.
- `Bsylvaticum.trees` that contains the bootstrap trees obtained with -wbtl in IQTree2.
- `Brachypodium_Config.txt` that contains the config file for the main TreePL analysis.
- `01.TreePL_Bootstrap.sh` that contains the main script to perform the bootstrap replicates of TreePL.


## BioGeoBEARS
For ancestral area reconstruction using BioGeoBEARS, we need: 
- `02.SplitTrees.R` that contains the script to split the dated tree in Eastern and Western Clades.
- `Biogeobears_Geo_Eastern.txt` and `Biogeobears_Geo_Western.txt` that contains the geographic areas per individual.
- `03.ModelSelection.R` that contains the vcf file obtained from SNP calling and after filtering for LD.
- `04.BioGeoBEARS_DIVAJ.R` that contains the main script to run BioGeoBEARS, as we select the same model, you just need to use the same script for both trees, changing the config input and names.
