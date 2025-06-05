# PACo, ParaFit and CoPhylogenies

In this section, you can find the scripts used for CoPhylo, and PACo/ParaFit analysis

> This repository contains data and instructions to run scripts used in different analyses of *Brachypodium sylvaticum* data included in the paper:
>
> "A Palearctic divide, niche conservatism and host-fungal endophyte interactions shaped the phylogeography of the grass *Brachypodium sylvaticum*" by Mª Ángeles Decena & Miguel Campos.
>
> And co-authored by Diana Calderón Pardo, Valeriia Shiposha, Marina Olonova, Ernesto Pérez-Collazos and Pilar Catalán. 

## Co-Phylogenies
To run the Co-Phylo script, we need:
- `01.CoPhylo.R` that contains the coordinates to run the niche modelling.
- `Nuclear.tree` and `Plastome.tree` that contains phylogenetic nuclear and plastomic trees.
- `TableNP.tsv` that contains the association table with Nuclear / Plastome identities.

## ParaFit
To run ParaFit we need: 
- `02.ParaFit.R` that contains the script to run ParaFit and asses coevolution in the branches with a 1000 bootstrap runs.
- `Bsylvaticum.tree` and `Esylvatica.tree` that contain the main nuclear and endophytic trees.
- `Matrix.tsv` that contains the associations matrix between host/parasite.
- `TableBE.tsv` that contains the association table with host/parasite identities.

## PACo
To run PACo we need: 
- `03.PACo.R` that contains the script to run PACo and asses coevolution in the branches using the main tree with 1000 permutations and using the 1000 bootstrap trees.
- `Bsylvaticum.tree` and `Esylvatica.tree` that contain the main nuclear and endophytic trees.
- `Bsylvaticum.ufboot` and `Esylvatica.ufboot` that contain the 1000 bootstrap nuclear and endophytic trees.
- `Matrix.tsv` that contains the associations matrix between host/parasite.
- `TableBE.tsv` that contains the association table with host/parasite identities.
