# Environmental Niche Modelling, Equivalence, Identity and Overlaps

In this section, you can find the scripts used for ENM, and compare niches

> This repository contains data and instructions to run scripts used in different analyses of *Brachypodium sylvaticum* data included in the paper:
>
> "A Palearctic divide, niche conservatism and host-fungal endophyte interactions shaped the phylogeography of the grass *Brachypodium sylvaticum*" by Mª Ángeles Decena & Miguel Campos.
>
> And co-authored by Diana Calderón Pardo, Valeriia Shiposha, Marina Olonova, Ernesto Pérez-Collazos and Pilar Catalán. 

## Environmental Niche Modelling (ENM)
For ENM, we need:
- `Coordinates.tsv` that contains the coordinates to run the niche modelling.
- `01.Bsylvaticum_ENM.R` that contains the main script to run ENM and the projections to LGM (change inputs for Eastern and Western sub-enms.

## Overlaps, niche breadth, similarity, and equivalence
For niche comparison, we need: 
- `02.Overlaps.R` that contains the script to compare all the niches.
- `Coordinates.tsv` that contains the coordinates used before.
- `Maxent_West.tsv` and `Maxent_East.tsv` that contain the presence/absence datasets for Western vs Eastern comparisons.
- `Model` and `Model_LGM` that contain the models obtained from the ENM before (Present and LGM).
- `02.Overlaps.R` that contains the script to compare all the niches.


