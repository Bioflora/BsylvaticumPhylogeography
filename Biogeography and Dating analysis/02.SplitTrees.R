###############################################
# Splitting a Dated Tree into Eastern and Western Clades
# Author: Miguel Campos
# Description:
# This script loads a dated phylogenetic tree and separates it
# into two subtrees corresponding to the Eastern and Western clades.
###############################################

library(ape)

# Load the dated phylogenetic tree
trees1 <- read.tree("Bsylvaticum_Dated.tree")

# Define taxon labels for each clade
Eastern <- c("Bbrev33H", "Bbrev34H", "Bkurilense", "Bmis66-3", "Bmis66-4", "Bmis67-2", "Bmis67-1")

Western <- c("Bglaucovirens", "Bspryginii", "Bsyl29H", "Bsyl30H", "Bsyl31H", "Bsyl32H", "Bsyl35H", "Bsyl36H",
  "Bsyl37H", "Bsyl38H", "Bsyl39H", "Bsyl41H", "Bsyl466-13", "Bsyl466-2", "Bsyl466-7b", "Bsyl467-10",
  "Bsyl467-2", "Bsyl467-7", "Bsyl467-9b", "Bsyl470-4", "Bsyl470-8b", "Bsyl470-9", "Bsyl476-10",
  "Bsyl476-11", "Bsyl476-14", "Bsyl476-9", "Bsyl477-10", "Bsyl477-11", "Bsyl477-1", "Bsyl477-3",
  "Bsyl500-2", "Bsyl500-3", "Bsyl500-6", "Bsyl501-1", "Bsyl501-5", "Bsyl501-6", "Bsyl501-7",
  "Bsyl502-2", "Bsyl502-3", "Bsyl502-4", "Bsyl502-5", "Bsyl505-2", "Bsyl505-4", "Bsyl505-6",
  "Bsyl506-2", "Bsyl506-3", "Bsyl506-6", "Bsyl508-2", "Bsyl508-3", "Bsyl508-6", "Bsyl54-3",
  "Bsyl54-9", "Bsyl550-10", "Bsyl550-1B", "Bsyl550-3B", "Bsyl550-8b", "Bsyl552-1", "Bsyl552-2",
  "Bsyl552-5", "Bsyl553-2B", "Bsyl553-3b", "Bsyl553-5b", "Bsyl554-1", "Bsyl554-3B", "Bsyl554-6",
  "Bsyl554-7", "Bsyl555-1", "Bsyl555-3B", "Bsyl555-5", "Bsyl555-8", "Bsyl557-2", "Bsyl557-7",
  "Bsyl59-1", "Bsyl59-2", "Bsyl59-4", "Bsyl59-5", "Bsyl62-8", "Bsyl63-2A", "Bsyl63-2B", "Bsyl63-4",
  "Bsyl63-5", "Bsyl72-1", "Bsyl72-3", "Bsyl72-4", "Bsyl73-1", "Bsyl73-3", "Bsyl73-4")

# Generate subtree for Eastern clade
tr1 <- drop.tip(trees1, Eastern)
write.tree(tr1, file = "Eastern_Dated.tree")

# Generate subtree for Western clade
tr2 <- drop.tip(trees1, Western)
write.tree(tr2, file = "Western_Dated.tree")
