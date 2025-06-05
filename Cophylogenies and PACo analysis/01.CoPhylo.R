############################################################
# Co-phylogenetic Tree Visualization: Nuclear vs Plastome
# Author: Miguel Campos
# Description:
# This script generates a tanglegram comparing nuclear and plastome trees
# using a predefined association table.
############################################################

# Load required libraries
library(phytools)
library(TreeTools)
library(ape)
library(phangorn)

# Set working directory
setwd("./CoPhylo")

# Load association table (specimen mappings between trees)
assoc_table <- read.table("TableNP.tsv", sep = '\t', header = FALSE)

# Load and format nuclear tree
Nuclear <- read.tree("Nuclear.tree")
Nuclear <- ladderize(Nuclear, right = TRUE)
plot(Nuclear, cex = 0.5, use.edge.length = FALSE)

# Load and format plastome tree
Plastome <- read.tree("Plastome.tree")
Plastome <- ladderize(Plastome, right = TRUE)
plot(Plastome, cex = 0.5, use.edge.length = FALSE)

# Create co-phylogenetic object
cophylo_obj <- cophylo(Nuclear, Plastome, assoc = assoc_table, use.edge.length = FALSE)

# Plot cophylogeny to PDF
pdf("CoPhylo_NuclearPlastome.pdf", width = 10, height = 10)
plot(
  cophylo_obj,
  link.type = "curved",
  link.lwd = 1.5,
  link.lty = "solid",
  link.col = make.transparent("darkgreen", 0.7),
  fsize = 0.5
)
title(main = "Nuclear vs Plastome", font.main = 1)
dev.off()
