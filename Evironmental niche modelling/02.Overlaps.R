#######################################################################################
# Niche Overlap Analysis – Dioscorea spp.
# Author: Miguel Campos
# Description:
# This script calculates ecological niche overlap, equivalency, and similarity tests 
# using PCA-based niche space with `ecospat` and performs temporal niche comparisons.
#######################################################################################

# Load required libraries
library(HH)         # Variance Inflation Factor
library(raster)     # Raster operations
library(dismo)      # Species distribution modeling
library(vegan)      # Ordinations
library(rgeos)      # Geometric operations
library(maxnet)
library(ecospat)    
library(devtools)
library(plyr)
library(car)
library(corrgram)
library(ENMTools)
library(terra)
library(sf)
library(ade4)

# Load custom functions
setwd("./NicheOverlaps")
source("functionsENM.R")

# Load occurrence data
data <- read.csv("W_Coordinates.tsv", sep = "\t", dec = ".")

# Load and crop bioclimatic variables
bioclim <- stack(list.files(path = "../bioclimaticas", pattern = '*.bil', full.names = TRUE))
bioclim <- crop(bioclim, extent(-50, 180, -20, 80))

# Select relevant variables
selected_vars <- c("bio7", "bio8", "bio15", "bio2", "bio19", "bio14")
bioclim <- brick(bioclim[[selected_vars]])

######################################
# Niche Overlap: Western vs Eastern
######################################

# Load predicted suitability values
pres.A <- na.omit(read.table("Maxent_West.tsv", header = TRUE, sep = "\t"))
pres.B <- na.omit(read.table("Maxent_East.tsv", header = TRUE, sep = "\t"))

# PCA of combined background
pca <- dudi.pca(rbind(pres.A, pres.B)[, selected_vars], scannf = FALSE, nf = 2)
scores <- pca$li

# Project presence and background into PCA space
scores.A.pres <- suprow(pca, pres.A[pres.A$presencia == 1, selected_vars])$li
scores.A.area <- suprow(pca, pres.A[, selected_vars])$li
scores.B.pres <- suprow(pca, pres.B[pres.B$presencia == 1, selected_vars])$li
scores.B.area <- suprow(pca, pres.B[, selected_vars])$li

# Create PCA environmental grids
grid.A <- ecospat.grid.clim.dyn(glob = scores, glob1 = scores.A.area, sp = scores.A.pres, R = 100, th.sp = 0)
grid.B <- ecospat.grid.clim.dyn(glob = scores, glob1 = scores.B.area, sp = scores.B.pres, R = 100, th.sp = 0)

# Overlap indices
D <- ecospat.niche.overlap(grid.A, grid.B, cor = TRUE)$D
I <- ecospat.niche.overlap(grid.A, grid.B, cor = TRUE)$I
print(paste("Schoener's D:", D))
print(paste("Hellinger’s I:", I))

# Plot overlap
ecospat.plot.niche.dyn(grid.A, grid.B, quant = 0.25, interest = 2,
                       title = "Niche Overlap West - East", name.axis1 = "PC1", name.axis2 = "PC2")

# Niche Equivalency Test
eq.test <- ecospat.niche.equivalency.test(grid.A, grid.B, rep = 1000, intersection = 0,
                                          overlap.alternative = "higher", expansion.alternative = "lower",
                                          stability.alternative = "higher", unfilling.alternative = "lower",
                                          ncores = 8)

pdf("Niche_Equivalency_West_East.pdf", width = 15, height = 10)
ecospat.plot.overlap.test(eq.test, "D", "Equivalency")
dev.off()

# Niche Similarity Test
sim.test <- ecospat.niche.similarity.test(grid.A, grid.B, rep = 1000, intersection = 0,
                                          overlap.alternative = "higher", expansion.alternative = "lower",
                                          stability.alternative = "higher", unfilling.alternative = "lower",
                                          ncores = 8)

pdf("Niche_Similarity_West_East.pdf", width = 15, height = 10)
ecospat.plot.overlap.test(sim.test, "D", "Similarity")
dev.off()

######################################
# Temporal Niche Overlap: Present vs LGM
######################################

# Load models for current and LGM
West <- raster("West_Model")
WestLGM <- raster("West_LGM")
West <- extend(West, extent(WestLGM))

East <- raster("East_Model")
EastLGM <- raster("East_LGM")
East <- extend(East, extent(EastLGM))

# Calculate Schoener’s D
D_West <- nicheOverlap(West, WestLGM, stat = "D")
D_East <- nicheOverlap(East, EastLGM, stat = "D")
print(paste("West Present vs LGM D:", D_West))
print(paste("East Present vs LGM D:", D_East))

######################################
# Niche Breadth
######################################

nicheBreadth(model.map = West)
nicheBreadth(model.map = WestLGM)
nicheBreadth(model.map = East)
nicheBreadth(model.map = EastLGM)
