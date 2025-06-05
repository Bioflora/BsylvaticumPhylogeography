#######################################################################################
# Ecological Niche Modeling (ENM) Pipeline using MaxEnt
# Author: Miguel Campos
# Description:
# This script performs general ENM using MaxEnt with bootstrapping,
# VIF-based variable selection, thresholding, and temporal projections to the LGM.
#######################################################################################

####################################
## Load required libraries
####################################
library(plotmo)
library(dismo) 
library(raster) 
library(rgeos)
library(HH) 
library(vegan)
library(rgdal)
library(ecospat) 
library(devtools)
library(ENMTools)
library(maxnet)
library(dplyr)
library(maptools)

setwd("./ENM/All")
source("functionsENM.R")  # Includes bootstrapMaxent()

####################################
## Load and clean occurrence data
####################################
data <- read.csv("Coordinates.tsv", sep = "\t", dec = ".")
duplicates <- duplicated(data[, c("Lat", "Lon")])
data <- data[!duplicates, ]
presences <- dplyr::select(data, -Ind)

####################################
## Load bioclimatic variables
####################################
bioclim <- stack(list.files(path = "../bioclimaticas", pattern = '*.bil', full.names = TRUE))
bioclim <- crop(bioclim, extent(-25, 190, -50, 80))

pdf("Presences.pdf", width = 15, height = 10)
plot(bioclim[[1]], main = "Presences on bio1")
points(presences)
dev.off()

####################################
## Variable selection
####################################
bioclim_df <- na.omit(as.data.frame(bioclim))
cor_matrix <- cor(bioclim_df)
abs_dist <- as.dist(abs(cor_matrix))
cluster <- hclust(1 - abs_dist)

pdf("Bioclimatic_correlation.pdf", width = 15, height = 10)
plot(cluster)
dev.off()

selected_vars <- c("bio7", "bio8", "bio15", "bio2", "bio12", "bio17")
vif_df <- na.omit(bioclim_df[, selected_vars])
vif(vif_df)  # Ensure all < 5
bioclim <- brick(bioclim[[selected_vars]])

####################################
## Create background points
####################################
names(presences) <- c("x", "y")
coordinates(presences) <- ~x + y
background <- randomPoints(bioclim[[1]], n = 20000, p = presences, excludep = TRUE)

####################################
## Prepare input matrix
####################################
sp <- list()

# Presence
sp$presence <- data.frame(
  coordinates(presences),
  extract(bioclim, presences),
  presencia = 1
)

# Background
sp$background <- data.frame(
  background,
  extract(bioclim, background),
  presencia = 0
)

pb <- rbind(sp$presence, sp$background)
pb <- pb[, c(selected_vars, "presencia")]
pb <- na.omit(pb)

write.csv(pb, "maxent.csv", row.names = FALSE)

####################################
## MaxEnt model fitting
####################################
pb <- read.csv("maxent.csv")
maxent_model <- maxnet(p = pb$presencia, data = pb[, 1:5])

# Predictor importance
data.frame(importance = sort(abs(maxent_model$betas)[1:9], decreasing = TRUE))

# Response curves
pdf("ResponseCurves.pdf", width = 15, height = 10)
plot(maxent_model, type = "logistic")
dev.off()

# Prediction map
prediction_map <- predict(bioclim, maxent_model, type = "logistic")
pdf("Model.pdf", width = 15, height = 10)
plot(prediction_map, col = viridis::turbo(100))
dev.off()
writeRaster(prediction_map, filename = "ModelMap", format = "EHdr")

# AUC and threshold
pres_vals <- extract(prediction_map, presences)
pres_vals <- pres_vals[!is.na(pres_vals)]
bg_vals <- extract(prediction_map, background)
eval <- dismo::evaluate(p = pres_vals, a = bg_vals)
dismo::threshold(eval, sensitivity = 0.9)

####################################
## Bootstrap evaluation
####################################
bootstrap <- bootstrapMaxent(
  data = pb,
  columna.presencia = "presencia",
  columnas.variables = names(bioclim),
  iteraciones = 100,
  porcentaje.evaluacion = 25,
  variables = bioclim
)

bootstrap$auc.mean
bootstrap$auc.sd
plot(bootstrap$models.mean)
plot(bootstrap$models.sd)
writeRaster(bootstrap$models.mean, "ModelMean", format = "EHdr")
writeRaster(bootstrap$models.sd, "ModelSD", format = "EHdr")

# Apply threshold
Model <- bootstrap$models.mean
Threshold <- Model
Threshold[Threshold < 0.2568011] <- 0
plot(Threshold, main = "Model with Threshold 10%", col = viridis::turbo(100))
writeRaster(Threshold, filename = "ModelThreshold", format = "EHdr")

# Export for niche breadth
writeRaster(bioclim, filename = "bioclimaticas_mcampos_dec22.tif", options = "INTERLEAVE=BAND", overwrite = TRUE)
write.table(pb, file = "presencia.background_mcampos_dec22.csv", sep = ",", row.names = FALSE)
write.table(presences, file = "presencia_limpia_mcampos_dec22.csv", sep = ",", row.names = FALSE)

####################################
## Temporal projection to LGM
####################################
lgm_vars <- stack(list.files(path = "../lgmbioblimaticas", pattern = '*.tif', full.names = TRUE))
lgm_vars <- crop(lgm_vars, extent(-25, 190, -50, 80))
lgm_vars <- resample(lgm_vars, bioclim, method = "bilinear")
lgm_vars <- lgm_vars[[selected_vars]]

lgm_map <- predict(lgm_vars, maxent_model, type = "logistic")

# Bootstrap LGM
bootstrap_lgm <- bootstrapMaxent(
  data = pb,
  columna.presencia = "presencia",
  columnas.variables = names(lgm_vars),
  iteraciones = 100,
  porcentaje.evaluacion = 25,
  variables = lgm_vars
)

bootstrap_lgm$auc.mean
bootstrap_lgm$auc.sd

# Thresholded LGM
pres_lgm <- extract(bootstrap_lgm$models.mean, presences)
pres_lgm <- pres_lgm[!is.na(pres_lgm)]
aus_lgm <- extract(bootstrap_lgm$models.mean, background)
eval_lgm <- dismo::evaluate(p = pres_lgm, a = aus_lgm)
dismo::threshold(eval_lgm, sensitivity = 0.9)

Threshold_lgm <- bootstrap_lgm$models.mean
Threshold_lgm[Threshold_lgm < 0.04313257] <- 0

pdf("ModeloLGM_Threshold.pdf", width = 15, height = 10)
plot(Threshold_lgm, main = "LGM Model with Threshold 10%", col = viridis::turbo(100))
dev.off()
writeRaster(Threshold_lgm, filename = "ModelLGM_Threshold", format = "EHdr")
