###############################################
## Script for IBD/IBE Analysis
## Author: Miguel Campos
## Date: June 2025
###############################################

## Load required libraries
library(vegan)
library(vcfR)
library(adegenet)
library(spaMM)
library(geosphere)
library(car)
library(dplyr)

###############################################
## Load and process genetic data
###############################################
# Set working directory and load VCF file
setwd("./ddRDA")
vcf_file <- "./Bsylvaticum_Filtered.vcf"
vcf <- read.vcfR(vcf_file)
gen_data <- vcfR2genlight(vcf)

# Compute genetic distance matrix (Manhattan)
gen_dist <- vegdist(gen_data, method = "manhattan", na.rm = TRUE)

###############################################
## Load bioclimatic and geographic data
###############################################
env_data <- read.csv("Climatic.tsv", sep = '\t')

# Geographic distance matrix (based on longitude and latitude)
geo_dist <- distm(env_data[, c("Longitude", "Latitude")], fun = distHaversine)
rownames(geo_dist) <- env_data$Pop
colnames(geo_dist) <- env_data$Pop

# Environmental distance matrix
env_dist <- dist(env_data[, grep("bio", colnames(env_data))])
rownames(env_dist) <- env_data$Pop
colnames(env_dist) <- env_data$Pop


###############################################
## Export distance matrices
###############################################
write.csv(as.matrix(gen_dist), "gen_dist.tsv", row.names = TRUE)
write.csv(as.matrix(geo_dist), "geo_dist.tsv", row.names = TRUE)
write.csv(as.matrix(env_dist), "env_dist.tsv", row.names = TRUE)

## MANUALLY SPLIT THE DISTANCE MATRICES WITH EASTERN / WESTERN SAMPLES

###############################################
## EASTERN CLADE ANALYSIS
###############################################
gen_dist <- as.dist(read.csv("E_gen_dist.tsv", row.names = 1, sep = "\t"))
env_dist <- as.dist(read.csv("E_env_dist.tsv", row.names = 1, sep = "\t"))
geo_dist <- as.dist(read.csv("E_geo_dist.tsv", row.names = 1, sep = "\t"))

env_pcoa <- as.data.frame(cmdscale(env_dist, k = 10))
geo_pcoa <- as.data.frame(cmdscale(geo_dist, k = 10))
colnames(env_pcoa) <- paste0("env_PC", 1:10)
colnames(geo_pcoa) <- paste0("geo_PC", 1:10)
predictors <- cbind(env_pcoa, geo_pcoa)

mod_full <- capscale(gen_dist ~ ., data = predictors)
mod_env  <- capscale(gen_dist ~ ., data = env_pcoa)
mod_geo  <- capscale(gen_dist ~ ., data = geo_pcoa)

VarFull <- anova(mod_full, permutations = 9999)
VarIBE  <- anova(mod_env,  permutations = 9999)
VarIBD  <- anova(mod_geo,  permutations = 9999)

anova(mod_full, by = "terms", permutations = 9999)
anova(mod_env,  by = "terms", permutations = 9999)
anova(mod_geo,  by = "terms", permutations = 9999)

summary(mod_full)
summary(mod_env)
summary(mod_geo)

VarFull
VarIBE
VarIBD

anova.cca(mod_full, mod_geo, permutations = 9999)
anova.cca(mod_full, mod_env, permutations = 9999)


###############################################
## WESTERN CLADE ANALYSIS
###############################################
gen_dist <- as.dist(read.csv("W_gen_dist.tsv", row.names = 1, sep = "\t"))
env_dist <- as.dist(read.csv("W_env_dist.tsv", row.names = 1, sep = "\t"))
geo_dist <- as.dist(read.csv("W_geo_dist.tsv", row.names = 1, sep = "\t"))

env_pcoa <- as.data.frame(cmdscale(env_dist, k = 10))
geo_pcoa <- as.data.frame(cmdscale(geo_dist, k = 10))
colnames(env_pcoa) <- paste0("env_PC", 1:10)
colnames(geo_pcoa) <- paste0("geo_PC", 1:10)
predictors <- cbind(env_pcoa, geo_pcoa)

mod_full <- capscale(gen_dist ~ ., data = predictors)
mod_env  <- capscale(gen_dist ~ ., data = env_pcoa)
mod_geo  <- capscale(gen_dist ~ ., data = geo_pcoa)

VarFull <- anova(mod_full, permutations = 9999)
VarIBE  <- anova(mod_env,  permutations = 9999)
VarIBD  <- anova(mod_geo,  permutations = 9999)

anova(mod_full, by = "terms", permutations = 9999)
anova(mod_env,  by = "terms", permutations = 9999)
anova(mod_geo,  by = "terms", permutations = 9999)

summary(mod_full)
summary(mod_env)
summary(mod_geo)

VarFull
VarIBE
VarIBD

anova.cca(mod_full, mod_geo, permutations = 9999)
anova.cca(mod_full, mod_env, permutations = 9999)


###############################################
## Notes:
## - To interpret the percentage of explained variance, check the "constrained" value in the summary output.
## - Significance is assessed via permutation tests ("VarianzaXXX" objects).
## - This script follows a framework for testing Isolation by Distance (IBD) vs. Isolation by Environment (IBE).
###############################################
