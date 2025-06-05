################################################################
# PACo Analysis – Procrustean Approach to Cophylogeny
# Author: Miguel Campos Cáceres
# Date: August 2024
# Description:
# This script performs PACo analysis using nuclear and endophytic trees,
# including bootstrap support, normalized residuals, clustering of outliers,
# and final visualization as a tanglegram.
################################################################

library(paco)
library(ape)
library(cluster)
library(gplots)
library(ggplot2)
library(phytools)
library(vegan)
library(PhylogeneticEM)
library(geiger)

setwd("./PACO")

############################################
# 1. Define helper functions
############################################

PACo.dV <- function(H.dist, P.dist, HP.bin) {
  HP.bin <- which(HP.bin > 0, arr.in = TRUE)
  H.PCo <- pcoa(sqrt(H.dist))$vectors
  P.PCo <- pcoa(sqrt(P.dist))$vectors
  list(
    H.PCo = H.PCo[HP.bin[, 1], ],
    P.PCo = P.PCo[HP.bin[, 2], ]
  )
}

D.wrapper <- function(n) {
  DH.add <- cophenetic(treeH[[n]])
  DP.add <- cophenetic(treeP[[n]])
  DH.top <- cophenetic(compute.brlen(treeH[[n]], 1))
  DP.top <- cophenetic(compute.brlen(treeP[[n]], 1))

  DH.add <- DH.add[rownames(NCP), rownames(NCP)]
  DP.add <- DP.add[rownames(NCP), rownames(NCP)]
  DH.top <- DH.top[rownames(NCP), rownames(NCP)]
  DP.top <- DP.top[rownames(NCP), rownames(NCP)]

  PACo.add <- PACo.dV(DH.add, DP.add, HP)
  PACo.top <- PACo.dV(DH.top, DP.top, HP)

  Proc.add <- procrustes(PACo.add$H.PCo, PACo.add$P.PCo)
  Proc.top <- procrustes(PACo.top$H.PCo, PACo.top$P.PCo)

  add.res <- residuals(Proc.add)
  top.res <- residuals(Proc.top)

  HostX <- Proc.add$X
  ParY <- Proc.add$Yrot
  colnamesPACo <- paste(rownames(HostX), rownames(ParY), sep = "_")

  write(add.res, file = "PACo_res_add.txt", ncolumns = NLinks, append = TRUE, sep = "\t")
  write(top.res, file = "PACo_res_top.txt", ncolumns = NLinks, append = TRUE, sep = "\t")
  write(colnamesPACo, "colnamesPACo.txt", ncolumns = NLinks, sep = "\t")
}

############################################
# 2. Load trees and association matrix
############################################

treeH <- read.newick("Bsylvaticum.ufboot")     # Host trees (bootstrap)
treeP <- read.newick("Esylvatica.ufboot")      # Parasite trees (bootstrap)
NCP <- as.matrix(read.table("Matrix.tsv", header = TRUE, sep = "\t", row.names = 1, check.names = FALSE))
HP <- NCP
NLinks <- sum(NCP)

# Run PACo on all replicates
lapply(1:length(treeH), D.wrapper)

############################################
# 3. Normalize residuals and plot summaries
############################################

colnamesPACo <- colnames(read.table("colnamesPACo.txt", header = TRUE))
pac.add <- read.table("PACo_res_add.txt", header = FALSE, col.names = colnamesPACo)
pac.top <- read.table("PACo_res_top.txt", header = FALSE, col.names = colnamesPACo)

m2A <- apply(pac.add, 1, sum)
pac.norm.add <- pac.add / m2A
m2T <- apply(pac.top, 1, sum)
pac.norm.top <- pac.top / m2T

write.table(pac.norm.add, "PACo_res_add_norm.txt")
mean_resid_add <- mean(pac.norm.add)
cat("Mean normalized residual (additive trees):", round(mean_resid_add, 4), "\n")

# Plot normalized residuals – Additive trees
pdf("PACo_Residuals_Additive.pdf")
mA <- apply(pac.norm.add, 2, median)
uCI <- apply(pac.norm.add, 2, quantile, probs = 0.975)
lCI <- apply(pac.norm.add, 2, quantile, probs = 0.025)
cols <- c("#A300FF", "#5BFF5B")[(mA > 1/NLinks) + 1]
barplot2(mA, main = "PACo residuals (additive trees)", xlab = "Association", ylab = "Normalized residual",
         col = cols, border = "lightgrey", names.arg = colnamesPACo, las = 2, cex.names = 0.5,
         plot.ci = TRUE, ci.l = lCI, ci.u = uCI, ci.color = "blue")
abline(h = 1/NLinks, col = "red")
dev.off()

# Plot normalized residuals – Topological trees
pdf("PACo_Residuals_Topological.pdf")
mT <- apply(pac.norm.top, 2, median)
uCI <- apply(pac.norm.top, 2, quantile, probs = 0.975)
lCI <- apply(pac.norm.top, 2, quantile, probs = 0.025)
cols <- c("#A300FF", "#5BFF5B")[(mT > 1/NLinks) + 1]
barplot2(mT, main = "PACo residuals (unit branch length)", xlab = "Association", ylab = "Normalized residual",
         col = cols, border = "lightgrey", names.arg = colnamesPACo, las = 2, cex.names = 0.5,
         plot.ci = TRUE, ci.l = lCI, ci.u = uCI, ci.color = "blue")
abline(h = 1/NLinks, col = "red")
dev.off()

############################################
# 4. Clustering outlier terminals (optional)
############################################

# Median-centered residuals
pac.add <- pac.add / apply(pac.add, 1, sum) - 1/NLinks
pac.top <- pac.top / apply(pac.top, 1, sum) - 1/NLinks

# Incongruence metric
im.add <- apply(pac.add, 2, median)
out <- ifelse(im.add > 1/NLinks, 1, 0)
names(out) <- unique(unlist(strsplit(colnamesPACo, "_")))

# Visualize on trees
host_tree <- ladderize(read.tree("Bsylvaticum.tree"), right = TRUE)
parasite_tree <- ladderize(read.tree("Esylvatica.tree"), right = TRUE)

pdf("PACo_Outlier_Trees.pdf", width = 10, height = 6)
par(mfrow = c(1, 2), mar = c(1, 1, 4, 4))
plotTree(host_tree, offset = 0.5, fsize = 0.5)
title("B. sylvaticum – PACo outliers", font.main = 1, cex.main = 0.8)
tiplabels(pie = to.matrix(out, sort(unique(out))), piecol = c("#A300FF", "#5BFF5B"), cex = 0.4)

plotTree(parasite_tree, offset = 0.5, fsize = 0.5)
title("E. sylvatica – PACo outliers", font.main = 1, cex.main = 0.8)
tiplabels(pie = to.matrix(out, sort(unique(out))), piecol = c("#A300FF", "#5BFF5B"), cex = 0.4)
dev.off()

############################################
# 5. Run global PACo test with statistics
############################################

host_tree <- read.tree("Esylvatica.tree")
parasite_tree <- read.tree("Bsylvaticum.tree")
NCP <- as.matrix(read.table("Matrix.tsv", header = TRUE, sep = "\t", row.names = 1))

H <- cophenetic(host_tree)
P <- cophenetic(parasite_tree)
D <- prepare_paco_data(H, P, NCP)
D <- add_pcoord(D, correction = "cailliez")
D <- PACo(D, nperm = 1000, seed = 13, method = "quasiswap", symmetric = TRUE)

cat("PACoGlobal:", D$gof$ss, "p.global:", D$gof$p, "\n")

# Extract residuals
res <- residuals_paco(D$proc)

# Format for plotting
assoc <- data.frame(
  pol = rownames(NCP)[which(NCP == 1, arr.ind = TRUE)[, "row"]],
  pla = colnames(NCP)[which(NCP == 1, arr.ind = TRUE)[, "col"]]
)
weight <- (res^-2) / 50
cols <- c("blue", "red")[(as.vector(out) == 1) + 1]

# Save final cophylogenetic plot
data <- read.table("TableBE.tsv", sep = "\t", header = FALSE)
cophyloplot <- cophylo(host_tree, parasite_tree, assoc = data, use.edge.length = FALSE)

pdf("PACo_CoPhylo_Final.pdf", width = 10, height = 10)
plot(cophyloplot, link.type = "curved", link.lwd = 1.5, link.lty = "solid", link.col = cols, fsize = 0.5)
title(main = "B. sylvaticum vs E. sylvatica – PACo Analysis", font.main = 1)
legend("topright", legend = c("Congruent (p < 0.05)", "Outlier (p ≥ 0.05)"), lty = 1, col = c("blue", "red"), bty = "n")
dev.off()
