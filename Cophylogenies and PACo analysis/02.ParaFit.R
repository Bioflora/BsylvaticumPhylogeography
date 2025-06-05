################################################################
# ParaFit Analysis – 1000 Replicates (Parallelized)
# Author: Miguel Campos Cáceres
# Description:
# This script performs a co-phylogenetic congruence test using ParaFit
# with 1000 replicates and parallel processing. It summarizes global and
# per-link statistics, and visualizes results on a cophylogenetic plot.
################################################################

# Load required libraries
library(ape)
library(ggplot2)
library(phytools)
library(parallel)

# Set working directory
setwd("./ParaFit/")

# 1. Load trees and association matrix
host_tree <- ladderize(read.tree("Bsylvaticum.tree"), right = TRUE)
parasite_tree <- ladderize(read.tree("Esylvatica.tree"), right = TRUE)
assoc_matrix <- as.matrix(read.table("Matrix.tsv", sep = "\t", header = TRUE, row.names = 1, check.names = FALSE))

# 2. Compute patristic distance matrices
host_dist <- cophenetic(host_tree)[rownames(assoc_matrix), rownames(assoc_matrix)]
parasite_dist <- cophenetic(parasite_tree)[colnames(assoc_matrix), colnames(assoc_matrix)]

# 3. Run ParaFit in parallel
n_reps <- 1000
n_cores <- detectCores() - 1

cat("Running", n_reps, "ParaFit replicates using", n_cores, "cores...\n")

run_parafit <- function(i) {
  cat("Run", i, "started\n")
  result <- parafit(
    host_dist,
    parasite_dist,
    assoc_matrix,
    nperm = 999,
    test.links = TRUE,
    silent = TRUE,
    correction = "cailliez"
  )
  list(
    global_stat = result$ParaFitGlobal,
    global_pval = result$p.global,
    link_table = as.matrix(result$link.table[, 1:5])
  )
}

results <- mclapply(1:n_reps, run_parafit, mc.cores = n_cores)

# 4. Summarize results
global_stats <- sapply(results, function(x) x$global_stat)
global_pvals <- sapply(results, function(x) x$global_pval)

link_sum <- Reduce("+", lapply(results, function(x) x$link_table))
link_avg <- link_sum / n_reps
colnames(link_avg) <- colnames(results[[1]]$link_table)
rownames(link_avg) <- rownames(results[[1]]$link_table)
link_avg_df <- as.data.frame(link_avg)

write.table(link_avg_df, file = "ParaFit_Links_Avg.tsv", sep = "\t", quote = FALSE)

cat("\n--- GLOBAL RESULTS (averaged across replicates) ---\n")
cat("Mean ParaFitGlobal statistic:", mean(global_stats), "\n")
cat("Mean global p-value:", mean(global_pvals), "\n")

# 5. Visualization (optional)
parafit_result <- read.table("ParaFit_Links_Avg.tsv", sep = '\t', header = TRUE)
link_colors <- ifelse(parafit_result$p.F1 <= 0.05, "blue", "red")

assoc <- read.table("TableBE.tsv", sep = '\t', header = FALSE)
BrachyEpi.cophylo <- cophylo(host_tree, parasite_tree, assoc = assoc, use.edge.length = FALSE)

pdf("ParaFit_CoPhylo.pdf", width = 10, height = 10)
plot(
  BrachyEpi.cophylo,
  link.type = "curved",
  link.lwd = 1.5,
  link.lty = "solid",
  link.col = link_colors,
  fsize = 0.5
)
title(main = "B. sylvaticum vs E. sylvatica (average of 1000 runs)", font.main = 1)
legend("topright", c("Congruent (p < 0.05)", "Incongruent (p ≥ 0.05)"), lty = 1, col = c("blue", "red"), bty = "n")
dev.off()
