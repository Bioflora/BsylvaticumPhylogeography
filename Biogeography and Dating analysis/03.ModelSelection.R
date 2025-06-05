############################################################
# BioGeoBEARS Model Selection and Comparative Analysis
# Author: Miguel Campos (adapted)
# Date: June 2025
# Description:
# This script runs multiple biogeographic models using BioGeoBEARS,
# compares them using AIC, and saves parameters and plots for evaluation.
# This script is made for Eastern clade, just repeat changing the input for Western
############################################################

# Load required libraries
library(BioGeoBEARS)
library(ape)
library(ggplot2)
library(reshape2)
library(RColorBrewer)
library(stringr)
library(phytools)

# Set working directory
setwd("./BioGeoBEARS/Eastern_ModelTest")

# Function to run and save a BioGeoBEARS model
run_model <- function(model_name, trfn, geogfn, max_range_size=9, include_null_range=TRUE, j_param=FALSE) {
  run_object <- define_BioGeoBEARS_run()
  run_object$trfn <- trfn
  run_object$geogfn <- geogfn
  run_object$max_range_size <- max_range_size
  run_object$include_null_range <- include_null_range
  run_object$use_optimx <- TRUE
  run_object$num_cores_to_use <- 100
  run_object$return_condlikes_table = TRUE
  run_object$calc_TTL_loglike_from_condlikes_table = TRUE
  run_object$calc_ancprobs = TRUE

  # Customize models
  if (model_name == "DIVALIKE") {
    run_object$BioGeoBEARS_model_object@params_table["s","type"] <- "fixed"
    run_object$BioGeoBEARS_model_object@params_table["s","init"] <- 0.0
    run_object$BioGeoBEARS_model_object@params_table["s","est"] <- 0.0
    run_object$BioGeoBEARS_model_object@params_table["ysv","type"] <- "2-j"
    run_object$BioGeoBEARS_model_object@params_table["ys","type"] <- "ysv*1/2"
    run_object$BioGeoBEARS_model_object@params_table["y","type"] <- "ysv*1/2"
    run_object$BioGeoBEARS_model_object@params_table["v","type"] <- "ysv*1/2"
    run_object$BioGeoBEARS_model_object@params_table["mx01v","type"] <- "fixed"
    run_object$BioGeoBEARS_model_object@params_table["mx01v","init"] <- 0.5
    run_object$BioGeoBEARS_model_object@params_table["mx01v","est"] <- 0.5
  }

  if (model_name == "BAYAREALIKE") {
    run_object$BioGeoBEARS_model_object@params_table["s","type"] <- "fixed"
    run_object$BioGeoBEARS_model_object@params_table["s","init"] <- 0.0
    run_object$BioGeoBEARS_model_object@params_table["v","type"] <- "fixed"
    run_object$BioGeoBEARS_model_object@params_table["v","init"] <- 0.0
    run_object$BioGeoBEARS_model_object@params_table["ysv","type"] <- "1-j"
    run_object$BioGeoBEARS_model_object@params_table["ys","type"] <- "ysv*1/1"
    run_object$BioGeoBEARS_model_object@params_table["y","type"] <- "1-j"
    run_object$BioGeoBEARS_model_object@params_table["mx01y","type"] <- "fixed"
    run_object$BioGeoBEARS_model_object@params_table["mx01y","init"] <- 0.9999
  }

  # Add +J parameter if required
  if (j_param) {
    run_object$BioGeoBEARS_model_object@params_table["j","type"] <- "free"
    run_object$BioGeoBEARS_model_object@params_table["j","init"] <- 0.0001
    run_object$BioGeoBEARS_model_object@params_table["j","est"] <- 0.0001
  }

  check_BioGeoBEARS_run(run_object, allow_huge_ranges = TRUE)
  res <- bears_optim_run(run_object)
  saveRDS(res, file = paste0("BioGeoBEARS_result_", model_name, if (j_param) "+J" else "", ".rds"))
  return(res)
}

# Define input files
trfn <- "Eastern_Dated.tree"
geogfn <- "Biogeobears_Geo_Eastern.txt"

# Run model set
resDEC         <- run_model("DEC", trfn, geogfn)
resDECJ        <- run_model("DEC", trfn, geogfn, j_param=TRUE)
resDIVALIKE    <- run_model("DIVALIKE", trfn, geogfn)
resDIVALIKEJ   <- run_model("DIVALIKE", trfn, geogfn, j_param=TRUE)
resBAYAREALIKE <- run_model("BAYAREALIKE", trfn, geogfn)
resBAYAREALIKEJ<- run_model("BAYAREALIKE", trfn, geogfn, j_param=TRUE)

# Collect results
model_results <- list(
  DEC = resDEC, DECJ = resDECJ,
  DIVALIKE = resDIVALIKE, DIVALIKEJ = resDIVALIKEJ,
  BAYAREALIKE = resBAYAREALIKE, BAYAREALIKEJ = resBAYAREALIKEJ
)

# Function to compare +J vs base model
get_stats <- function(name1, name2) {
  LnL1 <- get_LnL_from_BioGeoBEARS_results_object(model_results[[name1]])
  LnL2 <- get_LnL_from_BioGeoBEARS_results_object(model_results[[name2]])
  numparams1 <- nrow(model_results[[name1]]$output@params_table)
  numparams2 <- nrow(model_results[[name2]]$output@params_table)
  stats <- AICstats_2models(LnL1, LnL2, numparams1, numparams2)
  return(stats)
}

# Pairwise model comparisons
stats_DEC       <- get_stats("DECJ", "DEC")
stats_DIVALIKE  <- get_stats("DIVALIKEJ", "DIVALIKE")
stats_BAYAREA   <- get_stats("BAYAREALIKEJ", "BAYAREALIKE")

# Save comparisons
write.table(stats_DEC, "ModelComparison_DEC_vs_DECJ.tsv", sep="\t", quote=FALSE)
write.table(stats_DIVALIKE, "ModelComparison_DIVALIKE_vs_DIVALIKEJ.tsv", sep="\t", quote=FALSE)
write.table(stats_BAYAREA, "ModelComparison_BAYAREALIKE_vs_BAYAREALIKEJ.tsv", sep="\t", quote=FALSE)

# Extract and save all model parameters
extract_params <- function(name) {
  extract_params_from_BioGeoBEARS_results_object(
    model_results[[name]],
    returnwhat = "table",
    addl_params = c("j"),
    paramsstr_digits = 4
  )
}
restable <- do.call(rbind, lapply(names(model_results), extract_params))
row.names(restable) <- names(model_results)
write.table(restable, file="ModelParameters_All.tsv", sep="\t", quote=FALSE)

# Add AIC to summary table and plot
restable$LnL <- as.numeric(restable$LnL)
restable$numparams <- as.numeric(restable$numparams)
restable$AIC <- 2 * restable$numparams - 2 * restable$LnL

pdf("BioGeoBEARS_AIC_comparison.pdf", width = 8, height = 6)
ggplot(restable, aes(x = rownames(restable), y = AIC, fill = rownames(restable))) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  labs(title = "Model Comparison by AIC",
       x = "Model", y = "AIC") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none")
dev.off()
