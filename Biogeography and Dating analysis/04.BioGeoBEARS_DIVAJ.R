############################################################
# BioGeoBEARS Analysis: DIVALIKE +J Model
# Author: Miguel Campos
# Description:
# This script runs the DIVALIKE+J model in BioGeoBEARS,
# generates ancestral area reconstructions, classifies events, and produces visualizations.
# This is adapted to Eastern Clade, just change the input to Western to do the other one.
############################################################

# Load required libraries
library(BioGeoBEARS)
library(ape)
library(readr)
library(stringr)
library(phytools)

# Set working directory
setwd("./BioGeoBEARS/Eastern_DIVAJ")

# Define input files
trfn <- "East_Dated.tree"
geogfn <- "Biogeobears_Geo_Eastern.txt"
max_range_size <- 3

############################################################
# Configure and run DIVALIKE +J model
############################################################

BioGeoBEARS_run_object <- define_BioGeoBEARS_run()
BioGeoBEARS_run_object$trfn <- trfn
BioGeoBEARS_run_object$geogfn <- geogfn
BioGeoBEARS_run_object$max_range_size <- max_range_size
BioGeoBEARS_run_object$include_null_range <- TRUE
BioGeoBEARS_run_object$use_optimx <- TRUE
BioGeoBEARS_run_object$num_cores_to_use <- 4
BioGeoBEARS_run_object$return_condlikes_table <- TRUE
BioGeoBEARS_run_object$calc_TTL_loglike_from_condlikes_table <- TRUE
BioGeoBEARS_run_object$calc_ancprobs <- TRUE

# Set fixed subset sympatry to 0
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s", "type"] <- "fixed"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s", "init"] <- 0.0
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s", "est"] <- 0.0

# Set vicariance parameters
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["ysv", "type"] <- "2-j"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["ys", "type"] <- "ysv*1/2"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["y", "type"] <- "ysv*1/2"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["v", "type"] <- "ysv*1/2"

# Fix mx01v
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01v", "type"] <- "fixed"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01v", "init"] <- 0.5
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01v", "est"] <- 0.5

# Enable jump dispersal (+J)
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j", "type"] <- "free"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j", "init"] <- 0.01
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j", "est"] <- 0.01

# Run and save model
check_BioGeoBEARS_run(BioGeoBEARS_run_object)
res_DIVALIKE <- bears_optim_run(BioGeoBEARS_run_object)
saveRDS(res_DIVALIKE, file = "BioGeoBEARS_result_DIVALIKE.rds")

############################################################
# Plot ancestral range reconstruction
############################################################

pdf("BioGeoBEARS_DIVABASIC_plot.pdf", width = 10, height = 12)
plot_BioGeoBEARS_results(res_DIVALIKE, "DIVALIKE+J model", addl_params = list("j"))
dev.off()

pdf("BioGeoBEARS_DIVALIKE_plot.pdf", width = 10, height = 12)
plot_BioGeoBEARS_results(
  results_object = res_DIVALIKE,
  analysis_titletxt = "Ancestral area reconstruction (DIVALIKE+J)",
  plotwhat = "pie",
  label.offset = 0.05,
  tipcex = 0.6,
  statecex = 0.4,
  splitcex = 0.3,
  titlecex = 1.2,
  plotsplits = TRUE,
  plotlegend = TRUE,
  legend_ncol = 3,
  legend_cex = 0.8,
  include_null_range = TRUE
)
dev.off()

############################################################
# Event classification: Dispersal, Vicariance, Extinction
############################################################

tree <- read.tree(trfn)
tipranges_object <- getranges_from_LagrangePHYLIP(geogfn)

areas <- colnames(tipranges_object@df)
nareas <- length(areas)
combos_binarias <- as.matrix(expand.grid(rep(list(c(0, 1)), nareas)))
combos_binarias <- combos_binarias[rowSums(combos_binarias) > 0, ]
states_list <- apply(combos_binarias, 1, function(x) areas[x == 1])

edge <- tree$edge
parents <- edge[, 1]
children <- edge[, 2]
ML_states <- apply(res_DIVALIKE$ML_marginal_prob_each_state_at_branch_top_AT_node, 1, which.max)

state_to_areas <- function(state_index) {
  if (state_index > length(states_list)) return(NA)
  paste(states_list[[state_index]], collapse = "+")
}

events <- data.frame(
  parent = parents,
  child = children,
  parent_state = ML_states[parents],
  child_state = ML_states[children]
)
events$parent_state_areas <- sapply(events$parent_state, state_to_areas)
events$child_state_areas <- sapply(events$child_state, state_to_areas)

split_areas <- function(area_string) {
  if (is.na(area_string) || area_string == "") return(character(0))
  unlist(strsplit(area_string, split = "\\+"))
}

classify_event <- function(parent_areas, child_areas) {
  parent_set <- split_areas(parent_areas)
  child_set <- split_areas(child_areas)

  if (length(parent_set) == 0 || length(child_set) == 0) return("Unknown")
  if (all(child_set %in% parent_set) && length(child_set) < length(parent_set)) return("Extinction (loss of areas)")
  if (all(parent_set %in% child_set) && length(child_set) > length(parent_set)) return("Dispersal (gain of areas)")
  if (!all(child_set %in% parent_set) && !all(parent_set %in% child_set)) return("Vicariance (area split)")
  if (setequal(parent_set, child_set)) return("No change")
  return("Other")
}

events$event_type <- mapply(classify_event, events$parent_state_areas, events$child_state_areas)

write.table(events, file = "BioGeoBEARS_DIVALIKE_events.tsv", sep = "\t", row.names = FALSE, quote = FALSE)

############################################################
# Visualize events on the phylogeny
############################################################

event_colors <- c(
  "Dispersal (gain of areas)" = "blue",
  "Extinction (loss of areas)" = "red",
  "Vicariance (area split)"    = "green",
  "No change"                  = "grey",
  "Other"                      = "black",
  "Unknown"                    = "black"
)
events$color <- event_colors[events$event_type]

pdf("BioGeoBEARS_DIVALIKE_tree_events.pdf", width = 10, height = 12)
plot(tree, cex = 0.6, no.margin = TRUE)
title("Tree with classified biogeographical events", cex.main = 1.5)
for (i in 1:nrow(events)) {
  node_id <- events$parent[i]
  point_color <- events$color[i]
  nodelabels(node = node_id, pch = 21, bg = point_color, cex = 1.2)
}
legend("topright", legend = names(event_colors), fill = event_colors, cex = 0.8, bty = "n")
dev.off()
