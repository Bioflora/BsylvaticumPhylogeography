#!/bin/bash

############################################################
# TreePL Bootstrapping Pipeline
# Author: Miguel Campos
# Date: June 2025
# Description:
# This script runs TreePL on a set of 1000 bootstrap trees.
# For each tree, it creates a custom configuration file,
# replaces the tree and output paths, and runs TreePL via Docker.
############################################################

# Create output directory
mkdir -p bootstraps_results

# Initialize counter
i=0

# Loop through each bootstrap tree
for treefile in bootstraps_rooted/tree_*.tre; do
    echo "Processing $treefile..."

    # Get tree base name and define config/output filenames
    tree_base=$(basename "$treefile" .tre)
    conf_copy="conf_${i}.txt"
    dated_out="bootstraps_results/${tree_base}_dated.tree"

    # Copy the base TreePL config file
    cp conf.txt "$conf_copy"

    # Update treefile path in config
    sed -i "s|^treefile = .*|treefile = /TreePL/${treefile}|" "$conf_copy"

    # Update output file path in config
    sed -i "s|^outfile = .*|outfile = /TreePL/${dated_out}|" "$conf_copy"

    # Run TreePL via Docker
    docker run --rm -v "$PWD":/TreePL naturalis/docker-treepl "/TreePL/${conf_copy}"

    # Remove temporary config file
    rm "$conf_copy"

    ((i++))
done
