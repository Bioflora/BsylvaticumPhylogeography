#!/bin/bash

##############################################################################
# SNP Calling Pipeline for Brachypodium Mapping Data                         #
# Author: Miguel Campos Cáceres                                              #
# Description:                                                               #
# This script performs variant calling on Brachypodium BAM files mapped      #
# against a reference genome. It includes indexing, per-sample variant       #
# calling with bcftools, VCF merging, and filtering for biallelic SNPs       #
# with no missing data.                                                      #
##############################################################################

# Step 0: Create output directory
mkdir -p ./Brachypodium/VCF/

# Step 1: Index BAM files
echo "Indexing Brachypodium BAM mappings..."
for bam in ./*.brachypodium.bam; do
    samtools index "${bam}"
done

# Step 2: Convert BAM to VCF using bcftools
echo "Generating VCF files from BAMs..."
for bam in ./*.brachypodium.sort.bam; do
    base_name=$(basename "${bam}" .brachypodium.sort.bam)
    bcftools mpileup --threads 100 -f ../Complete.fa -a DP "${bam}" | \
    bcftools call --threads 100 -O v -m | \
    bcftools filter --threads 100 -e 'INFO/DP<10' -O z -o ./Brachypodium/VCF/"${base_name}".vcf.gz
    bcftools index ./Brachypodium/VCF/"${base_name}".vcf.gz
done

# Step 3: Merge all individual VCFs into one
echo "Merging all VCF files into a single file..."
bcftools merge --threads 100 -i - -m both -O z -o ./Brachypodium/VCF/Brachypodium.vcf.gz ./Brachypodium/VCF/*.vcf.gz

# Step 4: Filter multiallelic positions and missing data
echo "Filtering merged VCF to retain only biallelic and complete data..."
vcftools --gzvcf ./Brachypodium/VCF/Brachypodium.vcf.gz --max-alleles 2 --max-missing 1 --recode --stdout | \
bgzip > ./Brachypodium/VCF/Brachypodium_filtered.vcf.gz

echo "Filtered Brachypodium VCF file successfully generated."
