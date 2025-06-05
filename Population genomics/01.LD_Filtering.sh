###############################################
# SNP Filtering and LD Pruning
# Author: Miguel Campos
# Description:
# - Filters multiallelic sites and ambiguous bases
# - Converts to PLINK binary format
# - Performs LD pruning
# - Converts final dataset back to VCF
###############################################

# Step 1: Keep only biallelic SNPs and remove ambiguous nucleotides
bcftools view -i 'TYPE="snp" & (ALT="A" | ALT="C" | ALT="G" | ALT="T" | ALT="-")' Bsylvaticum.vcf -o Bsylvaticum_step1.vcf

# Step 2: Ensure only biallelic SNPs remain
vcftools --vcf Bsylvaticum_step1.vcf --max-alleles 2 --recode --out Bsylvaticum_step2

# Step 3: Convert VCF to PLINK2 binary format
plink2 --vcf Bsylvaticum_step2.recode.vcf --make-bed --out Bsylvaticum_plink --allow-extra-chr --set-missing-var-ids @:#

# Step 4: Remove SNPs with missing data and non-SNP variants
plink2 --bfile Bsylvaticum_plink --geno 0.0 --snps-only --make-bed --out Bsylvaticum_clean --allow-extra-chr

# Step 5: Perform LD pruning (500 bp window, step 5, r² < 0.2)
plink2 --bfile Bsylvaticum_clean --indep-pairwise 500 5 0.2 --out Bsylvaticum_pruned --allow-extra-chr

# Step 6: Extract pruned SNPs
plink2 --bfile Bsylvaticum_clean --extract Bsylvaticum_pruned.prune.in --make-bed --out Bsylvaticum_pruned_final --allow-extra-chr

# Step 7: Export final pruned dataset back to VCF format
plink2 --bfile Bsylvaticum_pruned_final --recode vcf --out Bsylvaticum_Filtered --allow-extra-chr
