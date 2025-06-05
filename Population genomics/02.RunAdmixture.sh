###############################################
# Script for Admixture Cross-Validation Analysis
# Author: Miguel Campos
# Description: Performs multiple ADMIXTURE runs across a range of K values,
#              storing output and extracting cross-validation (CV) errors.
###############################################

for run in {1..10}  # Number of independent ADMIXTURE runs
do
    mkdir "Run$run"
    cd "Run$run"

    for K in {1..10}  # Range of K (number of ancestral populations)
    do
        admixture -j50 -s $RANDOM --cv=10 ./Bsylvaticum_Filtered.bed $K > "log${K}.out"
    done

    # Extract CV error values and format output
    grep "CV" *out | awk '{print $3,$4}' | sed -e 's/(//;s/)//;s/://;s/K=//' > cv.error

    cd ..
done
