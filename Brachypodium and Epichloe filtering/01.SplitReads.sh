#!/bin/bash

# Full workflow for read mapping and sequence extraction
# MIGUEL CAMPOS CACERES

# Step 0: Create organized directories for results
mkdir -p ./Mappings
mkdir -p ./Unmapped
mkdir -p ./Separated/Epichloe_Mapping
mkdir -p ./Separated/Brachypodium_Mapping
mkdir -p ./Epichloe/{Mappings_CDS,FastQ,Stats,Fastas,Genes,Genes_Filtered}
mkdir -p ./Brachypodium/{Mappings_CDS,FastQ,Stats,Fastas,Genes,Genes_Filtered}

# Step 1: Map all reads against the concatenated reference
echo "Mapping reads against the concatenated reference..."
for i in ./*.fastq.gz; do
    base=$(basename "${i}" .fastq.gz)
    minimap2 -t 100 -ax sr ../Complete.fa "${i}" | \
    samtools view -bS | \
    samtools sort -@ 100 -o ./Mappings/"${base}".sort.bam
    samtools index ./Mappings/"${base}".sort.bam
done

# Step 2: Separate reads mapped to Brachypodium chromosomes (Chr01–Chr09)
echo "Separating reads mapped to Brachypodium genome..."
for bam in ./Mappings/*.sort.bam; do
    base=$(basename "${bam}" .sort.bam)
    samtools view -b "${bam}" Chr01 Chr02 Chr03 Chr04 Chr05 Chr06 Chr07 Chr08 Chr09 \
        > ./Separated/Brachypodium_Mapping/"${base}".brachypodium.bam
done

# Step 3: Convert BAMs mapped to E. typhina to FASTQ
echo "Converting Epichloe BAMs to FASTQ..."
for bam in ./Separated/Epichloe_Mapping/*.etyphina.bam; do
    base=$(basename "${bam}" .bam)
    samtools fastq "${bam}" > ./Epichloe/FastQ/"${base}".fastq
done

# Step 4: Remap unmapped reads to E. sylvatica genome
echo "Remapping unmapped reads to E. sylvatica..."
for bam in ./Unmapped/*.unmapped.bam; do
    base=$(basename "${bam}" .unmapped.bam)
    samtools fastq "${bam}" > ./Unmapped/"${base}".fastq
    minimap2 -t 100 -ax sr ./Esylvatica.fa ./Unmapped/"${base}".fastq | \
    samtools view -bS -F 4 | \
    samtools sort -@ 100 -o ./Separated/Epichloe_Mapping/"${base}".esylvatica.sort.bam
    samtools fastq ./Separated/Epichloe_Mapping/"${base}".esylvatica.sort.bam > \
        ./Epichloe/FastQ/"${base}".esylvatica.fastq
done

# Step 5: Combine reads from E. typhina and E. sylvatica into unified Epichloe reads
echo "Concatenating Epichloe reads..."
for prefix in $(ls ./Epichloe/FastQ/*.etyphina.fastq | sed 's/\.etyphina\.fastq//' | sort | uniq); do
    cat "${prefix}.etyphina.fastq" "${prefix}.esylvatica.fastq" > "${prefix}.epichloe.fastq"
    echo "Combined: ${prefix}.epichloe.fastq"
done

rm ./Epichloe/FastQ/*etyphina.fastq
rm ./Epichloe/FastQ/*esylvatica.fastq
echo "Concatenation complete."

# Step 6: Map Epichloe reads to E. festucae CDS
echo "Mapping Epichloe reads to E. festucae CDS..."
for fastq in ./Epichloe/FastQ/*.epichloe.fastq; do
    base=$(basename "${fastq}" .epichloe.fastq)
    minimap2 -t 100 -ax sr ./Efestucae.cds.fa "${fastq}" | \
    samtools view -bS -F 4 | \
    samtools sort -@ 100 -o ./Epichloe/Mappings_CDS/"${base}".sort.bam
done

# Step 7: Extract coding sequences from mapped BAMs
echo "Extracting coding sequences from Epichloe BAMs..."
for bam in ./Epichloe/Mappings_CDS/*.sort.bam; do
    base=$(basename "${bam}" .sort.bam)
    /SOFT/samtools-1.16.1/samtools consensus -a -d 3 --threads 100 "${bam}" \
        -o ./Epichloe/Fastas/"${base}".fasta
done

# Step 8: Generate gene stats using Python scripts
echo "Calculating gene stats and missing data for Epichloe..."
python3 02.Counts.py
python3 03.Summary.py
python3 04.NameGenes.py
python3 05.Extract_Genes.py
python3 06.FilterGenes.py

# Step 9: Convert Brachypodium BAMs to FASTQ
echo "Converting Brachypodium BAMs to FASTQ..."
for bam in ./Separated/Brachypodium_Mapping/*.brachypodium.bam; do
    base=$(basename "${bam}" .bam)
    samtools fastq "${bam}" > ./Brachypodium/FastQ/"${base}".fastq
done

# Step 10: Map Brachypodium reads to their CDS
echo "Mapping Brachypodium reads to B. sylvaticum CDS..."
for fastq in ./Brachypodium/FastQ/*.fastq; do
    base=$(basename "${fastq}" .fastq)
    minimap2 -t 100 -ax sr ./Bsylvaticum.cds.fa "${fastq}" | \
    samtools view -bS -F 4 | \
    samtools sort -@ 100 -o ./Brachypodium/Mappings_CDS/"${base}".sort.bam
done

# Step 11: Extract coding sequences
echo "Extracting coding sequences from Brachypodium BAMs..."
for bam in ./Brachypodium/Mappings_CDS/*.sort.bam; do
    base=$(basename "${bam}" .sort.bam)
    /SOFT/samtools-1.16.1/samtools consensus -a -d 10 --threads 100 "${bam}" \
        -o ./Brachypodium/Fastas/"${base}".fasta
done

# Step 12: Generate gene stats using Python scripts
echo "Calculating gene stats and missing data for Brachypodium..."
python3 07.Counts2.py
python3 08.Summary2.py
python3 09.NameGenes2.py
python3 10.Extract_Genes2.py
python3 FilterGenes2.py

# Final message
echo "Workflow completed. Results are organized in the following directories:"
echo "1. Mappings/: All full BAM mappings"
echo "2. Separated/Epichloe_Mapping/: Reads mapped to Epichloe"
echo "3. Separated/Brachypodium_Mapping/: Reads mapped to Brachypodium"
echo "4. Epichloe/: CDS sequences and stats for Epichloe"
echo "5. Brachypodium/: CDS sequences and stats for Brachypodium"
