#!/bin/bash

# Check if a gene name was provided
if [ -z "$1" ]; then
  echo "Usage: $0 <gene_name>"
  exit 1
fi

GENE=$1
GENOME_DIR="Ecoli_genome/GCA_000005845.2"
GFF_FILE="$GENOME_DIR/genomic.gff"
FASTA_FILE="$GENOME_DIR/GCA_000005845.2_ASM584v2_genomic.fna"

# Step 1: Extract gene annotation from GFF
grep -w "$GENE" "$GFF_FILE" | grep 'cds' > "${GENE}.gff"

# Step 2: Create BED file with 250 bp flanking regions
awk -v OFS='\t' '{
    start=$4-251; if (start<0) start=0;
    end=$5+250;
    print $1, start, end, $9
}' "${GENE}.gff" > "${GENE}_flanked.bed"

# Step 3: Extract flanked sequence using bedtools (assumed installed)
bedtools getfasta -fi "$FASTA_FILE" -bed "${GENE}_flanked.bed" -name > "${GENE}_flanked.fasta"
