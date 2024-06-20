#!/usr/bin/env bash

# Directory to store the genome
mkdir -p "$BASE_DIR/hg38"
GENOME_DIR="/data4/msc19104442/liam_alignment/hg38"
cd "$GENOME_DIR"

UNMAPPED_FILE="${SAMPLE_DIR}/${SAMPLE_NAME}_unmapped.fq"
HG38_MAPPED_FILE="${HG38_DIR}/${SAMPLE_NAME}_hg38_mapped.sam"

# Download the reference genome
wget ftp://ftp.ensembl.org/pub/release-100/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz

# Unzip the downloaded file
gunzip Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz

# Rename the file for convenience
mv Homo_sapiens.GRCh38.dna.primary_assembly.fa hg38.fa

# Index the reference genome with BWA
bwa index hg38.fa

echo "Download and preparation of the hg38 reference genome is complete."

# Map unmapped reads to human genome (hg38)
bowtie -x hg38 -q "$UNMAPPED_FILE" -k 50 -v 2 --best --strata --sam --trim5 2 --trim3 2 > "$HG38_MAPPED_FILE" 2> "${HG38_DIR}/${SAMPLE_NAME}_hg38.err"

echo "Unmapped files mapped to human genome."

