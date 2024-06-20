#!/usr/bin/env bash

# Set parameters
BASE_DIR="/data4/msc19104442/liam_alignment"
SAMPLES_FILE="$BASE_DIR/samples.txt"
HAIRPIN_FILE="$BASE_DIR/hairpin_miRBase.fa"
SPECIES="mmu"  # Adjusted for mouse-specific miRNA
COMMON_SEQUENCE="AACTGTAGGCACCATCAAT"

# Download miraligner from seqbuster
wget -O seqbuster-miraligner-3.5.zip "https://github.com/lpantano/seqbuster/archive/refs/tags/miraligner-3.5.zip"
unzip seqbuster-miraligner-3.5.zip -d "$BASE_DIR"

# Path to miraligner.jar
MIRALIGNERJAR="$BASE_DIR/seqbuster-miraligner-3.5/miraligner/miraligner-3.2/miraligner.jar"

# Ensure the miraligner JAR file exists
if [ ! -f "$MIRALIGNERJAR" ]; then
    echo "Error: miraligner.jar file not found at $MIRALIGNERJAR"
    exit 1
fi

# Create miraligner output directory
mkdir -p "$BASE_DIR/miraligner"

MIRALIGNER_DIR="$BASE_DIR/miraligner"

# Ensure MIRTOP directory exists
MIRTOP="$BASE_DIR/mirtop"
if [ ! -d "$MIRTOP" ]; then
    mkdir -p "$MIRTOP"
fi

# Download GFF file if not present
if [ ! -f "$MIRTOP/mmu.gff" ]; then
    wget --no-check-certificate https://www.mirbase.org/download/mmu.gff3 -O "$MIRTOP/mmu.gff"
fi

# Read sample list
if [ ! -f "$SAMPLES_FILE" ]; then
    echo "Sample file not found: $SAMPLES_FILE"
    exit 1
fi

SAMPLE_DIR=$SAMPLE_DIR/$SAMPLES

SAMPLES=$(cat "$SAMPLES_FILE")

# Run miRAligner for each sample
for SAMPLE in $SAMPLES; do
    SAMPLE_NAME=$(basename "$SAMPLE" .R1.fastq.gz)
    TRIMMED_FILE="$BASE_DIR/$SAMPLE_NAME/collapsed/${SAMPLE_NAME}_deduped_trimmed.fastq"
    if [ -f "$TRIMMED_FILE" ]; then
        echo "Processing sample: $SAMPLE_NAME"
        echo "Running miRAligner on $TRIMMED_FILE"
        
        java -jar "$MIRALIGNERJAR" -sub 1 -trim 3 -add 3 -s mmu -i "$TRIMMED_FILE" -db "$BASE_DIR/seqbuster-miraligner-3.5/miraligner/DB" -o "$BASE_DIR/miraligner/$SAMPLE_NAME" > "$BASE_DIR/miraligner/${SAMPLE_NAME}_log.txt" 2>&1

        if [ $? -ne 0 ]; then
            echo "Error: miRAligner failed for sample $SAMPLE_NAME. Check log file for details: $BASE_DIR/miraligner/${SAMPLE_NAME}_log.txt"
            continue
        fi
    else
        echo "Trimmed file not found for sample: $SAMPLE_NAME"
    fi
done

# Convert miRAligner output to GFF3-srna format
MIRNA_FILES=($MIRALIGNER_DIR/*.mirna)

if [ ${#MIRNA_FILES[@]} -gt 0 ]; then
    mirtop gff --format seqbuster --sps mmu --hairpin "$HAIRPIN_FILE" --gtf "$MIRTOP/mmu.gff" -o "$MIRALIGNER_DIR/miraligner_gff3" "${MIRNA_FILES[@]}"
else
    echo "No .mirna files found in $MIRALIGNER_DIR"
    exit 1
fi

# Quantify the *merged.bam files for each sample in their respective directories
for SAMPLE in $SAMPLES; do
    SAMPLE_NAME=$(basename "$SAMPLE" .R1.fastq.gz)
    SAMPLE_DIR="$BASE_DIR/$SAMPLE_NAME"
    MERGED_BAM_FILE="$SAMPLE_DIR/${SAMPLE_NAME}.merged.bam"
    
    if [ -f "$MERGED_BAM_FILE" ]; then
        echo "Processing merged BAM file for sample: $SAMPLE_NAME"
        mirtop gff --hairpin "$HAIRPIN_FILE" --gtf "$MIRTOP/mmu.gff" -o "$MIRTOP" --sps mmu "$MERGED_BAM_FILE"
    else
        echo "Merged BAM file not found for sample: $SAMPLE_NAME"
    fi
done

# Additional mirtop steps
mirtop counts --hairpin "$HAIRPIN_FILE" --gtf "$MIRTOP/mmu.gff" -o "$MIRTOP" --sps mmu --add-extra --gff "$MIRTOP/mirtop.gff"
mirtop export --format isomir --hairpin "$HAIRPIN_FILE" --gtf "$MIRTOP/mmu.gff" --sps mmu -o "$MIRTOP" "$MIRTOP/mirtop.gff"
mirtop stats "$MIRTOP/mirtop.gff" --out "$MIRTOP/stats"
Rscript collapse_mirtop.r "$MIRTOP/mirtop.tsv"

# Run MultiQC to aggregate FastQC reports
multiqc "$BASE_DIR" "$SAMPLE_DIR" -o "$BASE_DIR/multiqc"

# Notify user of completion
echo "Processing complete."

