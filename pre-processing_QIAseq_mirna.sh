#!/usr/bin/env bash

# Set parameters
BASE_DIR="/path/to/base_dir"
SAMPLES_FILE="$BASE_DIR/samples.txt"
SPECIES="Mouse"
MIN_LENGTH=15
MAX_LENGTH=55

# Create output directories if they don't exist
mkdir -p "$BASE_DIR/trimmed" "$BASE_DIR/fastqc" "$BASE_DIR/multiqc"

# Read sample list
if [ ! -f "$SAMPLES_FILE" ]; then
    echo "Sample file not found: $SAMPLES_FILE"
    exit 1
fi

SAMPLES=$(cat "$SAMPLES_FILE")

process_sample() {
    local BASE_DIR=$1
    local SAMPLE=$2
    local SAMPLE_NAME=$(basename "$SAMPLE" .R1.fastq.gz)
    local SAMPLE_DIR="$BASE_DIR/$SAMPLE_NAME"
    local TRIMMED=$BASE_DIR/trimmed
    local FASTQC=$BASE_DIR/fastqc
    local MULTIQC=$BASE_DIR/multiqc
    local LOG_FILE="$BASE_DIR/${SAMPLE_NAME}_processing.log"
    local BBDUK="/home/msc19104442/.conda/envs/mirna_env/bin/bbduk.sh"
    local ADAPTERS_FILE="/home/msc19104442/.conda/envs/mirna_env/opt/bbmap-39.06-1/resources/adapters.fa"
    local CONTAMINANTS_FILE="/home/msc19104442/.conda/envs/mirna_env/opt/bbmap-39.06-1/resources/sequencing_artifacts.fa"
    local BBDUK_TRIMMED_FILE="$TRIMMED/${SAMPLE_NAME}_bbduk_trimmed.fastq.gz"
    local UMI_EXTRACTED_FILE="$TRIMMED/${SAMPLE_NAME}_umi_extracted.fastq.gz"
    local DEDUPED_FILE="$TRIMMED/${SAMPLE_NAME}_deduped_trimmed.fastq"

    echo "Processing sample: $SAMPLE_NAME" | tee -a "$LOG_FILE"
    echo "Sample directory: $SAMPLE_DIR" | tee -a "$LOG_FILE"
    echo "Adapters file: $ADAPTERS_FILE" | tee -a "$LOG_FILE"
    echo "Contaminants file: $CONTAMINANTS_FILE" | tee -a "$LOG_FILE"

    # Create necessary directories
    mkdir -p "$SAMPLE_DIR" "$TRIMMED" "$FASTQC"

    # Step 1: Use umi_tools to extract UMIs
    umi_tools extract --extract-method=regex \
                      --bc-pattern="(?P<smallRNA>[ATCG]{$MIN_LENGTH,$MAX_LENGTH})(?P<common>AACTGTAGGCACCATCAAT)(?P<umi_1>[ATCG]{12})(?P<junk>.*)" \
                      --log="$TRIMMED/${SAMPLE_NAME}_umi_tools.log" \
                      -I "$SAMPLE" -S "$UMI_EXTRACTED_FILE" 2>&1 | tee -a "$LOG_FILE"

    if [ $? -ne 0 ]; then
        echo "UMI extraction failed for sample: $SAMPLE_NAME" | tee -a "$LOG_FILE"
        exit 1
    fi

    # Step 2: Use BBduk to detect and remove adapters and contaminants
    $BBDUK in="$UMI_EXTRACTED_FILE" out="$BBDUK_TRIMMED_FILE" ref="$ADAPTERS_FILE,$CONTAMINANTS_FILE" literal=AACTGTAGGCACCATCAAT ktrim=r k=19 mink=9 hdist=1 tpe tbo qtrim=rl trimq=25 maxlength=100 forbidn=t 2> "$TRIMMED/${SAMPLE_NAME}_bbduk.log" | tee -a "$LOG_FILE"
    
    if [ $? -ne 0 ]; then
        echo "BBduk failed for sample: $SAMPLE_NAME" | tee -a "$LOG_FILE"
        exit 1
    fi

    # Step 3: Filter reads by length using cutadapt
    cutadapt -m "$MIN_LENGTH" -M "$MAX_LENGTH" -o "$DEDUPED_FILE" "$BBDUK_TRIMMED_FILE" 2> "$TRIMMED/${SAMPLE_NAME}_cutadapt.log" | tee -a "$LOG_FILE"
    
    if [ $? -ne 0 ]; then
        echo "Length filtering failed for sample: $SAMPLE_NAME" | tee -a "$LOG_FILE"
        exit 1
    fi

    # Step 4: Run FastQC on the final deduplicated file
    fastqc -o "$FASTQC" "$DEDUPED_FILE" 2>&1 | tee -a "$LOG_FILE"
    
    if [ $? -ne 0 ]; then
        echo "FastQC failed for sample: $SAMPLE_NAME" | tee -a "$LOG_FILE"
        exit 1
    fi

    # Step 5: Run MULTIQC on the final deduplicated file
    multiqc -o "$MULTIQC" "$BASE_DIR" 2>&1 | tee -a "$LOG_FILE"
    
    if [ $? -ne 0 ]; then
        echo "MultiQC failed for sample: $SAMPLE_NAME" | tee -a "$LOG_FILE"
        exit 1
    fi

    echo "Processing completed for sample: $SAMPLE_NAME" | tee -a "$LOG_FILE"
}

# Process each sample in the list
for SAMPLE in $SAMPLES; do
    process_sample "$BASE_DIR" "$SAMPLE"
done

echo "All samples processed successfully."

