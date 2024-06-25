#!/usr/bin/env bash

# Set parameters
BASE_DIR="/path/to/base_dir"
SAMPLES_FILE="$BASE_DIR/samples.txt"
HAIRPIN_FILE="$BASE_DIR/hairpin_miRBase.fa"
SPECIES="Mouse"
MIN_LENGTH=15
MAX_LENGTH=55

# Process hairpin file
echo "Processing hairpin file..."
cp "$HAIRPIN_FILE" hairpin.fa
sed 's#^[^>]#N#g' hairpin.fa > hairpin_parse.fa
sed -i 's#\s.*##' hairpin_parse.fa
seqkit grep -r --pattern ".*mmu-.*" hairpin_parse.fa > hairpin_mmu.fa
seqkit seq --rna2dna hairpin_mmu.fa > tmp.fa
fasta_formatter -w 0 -i tmp.fa -o hairpin_processed.fa
rm hairpin.fa hairpin_mmu.fa hairpin_parse.fa tmp.fa
mv hairpin_processed.fa "$HAIRPIN_FILE"

# Index miRBase hairpin and mature files
echo "Building Bowtie index for hairpin file..."
bowtie-build "$HAIRPIN_FILE" hairpin

# Create output directories if they don't exist
mkdir -p "$BASE_DIR/mirtop" "$BASE_DIR/mirdeep2" "$BASE_DIR/trimmed" "$BASE_DIR/fastqc" "$BASE_DIR/multiqc"

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
    local MIRTOP=$BASE_DIR/mirtop
    local TRIMMED=$BASE_DIR/trimmed
    local FASTQC=$BASE_DIR/fastqc
    local MIRDEEP2=$BASE_DIR/mirdeep2
    local MULTIQC=$BASE_DIR/multiqc
    local BBDUK="/home/msc19104442/.conda/envs/mirna_env/bin/bbduk.sh"
    local ADAPTERS_FILE="/home/msc19104442/.conda/envs/mirna_env/opt/bbmap-39.06-1/resources/adapters.fa"
    local CONTAMINANTS_FILE="/home/msc19104442/.conda/envs/mirna_env/opt/bbmap-39.06-1/resources/sequencing_artifacts.fa"
    local BBDUK_TRIMMED_FILE="$TRIMMED/${SAMPLE_NAME}_bbduk_trimmed.fastq.gz"
    local UMI_EXTRACTED_FILE="$TRIMMED/${SAMPLE_NAME}_umi_extracted.fastq.gz"
    local COLLAPSED_DIR="$SAMPLE_DIR/collapsed"
    local COLLAPSED_FILE="$COLLAPSED_DIR/${SAMPLE_NAME}_deduped_trimmed.fastq"

    echo "Processing sample: $SAMPLE_NAME"
    echo "Sample directory: $SAMPLE_DIR"
    echo "Adapters file: $ADAPTERS_FILE"
    echo "Contaminants file: $CONTAMINANTS_FILE"

    # Create necessary directories
    mkdir -p "$SAMPLE_DIR" "$COLLAPSED_DIR" "$TRIMMED" "$FASTQC"

    # Step 1: Use umi_tools to extract UMIs
    umi_tools extract --extract-method=regex \
                      --bc-pattern="(?P<smallRNA>[ATCG]{15,55})(?P<common>AACTGTAGGCACCATCAAT)(?P<umi_1>[ATCG]{12})(?P<junk>.*)" \
                      --log="$TRIMMED/${SAMPLE_NAME}_umi_tools.log" \
                      -I "$SAMPLE" -S "$UMI_EXTRACTED_FILE"

    if [ $? -ne 0 ]; then
        echo "UMI extraction failed for sample: $SAMPLE_NAME"
        exit 1
    fi

    # Step 2: Use BBduk to detect and remove adapters and contaminants
    if [ ! -f "$BBDUK_TRIMMED_FILE" ]; then
        $BBDUK in="$UMI_EXTRACTED_FILE" out="$BBDUK_TRIMMED_FILE" ref="$ADAPTERS_FILE,$CONTAMINANTS_FILE" literal=AACTGTAGGCACCATCAAT ktrim=r k=19 mink=9 hdist=1 tpe tbo qtrim=rl trimq=25 maxlength=100 forbidn=t 2> "$TRIMMED/${SAMPLE_NAME}_bbduk.log"
        if [ $? -ne 0 ]; then
            echo "BBduk failed for sample: $SAMPLE_NAME"
            exit 1
        fi
    fi

    # Step 3: Filter reads by length using cutadapt
    cutadapt -m 15 -M 55 -o "$DEDUPED_FILE" "$BBDUK_TRIMMED_FILE" 2> "$TRIMMED/${SAMPLE_NAME}_cutadapt.log"
    if [ $? -ne 0 ]; then
        echo "Length filtering failed for sample: $SAMPLE_NAME"
        exit 1
    fi

    # Step 4: Run FastQC on the final deduplicated file
    fastqc -o "$FASTQC" "$DEDUPED_FILE"
    if [ $? -ne 0 ]; then
        echo "FastQC failed for sample: $SAMPLE_NAME"
        exit 1
    fi

    # Collapse sequences
    if [ ! -f "$COLLAPSED_FILE" ]; then
        seqcluster collapse -f "$DEDUPED_FILE" -m 1 --min_size 15 -o "$COLLAPSED_DIR" > "$TRIMMED/${SAMPLE_NAME}_seqcluster.log" 2>&1
        if [ $? -ne 0 ]; then
            echo "Seqcluster collapse failed for sample: $SAMPLE"
            exit 1
        fi
        COLLAPSED_FILE=$(ls "$COLLAPSED_DIR/${SAMPLE_NAME}"*_deduped_trimmed.fastq | head -n 1)
        if [ ! -f "$COLLAPSED_FILE" ]; then
            echo "ERROR: Collapsed file not created: $COLLAPSED_FILE"
            exit 1
        fi
    else
        echo "Collapsed file already exists: $COLLAPSED_FILE. Skipping collapse step."
    fi

    (bowtie -x hairpin -q "$COLLAPSED_FILE" -k 50 -e 99999 --best --strata --sam > "${SAMPLE_DIR}/${SAMPLE_NAME}_run1.sam") 2> "${SAMPLE_DIR}/${SAMPLE_NAME}_run1.err"

    local processed=$(grep -o 'reads processed: [0-9]*' "${SAMPLE_DIR}/${SAMPLE_NAME}_run1.err" | awk '{print $3}')
    local alignment=$(grep -o 'reads with at least one alignment: [0-9]*' "${SAMPLE_DIR}/${SAMPLE_NAME}_run1.err" | awk '{print $7}')
    local failed=$(grep -o 'reads that failed to align: [0-9]*' "${SAMPLE_DIR}/${SAMPLE_NAME}_run1.err" | awk '{print $6}')
    local reported=$(grep -o 'Reported [0-9]* alignments' "${SAMPLE_DIR}/${SAMPLE_NAME}_run1.err" | awk '{print $2}')

    rm "${SAMPLE_DIR}/${SAMPLE_NAME}_run1.err"
    printf "%s %s %s %s %s %s\n" "$SAMPLE_NAME" "1" "$processed" "$alignment" "$failed" "$reported" > "${SAMPLE_DIR}/${SAMPLE_NAME}_run1_summary.err"

    for ((i=1; i<=20; i++)); do
        local SAM_FILE="${SAMPLE_DIR}/${SAMPLE_NAME}_run${i}.sam"
        local BAM_FILE="${SAMPLE_DIR}/${SAMPLE_NAME}_run${i}.bam"
        local SRT_BAM_FILE="${SAMPLE_DIR}/${SAMPLE_NAME}_run${i}_srt.bam"
        local UNMAPPED_FILE="${SAMPLE_DIR}/${SAMPLE_NAME}_run${i}_unmapped.fq"
        
        if [ ! -f "$SAM_FILE" ]; then
            echo "SAM file not found: $SAM_FILE. Skipping subsequent steps."
            break
        fi
        
        samtools view -F 4 -h -b "$SAM_FILE" > "$BAM_FILE"
        samtools sort -o "$SRT_BAM_FILE" "$BAM_FILE"
        samtools view -f 4 "$SAM_FILE" | samtools fastq - > "$UNMAPPED_FILE"
        local output=$(samtools view -f 4 "$SAM_FILE" | samtools fastq - 2>&1)

        if [[ $output == *"processed 0 reads"* ]]; then
            echo "Finished"
            break
        else
            local j=$((i + 1))
            local COLLAPSED_RUN_FILE=$(ls "$COLLAPSED_DIR/${SAMPLE_NAME}_run${i}"*_unmapped_trimmed.fastq | head -n 1)
            if [ ! -f "$COLLAPSED_RUN_FILE" ]; then
                seqcluster collapse -f "$UNMAPPED_FILE" -m 1 --min_size 15 -o "$COLLAPSED_DIR"
                if [ $? -ne 0 ]; then
                    echo "Seqcluster collapse failed for sample: $SAMPLE"
                    exit 1
                fi
                COLLAPSED_RUN_FILE=$(ls "$COLLAPSED_DIR/${SAMPLE_NAME}_run${i}"*_unmapped_trimmed.fastq | head -n 1)
                if [ ! -f "$COLLAPSED_RUN_FILE" ]; then
                    echo "ERROR: Collapsed run file not created: $COLLAPSED_RUN_FILE"
                    exit 1
                fi
            else
                echo "Collapsed run file already exists: $COLLAPSED_RUN_FILE. Skipping collapse step."
            fi
            
            (bowtie -x hairpin -q "$COLLAPSED_RUN_FILE" -k 50 -e 99999 --best --strata --sam --trim5 1 --trim3 1 > "${SAMPLE_DIR}/${SAMPLE_NAME}_run${j}.sam") 2> "${SAMPLE_DIR}/${SAMPLE_NAME}_run${j}.err"

            local processed=$(grep -o 'reads processed: [0-9]*' "${SAMPLE_DIR}/${SAMPLE_NAME}_run${j}.err" | awk '{print $3}')
            local alignment=$(grep -o 'reads with at least one alignment: [0-9]*' "${SAMPLE_DIR}/${SAMPLE_NAME}_run${j}.err" | awk '{print $7}')
            local failed=$(grep -o 'reads that failed to align: [0-9]*' "${SAMPLE_DIR}/${SAMPLE_NAME}_run${j}.err" | awk '{print $6}')
            local reported=$(grep -o 'Reported [0-9]* alignments' "${SAMPLE_DIR}/${SAMPLE_NAME}_run${j}.err" | awk '{print $2}')

            rm "${SAMPLE_DIR}/${SAMPLE_NAME}_run${j}.err"
            printf "%s %s %s %s %s %s\n" "$SAMPLE_NAME" "$j" "$processed" "$alignment" "$failed" "$reported" > "${SAMPLE_DIR}/${SAMPLE_NAME}_run${j}_summary.err"
        fi
    done

    # Collect all sorted BAM files for merging
    local SORTED_BAM_FILES=$(ls "${SAMPLE_DIR}/${SAMPLE_NAME}_run"*_srt.bam)
    if [ -z "$SORTED_BAM_FILES" ]; then
        echo "ERROR: No sorted BAM files found for merging in $SAMPLE_DIR"
        exit 1
    fi

    local MERGED_BAM_FILE="${SAMPLE_DIR}/${SAMPLE_NAME}.merged.bam"
    if [ ! -f "$MERGED_BAM_FILE" ]; then  # Check if merged BAM file exists
        samtools merge "$MERGED_BAM_FILE" $SORTED_BAM_FILES
        samtools index "$MERGED_BAM_FILE"
    else
        echo "Merged BAM file already exists: $MERGED_BAM_FILE. Skipping merging step."
    fi

    cat "${SAMPLE_DIR}/${SAMPLE_NAME}_run"*"_summary.err" > "${SAMPLE_DIR}/${SAMPLE_NAME}.err" && rm "${SAMPLE_DIR}/${SAMPLE_NAME}_run"*"_summary.err"
}

# Process each sample sequentially
for SAMPLE in $SAMPLES; do
    process_sample "$BASE_DIR" "$SAMPLE"
done

# Notify user of completion
echo "Processing complete."




