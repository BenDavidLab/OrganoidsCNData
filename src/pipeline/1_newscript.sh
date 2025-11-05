#!/bin/bash


################################################################################
# WES Alignment Pipeline - lane aware, fully automated
# Author: Haia K. (adapted & enhanced)
# Date: 2025-08-01
#
# Description:
#   This script automates alignment of all paired-end WES FASTQ files found in a
#   directory. It processes each lane separately, adds detailed read groups,
#   marks duplicates, and produces sorted, deduplicated BAM files per lane.
#
# Inputs:
#   - FASTQ_DIR: directory with all FASTQ files (with lanes in filenames)
#   - REFERENCE_GENOME: indexed reference genome fasta file
#   - OUT_DIR: output directory for BAMs and metrics
#   - THREADS: number of CPU threads for alignment and sorting
#   - LIB: library name for RG tags
#   - PLATFORM: sequencing platform for RG tags
#
# Outputs:
#   - Per-lane sorted, deduplicated BAMs and their indexes
#   - Deduplication metrics files
#   - Log file with timestamps and sample processing info
#
# Dependencies:
#   - bwa-mem2
#   - samtools 
# 
# Usage:
#   nohup bash newscript.sh > pipeline.log 2>&1 &
################################################################################


# Path to the alignment script
ALIGN_SCRIPT="./align_WES_by_lane.sh"

# Reference genome
REFERENCE_GENOME="/workspace/RefGen/Homo_sapiens_assembly19.fasta"

# Output BAMs and logs
OUTPUT_DIR="/workspace/bams"
mkdir -p "$OUTPUT_DIR"
mkdir -p ./logs

# FASTQ directory
FASTQ_DIR="/workspace/EGAD00001010134"

# Metadata
LIB="SureSelectV8"
PLATFORM="ILLUMINA"
THREADS=40

# Sample-to-FASTQ mapping
declare -A SAMPLES=(
    [46TOrg]="46TOrg_XXXXXX_s1_L001_R1_001.fastq.gz,46TOrg_XXXXXX_s1_L001_R2_001.fastq.gz"
    [46TTtis]="46TTtis_XXXXXX_s1_L001_R1_001.fastq.gz,46TTtis_XXXXXX_s1_L001_R2_001.fastq.gz"
    [M19-00991Org]="M19-00991Org_XXXXX_S1_L001_R1_001.fastq.gz,M19-00991Org_XXXXX_S1_L001_R2_001.fastq.gz"
    [M19-00991Ttis]="M19-00991Ttis_XXXXXXX_S1_L001_R1_001.fastq.gz,M19-00991Ttis_XXXXXXX_S1_L001_R2_001.fastq.gz"
    [M19-258Org]="M19-258Org_XXXXXX_s1_L001_R1_001.fastq.gz,M19-258Org_XXXXXX_s1_L001_R2_001.fastq.gz"
    [M19-258Ttis]="M19-258Ttis_XXXXXX_s1_L001_R1_001.fastq.gz,M19-258Ttis_XXXXXX_s1_L001_R2_001.fastq.gz"
    [M19-432Org]="M19-432Org_XXXXXX_s1_L001_R1_001.fastq.gz,M19-432Org_XXXXXX_s1_L001_R2_001.fastq.gz"
    [M19-432Ttis]="M19-432Ttis_XXXXXX_s1_L001_R1_001.fastq.gz,M19-432Ttis_XXXXXX_s1_L001_R2_001.fastq.gz"
    [M20-00178Org]="M20-00178Org_XXXXXX_S1_L001_R1_001.fastq.gz,M20-00178Org_XXXXXX_S1_L001_R2_001.fastq.gz"
    [M20-00178Ttis]="M20-00178Ttis_XXXXXX_S1_L001_R1_001.fastq.gz,M20-00178Ttis_XXXXXX_S1_L001_R2_001.fastq.gz"
    [M20-00383Org]="M20-00383Org_XXXXXX_S1_L001_R1_001.fastq.gz,M20-00383Org_XXXXXX_S1_L001_R2_001.fastq.gz"
    [M20-00383Ttis]="M20-00383Ttis_XXXXXX_S1_L001_R1_001.fastq.gz,M20-00383Ttis_XXXXXX_S1_L001_R2_001.fastq.gz"
    [M22-32Org]="M22-32Org_XXXXXX_s1_L001_R1_001.fastq.gz,M22-32Org_XXXXXX_s1_L001_R2_001.fastq.gz"
    [M22-32Ttis]="M22-32Ttis_XXXXXX_s1_L001_R1_001.fastq.gz,M22-32Ttis_XXXXXX_s1_L001_R2_001.fastq.gz"
)

# Run each sample alignment in parallel
for SAMPLE_NAME in "${!SAMPLES[@]}"; do
    echo "Launching alignment for sample: $SAMPLE_NAME"

    FASTQ_FILES="${SAMPLES[$SAMPLE_NAME]}"
    IFS=',' read -r -a FASTQ_ARRAY <<< "$FASTQ_FILES"

    R1_FILES=()
    R2_FILES=()
    for FASTQ_FILE in "${FASTQ_ARRAY[@]}"; do
        if [[ "$FASTQ_FILE" == *_R1_* ]]; then
            R1_FILES+=("${FASTQ_DIR}/${FASTQ_FILE}")
        elif [[ "$FASTQ_FILE" == *_R2_* ]]; then
            R2_FILES+=("${FASTQ_DIR}/${FASTQ_FILE}")
        fi
    done

    R1_FILES_JOINED=$(IFS=','; echo "${R1_FILES[*]}")
    R2_FILES_JOINED=$(IFS=','; echo "${R2_FILES[*]}")

    # Launch alignment in background
    bash "$ALIGN_SCRIPT" \
        -s "$SAMPLE_NAME" \
        -r "$REFERENCE_GENOME" \
        -o "$OUTPUT_DIR" \
        -1 "$R1_FILES_JOINED" \
        -2 "$R2_FILES_JOINED" \
        -l "$LIB" \
        -p "$PLATFORM" \
        -t "$THREADS" > ./logs/$SAMPLE_NAME.log 2>&1 
done

# Wait for all background processes
wait
echo "🎉 All parallel WES alignments completed."


