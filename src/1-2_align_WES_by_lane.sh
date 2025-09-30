#!/bin/bash

# Help message
function usage() {
    echo "Usage: $0 -s SAMPLE_NAME -r REFERENCE_GENOME -o BAM_DIR -1 LANE_R1_FILES -2 LANE_R2_FILES -l LIB -p PLATFORM -t THREADS"
    echo
    echo "Parameters:"
    echo "  -s SAMPLE_NAME       Sample name (e.g., SAMPLE1)"
    echo "  -r REFERENCE_GENOME  Path to the reference genome FASTA file"
    echo "  -o BAM_DIR           Output directory for BAM files"
    echo "  -1 LANE_R1_FILES     Comma-separated list of R1 FASTQ files for each lane"
    echo "  -2 LANE_R2_FILES     Comma-separated list of R2 FASTQ files for each lane"
    echo "  -l LIB               Library name (e.g., ExomeTwist)"
    echo "  -p PLATFORM          Sequencing platform (e.g., ILLUMINA)"
    echo "  -t THREADS           Number of threads (default: 8)"
    exit 1
}

# Default number of threads
THREADS=40 

# Parse named parameters
while getopts "s:r:o:1:2:l:p:t:" opt; do
    case $opt in
        s) SAMPLE_NAME="$OPTARG" ;;
        r) REFERENCE_GENOME="$OPTARG" ;;
        o) BAM_DIR="$OPTARG" ;;
        1) LANE_R1_FILES="$OPTARG" ;;
        2) LANE_R2_FILES="$OPTARG" ;;
        l) LIB="$OPTARG" ;;
        p) PLATFORM="$OPTARG" ;;
        t) THREADS="$OPTARG" ;;
        *) usage ;;
    esac
done

# Ensure all required parameters are provided
if [[ -z "$SAMPLE_NAME" || -z "$REFERENCE_GENOME" || -z "$BAM_DIR" || -z "$LANE_R1_FILES" || -z "$LANE_R2_FILES" || -z "$LIB" || -z "$PLATFORM" ]]; then
    usage
fi

# Create output directory if it doesn't exist
if [[ ! -d "$BAM_DIR" ]]; then
    echo "Creating output directory: $BAM_DIR"
    mkdir -p "$BAM_DIR"
fi

# Convert comma-separated lists to arrays
IFS=',' read -r -a R1_FILES <<< "$LANE_R1_FILES"
IFS=',' read -r -a R2_FILES <<< "$LANE_R2_FILES"

# Ensure the number of R1 and R2 files match
if [[ ${#R1_FILES[@]} -ne ${#R2_FILES[@]} ]]; then
    echo "Error: Mismatch in the number of R1 and R2 files."
    exit 1
fi


START_TIME=$(date +%s)
echo "[$(date)] Starting pipeline for sample: $SAMPLE_NAME"

# Sample and lane-specific BAM outputs
LANE_BAMS=()

# Process each lane
LANE_BAMS=()
for i in "${!R1_FILES[@]}"; do
    R1="${R1_FILES[i]}"
    R2="${R2_FILES[i]}"

    # Extract read group details
    first_header=$(zcat "$R1" | head -n 1)
    if [[ -z "$first_header" ]]; then
        echo "Error: Unable to extract header from FASTQ file $R1."
        exit 1
    fi

    flowcell=$(echo "$first_header" | cut -d':' -f3)
    lane_id=$(echo "$first_header" | cut -d':' -f4)
    
    read_group="@RG\tID:${flowcell}_${lane_id}\tSM:${SAMPLE_NAME}\tLB:${LIB}\tPL:${PLATFORM}\tPU:${flowcell}_${lane_id}"

    echo "Processing lane $lane_id with read group: $read_group [$(date '+%Y-%m-%d %H:%M:%S')]"

    # Define lane-specific BAM output
    LANE_BAM="$BAM_DIR/${SAMPLE_NAME}_${lane_id}.bam"
    LANE_BAMS+=("$LANE_BAM")

    # Calculate mean read length
    mean_read_length=$(zcat "$R1" | awk 'NR % 4 == 2 { total += length($0); count++ } END { if (count > 0) print total / count; else print 0 }')
    echo "Mean read length for lane $lane_id: $mean_read_length"

    # Align reads based on mean read length
#    if (( $(echo "$mean_read_length >= 70" | bc -l) )); then
        echo "Running BWA-MEM for lane $lane_id ... [$(date '+%Y-%m-%d %H:%M:%S')] "
        bwa mem -t "$THREADS" -T 0 -R "$read_group" "$REFERENCE_GENOME" "$R1" "$R2" | samtools view -Shb -o "$LANE_BAM" -
#    else
#        echo "Running BWA-ALN for lane $lane_id..."
#        bwa aln -t "$THREADS" "$REFERENCE_GENOME" "$R1" > sai_1.sai
#        bwa aln -t "$THREADS" "$REFERENCE_GENOME" "$R2" > sai_2.sai
#        bwa sampe -r "$read_group" "$REFERENCE_GENOME" sai_1.sai sai_2.sai "$R1" "$R2" | samtools view -Shb -o "$LANE_BAM" -
#        rm sai_1.sai sai_2.sai
#    fi
done

# Merge lane BAMs
MERGED_BAM="$BAM_DIR/${SAMPLE_NAME}_merged.bam"
if [ "${#LANE_BAMS[@]}" -gt 1 ]; then
    samtools merge -@ "$THREADS" -f "$MERGED_BAM" "${LANE_BAMS[@]}"
else
    cp "${LANE_BAMS[0]}" "$MERGED_BAM"
fi

# Fix mates
FIXMATE_BAM="$BAM_DIR/${SAMPLE_NAME}_fixmate.bam"
samtools fixmate -@ "$THREADS" -m "$MERGED_BAM" "$FIXMATE_BAM"

# Sort BAM
SORTED_BAM="$BAM_DIR/${SAMPLE_NAME}_sorted.bam"
samtools sort -@ "$THREADS" -o "$SORTED_BAM" "$FIXMATE_BAM"

# Mark duplicates
FINAL_BAM="$BAM_DIR/${SAMPLE_NAME}_final.bam"
samtools markdup -@ "$THREADS" "$SORTED_BAM" "$FINAL_BAM"

# BAM
samtools index "$FINAL_BAM" 

END_TIME=$(date +%s)
DURATION=$((END_TIME - START_TIME))
echo "[$(date)] Finished sample: $SAMPLE_NAME in $((DURATION / 60)) min $((DURATION % 60)) sec" | tee -a pipeline_times.log
echo "Pipeline completed successfully. Final BAM: $FINAL_BAM"

rm -f "$MERGED_BAM"
rm -f "$FIXMATE_BAM"
rm -f "$SORTED_BAM"