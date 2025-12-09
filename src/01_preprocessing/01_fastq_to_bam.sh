#!/bin/bash


################################################################################
# WES Alignment Pipeline - lane aware, fully automated
# Author: Haia K. (adapted & enhanced)
# Date: 2025-08-01
#
# Description:
#   Automated alignment pipeline for whole exome sequencing data. Processes
#   all FASTQ files in the input directory with lane-aware processing. 
#   For each sample:
#   - Aligns reads with BWA-MEM2
#   - Adds detailed read groups per lane
#   - Merges lanes if multiple
#   - Fixes mate information
#   - Sorts and marks duplicates
#   - Indexes final BAMs
#
# Inputs:
#   - FASTQ files (paired-end, gzipped)
#     Location: data/fastq/ or user-specified
#     Format: *_R1_*.fastq.gz, *_R2_*.fastq.gz
#   - Reference genome (indexed with BWA)
#     Location: data/reference_files/
#   - Sample list (to be configured via cohort config - see Point #2)
#
# Outputs:
#   - results/bams/{sample}_final.bam - Sorted, deduplicated BAM per sample
#   - results/bams/{sample}_final.bam.bai - BAM indexes
#   - logs/{sample}.log - Per-sample processing logs
#
# Dependencies:
#   - bwa-mem (or bwa)
#   - samtools (≥1.10)
#
# Configuration:
#   Edit variables in "User Configuration" section below
#   Sample-specific mapping will be moved to cohort_config.yaml (Point #2)
#
# Usage:
#   # From repository root:
#   bash src/preprocessing/01_fastq_to_bam.sh
#   
#   # Or run in background with logging:
#   nohup bash src/preprocessing/01_fastq_to_bam.sh > alignment.log 2>&1 &
#
################################################################################
set -euo pipefail

# ============================================================================
# AUTO-DETECT REPOSITORY STRUCTURE
# ============================================================================

# Get script location (src/preprocessing/)
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
# Go up 2 levels to repo root
REPO_ROOT="$(dirname "$(dirname "$SCRIPT_DIR")")"

# Define standard paths relative to repo root
DATA_DIR="$REPO_ROOT/data"
REFERENCE_DIR="$DATA_DIR/reference_files"
RESULTS_DIR="$REPO_ROOT/results"
BAMS_DIR="$RESULTS_DIR/bams"
LOGS_DIR="$REPO_ROOT/logs"

# Create output directories if they don't exist
mkdir -p "$BAMS_DIR"
mkdir -p "$LOGS_DIR"

# ============================================================================
# USER CONFIGURATION
# ============================================================================

# Reference genome (place in data/reference_files/)
REFERENCE_GENOME="$REFERENCE_DIR/Homo_sapiens_assembly19.fasta"

# FASTQ directory (adjust if needed)
FASTQ_DIR="$DATA_DIR/fastq"

# Metadata for read groups
LIB="SureSelectV8"          # Library preparation kit
PLATFORM="ILLUMINA"         # Sequencing platform
THREADS=40                  # CPU threads for alignment

# ============================================================================
# SAMPLE CONFIGURATION (TO BE REPLACED WITH cohort_config.yaml)
# ============================================================================
# NOTE: This sample mapping will be moved to a cohort configuration file
# in Point #2. For now, configure your samples here OR use auto-discovery.

# Sample discovery mode: "auto" or "manual"
SAMPLE_DISCOVERY="auto"  # Set to "manual" to use SAMPLES array below

# Manual sample definition (used only if SAMPLE_DISCOVERY="manual")
declare -A SAMPLES=(
    # Format: [SampleName]="R1_file1.fastq.gz,R2_file1.fastq.gz,R1_file2.fastq.gz,R2_file2.fastq.gz"
    # Add your samples below:
    
    # Example:
    # [Sample1]="Sample1_L001_R1_001.fastq.gz,Sample1_L001_R2_001.fastq.gz"
    # [Sample2]="Sample2_L001_R1_001.fastq.gz,Sample2_L001_R2_001.fastq.gz,Sample2_L002_R1_001.fastq.gz,Sample2_L002_R2_001.fastq.gz"
)

# ============================================================================
# SAMPLE DISCOVERY FUNCTION
# ============================================================================

function discover_samples() {
    echo "[$(date)] Auto-discovering samples from FASTQ directory..."
    
    # Find all unique sample names from R1 files
    # Assumes naming pattern: {SAMPLE}_*_R1_*.fastq.gz
    local sample_names=()
    
    for r1_file in "$FASTQ_DIR"/*_R1_*.fastq.gz; do
        [[ -f "$r1_file" ]] || continue
        
        # Extract sample name (everything before first underscore or common pattern)
        local basename=$(basename "$r1_file")
        
        # Try to extract sample name (adjust pattern as needed for your naming)
        # Common patterns:
        # - SampleName_S1_L001_R1_001.fastq.gz → SampleName
        # - SampleName_L001_R1_001.fastq.gz → SampleName
        # - SampleName_R1.fastq.gz → SampleName
        
        local sample_name
        if [[ "$basename" =~ ^([^_]+)_.*_R1_ ]]; then
            sample_name="${BASH_REMATCH[1]}"
        elif [[ "$basename" =~ ^([^_]+)_R1 ]]; then
            sample_name="${BASH_REMATCH[1]}"
        else
            # Fallback: use everything before _R1
            sample_name="${basename%%_R1*}"
        fi
        
        # Add to list if not already present
        if [[ ! " ${sample_names[@]} " =~ " ${sample_name} " ]]; then
            sample_names+=("$sample_name")
        fi
    done
    
    if [[ ${#sample_names[@]} -eq 0 ]]; then
        echo "ERROR: No samples discovered in $FASTQ_DIR"
        echo "       No files matching pattern *_R1_*.fastq.gz found"
        echo "       Either add FASTQ files or use SAMPLE_DISCOVERY=\"manual\""
        exit 1
    fi
    
    echo "  Discovered ${#sample_names[@]} samples:"
    
    # Build SAMPLES array automatically
    declare -gA SAMPLES
    for sample_name in "${sample_names[@]}"; do
        # Find all R1 and R2 files for this sample
        local r1_files=()
        local r2_files=()
        
        for file in "$FASTQ_DIR"/${sample_name}*_R1_*.fastq.gz; do
            [[ -f "$file" ]] && r1_files+=($(basename "$file"))
        done
        
        for file in "$FASTQ_DIR"/${sample_name}*_R2_*.fastq.gz; do
            [[ -f "$file" ]] && r2_files+=($(basename "$file"))
        done
        
        # Combine into comma-separated list
        local all_files=()
        all_files+=("${r1_files[@]}")
        all_files+=("${r2_files[@]}")
        
        if [[ ${#all_files[@]} -gt 0 ]]; then
            local files_string=$(IFS=','; echo "${all_files[*]}")
            SAMPLES["$sample_name"]="$files_string"
            echo "    - $sample_name (${#r1_files[@]} lane(s))"
        fi
    done
    
    echo ""
}

# ============================================================================
# VALIDATION
# ============================================================================

# Check dependencies
command -v bwa >/dev/null 2>&1 || { echo "ERROR: bwa not found. Install bwa-mem2 or bwa."; exit 1; }
command -v samtools >/dev/null 2>&1 || { echo "ERROR: samtools not found."; exit 1; }

# Check reference genome
if [[ ! -f "$REFERENCE_GENOME" ]]; then
    echo "ERROR: Reference genome not found: $REFERENCE_GENOME"
    echo "Please place reference genome in: $REFERENCE_DIR/"
    exit 1
fi

# Check FASTQ directory
if [[ ! -d "$FASTQ_DIR" ]]; then
    echo "ERROR: FASTQ directory not found: $FASTQ_DIR"
    exit 1
fi

# ============================================================================
# ALIGNMENT FUNCTION (merged from 1-2_align_WES_by_lane.sh)
# ============================================================================

function align_sample() {
    local SAMPLE_NAME="$1"
    local R1_FILES="$2"
    local R2_FILES="$3"
    
    local START_TIME=$(date +%s)
    echo "[$(date)] ===== Starting alignment for: $SAMPLE_NAME ====="
    
    # Convert comma-separated strings to arrays
    IFS=',' read -r -a R1_ARRAY <<< "$R1_FILES"
    IFS=',' read -r -a R2_ARRAY <<< "$R2_FILES"
    
    # Validate R1/R2 pairing
    if [[ ${#R1_ARRAY[@]} -ne ${#R2_ARRAY[@]} ]]; then
        echo "ERROR: Mismatch in R1/R2 files for $SAMPLE_NAME"
        return 1
    fi
    
    # Process each lane
    local LANE_BAMS=()
    for i in "${!R1_ARRAY[@]}"; do
        local R1="${R1_ARRAY[i]}"
        local R2="${R2_ARRAY[i]}"
        
        # Extract read group info from FASTQ header
        local first_header=$(zcat "$R1" | head -n 1)
        if [[ -z "$first_header" ]]; then
            echo "ERROR: Cannot read FASTQ header from $R1"
            return 1
        fi
        
        local flowcell=$(echo "$first_header" | cut -d':' -f3)
        local lane_id=$(echo "$first_header" | cut -d':' -f4)
        local read_group="@RG\tID:${flowcell}_${lane_id}\tSM:${SAMPLE_NAME}\tLB:${LIB}\tPL:${PLATFORM}\tPU:${flowcell}_${lane_id}"
        
        echo "  [$(date '+%H:%M:%S')] Processing lane $lane_id"
        
        # Align with BWA-MEM
        local LANE_BAM="$BAMS_DIR/${SAMPLE_NAME}_${lane_id}.bam"
        LANE_BAMS+=("$LANE_BAM")
        
        bwa mem -t "$THREADS" -T 0 -R "$read_group" "$REFERENCE_GENOME" "$R1" "$R2" | \
            samtools view -Shb -o "$LANE_BAM" -
    done
    
    # Merge lanes if multiple
    local MERGED_BAM="$BAMS_DIR/${SAMPLE_NAME}_merged.bam"
    if [[ ${#LANE_BAMS[@]} -gt 1 ]]; then
        echo "  [$(date '+%H:%M:%S')] Merging ${#LANE_BAMS[@]} lanes"
        samtools merge -@ "$THREADS" -f "$MERGED_BAM" "${LANE_BAMS[@]}"
        rm -f "${LANE_BAMS[@]}"  # Clean up lane BAMs
    else
        mv "${LANE_BAMS[0]}" "$MERGED_BAM"
    fi
    
    # Fix mate information
    local FIXMATE_BAM="$BAMS_DIR/${SAMPLE_NAME}_fixmate.bam"
    echo "  [$(date '+%H:%M:%S')] Fixing mate pairs"
    samtools fixmate -@ "$THREADS" -m "$MERGED_BAM" "$FIXMATE_BAM"
    rm -f "$MERGED_BAM"
    
    # Sort
    local SORTED_BAM="$BAMS_DIR/${SAMPLE_NAME}_sorted.bam"
    echo "  [$(date '+%H:%M:%S')] Sorting"
    samtools sort -@ "$THREADS" -o "$SORTED_BAM" "$FIXMATE_BAM"
    rm -f "$FIXMATE_BAM"
    
    # Mark duplicates
    local FINAL_BAM="$BAMS_DIR/${SAMPLE_NAME}_final.bam"
    echo "  [$(date '+%H:%M:%S')] Marking duplicates"
    samtools markdup -@ "$THREADS" "$SORTED_BAM" "$FINAL_BAM"
    rm -f "$SORTED_BAM"
    
    # Index
    echo "  [$(date '+%H:%M:%S')] Indexing"
    samtools index "$FINAL_BAM"
    
    local END_TIME=$(date +%s)
    local DURATION=$((END_TIME - START_TIME))
    echo "[$(date)] ===== Completed: $SAMPLE_NAME in $((DURATION / 60))m $((DURATION % 60))s ====="
    echo ""
}

# ============================================================================
# MAIN EXECUTION
# ============================================================================

echo "################################################################################"
echo "# WES Alignment Pipeline"
echo "# Repository: $REPO_ROOT"
echo "# Output: $BAMS_DIR"
echo "# Reference: $REFERENCE_GENOME"
echo "# Threads: $THREADS"
echo "################################################################################"
echo ""

# ============================================================================
# Sample Discovery or Manual Definition
# ============================================================================

if [[ "$SAMPLE_DISCOVERY" == "auto" ]]; then
    discover_samples
else
    if [[ ${#SAMPLES[@]} -eq 0 ]]; then
        echo "ERROR: No samples configured!"
        echo "Either:"
        echo "  - Set SAMPLE_DISCOVERY=\"auto\" for automatic discovery"
        echo "  - Or define samples manually in SAMPLES array"
        echo "  - Or wait for cohort_config.yaml (Point #2)"
        exit 1
    fi
    echo "[$(date)] Using ${#SAMPLES[@]} manually defined samples"
    echo ""
fi

# Process each sample
for SAMPLE_NAME in "${!SAMPLES[@]}"; do
    FASTQ_FILES="${SAMPLES[$SAMPLE_NAME]}"
    IFS=',' read -r -a FASTQ_ARRAY <<< "$FASTQ_FILES"
    
    # Separate R1 and R2 files
    R1_FILES=()
    R2_FILES=()
    for FASTQ_FILE in "${FASTQ_ARRAY[@]}"; do
        FULL_PATH="$FASTQ_DIR/$FASTQ_FILE"
        
        if [[ ! -f "$FULL_PATH" ]]; then
            echo "WARNING: File not found: $FULL_PATH (skipping)"
            continue
        fi
        
        if [[ "$FASTQ_FILE" == *_R1_* ]]; then
            R1_FILES+=("$FULL_PATH")
        elif [[ "$FASTQ_FILE" == *_R2_* ]]; then
            R2_FILES+=("$FULL_PATH")
        fi
    done
    
    # Validate we have files
    if [[ ${#R1_FILES[@]} -eq 0 ]]; then
        echo "WARNING: No R1 files found for $SAMPLE_NAME (skipping)"
        continue
    fi
    
    # Join arrays back to comma-separated strings
    R1_FILES_JOINED=$(IFS=','; echo "${R1_FILES[*]}")
    R2_FILES_JOINED=$(IFS=','; echo "${R2_FILES[*]}")
    
    # Run alignment (can be parallelized by adding & and wait at end)
    align_sample "$SAMPLE_NAME" "$R1_FILES_JOINED" "$R2_FILES_JOINED" 2>&1 | tee "$LOGS_DIR/${SAMPLE_NAME}.log"
done

echo "################################################################################"
echo "# All alignments completed!"
echo "# Output BAMs: $BAMS_DIR/*.bam"
echo "# Logs: $LOGS_DIR/*.log"
echo "################################################################################"