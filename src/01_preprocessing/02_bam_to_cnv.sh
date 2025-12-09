#!/bin/bash


################################################################################
# WES CNVkit Pipeline - multi-sample, fully automated
# Author: Haia K. (adapted & enhanced)
# Date: 2025-08-18
#
# Description:
#   Automated CNV analysis pipeline using CNVkit. For each cohort:
#   - Generates target and antitarget bins from BED file
#   - Calculates coverage for all BAMs (target + antitarget)
#   - Builds pooled reference from normal samples
#   - Normalizes tumor/organoid samples against reference
#   - Segments CNV regions
#   - Calls absolute copy numbers
#   - Optional: Filters extreme log2 values
#   - Exports to .seg format
#   - Generates visualization plots
#
# Inputs:
#   - BAM files from script 01
#     Location: results/bams/*_final.bam
#   - Target regions BED file
#     Location: data/reference_files/
#   - Reference genome (indexed)
#     Location: data/reference_files/
#
# Outputs:
#   - results/cnv/{cohort}/*.cnn - Coverage files
#   - results/cnv/{cohort}/cnv_reference.cnn - Reference
#   - results/cnv/{cohort}/*.cnr - Normalized log2 ratios
#   - results/cnv/{cohort}/*.cns - Segmented CNVs
#   - results/cnv/{cohort}/*.called.cns - CNV calls with absolute CN
#   - results/cnv/{cohort}/*.seg - Exported segments
#   - results/cnv/{cohort}/*_scatter.pdf - Visualization
#   - results/cnv/{cohort}/*_diagram.pdf - Genome-wide plots
#
# Dependencies:
#   - CNVkit (>= 0.9.9)
#   - Python 3 with pandas
#   - samtools
#   - R with DNAcopy (for CBS segmentation method)
#
# Configuration:
#   Edit variables in "User Configuration" section below
#   Normal sample patterns will be moved to cohort_config.yaml (Point #2)
#
# Usage:
#   # From repository root:
#   bash src/preprocessing/02_bam_to_cnv.sh
#   
#   # Or with specific cohort:
#   bash src/preprocessing/02_bam_to_cnv.sh COHORT_NAME
################################################################################

set -euo pipefail
shopt -s nullglob

# ============================================================================
# AUTO-DETECT REPOSITORY STRUCTURE
# ============================================================================

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
REPO_ROOT="$(dirname "$(dirname "$SCRIPT_DIR")")"

# Define standard paths
DATA_DIR="$REPO_ROOT/data"
REFERENCE_DIR="$DATA_DIR/reference_files"
RESULTS_DIR="$REPO_ROOT/results"
BAMS_DIR="$RESULTS_DIR/bams"
CNV_DIR="$RESULTS_DIR/cnv"

# ============================================================================
# USER CONFIGURATION
# ============================================================================

# Cohort name (can be passed as argument or set here)
COHORT_NAME="${1:-default_cohort}"
COHORT_CNV_DIR="$CNV_DIR/$COHORT_NAME"
mkdir -p "$COHORT_CNV_DIR"

# Reference files (place in data/reference_files/)
REFERENCE_GENOME="$REFERENCE_DIR/Homo_sapiens_assembly19.fasta"
TARGETS_BED="$REFERENCE_DIR/exome_targets.bed"  # Adjust filename

# CNVkit parameters
THREADS=32
SEG_METHOD="cbs"  # Options: cbs, hmm, haar, flasso

# Normal sample pattern (will be moved to cohort_config.yaml in Point #2)
# Adjust this glob pattern to match your normal/control samples
NORMAL_GLOB="*T_final.bam"  # e.g., "*Normal*_final.bam", "*T_final.bam"

# Optional CNR filtering (integrated from 2-2_filter_cnr.py)
ENABLE_CNR_FILTERING=false
CNR_LOG2_THRESHOLD=-2.0  # Remove bins with log2 < this value

# ============================================================================
# HELPER FUNCTIONS
# ============================================================================

timestamp() { date "+%Y-%m-%d %H:%M:%S"; }

# Integrated CNR filtering function (from 2-2_filter_cnr.py)
filter_cnr() {
    local input_cnr="$1"
    local output_cnr="$2"
    
    if [[ "$ENABLE_CNR_FILTERING" != "true" ]]; then
        return 0
    fi
    
    echo "  [$(timestamp)] Filtering CNR (removing log2 < $CNR_LOG2_THRESHOLD)"
    
    python3 - "$input_cnr" "$output_cnr" "$CNR_LOG2_THRESHOLD" << 'PYTHON_FILTER'
import sys
import pandas as pd

cnr_file = sys.argv[1]
output_file = sys.argv[2]
threshold = float(sys.argv[3])

# Load CNR file
cnr_data = pd.read_csv(cnr_file, sep="\t")

# Filter rows where log2 >= threshold
filtered_cnr = cnr_data[cnr_data["log2"] >= threshold]

# Save filtered file
filtered_cnr.to_csv(output_file, sep="\t", index=False)
print(f"  Filtered: {len(cnr_data) - len(filtered_cnr)} bins removed")
PYTHON_FILTER
}

# ============================================================================
# VALIDATION
# ============================================================================

# Check dependencies
command -v cnvkit.py >/dev/null 2>&1 || { echo "ERROR: CNVkit not found. Install: pip install cnvkit"; exit 1; }
command -v python3 >/dev/null 2>&1 || { echo "ERROR: python3 not found."; exit 1; }
command -v samtools >/dev/null 2>&1 || { echo "ERROR: samtools not found."; exit 1; }

# Check reference genome
if [[ ! -f "$REFERENCE_GENOME" ]]; then
    echo "ERROR: Reference genome not found: $REFERENCE_GENOME"
    exit 1
fi

# Check targets BED
if [[ ! -f "$TARGETS_BED" ]]; then
    echo "ERROR: Targets BED not found: $TARGETS_BED"
    echo "Please place your exome targets BED file in: $REFERENCE_DIR/"
    exit 1
fi

# Check BAMs directory
if [[ ! -d "$BAMS_DIR" ]]; then
    echo "ERROR: BAMs directory not found: $BAMS_DIR"
    echo "Run 01_fastq_to_bam.sh first!"
    exit 1
fi

# ============================================================================
# SETUP TARGET/ANTITARGET REGIONS
# ============================================================================

echo "################################################################################"
echo "# CNVkit Pipeline - Cohort: $COHORT_NAME"
echo "# Working directory: $COHORT_CNV_DIR"
echo "# BAMs: $BAMS_DIR"
echo "# Reference: $REFERENCE_GENOME"
echo "################################################################################"
echo ""

cd "$COHORT_CNV_DIR"

TARGETS="${TARGETS_BED%.bed}.split.bed"
ANTITARGETS="${TARGETS_BED%.bed}.antitargets.bed"

# Generate target regions if not present
if [[ ! -s "$TARGETS" ]]; then
    echo "[$(timestamp)] Generating target regions"
    cnvkit.py target "$TARGETS_BED" --split -o "$TARGETS"
fi

# Generate antitarget regions if not present
if [[ ! -s "$ANTITARGETS" ]]; then
    echo "[$(timestamp)] Generating antitarget regions"
    cnvkit.py antitarget "$TARGETS_BED" -o "$ANTITARGETS"
fi

# ============================================================================
# BUILD POOLED REFERENCE FROM NORMALS
# ============================================================================

REFERENCE_CNN="cnv_reference.cnn"

if [[ ! -s "$REFERENCE_CNN" ]]; then
    echo "[$(timestamp)] Building CNVkit reference from normals"
    echo "  Normal pattern: $BAMS_DIR/$NORMAL_GLOB"
    
    declare -a NORMAL_TARGET_CNNS=()
    declare -a NORMAL_ANTITARGET_CNNS=()
    
    for NORMAL_BAM in "$BAMS_DIR"/$NORMAL_GLOB; do
        [[ -f "$NORMAL_BAM" ]] || continue
        
        BASENAME=$(basename "$NORMAL_BAM" _final.bam)
        T_CNN="${BASENAME}.targetcoverage.cnn"
        A_CNN="${BASENAME}.antitargetcoverage.cnn"
        
        echo "  Processing normal: $BASENAME"
        
        if [[ ! -s "$T_CNN" ]]; then
            cnvkit.py coverage "$NORMAL_BAM" "$TARGETS" -o "$T_CNN" --processes "$THREADS"
        fi
        if [[ ! -s "$A_CNN" ]]; then
            cnvkit.py coverage "$NORMAL_BAM" "$ANTITARGETS" -o "$A_CNN" --processes "$THREADS"
        fi
        
        NORMAL_TARGET_CNNS+=("$T_CNN")
        NORMAL_ANTITARGET_CNNS+=("$A_CNN")
    done
    
    if [[ ${#NORMAL_TARGET_CNNS[@]} -eq 0 ]]; then
        echo "ERROR: No normal samples matched pattern '$NORMAL_GLOB' in $BAMS_DIR"
        echo "Please adjust NORMAL_GLOB variable or check your BAM file naming"
        exit 1
    fi
    
    echo "  Creating reference from ${#NORMAL_TARGET_CNNS[@]} normal samples"
    cnvkit.py reference "${NORMAL_TARGET_CNNS[@]}" "${NORMAL_ANTITARGET_CNNS[@]}" \
        -f "$REFERENCE_GENOME" -o "$REFERENCE_CNN"
    
    echo "[$(timestamp)] Reference created: $REFERENCE_CNN"
else
    echo "[$(timestamp)] Using existing reference: $REFERENCE_CNN"
fi

# ============================================================================
# PROCESS MATCHED PAIRS (Tumor/Organoid ↔ Normal Tissue)
# ============================================================================

echo ""
echo "[$(timestamp)] Processing matched sample pairs"
echo "------------------------------------------------"

for ORG_BAM in "$BAMS_DIR"/$NORMAL_GLOB; do
    [[ -f "$ORG_BAM" ]] || continue
    
    base=$(basename "$ORG_BAM")
    # Pattern matching: change T to NT (adjust as needed for your naming)
    mate=${base/T_final/NT_final}
    TISSUE_BAM="$BAMS_DIR/$mate"
    
    if [[ ! -f "$TISSUE_BAM" ]]; then
        echo "[WARN $(timestamp)] No matching tissue BAM for $(basename $ORG_BAM)"
        continue
    fi
    
    echo ""
    echo "[$(timestamp)] Processing pair:"
    echo "  Sample 1: $(basename $ORG_BAM)"
    echo "  Sample 2: $(basename $TISSUE_BAM)"
    
    # Process both samples in the pair
    for BAM in "$ORG_BAM" "$TISSUE_BAM"; do
        SAMPLE=$(basename "$BAM" .bam)
        
        # Skip if already processed (check for final diagram)
        if [[ -s "${SAMPLE}_diagram.pdf" ]]; then
            echo "  [$(timestamp)] $SAMPLE already processed, skipping"
            continue
        fi
        
        echo "  [$(timestamp)] Processing: $SAMPLE"
        
        # Step 1: Calculate coverage
        T_CNN="${SAMPLE}.targetcoverage.cnn"
        A_CNN="${SAMPLE}.antitargetcoverage.cnn"
        
        if [[ ! -s "$T_CNN" ]]; then
            cnvkit.py coverage "$BAM" "$TARGETS" -o "$T_CNN" --processes "$THREADS"
        fi
        if [[ ! -s "$A_CNN" ]]; then
            cnvkit.py coverage "$BAM" "$ANTITARGETS" -o "$A_CNN" --processes "$THREADS"
        fi
        
        # Step 2: Normalize (fix)
        CNR="${SAMPLE}.cnr"
        cnvkit.py fix "$T_CNN" "$A_CNN" "$REFERENCE_CNN" -o "$CNR"
        
        # Step 3: Optional filtering (integrated from 2-2_filter_cnr.py)
        if [[ "$ENABLE_CNR_FILTERING" == "true" ]]; then
            FILTERED_CNR="${SAMPLE}.filtered.cnr"
            filter_cnr "$CNR" "$FILTERED_CNR"
            CNR="$FILTERED_CNR"
        fi
        
        # Step 4: Segment
        CNS="${SAMPLE}.cns"
        if [[ "$SEG_METHOD" == "hmm" ]]; then
            cnvkit.py segment "$CNR" -o "$CNS"
        else
            cnvkit.py segment "$CNR" --method "$SEG_METHOD" -o "$CNS"
        fi
        
        # Step 5: Call absolute copy numbers
        CALLED="${SAMPLE}.called.cns"
        cnvkit.py call "$CNS" -o "$CALLED"
        
        # Step 6: Export to SEG format
        cnvkit.py export seg "$CALLED" -o "${SAMPLE}.seg"
        
        # Step 7: Generate plots
        cnvkit.py scatter "$CNR" -s "$CALLED" -o "${SAMPLE}_scatter.pdf"
        cnvkit.py diagram "$CNR" -s "$CALLED" -o "${SAMPLE}_diagram.pdf"
        
        echo "  [$(timestamp)] Completed: $SAMPLE"
    done
done

echo ""
echo "################################################################################"
echo "# CNV calling completed for cohort: $COHORT_NAME"
echo "# Output directory: $COHORT_CNV_DIR"
echo "# Reference: $REFERENCE_CNN"
echo "# Next step: Run 03_aggregate_segments.py"
echo "################################################################################"