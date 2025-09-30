#!/bin/bash


################################################################################
# WES CNVkit Pipeline - multi-sample, fully automated
# Author: Haia K. (adapted & enhanced)
# Date: 2025-08-18
#
# Description:
#   This script automates CNV analysis for multiple WES samples (tumor, organoid, normal).
#   It generates target and antitarget bins, calculates coverage per BAM, builds a reference
#   from normal tissue, normalizes tumor/organoid samples, segments CNVs, calls absolute copy numbers,
#   exports results, and generates visualization plots.
#
# Inputs:
#   - BAM_DIR: directory containing all sorted and indexed BAM files
#   - COMMON_BED: merged/intersected BED file of all target regions
#   - REFERENCE_GENOME: indexed reference genome fasta file
#   - THREADS: number of CPU threads to use
#
# Outputs:
#   - .cnn files: per-sample target and antitarget coverage
#   - cnv_reference.cnn: reference built from normal sample(s)
#   - .cnr files: normalized log2 copy ratios per sample
#   - .cns files: segmented CNV regions per sample
#   - .called.cns files: CNV calls with absolute copy numbers
#   - .seg files: exported CNVs for downstream analysis
#   - scatter and diagram PDFs: visualization of CNVs per sample
#
# Dependencies:
#   - CNVkit
#   - Python 3 + pandas (optional for filtering extreme log2 values)
#   - samtools (for BAM indexing/validation)
#
# Usage:
#   bash 2_newCNV.sh
#   (optionally filter .cnr files with filter_cnr.py)
################################################################################

set -euo pipefail
shopt -s nullglob

# ======== User-configurable variables ========
BAM_DIR="/workspace/final_bams/EGAD00001004205"
COMMON_BED="/workspace/EGAD00001004205/S30409818_Regions.bed"
REFERENCE_GENOME="/workspace/RefGen/Homo_sapiens_assembly19.fasta"
THREADS=32

# Which BAMs should be treated as "normals" to build the pooled reference?
# Adjust this pattern to your true normals (e.g., "*NORMAL*_final.bam").
NORMAL_GLOB="*T.bam"

# Optional filtering of CNRs (your Python script)
FILTER_CNR=false
FILTER_SCRIPT="/workspace/filter_cnr.py"

# Segmentation method: cbs (needs R), hmm, haar, flasso
SEG_METHOD="${SEG_METHOD:-cbs}"

# ======== Derived file paths ========
TARGETS="${COMMON_BED%.bed}.split.bed"
ANTITARGETS="${COMMON_BED%.bed}.antitargets.bed"
REFERENCE_CNN="cnv_reference.cnn"

# ======== Helpers ========
timestamp() { date "+%Y-%m-%d %H:%M:%S"; }

echo "[$(timestamp)] CNVkit pipeline starting…"

# Optional: generate target/antitarget beds if not present
if [[ ! -s "$TARGETS" ]]; then
  echo "[$(timestamp)] Generating TARGETS from $COMMON_BED → $TARGETS"
  cnvkit.py target "$COMMON_BED" --split -o "$TARGETS"
fi
if [[ ! -s "$ANTITARGETS" ]]; then
  echo "[$(timestamp)] Generating ANTITARGETS → $ANTITARGETS"
  cnvkit.py antitarget "$COMMON_BED" -o "$ANTITARGETS"
fi

# ======== Build pooled reference from normals ========
if [[ ! -s "$REFERENCE_CNN" ]]; then
  declare -a NORMAL_TARGET_CNNS=()
  declare -a NORMAL_ANTITARGET_CNNS=()
  
  echo "[$(timestamp)] Building reference from normals in: $BAM_DIR/$NORMAL_GLOB"
  for NORMAL_BAM in "$BAM_DIR"/$NORMAL_GLOB; do
    [[ -f "$NORMAL_BAM" ]] || continue
    BASENAME=$(basename "$NORMAL_BAM" _final.bam)
    T_CNN="${BASENAME}.targetcoverage.cnn"
    A_CNN="${BASENAME}.antitargetcoverage.cnn"
  
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
    echo "ERROR: No normals matched pattern '$NORMAL_GLOB' in $BAM_DIR" >&2
    exit 1
  fi
  
  echo "[$(timestamp)] Creating CNVkit reference: $REFERENCE_CNN"
  cnvkit.py reference "${NORMAL_TARGET_CNNS[@]}" "${NORMAL_ANTITARGET_CNNS[@]}" \
    -f "$REFERENCE_GENOME" -o "$REFERENCE_CNN"
fi

# ======== Process matched pairs (ORG ↔ tissue) ========
echo "[$(timestamp)] Processing matched pairs…"
for ORG_BAM in "$BAM_DIR"/$NORMAL_GLOB; do
  [[ -f "$ORG_BAM" ]] || continue

  base=$(basename "$ORG_BAM")
  mate=${base/T/NT} #automate to accept t and T
  TISSUE_BAM="$BAM_DIR/$mate"

  if [[ ! -f "$TISSUE_BAM" ]]; then
    echo "[WARN $(timestamp)] No matching tissue BAM for $ORG_BAM"
    continue
  fi

  echo "[$(timestamp)] Pair found:"
  echo "  ORG:    $ORG_BAM"
  echo "  tissue: $TISSUE_BAM"
  echo "------------------------------------------------"
  echo "------------------------------------------------"
  echo "------------------------------------------------"
    
  # Process both samples in the pair
  for BAM in "$ORG_BAM" "$TISSUE_BAM"; do
    SAMPLE=$(basename "$BAM" .bam)
    if [[ ! -s "${SAMPLE}_diagram.pdf" ]]; then
    
      # Coverage
      T_CNN="${SAMPLE}.targetcoverage.cnn"
      A_CNN="${SAMPLE}.antitargetcoverage.cnn"
      if [[ ! -s "$T_CNN" ]]; then
        cnvkit.py coverage "$BAM" "$TARGETS" -o "$T_CNN" --processes "$THREADS"
      fi
      if [[ ! -s "$A_CNN" ]]; then
        cnvkit.py coverage "$BAM" "$ANTITARGETS" -o "$A_CNN" --processes "$THREADS"
      fi
  
      # Fix / normalize
      CNR="${SAMPLE}.cnr"
      cnvkit.py fix "$T_CNN" "$A_CNN" "$REFERENCE_CNN" -o "$CNR"
  
      # Optional filtering
      if [[ "$FILTER_CNR" == true ]]; then
        FILTERED="${SAMPLE}.filtered.cnr"
        python3 "$FILTER_SCRIPT" "$CNR" "$FILTERED"
        CNR="$FILTERED"
      fi
  
      # Segment (choose method)
      CNS="${SAMPLE}.cns"
      if [[ "$SEG_METHOD" == "hmm" ]]; then
        cnvkit.py segment "$CNR" -o "$CNS"
      else
        cnvkit.py segment "$CNR" --method "$SEG_METHOD" -o "$CNS"
      fi
  
      # Call absolute CN
      CALLED="${SAMPLE}.called.cns"
      cnvkit.py call "$CNS" -o "$CALLED"
  
      # Export + Plots
      cnvkit.py export seg "$CALLED" -o "${SAMPLE}.seg"
      cnvkit.py scatter "$CNR" -s "$CALLED" -o "${SAMPLE}_scatter.pdf"
      cnvkit.py diagram "$CNR" -s "$CALLED" -o "${SAMPLE}_diagram.pdf"
    fi
    echo "[$(timestamp)] Finished sample: $SAMPLE"
  done
done

echo "[$(timestamp)] All done. Reference: $REFERENCE_CNN"