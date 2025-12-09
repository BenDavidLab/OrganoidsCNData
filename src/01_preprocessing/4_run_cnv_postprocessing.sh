#!/bin/bash

# -----------------------------
# Directories (adjust if needed)
# -----------------------------
OUTPUT_DIR="/workspace/final_results"
DATASET="4205TNT" #edit later 
SEG_FILE="/workspace/SEGs/segments_of_all_samples_4205TNT.seg"

# Output files
MB_BINS="$OUTPUT_DIR/1mb_bins_normalized_4205TNT.csv"
ARM_BINS="$OUTPUT_DIR/arm_level_normalized_4205TNT.csv"

# -----------------------------
# Step 1: 1Mb binning
# -----------------------------
echo "▶ Running 1Mb binning..."
python3 4-2_mb_binning.py \
    --input $SEG_FILE \
    --output $MB_BINS \
    --bin_size 1000000 \
    --normalize

if [ $? -ne 0 ]; then
    echo "❌ Error in mb_binning.py"
    exit 1
fi
echo "✅ 1Mb bins saved to $MB_BINS"

# -----------------------------
# Step 2: Arm-level normalization
# -----------------------------
echo "▶ Running arm-level binning..."
python3 4-3_arm_level_binning.py \
    --input $SEG_FILE \
    --output $ARM_BINS

if [ $? -ne 0 ]; then
    echo "❌ Error in arm_level_binning.py"
    exit 1
fi
echo "✅ Arm-level bins saved to $ARM_BINS"

# -----------------------------
# Done
# -----------------------------
echo "🎉 CNV postprocessing finished!"
echo "   - 1Mb bins: $MB_BINS"
echo "   - Arm-level bins: $ARM_BINS"
