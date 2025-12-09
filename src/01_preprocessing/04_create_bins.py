#!/usr/bin/env python3

"""
================================================================================
# WES CNVkit Pipeline - multi-sample, fully automated
# Author: Haia K. (adapted & enhanced)
# Date: 2025-08-18

Description:
    Final step of pre-processing pipeline. Takes aggregated segments file and
    generates two normalized output formats:
    1. 1MB bins: Genome-wide fixed bins at 1 megabase resolution
    2. Arm-level: 46 chromosome arm values (1p, 1q, ..., 22p, 22q, Xp, Xq)
    
    Both outputs use weighted averaging of segment log2 values based on
    overlap length, producing normalized copy number calls ready for
    discordance analysis.

Inputs:
    - Aggregated segments file from script 03
      Location: results/segments/segments_of_all_samples_{cohort}.seg
    - Chromosome arm annotations (hg19/GRCh37)
      Location: data/reference_files/chr_arm_annotations_hg19.xlsx

Outputs:
    - results/Cohort_data/{cohort}/{cohort}_mb_bins_normalized.csv
      ~3000 columns (1MB bins across genome)
    - results/Cohort_data/{cohort}/{cohort}_arm_level_normalized.csv
      46 columns (chromosome arms)

Output Format:
    - First column: Sample ID
    - Remaining columns: Normalized log2 copy number ratios
    - Missing values: "NA" for low-quality or missing regions

Dependencies:
    - Python 3
    - pandas
    - openpyxl (for reading Excel arm annotations)

Usage:
    # From repository root:
    python3 src/preprocessing/04_create_bins.py --cohort COHORT_NAME
    
    # Or with custom paths:
    python3 src/preprocessing/04_create_bins.py \\
        --input results/segments/cohort1_segments.seg \\
        --output results/Cohort_data/cohort1/ \\
        --arms data/reference_files/chr_arm_annotations_hg19.xlsx

Note: This script combines functionality from:
    - 4_run_cnv_postprocessing.sh (master script)
    - 4-2_mb_binning.py (1MB binning)
    - 4-3_arm_level_binning.py (arm-level)
================================================================================
"""

import os
import sys
import pandas as pd
import warnings
import argparse
from pathlib import Path

warnings.filterwarnings('ignore')

# ============================================================================
# AUTO-DETECT REPOSITORY STRUCTURE
# ============================================================================

def get_repo_root():
    """Auto-detect repository root from script location"""
    script_path = Path(__file__).resolve()
    repo_root = script_path.parent.parent.parent
    return repo_root

# ============================================================================
# ARGUMENT PARSING
# ============================================================================

def parse_arguments():
    parser = argparse.ArgumentParser(
        description="Generate normalized 1MB bins and arm-level CNV data",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
    # Using cohort name (auto-detects paths):
    python3 04_create_bins.py --cohort COHORT_NAME
    
    # Using custom paths:
    python3 04_create_bins.py \\
        --input results/segments/cohort1_segments.seg \\
        --output results/Cohort_data/cohort1/ \\
        --arms data/reference_files/chr_arm_annotations_hg19.xlsx
        """
    )
    
    parser.add_argument(
        "--cohort",
        type=str,
        help="Cohort name (auto-detects paths)"
    )
    
    parser.add_argument(
        "--input",
        type=str,
        help="Input aggregated segments file (overrides --cohort)"
    )
    
    parser.add_argument(
        "--output",
        type=str,
        help="Output directory (overrides --cohort)"
    )
    
    parser.add_argument(
        "--arms",
        type=str,
        help="Chromosome arm annotations file (Excel or CSV)"
    )
    
    parser.add_argument(
        "--bin-size",
        type=int,
        default=1000000,
        help="Bin size in bp (default: 1000000 = 1MB)"
    )
    
    args = parser.parse_args()
    
    if not args.cohort and not args.input:
        parser.error("Either --cohort or --input must be specified")
    
    return args

# ============================================================================
# LOAD CHROMOSOME ARM ANNOTATIONS
# ============================================================================

def load_arm_annotations(arms_file):
    """
    Load chromosome arm start/end positions
    
    Returns:
        arms_borders_dict: {arm_name: (start, end), ...}
        chromosomes_borders_dict: {chrom: (start, end), ...}
    """
    
    print(f"Loading chromosome arm annotations: {arms_file}")
    
    # Read file (Excel or CSV)
    if arms_file.endswith('.xlsx') or arms_file.endswith('.xls'):
        arms_df = pd.read_excel(arms_file)
    else:
        arms_df = pd.read_csv(arms_file)
    
    # Expected columns: Arm, Start, End (or similar)
    if 'Arm' not in arms_df.columns:
        # Try to find arm column
        possible_names = ['arm', 'Chromosome', 'Chr']
        for col in possible_names:
            if col in arms_df.columns:
                arms_df = arms_df.rename(columns={col: 'Arm'})
                break
    
    # Create arm list
    chromosome_arm_list = sorted(
        [str(i) + "p" for i in range(1, 23)] + [str(i) + "q" for i in range(1, 23)],
        key=lambda x: int(x[:-1])
    ) + ["Xp", "Xq"]
    
    # Create dictionary: arm -> (start, end)
    arms_borders_dict = dict(
        zip(
            chromosome_arm_list,
            zip(arms_df["Start"].to_list(), arms_df["End"].to_list())
        )
    )
    
    # Create chromosome borders
    chromosomes_list = [str(i) for i in range(1, 23)] + ["X"]
    chromosome_borders_list = [
        (arms_borders_dict[str(chrom) + "p"][0], arms_borders_dict[str(chrom) + "q"][1])
        for chrom in chromosomes_list
    ]
    chromosomes_borders_dict = dict(zip(chromosomes_list, chromosome_borders_list))
    
    print(f"  Loaded {len(arms_borders_dict)} chromosome arms")
    return arms_borders_dict, chromosomes_borders_dict, chromosome_arm_list

# ============================================================================
# 1MB BINNING
# ============================================================================

def generate_1mb_bins(segments_df, chromosomes_borders_dict, bin_size, output_file):
    """
    Generate 1MB fixed bins across the genome
    
    Args:
        segments_df: DataFrame with columns [ID, chrom, loc.start, loc.end, seg.mean]
        chromosomes_borders_dict: {chrom: (start, end)}
        bin_size: Bin size in bp (default 1000000)
        output_file: Output CSV file path
    """
    
    print(f"\n{'='*80}")
    print(f"Generating 1MB bins (bin size: {bin_size:,} bp)")
    print(f"{'='*80}")
    
    chromosomes_list = [str(i) for i in range(1, 23)] + ["X"]
    
    # Generate bin list
    bins_list = []
    for chromosome in chromosomes_list:
        start, end = chromosomes_borders_dict[chromosome]
        for i in range(start, end + 1, bin_size):
            bins_list.append(f"{chromosome}:{i}-{i + bin_size}")
    
    print(f"Total bins: {len(bins_list):,}")
    
    # Get unique samples
    cells_set = list(segments_df["ID"].unique())
    print(f"Total samples: {len(cells_set)}")
    print("")
    
    # Initialize output data
    data = [["ID"] + bins_list]
    
    # Process each sample
    for idx, cell in enumerate(cells_set, 1):
        print(f"Processing {idx}/{len(cells_set)}: {cell}")
        
        segments_curr_cell = segments_df[segments_df["ID"] == cell]
        sample_bins = [cell]
        
        for bin_name in bins_list:
            # Parse bin coordinates
            chrom, coords = bin_name.split(":")
            start, end = map(int, coords.split("-"))
            
            # Find overlapping segments
            segs = segments_curr_cell[
                (segments_curr_cell["chrom"] == chrom) &
                (segments_curr_cell["loc.start"] < end) &
                (segments_curr_cell["loc.end"] > start)
            ].reset_index(drop=True)
            
            if len(segs) == 0:
                sample_bins.append("NA")
                continue
            
            # Weighted average based on segment length
            segs["segment_length"] = segs["loc.end"] - segs["loc.start"]
            segs["weighted_seg"] = segs["segment_length"] * segs["seg.mean"]
            bin_cn = segs["weighted_seg"].sum() / segs["segment_length"].sum()
            
            sample_bins.append(round(bin_cn, 3))
        
        data.append(sample_bins)
    
    # Save to CSV
    df = pd.DataFrame(columns=data[0], data=data[1:])
    df.fillna("NA").to_csv(output_file, index=False)
    
    print("")
    print(f"✅ 1MB bins saved: {output_file}")
    print(f"   Dimensions: {len(df)} samples × {len(bins_list)} bins")
    print(f"{'='*80}")

# ============================================================================
# ARM-LEVEL BINNING
# ============================================================================

def generate_arm_level(segments_df, arms_borders_dict, chromosome_arm_list, output_file):
    """
    Generate arm-level normalized values (46 arms)
    
    Args:
        segments_df: DataFrame with columns [ID, chrom, loc.start, loc.end, seg.mean]
        arms_borders_dict: {arm_name: (start, end)}
        chromosome_arm_list: List of 46 arm names
        output_file: Output CSV file path
    """
    
    print(f"\n{'='*80}")
    print(f"Generating arm-level values (46 arms)")
    print(f"{'='*80}")
    
    # Get unique samples
    cells_set = list(segments_df["ID"].unique())
    print(f"Total samples: {len(cells_set)}")
    print("")
    
    # Initialize output data
    data = [["ID"] + chromosome_arm_list]
    
    # Process each sample
    for idx, cell in enumerate(cells_set, 1):
        print(f"Processing {idx}/{len(cells_set)}: {cell}")
        
        segments_curr_cell = segments_df[segments_df["ID"] == cell]
        sample_values = [cell]
        
        for arm in chromosome_arm_list:
            # Parse arm (e.g., "1p" -> chrom="1", arm="p")
            chrom = arm[:-1]
            start, end = arms_borders_dict[arm]
            
            # Find overlapping segments
            segs = segments_curr_cell[
                (segments_curr_cell["chrom"] == chrom) &
                (segments_curr_cell["loc.start"] < end) &
                (segments_curr_cell["loc.end"] > start)
            ].reset_index(drop=True)
            
            if len(segs) == 0:
                sample_values.append("NA")
                continue
            
            # Weighted average based on segment length
            segs["segment_length"] = segs["loc.end"] - segs["loc.start"]
            segs["weighted_seg"] = segs["segment_length"] * segs["seg.mean"]
            arm_cn = segs["weighted_seg"].sum() / segs["segment_length"].sum()
            
            sample_values.append(round(arm_cn, 3))
        
        data.append(sample_values)
    
    # Save to CSV
    df = pd.DataFrame(columns=data[0], data=data[1:])
    df.fillna("NA").to_csv(output_file, index=False)
    
    print("")
    print(f"✅ Arm-level data saved: {output_file}")
    print(f"   Dimensions: {len(df)} samples × 46 arms")
    print(f"{'='*80}")

# ============================================================================
# MAIN EXECUTION
# ============================================================================

def main():
    args = parse_arguments()
    repo_root = get_repo_root()
    
    # Determine paths
    if args.input:
        input_file = args.input
    else:
        input_file = os.path.join(
            repo_root, "results", "segments", 
            f"segments_of_all_samples_{args.cohort}.seg"
        )
    
    if args.output:
        output_dir = args.output
    else:
        output_dir = os.path.join(
            repo_root, "results", "Cohort_data", args.cohort
        )
    
    if args.arms:
        arms_file = args.arms
    else:
        arms_file = os.path.join(
            repo_root, "data", "reference_files", 
            "chr_arm_annotations_hg19.xlsx"
        )
    
    # Create output directory
    os.makedirs(output_dir, exist_ok=True)
    
    # Define output files
    cohort_name = args.cohort if args.cohort else "output"
    mb_bins_file = os.path.join(output_dir, f"{cohort_name}_mb_bins_normalized.csv")
    arm_level_file = os.path.join(output_dir, f"{cohort_name}_arm_level_normalized.csv")
    
    print(f"{'='*80}")
    print(f"CNV Binning Pipeline - Cohort: {cohort_name}")
    print(f"{'='*80}")
    print(f"Input segments: {input_file}")
    print(f"Output directory: {output_dir}")
    print(f"Arm annotations: {arms_file}")
    print(f"{'='*80}")
    
    # Validate input file
    if not os.path.exists(input_file):
        print(f"ERROR: Input file not found: {input_file}")
        sys.exit(1)
    
    # Validate arms file
    if not os.path.exists(arms_file):
        print(f"ERROR: Chromosome arm annotations not found: {arms_file}")
        print(f"Please place arm annotations file in: {os.path.dirname(arms_file)}/")
        sys.exit(1)
    
    # Load segments
    print(f"\nLoading segments...")
    try:
        segments_df = pd.read_csv(input_file, sep="\t")
    except:
        # Try space-separated
        segments_df = pd.read_csv(input_file, sep=r"\s+", engine="python")
    
    print(f"  Loaded {len(segments_df):,} segments from {segments_df['ID'].nunique()} samples")
    
    # Load chromosome arm annotations
    arms_borders_dict, chromosomes_borders_dict, chromosome_arm_list = load_arm_annotations(arms_file)
    
    # Generate 1MB bins
    generate_1mb_bins(segments_df, chromosomes_borders_dict, args.bin_size, mb_bins_file)
    
    # Generate arm-level values
    generate_arm_level(segments_df, arms_borders_dict, chromosome_arm_list, arm_level_file)
    
    print(f"\n{'='*80}")
    print(f"CNV binning completed successfully!")
    print(f"{'='*80}")
    print(f"Output files:")
    print(f"  1MB bins: {mb_bins_file}")
    print(f"  Arm-level: {arm_level_file}")
    print(f"\nNext step: Run post-processing scripts in src/postprocessing/")
    print(f"{'='*80}")

if __name__ == "__main__":
    main()
