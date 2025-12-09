#!/usr/bin/env python3

"""
================================================================================
# WES CNVkit Pipeline - multi-sample, fully automated
# Author: Haia K. (adapted & enhanced)
# Date: 2025-08-18

Description:
    Combines all individual .seg files from CNVkit into a single aggregated
    file for downstream analysis. Processes all .seg files in the input
    directory, adds sample IDs, and merges into one unified segments file.

Inputs:
    - Individual .seg files from script 02
      Location: results/cnv/{cohort}/*.seg
    
Outputs:
    - results/segments/segments_of_all_samples_{cohort}.seg
      Unified tab-delimited segments file
    - results/segments/segments_of_all_samples_{cohort}.csv
      Same data in CSV format

Output Columns:
    - ID: Sample identifier
    - chrom: Chromosome
    - loc.start: Segment start position
    - loc.end: Segment end position
    - num.mark: Number of markers in segment
    - seg.mean: Log2 copy ratio of segment

Dependencies:
    - Python 3
    - pandas

Usage:
    # From repository root:
    python3 src/preprocessing/03_aggregate_segments.py --cohort COHORT_NAME
    
    # Or specify custom paths:
    python3 src/preprocessing/03_aggregate_segments.py \\
        --input results/cnv/cohort1 \\
        --output results/segments/cohort1_segments.seg
================================================================================
"""

import os
import sys
import pandas as pd
import argparse
from pathlib import Path

# ============================================================================
# AUTO-DETECT REPOSITORY STRUCTURE
# ============================================================================

def get_repo_root():
    """Auto-detect repository root from script location"""
    script_path = Path(__file__).resolve()
    # Go up from src/preprocessing/ to repo root
    repo_root = script_path.parent.parent.parent
    return repo_root

# ============================================================================
# ARGUMENT PARSING
# ============================================================================

def parse_arguments():
    parser = argparse.ArgumentParser(
        description="Aggregate individual .seg files into unified dataset",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
    # Using cohort name (auto-detects paths):
    python3 03_aggregate_segments.py --cohort COHORT_NAME
    
    # Using custom paths:
    python3 03_aggregate_segments.py \\
        --input results/cnv/cohort1 \\
        --output results/segments/cohort1_segments.seg
        """
    )
    
    parser.add_argument(
        "--cohort",
        type=str,
        help="Cohort name (auto-detects input/output paths)"
    )
    
    parser.add_argument(
        "--input",
        type=str,
        help="Input directory containing .seg files (overrides --cohort)"
    )
    
    parser.add_argument(
        "--output",
        type=str,
        help="Output file path (overrides --cohort)"
    )
    
    args = parser.parse_args()
    
    # Validate arguments
    if not args.cohort and not args.input:
        parser.error("Either --cohort or --input must be specified")
    
    return args

# ============================================================================
# MAIN PROCESSING
# ============================================================================

def aggregate_segments(input_dir, output_file, cohort_name=None):
    """
    Aggregate all .seg files in input directory into single file
    
    Args:
        input_dir: Directory containing .seg files
        output_file: Path for output aggregated file
        cohort_name: Optional cohort name for display
    """
    
    print(f"{'='*80}")
    print(f"Aggregating Segments{f' - Cohort: {cohort_name}' if cohort_name else ''}")
    print(f"{'='*80}")
    print(f"Input directory: {input_dir}")
    print(f"Output file: {output_file}")
    print("")
    
    # Check input directory exists
    if not os.path.exists(input_dir):
        print(f"ERROR: Input directory not found: {input_dir}")
        sys.exit(1)
    
    # Find all .seg files
    seg_files = [f for f in os.listdir(input_dir) if f.endswith(".seg")]
    
    if not seg_files:
        print(f"ERROR: No .seg files found in {input_dir}")
        sys.exit(1)
    
    print(f"Found {len(seg_files)} .seg files")
    print("")
    
    # Process each .seg file
    df_samples_list = []
    for i, filename in enumerate(seg_files, 1):
        print(f"Processing {i}/{len(seg_files)}: {filename}")
        
        file_path = os.path.join(input_dir, filename)
        
        try:
            # Read segment file
            curr_seg = pd.read_csv(file_path, sep="\t")
            
            # Add sample ID from filename (remove extension)
            sample_id = filename.replace(".seg", "")
            curr_seg["ID"] = sample_id
            
            df_samples_list.append(curr_seg)
            
        except Exception as e:
            print(f"  WARNING: Error processing {filename}: {e}")
            continue
    
    if not df_samples_list:
        print("ERROR: No valid .seg files were processed")
        sys.exit(1)
    
    # Merge all samples into one dataframe
    print("")
    print("Merging all samples...")
    all_samples_df = pd.concat(df_samples_list, ignore_index=True)
    
    # Define expected columns and types
    wanted_cols = {
        "ID": "object",
        "chrom": "object",
        "loc.start": "Int64",
        "loc.end": "Int64",
        "num.mark": "float64",
        "seg.mean": "float64"
    }
    
    # Clean whitespace from string columns
    all_samples_df = all_samples_df.apply(
        lambda col: col.str.strip() if col.dtype == "object" else col
    )
    
    # Keep only wanted columns that exist
    available_cols = [c for c in wanted_cols if c in all_samples_df.columns]
    all_samples_df = all_samples_df[available_cols]
    
    # Convert to proper types
    for col, dtype in wanted_cols.items():
        if col in all_samples_df.columns:
            try:
                all_samples_df[col] = all_samples_df[col].astype(dtype)
            except:
                pass  # Keep original type if conversion fails
    
    # Create output directory if needed
    output_dir = os.path.dirname(output_file)
    if output_dir and not os.path.exists(output_dir):
        os.makedirs(output_dir)
        print(f"Created output directory: {output_dir}")
    
    # Save as tab-delimited .seg file
    all_samples_df.to_csv(output_file, sep="\t", index=False)
    
    # Also save as CSV
    csv_file = output_file.replace(".seg", ".csv")
    all_samples_df.to_csv(csv_file, index=False)
    
    # Report statistics
    print("")
    print(f"{'='*80}")
    print(f"Aggregation complete!")
    print(f"{'='*80}")
    print(f"Total samples: {len(df_samples_list)}")
    print(f"Total segments: {len(all_samples_df):,}")
    print(f"Chromosomes: {sorted(all_samples_df['chrom'].unique())}")
    print("")
    print(f"Output files:")
    print(f"  SEG: {output_file}")
    print(f"  CSV: {csv_file}")
    print(f"{'='*80}")

# ============================================================================
# MAIN EXECUTION
# ============================================================================

def main():
    args = parse_arguments()
    repo_root = get_repo_root()
    
    # Determine input and output paths
    if args.input:
        # User provided explicit input path
        input_dir = args.input
    else:
        # Auto-detect from cohort name
        input_dir = os.path.join(repo_root, "results", "cnv", args.cohort)
    
    if args.output:
        # User provided explicit output path
        output_file = args.output
    else:
        # Auto-detect from cohort name
        segments_dir = os.path.join(repo_root, "results", "segments")
        os.makedirs(segments_dir, exist_ok=True)
        output_file = os.path.join(
            segments_dir, 
            f"segments_of_all_samples_{args.cohort}.seg"
        )
    
    # Run aggregation
    aggregate_segments(input_dir, output_file, args.cohort)

if __name__ == "__main__":
    main()
