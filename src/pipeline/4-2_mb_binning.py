import pandas as pd
import warnings
import argparse

# Suppress warnings
warnings.filterwarnings('ignore')

# ---------------------------
# Parse arguments
# ---------------------------
parser = argparse.ArgumentParser(description="Generate 1Mb normalized bins from SEG files")
parser.add_argument("--input", required=True, help="Input aggregated SEG/CSV file")
parser.add_argument("--output", required=True, help="Output CSV file for 1Mb bins")
parser.add_argument("--bin_size", type=int, default=1000000, help="Bin size (default: 1Mb)")
parser.add_argument("--normalize", action="store_true", help="Apply normalization (currently placeholder)")
args = parser.parse_args()

BIN_STEP = args.bin_size
BUILD_FILE = r"/workspace/chr_arm_annotations_modified_xy_hg19.xlsx"  # hg19

# ---------------------------
# Load chromosome arm borders
# ---------------------------
arms_borders_df = pd.read_excel(BUILD_FILE)
chromosome_arm_list = sorted(
    [str(i) + "p" for i in range(1, 23)] + [str(i) + "q" for i in range(1, 23)],
    key=lambda x: int(x[:-1])
) + ["Xp", "Xq"]

arms_borders_dict = dict(
    zip(
        chromosome_arm_list,
        zip(arms_borders_df["Start"].to_list(), arms_borders_df["End"].to_list())
    )
)

chromosomes_list = [str(i) for i in range(1, 23)] + ["X"]
chromosome_borders_list = [
    (arms_borders_dict[str(chrom) + "p"][0], arms_borders_dict[str(chrom) + "q"][1])
    for chrom in chromosomes_list
]

chromosomes_borders_dict = dict(zip(chromosomes_list, chromosome_borders_list))

# ---------------------------
# Load segments
# ---------------------------
segments_df = pd.read_csv(args.input, sep="\t")
cells_set = list(segments_df["ID"].unique())

# ---------------------------
# Generate bins
# ---------------------------
bins_list = []
for chromosome in chromosomes_list:
    start, end = chromosomes_borders_dict[chromosome]
    for i in range(start, end + 1, BIN_STEP):
        bins_list.append(f"{chromosome}:{i}-{i + BIN_STEP}")

data = [["ID"] + bins_list]

# ---------------------------
# Process each sample
# ---------------------------
for i, cell in enumerate(cells_set):
    print(f"Processing {i+1}/{len(cells_set)}: {cell}")
    segments_curr_cell = segments_df[segments_df["ID"] == cell]
    sample_bins = [cell]
    for bin in bins_list:
        start, end = map(int, bin.split(":")[1].split("-"))
        chrom = bin.split(":")[0]
        segs = segments_curr_cell[
            (segments_curr_cell["chrom"] == chrom) &
            (segments_curr_cell["loc.start"] < end) &
            (segments_curr_cell["loc.end"] > start)
        ].reset_index(drop=True)

        if len(segs) == 0:
            sample_bins.append("NA")
            continue

        segs["segment_length"] = segs["loc.end"] - segs["loc.start"]
        segs["weighted_seg"] = segs["segment_length"] * segs["seg.mean"]
        arm_cn = segs["weighted_seg"].sum() / segs["segment_length"].sum()
        sample_bins.append(round(arm_cn, 3))

    data.append(sample_bins)

# ---------------------------
# Save output
# ---------------------------
x = pd.DataFrame(columns=data[0], data=data[1:])
x.fillna("NA").to_csv(args.output, index=False)

print(f"✅ 1Mb normalized bins saved to {args.output}")