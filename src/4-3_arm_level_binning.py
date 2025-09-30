import pandas as pd
import warnings
import argparse

warnings.filterwarnings('ignore')

# ---------------------------
# Parse arguments
# ---------------------------
parser = argparse.ArgumentParser(description="Generate arm-level normalized CN values from SEG files")
parser.add_argument("--input", required=True, help="Input aggregated SEG/CSV file")
parser.add_argument("--output", required=True, help="Output CSV file for arm-level values")
args = parser.parse_args()

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

# ---------------------------
# Load segments
# ---------------------------
#segments_df = pd.read_csv(args.input)
segments_df = pd.read_csv(args.input, sep="\t", engine="python")
if "ID" not in segments_df.columns:
    # fallback: try space separation
    segments_df = pd.read_csv(args.input, sep=r"\s+", engine="python")

# Standardize columns just in case
segments_df = segments_df.rename(columns={
    "Sample": "ID",
    "sample": "ID"
})
cells_set = list(segments_df["ID"].unique())

# ---------------------------
# Process arm-level values
# ---------------------------
data = [["ID"] + chromosome_arm_list]

for i, cell in enumerate(cells_set):
    print(f"Processing {i+1}/{len(cells_set)}: {cell}")
    segments_curr_cell = segments_df[segments_df["ID"] == cell]
    sample_values = [cell]

    for arm in chromosome_arm_list:
        chrom = arm[:-1]
        start, end = arms_borders_dict[arm]

        segs = segments_curr_cell[
            (segments_curr_cell["chrom"] == chrom) &
            (segments_curr_cell["loc.start"] < end) &
            (segments_curr_cell["loc.end"] > start)
        ].reset_index(drop=True)

        if len(segs) == 0:
            sample_values.append("NA")
            continue

        segs["segment_length"] = segs["loc.end"] - segs["loc.start"]
        segs["weighted_seg"] = segs["segment_length"] * segs["seg.mean"]
        arm_cn = segs["weighted_seg"].sum() / segs["segment_length"].sum()
        sample_values.append(round(arm_cn, 3))

    data.append(sample_values)

# ---------------------------
# Save output
# ---------------------------
x = pd.DataFrame(columns=data[0], data=data[1:])
x.fillna("NA").to_csv(args.output, index=False)

print(f"✅ Arm-level normalized file saved to {args.output}")
