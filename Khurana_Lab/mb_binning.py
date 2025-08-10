import pandas as pd
import warnings


ACROCENTRIC_LIST = [13, 14, 15, 21, 22, "13", "14", "15", "21", "22"]
warnings.filterwarnings('ignore')
BIN_STEP = pow(10, 6)
BUILD_FILE = r".\chr_arm_annotations_modified_xy_hg19.xlsx" ## hg19

arms_borders_df = pd.read_excel(BUILD_FILE)
chromosome_arm_list = sorted([str(i) + "p" for i in range(1, 23)] +
                        [str(i) + "q" for i in range(1, 23)], key=lambda x: int(x[:-1])) + ["Xp", "Xq"]
arms_borders_dict = dict(zip(chromosome_arm_list,
                             zip(arms_borders_df["Start"].to_list(), arms_borders_df["End"].to_list())))
chromosomes_list = [str(i) for i in range(1, 23)] + ["X"]
chromosome_borders_list = [(arms_borders_dict[str(chrom) + "p"][0], arms_borders_dict[str(chrom) + "q"][1])
                           for chrom in chromosomes_list]

chromosomes_borders_dict = dict(zip(chromosomes_list, chromosome_borders_list))
segments_df = pd.read_csv(r".\segments_of_all_samples.csv")


cells_set = list(segments_df["Sample_ID"].unique())

bins_list = []
for chromosome in chromosomes_list:
    start, end = chromosomes_borders_dict[chromosome]
    for i in range(start, end + 1, BIN_STEP):
        bins_list.append(str(chromosome) + ":" + str(i) + "-" + str(i + BIN_STEP))
data = [["Sample_ID"] + bins_list]
for i, cell in enumerate(cells_set):
    print(i, cell)
    segments_curr_cell = segments_df[segments_df["Sample_ID"] == cell]
    sample_bins = [cell]
    for bin in bins_list:
        start, end = int(bin.split(":")[1].split("-")[0]), int(bin.split(":")[1].split("-")[1])
        segments_curr_cell_arm = segments_curr_cell[(segments_curr_cell["Chromosome"] == bin.split(":")[0]) &
                                                    (segments_curr_cell["Start"] < end) &
                                                    (segments_curr_cell["End"] > start)].reset_index(drop=True)
        segments_curr_cell_arm["segment_length"] = segments_curr_cell_arm["End"] - segments_curr_cell_arm["Start"]
        segments_curr_cell_arm["weighted_seg"] = segments_curr_cell_arm["segment_length"] * segments_curr_cell_arm["SegmentMean"]
        arm_cn = segments_curr_cell_arm["weighted_seg"].sum() / segments_curr_cell_arm["segment_length"].sum()
        sample_bins.append(round(arm_cn, 3))
    data.append(sample_bins)
x = pd.DataFrame(columns=data[0], data=data[1:])
print(x.to_string())
x.fillna("NA").to_csv("mb_bins_normalized.csv", index=False)