import warnings

import pandas as pd

warnings.filterwarnings('ignore')
# ACROCENTRIC_LIST = [13, 14, 15, 21, 22, "13", "14", "15", "21", "22"]
BUILD_FILE = r"..\..\chr_arm_annotations_modified_xy_hg19.xlsx" #hg19

arms_borders_df = pd.read_excel(BUILD_FILE)
# arms_borders_df = arms_borders_df[~arms_borders_df["Arm"].isin(["13p", "14p", "15p", "21p", "22p"])]
chromosome_arm_list = sorted([str(i) + "p" for i in range(1, 23)] +
                        [str(i) + "q" for i in range(1, 23)], key=lambda x: int(x[:-1])) + ["Xp", "Xq"]
arms_borders_dict = dict(zip(chromosome_arm_list,
                             zip(arms_borders_df["Start"].to_list(), arms_borders_df["End"].to_list())))
segments_df = pd.read_csv(r".\segments_of_all_samples.csv")

cells_set = list(segments_df["Sample_ID"].unique())

data = [["Sample_ID"] + chromosome_arm_list]
for cell in cells_set:
    print(cell)
    segments_curr_cell = segments_df[segments_df["Sample_ID"] == cell]
    sample_karyotype = [cell]
    for chromosome_arm in chromosome_arm_list:
        start, end = arms_borders_dict[chromosome_arm][0], arms_borders_dict[chromosome_arm][1]
        segments_curr_cell_arm = segments_curr_cell[(segments_curr_cell["Chromosome"] == chromosome_arm[:-1]) &
                                                    (segments_curr_cell["Start"] < end) &
                                                    (segments_curr_cell["End"] > start)].reset_index(drop=True)
        if chromosome_arm[-1] == "q":
            if len(segments_curr_cell_arm) != 0:
                segments_curr_cell_arm.loc[0, "Start"] = max(start, segments_curr_cell_arm.loc[0, "Start"])
        else:
            if len(segments_curr_cell_arm) != 0:
                segments_curr_cell_arm.loc[len(segments_curr_cell_arm) - 1, "End"] =\
                    min(end, segments_curr_cell_arm.loc[len(segments_curr_cell_arm) - 1, "End"])
        segments_curr_cell_arm["segment_length"] = segments_curr_cell_arm["End"] - segments_curr_cell_arm["Start"]
        segments_curr_cell_arm["weighted_seg"] = segments_curr_cell_arm["segment_length"] * segments_curr_cell_arm["SegmentMean"]
        arm_cn = segments_curr_cell_arm["weighted_seg"].sum() / segments_curr_cell_arm["segment_length"].sum()
        sample_karyotype.append(round(arm_cn, 3))
    data.append(sample_karyotype)

x = pd.DataFrame(columns=data[0], data=data[1:])
print(x.to_string())
x.fillna("NA").to_csv("arm_level_normalized.csv", index=False)
