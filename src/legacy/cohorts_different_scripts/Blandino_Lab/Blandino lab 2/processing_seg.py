import os

import pandas as pd

input_path = r".\CNS_files"
df_samples_list = []
for f in os.listdir(input_path):
    curr_seg = pd.read_csv(input_path + "/" + f, sep="\t")
    curr_seg["Sample_ID"] = [f.split(".")[0][:-2] for i in range(len(curr_seg))]
    curr_seg = curr_seg.loc[:, ["Sample_ID", "chromosome", "start", "end", "log2"]]
    curr_seg["temp_sort_chromosome"] = curr_seg["chromosome"].apply(lambda x: int(x) if x not in ["X", "Y"] else 23 if x == "X" else 24)
    curr_seg = curr_seg.sort_values(["temp_sort_chromosome", "start", "end"]).reset_index(drop=True)
    curr_seg = curr_seg.drop(columns=["temp_sort_chromosome"])
    curr_seg = curr_seg.rename(columns={"chromosome": "Chromosome", "start": "Start", "end": "End", "log2": "SegmentMean"})
    df_samples_list.append(curr_seg)

all_samples_df = pd.concat(df_samples_list).reset_index(drop=True)
# Capping: -3 <= SegmentMean
all_samples_df["SegmentMean"] = [max(-3, x) for x in all_samples_df["SegmentMean"].to_list()]
all_samples_df.to_csv("segments_of_all_samples.csv", index=False)
