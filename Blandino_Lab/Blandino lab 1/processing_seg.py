import os

import pandas as pd

input_path = r".\Blandino_Samples.seg"
all_samples_df = pd.read_csv(input_path, sep="\t")
all_samples_df = all_samples_df.rename(columns={"ID": "Sample_ID", "chrom": "Chromosome", "loc.start": "Start", "loc.end": "End", "seg.mean": "SegmentMean"})
all_samples_df = all_samples_df.drop(columns=["num.mark"])

all_samples_df.to_csv("segments_of_all_samples.csv", index=False)
