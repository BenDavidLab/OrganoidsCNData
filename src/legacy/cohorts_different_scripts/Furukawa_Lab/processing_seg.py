import os
import pandas as pd

all_segs = []
for f in os.listdir(r"."):
    if f.split(".")[1] == "seg":
        all_segs.append(pd.read_csv(f, sep="\t"))
# Sample_ID,Chromosome,Start,End,SegmentMean
all_segs = pd.concat(all_segs).drop(columns=["C", "num.mark"])
all_segs = all_segs.rename(columns={"ID": "Sample_ID", "chrom": "Chromosome", "loc.start": "Start", "loc.end": "End",
                            "seg.mean": "SegmentMean"})
all_segs.to_csv("segments_of_all_samples.csv", index=False)
