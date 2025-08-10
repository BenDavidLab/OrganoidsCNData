import pandas as pd

all_segments = pd.read_csv(r".\bladder_columbia_msk_2018_segments.seg", sep="\t")
print(all_segments.head().to_string())
all_segments = all_segments.drop(columns=["num.mark"])
all_segments = all_segments.rename(columns={"ID": "Sample_ID", "loc.start": "Start", "chrom": "Chromosome",
                                            "loc.end": "End", "seg.mean": "SegmentMean"})

all_segments.to_csv("segments_of_all_samples.csv", index=False)

