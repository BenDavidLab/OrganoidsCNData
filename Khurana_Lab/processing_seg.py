import pandas as pd

segments1 = pd.read_csv(r".\prad_mskcc_cheny1_organoids_2014_segments.seg", sep="\t")
segments2 = pd.read_csv(r".\prad_organoids_msk_2022_segments.seg", sep="\t")
all_segments = pd.concat([segments1, segments2]).reset_index(drop=True)
all_segments = all_segments.drop(columns=["num.mark"])
all_segments = all_segments.rename(columns={"ID": "Sample_ID", "chrom":"Chromosome", "loc.start": "Start",
                                            "loc.end": "End", "seg.mean": "SegmentMean"})

all_segments.to_csv("segments_of_all_samples.csv", index=False)

