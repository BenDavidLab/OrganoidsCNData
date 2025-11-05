import os
import pandas as pd

seg_files_list = os.listdir(r".\seg_files")
seg_files_list = [pd.read_csv(r".\seg_files" +
                              "/" + x, sep="\t") for x in seg_files_list]
all_segments = pd.concat(seg_files_list)
all_segments = all_segments[all_segments["chrom"].isin([str(i) for i in range(1, 23)] + ["X", "Y"])]
all_segments.to_csv("segments_of_all_samples.seg", sep="\t", index=False)
print(all_segments)

all_segments = all_segments.rename(columns={"ID": "Sample_ID", "chrom": "Chromosome", "loc.start": "Start",
                                            "loc.end": "End", "seg.mean": "SegmentMean"})

all_segments.to_csv("segments_of_all_samples.csv", index=False)

