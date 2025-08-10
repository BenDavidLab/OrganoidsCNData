import pandas as pd

seg1 = pd.read_excel(r".\cell_stem_cell_2020_leung.xlsx")
seg2 = pd.read_csv(r".\gut_2018_leung_seg.csv")
all_segments = pd.concat([seg1, seg2]).loc[:,["sample", "chrom", "start", "end", "cnlr.median"]].reset_index(drop=True)
all_segments["chrom"] = all_segments["chrom"].apply(lambda x: x[3:])
print(all_segments.head().to_string())

all_segments = all_segments.rename(columns={"sample": "Sample_ID", "start": "Start",
                                             "end": "End", "cnlr.median": "SegmentMean", "chrom": "Chromosome"})

all_segments.to_csv("segments_of_all_samples.csv", index=False)

