import pandas as pd

all_segments = pd.read_csv(r".\samples_copy_number_segments.txt", sep="\t")

all_segments = all_segments.loc[:, ["sample_id", "chr", "start", "end", "log2"]]
all_segments = all_segments.rename(columns={"sample_id": "Sample_ID", "chr": "Chromosome", "start": "Start",
                                            "end": "End", "log2": "SegmentMean"})
all_segments = all_segments[all_segments["Chromosome"].isin(["chr" + str(i) for i in range(1, 23)] + ["chrX", "chrY"])]
all_segments["Chromosome"] = all_segments["Chromosome"].apply(lambda x: "chr23" if x == "chrX" else "chr24" if x == "chrY" else x)
all_segments["Chromosome"] = all_segments["Chromosome"].apply(lambda x: int(x[3:]))
all_segments = all_segments.sort_values(["Sample_ID", "Chromosome", "Start"])
all_segments["Chromosome"] = all_segments["Chromosome"].apply(lambda x: "X" if x == 23 else "Y" if x == 24 else str(x))
all_segments["Chromosome"] = all_segments["Chromosome"].astype("str")
all_segments.to_csv("segments_of_all_samples.csv", index=False)

