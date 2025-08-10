import pandas as pd

all_segments = pd.read_excel(r".\Segmented_Data_Sheet.xlsx")
print(all_segments.head().to_string())

all_segments = all_segments.drop(columns=["Patient ID", "Log2 Ratio"])
all_segments = all_segments.rename(columns={"Sample ID": "Sample_ID", "Start position": "Start",
                                            "End position": "End", "Log2 Ratio - CLONET adjusted": "SegmentMean"})

all_segments.to_csv("segments_of_all_samples.csv", index=False)

