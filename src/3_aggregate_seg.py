import os
import pandas as pd

input_path = "/workspace/SEGs"
dataset = "4205TNT"
df_samples_list = []
for f in os.listdir(input_path):
  if f.endswith(".seg"):  # only process .seg files
    curr_seg = pd.read_csv(os.path.join(input_path, f), sep="\t")
    
    # Add sample ID (remove extension and last 2 chars if needed)
    curr_seg["Sample_ID"] = [f.split(".")[0][:21] for i in range(len(curr_seg))]

    # Keep only relevant columns
    #curr_seg = curr_seg.loc[:, ["ID", "chrom", "loc.start", "loc.end", "seg.mean"]]
    
    # Sort chromosomes numerically
    #curr_seg["temp_sort_chromosome"] = curr_seg["chrom"].apply(lambda x: int(x) if x not in ["X", "Y"] else 23 if x == "X" else 24)
    #curr_seg = curr_seg.sort_values(["temp_sort_chromosome", "loc.start", "loc.end"]).reset_index(drop=True)
    #curr_seg = curr_seg.drop(columns=["temp_sort_chromosome"])
    
    # Rename columns for IGV style
    #curr_seg = curr_seg.rename(columns={"chrom": "Chromosome", "loc.start": "Start", "loc.end": "End", "seg.mean": "SegmentMean"})
    df_samples_list.append(curr_seg)
  
  
# Merge all samples into one
all_samples_df = pd.concat(df_samples_list).reset_index(drop=True)

# Cap extreme values at -3
#all_samples_df["seg.mean"] = [
#    max(-3, x) for x in all_samples_df["seg.mean"].to_list()
#]

# Define the wanted columns and types
wanted_cols = {
    "ID": "object",
    "chrom": "object",
    "loc.start": "Int64",
    "loc.end": "Int64",
    "num.mark": "float64",
    "seg.mean": "float64"
}

# Clean whitespace from all string/object columns
all_samples_df = all_samples_df.apply(lambda col: col.str.strip() if col.dtype == "object" else col)

# Keep only the wanted columns that exist in the file
all_samples_df = all_samples_df[[c for c in wanted_cols if c in all_samples_df.columns]]

# Convert types
for col, dtype in wanted_cols.items():
    if col in all_samples_df.columns:
        all_samples_df[col] = all_samples_df[col].astype(dtype, errors="ignore")

# Save output
all_samples_df.to_csv("segments_of_all_samples_{dataset}.csv", index=False)

# Save as tab-delimited .seg file
output_file = os.path.join(input_path, f"segments_of_all_samples_{dataset}.seg")
all_samples_df.to_csv(output_file, sep="\t", index=False)

print(f"✅ Unified SEG file written: {output_file}")