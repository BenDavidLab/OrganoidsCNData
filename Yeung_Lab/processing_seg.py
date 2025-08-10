import pandas as pd

input_path = ".\WES_seg_files_Yeung_Lab"
df_annotation = pd.read_csv(".\WES_seg_files_Yeung_Lab\Data_Info_Table.tsv", sep="\t")
df_annotation = df_annotation.sort_values(["Patient_Number", "Sample_Short_Form"]).reset_index(drop=True)
df_samples_list = []
for f in df_annotation["File Name"].to_list():
    with open(input_path + "/" + f) as pre_file:
        counter = 0
        lines = pre_file.readlines()
        for line in lines:
            if line[:3] == "@RG":
                break
            counter += 1
        with open(input_path + "/" + f[:-4] + "_no_header.txt", "w") as post_file:
            post_file.writelines(lines[counter + 1:])
    curr_seg = pd.read_csv(input_path + "/" + f[:-4] + "_no_header.txt", sep="\t")
    curr_seg["Sample_ID"] = [f.split(".")[0] for i in range(len(curr_seg))]
    curr_seg = curr_seg.loc[:, ["Sample_ID", "CONTIG", "START", "END", "LOG2_COPY_RATIO_POSTERIOR_50"]]
    curr_seg["temp_sort_chromosome"] = curr_seg["CONTIG"].apply(lambda x: int(x) if x not in ["X", "Y"] else 23 if x == "X" else 24)
    curr_seg = curr_seg.sort_values(["temp_sort_chromosome", "START", "END"]).reset_index(drop=True)
    curr_seg = curr_seg.drop(columns=["temp_sort_chromosome"])
    curr_seg = curr_seg.rename(columns={"CONTIG": "Chromosome", "START": "Start", "END": "End", "LOG2_COPY_RATIO_POSTERIOR_50": "SegmentMean"})
    df_samples_list.append(curr_seg)

all_samples_df = pd.concat(df_samples_list).reset_index(drop=True)
all_samples_df.to_csv("segments_of_all_samples.csv", index=False)
