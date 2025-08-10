import numpy as np
import pandas as pd

ACCROCENTRIC_CHROMOSOMES = [13, 14, 15, 21, 22, "13", "14", "15", "21", "22"]
seg_name_to_path = "segments_of_all_samples"
gene_matrix_path = r".\M_vals_Engel_etal_2022.xlsx"
BUILD_FILE = r"..\chr_arm_annotations_modified_xy_hg38.xlsx" #hg38

genes_location_matrix = pd.read_csv(r".\Genes_Location_HG38.csv")
genes_matrix = pd.read_excel(gene_matrix_path).fillna(0)

arms_borders_df = pd.read_excel(BUILD_FILE)
chromosome_arm_list = sorted([str(i) + "p" for i in range(1, 23)] +
                        [str(i) + "q" for i in range(1, 23)], key=lambda x: int(x[:-1])) + ["Xp", "Xq"]
arms_borders_dict = dict(zip(chromosome_arm_list,
                             zip(arms_borders_df["Start"].to_list(), arms_borders_df["End"].to_list())))
# gene_df = genes_matrix.reset_index().rename(columns={"index": "Gene"}).iloc[:,0:4]
gene_df = genes_matrix.rename(columns={"index": "Gene", "gene": "Gene"})
gene_df = genes_location_matrix.merge(gene_df, on="Gene")
count_gene_chromosomes = gene_df.groupby(["Gene"])["Chromosome"].nunique().reset_index().rename(columns={"Chromosome": "Count"})
gene_df = gene_df[~gene_df["Gene"].isin(count_gene_chromosomes[count_gene_chromosomes["Count"]>1]["Gene"].to_list())]

cell_list = []
print(arms_borders_dict)
for cell in gene_df.columns.to_list()[6:]:
    curr_cell_df = gene_df.loc[:, gene_df.columns.to_list()[:6] + [cell]]
    curr_cell_df[cell] = curr_cell_df.apply(lambda row: np.nan if ((row["Chromosome"] in ACCROCENTRIC_CHROMOSOMES) and (row["Start"]<arms_borders_dict[str(row["Chromosome"])+"p"][1])) else row[cell], axis=1)
    if len(curr_cell_df) == 0:
        continue
    curr_cell_df["Sample_ID"] = [cell for i in range(len(curr_cell_df))]
    curr_cell_df = curr_cell_df.loc[:, ["Sample_ID"] + curr_cell_df.columns.to_list()[:-1]]
    chromosome_arm_list = []
    for chromosome in gene_df["Chromosome"].unique().tolist():
        for arm in ["p", "q"]:
            curr_df = curr_cell_df[(gene_df["Chromosome"] == chromosome) & (curr_cell_df["Arm"] == arm)].copy()
            if len(curr_df) == 0:
                continue
            arm_start, arm_end = arms_borders_dict[str(chromosome) + arm]
            gaps = []
            genes_intervals = list(zip(curr_df["Start"].to_list(), curr_df["End"].to_list()))
            if genes_intervals[0][0] != 1:
                gaps.append((cell, "", chromosome, arm, arm_start, genes_intervals[0][0] - 1, genes_intervals[0][0] - 2, np.nan))
            for i in range(len(genes_intervals) - 1):
                if genes_intervals[i+1][0] - genes_intervals[i][1] > 1:
                    gaps.append((cell, "", chromosome, arm, genes_intervals[i][1] + 1, genes_intervals[i+1][0] - 1, genes_intervals[i+1][0] - 1 - (genes_intervals[i][1] + 1), np.nan))
            if genes_intervals[-1][0] != arm_end:
                gaps.append((cell, "", chromosome, arm, genes_intervals[-1][0] + 1, arm_end, arm_end - (genes_intervals[-1][0] + 1), np.nan))
            imputed_curr_df = pd.concat([curr_df, pd.DataFrame(columns=curr_df.columns, data=gaps)]).sort_values(["Start"]).reset_index(drop=True)
            imputed_curr_df[cell] = imputed_curr_df[cell].fillna((imputed_curr_df[cell].shift() + imputed_curr_df[cell].shift(-1)) / 2)
            imputed_curr_df.loc[0, cell] = imputed_curr_df.loc[1, cell]
            imputed_curr_df.loc[len(imputed_curr_df) - 1, cell] = imputed_curr_df.loc[len(imputed_curr_df) - 2, cell]
            print(cell, chromosome, arm)
            print(imputed_curr_df.isna().values.flatten().sum())
            chromosome_arm_list.append(imputed_curr_df)
    cell_df = pd.concat(chromosome_arm_list).reset_index(drop=True)
    cell_df = cell_df.drop(columns=["Gene", "Arm", "Gene_Length"])
    cell_df = cell_df.rename(columns={cell: "SegmentMean"})
    cell_list.append(cell_df)
all_cells_segment = pd.concat(cell_list).reset_index(drop=True)
print(all_cells_segment)
all_cells_segment.to_csv(seg_name_to_path + ".csv", index=False, na_rep="NA")
