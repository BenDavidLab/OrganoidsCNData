import pandas as pd

expression_matrix = pd.read_csv(r".\rsem.gene.estimated.counts_combined.csv")
expression_matrix.iloc[:, 0] = expression_matrix.iloc[:, 0].apply(lambda x: x.split(".")[0])

expression_matrix = expression_matrix.rename(columns={"Unnamed: 0": "ensembl_id"})
expression_matrix = expression_matrix.groupby(["ensembl_id"]).sum().reset_index()
expression_matrix = expression_matrix.set_index(expression_matrix.columns[0])
expression_matrix = expression_matrix.round(0).astype(int)

for sample in expression_matrix.columns.to_list():
    curr_df = expression_matrix.loc[:, [sample]]
    curr_df.to_csv(r".\per_sample/" + sample + ".counts", header=False, sep="\t")
