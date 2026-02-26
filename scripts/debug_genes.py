import pandas as pd

df = pd.read_csv(
    "results/paper2_bcells_scvi/tables/bcells_cluster0_vs_others_wilcoxon.tsv",
    sep="\t"
)

genes_of_interest = [
    "CD74", "HLA-DRA", "HLA-DRB1", "HLA-DPA1", "HLA-DPB1",
    "HLA-DQB1", "CD83", "IFI27", "NFKBIA", "RELA"
]

print(df[df["names"].isin(genes_of_interest)].sort_values("logfoldchanges", ascending=False))

