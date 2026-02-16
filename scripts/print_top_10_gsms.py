import pandas as pd
import scanpy as sc

adata = sc.read_h5ad("data/processed/GSE195452.scvi_scored_umap.h5ad")

ct = (
    adata.obs
    .groupby(["gsm","leiden_scvi"])
    .size()
    .unstack(fill_value=0)
)

cluster7 = ct["7"]
print("Number of GSMs with cluster 7 cells:", (cluster7 > 0).sum())
print("Top 10 GSMs by cluster7 count:")
print(cluster7.sort_values(ascending=False).head(10))

