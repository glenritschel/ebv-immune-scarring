import scanpy as sc
import pandas as pd

adata = sc.read_h5ad("data/processed/GSE195452.bcells_scvi.h5ad")

# Full score table
score_cols = ['score_ISG', 'score_HLA', 'score_NFkB', 'score_scarring']
full_table = adata.obs.groupby('leiden_scvi_bcell')[score_cols].mean().round(4)
full_table['n_cells'] = adata.obs['leiden_scvi_bcell'].value_counts().sort_index()
print(full_table)
