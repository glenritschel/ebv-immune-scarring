import scanpy as sc
import pandas as pd

adata = sc.read_h5ad("data/processed/GSE195452.bcells_scvi.h5ad")

# Cluster sizes
print("=== CLUSTER SIZES ===")
print(adata.obs["leiden_scvi_bcell"].value_counts().sort_index())

# Check for signature scores
score_cols = [c for c in adata.obs.columns if any(x in c.lower() 
              for x in ['score', 'isg', 'hla', 'nfkb', 'scarr'])]
print("\n=== SIGNATURE SCORE COLUMNS ===")
print(score_cols)

if score_cols:
    print("\n=== MEAN SCORES BY CLUSTER ===")
    print(adata.obs.groupby('leiden_scvi_bcell')[score_cols].mean().round(4))

# Check all obs columns
print("\n=== ALL OBS COLUMNS ===")
print(list(adata.obs.columns))
