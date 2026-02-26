# Pseudobulk DE: Cluster 0 vs Cluster 1 (Cluster 2 excluded as too small)
# Aggregate by GSM (sample identifier)

import scanpy as sc
import pandas as pd
import numpy as np
from scipy.sparse import issparse

adata = sc.read_h5ad("data/processed/GSE195452.bcells_scvi.h5ad")

# Keep only Clusters 0 and 1
adata_pb = adata[adata.obs['leiden_scvi_bcell'].isin(['0', '1'])].copy()

# Pseudobulk: sum raw counts per sample per cluster
records = []
for (gsm, cluster), grp in adata_pb.obs.groupby(['gsm', 'leiden_scvi_bcell']):
    if len(grp) < 5:  # skip samples with too few cells
        continue
    idx = grp.index
    X = adata_pb[idx].X
    if issparse(X):
        X = X.toarray()
    counts = X.sum(axis=0)
    records.append({
        'gsm': gsm,
        'cluster': cluster,
        'n_cells': len(grp),
        **dict(zip(adata_pb.var_names, counts))
    })

pb = pd.DataFrame(records)
print(f"Pseudobulk samples: {len(pb)}")
print(pb.groupby('cluster').size())

# Save for DESeq2 input
pb.to_csv("results/paper2_bcells_scvi/pseudobulk_cluster0_vs_1.csv", index=False)
print("Saved pseudobulk matrix.")
