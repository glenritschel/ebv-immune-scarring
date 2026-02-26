import scanpy as sc
import pandas as pd
import numpy as np
from scipy.sparse import issparse
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats

# Load
adata = sc.read_h5ad("data/processed/GSE195452.bcells_scvi.h5ad")

# Keep Clusters 0 and 1 only (Cluster 2 too small)
adata_pb = adata[adata.obs['leiden_scvi_bcell'].isin(['0', '1'])].copy()

# Pseudobulk: sum raw counts per (gsm, cluster)
records = []
for (gsm, cluster), grp in adata_pb.obs.groupby(['gsm', 'leiden_scvi_bcell']):
    if len(grp) < 5:
        continue
    X = adata_pb[grp.index].X
    if issparse(X):
        X = X.toarray()
    counts = X.sum(axis=0).astype(int)
    records.append({'gsm': gsm, 'cluster': cluster, 'n_cells': len(grp),
                    **dict(zip(adata_pb.var_names, counts))})

pb = pd.DataFrame(records).set_index('gsm')
print(f"Pseudobulk samples: {len(pb)}")
print(pb.groupby('cluster').size())

# Separate counts and metadata
meta = pb[['cluster', 'n_cells']].copy()
counts = pb.drop(columns=['cluster', 'n_cells']).T  # genes x samples

# Filter lowly expressed genes (sum >= 10 across all samples)
counts = counts[counts.sum(axis=1) >= 10]
print(f"Genes after filtering: {len(counts)}")

# PyDESeq2
dds = DeseqDataSet(
    counts=counts.T,          # samples x genes
    metadata=meta,
    design_factors="cluster"
)
dds.deseq2()

stats = DeseqStats(dds, contrast=["cluster", "0", "1"])
stats.summary()

results = stats.results_df.sort_values("padj")
print("\nTop 20 DE genes:")
print(results[['baseMean','log2FoldChange','pvalue','padj']].head(20))

# Save
results.to_csv(
    "results/paper2_bcells_scvi/pseudobulk_pydeseq2_cluster0_vs_1.csv"
)
print("\nSaved.")
