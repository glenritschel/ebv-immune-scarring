import scanpy as sc
adata = sc.read_h5ad("data/processed/GSE195452.bcells_scvi.h5ad")

# Check what clustering columns exist
print("Available obs columns:")
print([col for col in adata.obs.columns if 'leiden' in col.lower()])

# Check signature scores
score_cols = [col for col in adata.obs.columns if 'score' in col.lower()]
print("\nScore columns:", score_cols)

# Print signature scores by cluster
if score_cols:
    print("\nMean scores by cluster:")
    print(adata.obs.groupby('leiden_scvi_bcell')[score_cols].mean())
