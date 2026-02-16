import scanpy as sc
adata = sc.read_h5ad("data/processed/GSE195452.scvi_scored_umap.h5ad")

print("OBS columns:")
print(list(adata.obs.columns))

print("\nUnique gsm count:", adata.obs["gsm"].nunique())
print("\nSample of obs rows:")
print(adata.obs.head())

