import scanpy as sc

adata = sc.read_h5ad("data/processed/GSE195452.bcells_scvi.h5ad")

print("Resolution used?")
print(adata.uns.get("leiden_scvi_bcell"))


print(adata.uns.keys())

