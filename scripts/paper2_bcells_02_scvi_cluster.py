#!/usr/bin/env python3
import os
import numpy as np
import scanpy as sc

import scvi

IN_H5AD  = "data/processed/GSE195452.bcells.h5ad"
OUT_H5AD = "data/processed/GSE195452.bcells_scvi.h5ad"
OUTDIR   = "results/paper2_bcells_scvi"

BATCH_KEY_CANDIDATES = ["sample_id", "gsm"]  # use whichever exists

def pick_batch_key(adata):
    for k in BATCH_KEY_CANDIDATES:
        if k in adata.obs.columns:
            return k
    raise SystemExit(f"No batch key found. Tried: {BATCH_KEY_CANDIDATES}")

def ensure_counts_layer(adata):
    if "counts" in adata.layers:
        return
    # If X looks like counts, keep it as counts. Otherwise, you need to regenerate counts earlier in pipeline.
    adata.layers["counts"] = adata.X.copy()

def main():
    os.makedirs(OUTDIR, exist_ok=True)
    print("Loading:", IN_H5AD)
    adata = sc.read_h5ad(IN_H5AD)
    ensure_counts_layer(adata)

    batch_key = pick_batch_key(adata)
    print("Using batch_key:", batch_key)

    scvi.model.SCVI.setup_anndata(adata, layer="counts", batch_key=batch_key)
    model = scvi.model.SCVI(adata, n_latent=30)
    model.train(max_epochs=200)

    adata.obsm["X_scVI_bcell"] = model.get_latent_representation()

    # neighbors/umap/leiden on latent
    sc.pp.neighbors(adata, use_rep="X_scVI_bcell", n_neighbors=15)
    sc.tl.umap(adata)
    sc.tl.leiden(adata, key_added="leiden_scvi_bcell", resolution=0.6)

    adata.write(OUT_H5AD)
    print("Wrote:", OUT_H5AD)

if __name__ == "__main__":
    main()

