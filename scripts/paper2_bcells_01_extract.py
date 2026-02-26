#!/usr/bin/env python3
import os
import numpy as np
import scanpy as sc

IN_H5AD = "data/processed/GSE195452.scvi_scored_umap.h5ad"
OUT_H5AD = "data/processed/GSE195452.bcells.h5ad"

# marker genes
MARKERS = ["CD19", "MS4A1", "CD79A"]
EXCLUDE_PLASMA = ["MZB1", "XBP1", "SDC1"]  # optional

def ensure_log1p_layer(adata):
    if "log1p" in adata.layers:
        return
    if "counts" in adata.layers:
        X = adata.layers["counts"]
        adata.layers["log1p"] = X.copy()
        adata.layers["log1p"] = adata.layers["log1p"].astype(np.float32)
        if hasattr(adata.layers["log1p"], "data"):  # sparse
            adata.layers["log1p"].data = np.log1p(adata.layers["log1p"].data)
        else:
            adata.layers["log1p"] = np.log1p(adata.layers["log1p"])
        return
    # fallback: assume adata.X is already log-like
    adata.layers["log1p"] = adata.X.copy()

def expr_gt0(adata, gene, layer="log1p"):
    if gene not in adata.var_names:
        return None
    x = adata[:, gene].layers[layer]
    if hasattr(x, "toarray"):
        x = x.toarray()
    x = np.asarray(x).reshape(-1)
    return x > 0

def main():
    print("Loading:", IN_H5AD)
    adata = sc.read_h5ad(IN_H5AD)
    ensure_log1p_layer(adata)

    present = [g for g in MARKERS if g in adata.var_names]
    if len(present) == 0:
        raise SystemExit(f"No B cell markers found in var_names: {MARKERS}")

    # Tier 1: broad capture (OR)
    masks = [expr_gt0(adata, g) for g in present]
    m_or = masks[0].copy()
    for m in masks[1:]:
        m_or |= m

    # Tier 2: higher specificity (CD79A AND (CD19 OR MS4A1)) if available
    m_cd79a = expr_gt0(adata, "CD79A")
    m_cd19 = expr_gt0(adata, "CD19")
    m_ms4a1 = expr_gt0(adata, "MS4A1")
    if m_cd79a is not None and (m_cd19 is not None or m_ms4a1 is not None):
        m_or2 = np.zeros(adata.n_obs, dtype=bool)
        if m_cd19 is not None: m_or2 |= m_cd19
        if m_ms4a1 is not None: m_or2 |= m_ms4a1
        m_tier2 = m_cd79a & m_or2
    else:
        m_tier2 = m_or  # fallback

    adata.obs["is_bcell_tier1"] = m_or.astype(int)
    adata.obs["is_bcell_tier2"] = m_tier2.astype(int)

    # Use Tier 2 by default
    b = adata[m_tier2].copy()
    print("B cells selected (tier2):", b.n_obs)

    # OPTIONAL: exclude plasma-like cells (to keep “B cell compartment” clean)
    # Uncomment if you want
    # ex_present = [g for g in EXCLUDE_PLASMA if g in b.var_names]
    # if ex_present:
    #     ex_mask = np.zeros(b.n_obs, dtype=bool)
    #     for g in ex_present:
    #         ex_mask |= expr_gt0(b, g)
    #     b = b[~ex_mask].copy()
    #     print("After plasma exclusion:", b.n_obs)

    os.makedirs(os.path.dirname(OUT_H5AD), exist_ok=True)
    b.write(OUT_H5AD)
    print("Wrote:", OUT_H5AD)

if __name__ == "__main__":
    main()

