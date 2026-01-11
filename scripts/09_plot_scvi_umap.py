#!/usr/bin/env python3
"""
Plot scVI UMAPs (cluster and metadata).

Loads an AnnData with:
- scVI latent in .obsm["X_scVI"] OR existing UMAP in .obsm["X_umap"]
- cluster labels in .obs["leiden"] (or user-provided)

If UMAP not present, computes neighbors/UMAP from X_scVI.
"""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Optional

import matplotlib
matplotlib.use("Agg")

import matplotlib.pyplot as plt
import scanpy as sc


def ensure_umap(adata, use_rep: str = "X_scVI", neighbors_k: int = 15) -> None:
    if "X_umap" in adata.obsm:
        return
    if use_rep not in adata.obsm:
        raise ValueError(f"Expected latent rep in adata.obsm['{use_rep}'] but it is missing.")
    sc.pp.neighbors(adata, use_rep=use_rep, n_neighbors=neighbors_k)
    sc.tl.umap(adata)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--h5ad", required=True, help="Input .h5ad (typically *.scvi_scored.h5ad or *.scvi.h5ad)")
    ap.add_argument("--out", default="figures/Fig1_scvi_umap.pdf", help="Output PDF path")
    ap.add_argument("--cluster-key", default="leiden", help="Cluster column in .obs (default: leiden)")
    ap.add_argument("--color2", default="", help="Optional second color column (e.g., gsm or dataset)")
    ap.add_argument("--neighbors-k", type=int, default=15, help="Neighbors for UMAP if needed")
    args = ap.parse_args()

    adata = sc.read_h5ad(args.h5ad)
    ensure_umap(adata, neighbors_k=args.neighbors_k)

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)

    cluster_key = args.cluster_key
    if cluster_key not in adata.obs.columns:
        # fallback: first obs col containing 'leiden'
        fallback = [c for c in adata.obs.columns if "leiden" in c.lower()]
        if fallback:
            cluster_key = fallback[0]
            print(f"[WARN] cluster key '{args.cluster_key}' not found; using '{cluster_key}'")
        else:
            raise KeyError(f"Could not find cluster key '{args.cluster_key}' and no leiden-like columns exist.")
    color_list = [cluster_key]

    if args.color2 and args.color2 in adata.obs.columns:
        color_list.append(args.color2)
    elif args.color2:
        print(f"[WARN] color2 '{args.color2}' not found in .obs; skipping")

    sc.pl.umap(
        adata,
        color=color_list,
        wspace=0.35,
        show=False,
    )
    plt.savefig(out, bbox_inches="tight")
    plt.close("all")
    print(f"[OK] Wrote: {out}")


if __name__ == "__main__":
    main()

