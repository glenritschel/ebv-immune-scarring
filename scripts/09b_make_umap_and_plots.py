#!/usr/bin/env python3
"""
Create UMAP from scVI embedding (X_scVI) if missing, save new h5ad,
and write standard UMAP figures:
- umap_leiden_scvi.png
- umap_scarred_cluster7.png
- umap_signature_scores.png (all obs cols that start with score_)
"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")

import matplotlib.pyplot as plt
import scanpy as sc


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--in-h5ad", required=True)
    ap.add_argument("--out-h5ad", required=True)
    ap.add_argument("--figdir", required=True)
    ap.add_argument("--cluster", default="7")
    ap.add_argument("--groupby", default="leiden_scvi")
    ap.add_argument("--use-rep", default="X_scVI")
    ap.add_argument("--neighbors-k", type=int, default=15)
    args = ap.parse_args()

    adata = sc.read_h5ad(args.in_h5ad)
    if args.use_rep not in adata.obsm:
        raise SystemExit(f"ERROR: adata.obsm['{args.use_rep}'] not found")
    if args.groupby not in adata.obs.columns:
        raise SystemExit(f"ERROR: adata.obs['{args.groupby}'] not found")

    adata.obs[args.groupby] = adata.obs[args.groupby].astype(str)
    adata.obs["is_cluster7"] = (adata.obs[args.groupby] == str(args.cluster))

    # UMAP missing?
    has_umap = "X_umap" in adata.obsm
    if not has_umap:
        # Ensure neighbors exist / recompute neighbors from scVI embedding
        sc.pp.neighbors(adata, use_rep=args.use_rep, n_neighbors=args.neighbors_k)
        sc.tl.umap(adata)
        print("Computed UMAP from", args.use_rep)
    else:
        print("UMAP already present in adata.obsm['X_umap']")

    out_h5ad = Path(args.out_h5ad)
    out_h5ad.parent.mkdir(parents=True, exist_ok=True)
    adata.write_h5ad(out_h5ad)
    print("Wrote:", out_h5ad)

    figdir = Path(args.figdir)
    figdir.mkdir(parents=True, exist_ok=True)

    # UMAP colored by cluster labels
    out1 = figdir / "umap_leiden_scvi.png"
    sc.pl.umap(adata, color=args.groupby, show=False)
    plt.tight_layout()
    plt.savefig(out1, dpi=200)
    plt.close()
    print("Wrote:", out1)

    # UMAP colored by boolean cluster7
    out2 = figdir / "umap_scarred_cluster7.png"
    sc.pl.umap(adata, color="is_cluster7", show=False)
    plt.tight_layout()
    plt.savefig(out2, dpi=200)
    plt.close()
    print("Wrote:", out2)

    # Signature overlays: all obs columns starting with "score_"
    score_cols = [c for c in adata.obs.columns if c.startswith("score_")]
    if len(score_cols) > 0:
        out3 = figdir / "umap_signature_scores.png"
        sc.pl.umap(adata, color=score_cols, ncols=4, show=False)
        plt.tight_layout()
        plt.savefig(out3, dpi=200)
        plt.close()
        print("Wrote:", out3)
    else:
        print("No score_* columns found; skipping umap_signature_scores.png")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())

