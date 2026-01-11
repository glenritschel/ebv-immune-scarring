#!/usr/bin/env python3
"""
Plot signature overlays on UMAP for scored scVI AnnData.

Assumes:
- UMAP exists or computed from X_scVI
- signature scores are in .obs columns (e.g., 'sig_ISG', 'sig_AntigenPresentation', etc.)

You can either pass explicit columns, or let it auto-detect by prefix.
"""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import List

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


def autodetect_signatures(adata, prefix: str) -> List[str]:
    cols = [c for c in adata.obs.columns if c.startswith(prefix)]
    return sorted(cols)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--h5ad", required=True, help="Input scored .h5ad (e.g., *.scvi_scored.h5ad)")
    ap.add_argument("--out", default="figures/Fig2_signature_overlays.pdf", help="Output PDF path")
    ap.add_argument("--sig-prefix", default="sig_", help="Auto-detect signatures by prefix (default sig_)")
    ap.add_argument("--sigs", nargs="*", default=[], help="Explicit signature columns to plot (overrides prefix)")
    ap.add_argument("--neighbors-k", type=int, default=15, help="Neighbors for UMAP if needed")
    args = ap.parse_args()

    adata = sc.read_h5ad(args.h5ad)
    ensure_umap(adata, neighbors_k=args.neighbors_k)

    sigs = args.sigs if args.sigs else autodetect_signatures(adata, args.sig_prefix)
    if not sigs:
        raise SystemExit(
            f"No signature columns found. Provide --sigs or ensure .obs columns start with '{args.sig_prefix}'."
        )

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)

    sc.pl.umap(
        adata,
        color=sigs,
        ncols=3,
        wspace=0.4,
        show=False,
    )
    plt.savefig(out, bbox_inches="tight")
    plt.close("all")
    print(f"[OK] Wrote: {out}")


if __name__ == "__main__":
    main()

