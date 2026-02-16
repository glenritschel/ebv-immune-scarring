#!/usr/bin/env python3
"""
Immune scarring "heatmap-like" matrixplot across clusters (groupby=leiden_scvi).
Uses Agg backend to avoid Qt/Wayland issues.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import scanpy as sc


DEFAULT_GENES = [
    "ISG15","MX1","OAS1","IFIT1","IFIT3","IFI44L","IFI27",
    "HLA-DRA","HLA-DRB1","HLA-DPA1","HLA-DPB1","HLA-DQA1","HLA-DQB1",
    "CD74","CD83","NFKBIA","RELA","NFKB1","NFKB2","REL",
]


def ensure_log1p_layer(adata: sc.AnnData) -> str:
    if "log1p" in adata.layers:
        return "log1p"
    if "counts" not in adata.layers:
        return None
    X = adata.layers["counts"]
    adata.layers["log1p"] = X.copy()
    adata.layers["log1p"] = adata.layers["log1p"].astype(np.float32)
    if hasattr(adata.layers["log1p"], "data"):
        adata.layers["log1p"].data = np.log1p(adata.layers["log1p"].data)
    else:
        adata.layers["log1p"] = np.log1p(adata.layers["log1p"])
    return "log1p"


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--in-h5ad", required=True)
    ap.add_argument("--out-png", required=True)
    ap.add_argument("--groupby", default="leiden_scvi")
    ap.add_argument("--layer", default="log1p")
    ap.add_argument("--standard-scale", default="var", choices=["var", "obs", "none"])
    ap.add_argument("--width", type=float, default=12)
    ap.add_argument("--height", type=float, default=6)
    ap.add_argument("--genes", nargs="*", default=DEFAULT_GENES)
    args = ap.parse_args()

    adata = sc.read_h5ad(args.in_h5ad)
    if args.groupby not in adata.obs.columns:
        raise SystemExit(f"ERROR: adata.obs['{args.groupby}'] not found")

    layer = args.layer
    if layer == "log1p":
        layer = ensure_log1p_layer(adata)

    if layer is None:
        print("WARN: no log1p/counts; using adata.X for plotting")

    genes_present = [g for g in args.genes if g in adata.var_names]
    print(f"Immune-scarring genes present: {len(genes_present)}/{len(args.genes)}")
    if len(genes_present) == 0:
        raise SystemExit("ERROR: none of the requested genes found in adata.var_names")

    out_png = Path(args.out_png)
    out_png.parent.mkdir(parents=True, exist_ok=True)

    standard_scale = None if args.standard_scale == "none" else args.standard_scale

    sc.pl.matrixplot(
        adata,
        var_names=genes_present,
        groupby=args.groupby,
        layer=layer,
        standard_scale=standard_scale,
        dendrogram=False,
        show=False,
    )
    plt.gcf().set_size_inches(args.width, args.height)
    plt.tight_layout()
    plt.savefig(out_png, dpi=200)
    plt.close()

    print("Wrote:", out_png)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

