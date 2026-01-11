#!/usr/bin/env python3
"""
Plot cluster-level signature heatmap (publication-friendly; avoids raw gene heatmaps).

- Aggregates mean signature score per cluster
- Z-scores signatures across clusters
- Saves PDF + TSV table

Assumes signature scores live in .obs columns (prefix sig_ by default).
"""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import List

import numpy as np
import pandas as pd

import matplotlib
matplotlib.use("Agg")

import matplotlib.pyplot as plt
import scanpy as sc


def autodetect_signatures(adata, prefix: str) -> List[str]:
    cols = [c for c in adata.obs.columns if c.startswith(prefix)]
    return sorted(cols)


def zscore_cols(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    for c in out.columns:
        x = out[c].astype(float).values
        mu = np.nanmean(x)
        sd = np.nanstd(x)
        if sd == 0 or np.isnan(sd):
            out[c] = 0.0
        else:
            out[c] = (x - mu) / sd
    return out


def plot_heatmap(matrix: pd.DataFrame, out_pdf: Path, title: str = "") -> None:
    # simple matplotlib heatmap (no seaborn)
    fig = plt.figure(figsize=(max(8, 0.6 * matrix.shape[1] + 3), max(6, 0.35 * matrix.shape[0] + 2)))
    ax = fig.add_subplot(111)

    im = ax.imshow(matrix.values, aspect="auto")
    ax.set_xticks(np.arange(matrix.shape[1]))
    ax.set_xticklabels(list(matrix.columns), rotation=45, ha="right")
    ax.set_yticks(np.arange(matrix.shape[0]))
    ax.set_yticklabels([str(i) for i in matrix.index])

    if title:
        ax.set_title(title)

    fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    fig.tight_layout()
    out_pdf.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(out_pdf, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--h5ad", required=True, help="Input scored .h5ad (e.g., *.scvi_scored.h5ad)")
    ap.add_argument("--out-pdf", default="figures/Fig3_cluster_signature_heatmap.pdf", help="Output heatmap PDF")
    ap.add_argument("--out-tsv", default="results/tables/cluster_signature_means.tsv", help="Output table TSV")
    ap.add_argument("--cluster-key", default="leiden", help="Cluster key in .obs (default: leiden)")
    ap.add_argument("--sig-prefix", default="sig_", help="Signature prefix in .obs (default: sig_)")
    ap.add_argument("--sigs", nargs="*", default=[], help="Explicit signature columns to include (optional)")
    ap.add_argument("--zscore", action="store_true", help="Z-score signatures across clusters (recommended)")
    ap.add_argument("--title", default="Cluster-level signature means", help="Plot title")
    args = ap.parse_args()

    adata = sc.read_h5ad(args.h5ad)

    if args.cluster_key not in adata.obs.columns:
        raise SystemExit(f"Missing cluster key in .obs: {args.cluster_key}")

    sigs = args.sigs if args.sigs else autodetect_signatures(adata, args.sig_prefix)
    if not sigs:
        raise SystemExit(f"No signature columns found in .obs with prefix '{args.sig_prefix}' (or pass --sigs).")

    df = adata.obs[[args.cluster_key] + sigs].copy()
    for c in sigs:
        df[c] = pd.to_numeric(df[c], errors="coerce")

    means = df.groupby(args.cluster_key, observed=False)[sigs].mean()
    means = means.sort_index(key=lambda s: s.astype(str))

    out_tsv = Path(args.out_tsv)
    out_tsv.parent.mkdir(parents=True, exist_ok=True)
    means.to_csv(out_tsv, sep="\t")

    matrix = zscore_cols(means) if args.zscore else means

    out_pdf = Path(args.out_pdf)
    plot_heatmap(matrix, out_pdf=out_pdf, title=args.title)

    print(f"[OK] Wrote: {out_pdf}")
    print(f"[OK] Wrote: {out_tsv}")


if __name__ == "__main__":
    main()

