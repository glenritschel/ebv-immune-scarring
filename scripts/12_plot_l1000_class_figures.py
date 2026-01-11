#!/usr/bin/env python3
"""
12_plot_l1000_class_figures.py

Creates paper-ready figures from 08d outputs:
- Bar plot of top classes per signature
- Heatmap of class x signature (score_mean or n_hits)

Inputs:
  results/l1000cds2/l1000_class_enrichment.tsv
  results/l1000cds2/l1000_class_cross_signature.tsv  (optional)

Example:
  python scripts/12_plot_l1000_class_figures.py \
    --enrichment results/l1000cds2/l1000_class_enrichment.tsv \
    --out-dir figures \
    --metric score_mean \
    --topn 8
"""

from __future__ import annotations
import argparse
from pathlib import Path

import pandas as pd
import matplotlib.pyplot as plt


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--enrichment", required=True, help="l1000_class_enrichment.tsv from 08d")
    ap.add_argument("--out-dir", required=True, help="Output directory for PDFs")
    ap.add_argument("--metric", default="score_mean",
                    choices=["score_mean", "score_median", "score_max", "n_hits", "n_unique_compounds", "frac_hits"],
                    help="Metric to plot")
    ap.add_argument("--topn", type=int, default=8, help="Top N classes per signature for bar plots")
    ap.add_argument("--min-signatures", type=int, default=1,
                    help="Keep only classes appearing in >= this many signatures in heatmap")
    args = ap.parse_args()

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    df = pd.read_csv(args.enrichment, sep="\t", dtype=str).fillna("")

    # numeric coercion for metrics we might plot
    numeric_cols = {"score_mean","score_median","score_max","n_hits","n_unique_compounds","frac_hits"}
    for c in numeric_cols & set(df.columns):
        df[c] = pd.to_numeric(df[c], errors="coerce")

    metric = args.metric
    if metric not in df.columns:
        raise SystemExit(f"[FATAL] Metric '{metric}' not found in enrichment columns: {df.columns.tolist()}")

    # ---------------------------------------
    # Fig 4A: Top classes per signature (bars)
    # ---------------------------------------
    sigs = sorted(df["signature_name"].unique().tolist())

    # One PDF per signature (cleanest for paper assembly) + an all-in-one page
    for sig in sigs:
        sub = df[df["signature_name"] == sig].copy()
        sub = sub.sort_values(metric, ascending=False).head(args.topn)

        if len(sub) == 0:
            continue

        plt.figure()
        plt.barh(sub["class"], sub[metric])
        plt.gca().invert_yaxis()
        plt.xlabel(metric)
        plt.title(f"L1000 reversing classes: {sig} (top {args.topn} by {metric})")
        plt.tight_layout()

        out = out_dir / f"Fig4A_L1000_top_classes_{sig}.{metric}.pdf"
        plt.savefig(out)
        plt.close()
        print(f"[OK] Wrote: {out}")

    # Combined (signature facets in one plot) – small but useful
    # We’ll just concatenate and show signature:class labels
    top_all = (df.sort_values([ "signature_name", metric], ascending=[True, False])
                 .groupby("signature_name")
                 .head(args.topn)
                 .copy())
    if len(top_all):
        top_all["label"] = top_all["signature_name"] + " | " + top_all["class"]
        top_all = top_all.sort_values(metric, ascending=True)

        plt.figure(figsize=(8, max(4, 0.25 * len(top_all))))
        plt.barh(top_all["label"], top_all[metric])
        plt.xlabel(metric)
        plt.title(f"L1000 reversing classes (top {args.topn} per signature) by {metric}")
        plt.tight_layout()

        out = out_dir / f"Fig4A_L1000_top_classes_ALL.{metric}.pdf"
        plt.savefig(out)
        plt.close()
        print(f"[OK] Wrote: {out}")

    # ---------------------------------------
    # Fig 4B: Heatmap class x signature
    # ---------------------------------------
    # Filter classes by presence across signatures
    cls_counts = df.groupby("class")["signature_name"].nunique().rename("n_signatures").reset_index()
    keep = cls_counts[cls_counts["n_signatures"] >= args.min_signatures]["class"].tolist()
    dfh = df[df["class"].isin(keep)].copy()

    pivot = dfh.pivot_table(index="class", columns="signature_name", values=metric, aggfunc="mean")
    # Sort rows by overall mean
    pivot["__mean__"] = pivot.mean(axis=1, skipna=True)
    pivot = pivot.sort_values("__mean__", ascending=False).drop(columns="__mean__")

    plt.figure(figsize=(max(6, 1.8 * len(pivot.columns)), max(4, 0.35 * len(pivot.index))))
    plt.imshow(pivot.fillna(0.0).values, aspect="auto")
    plt.yticks(range(len(pivot.index)), pivot.index)
    plt.xticks(range(len(pivot.columns)), pivot.columns, rotation=45, ha="right")
    plt.colorbar(label=metric)
    plt.title(f"L1000 class heatmap (class × signature) | metric={metric}")
    plt.tight_layout()

    out = out_dir / f"Fig4B_L1000_class_heatmap.{metric}.pdf"
    plt.savefig(out)
    plt.close()
    print(f"[OK] Wrote: {out}")


if __name__ == "__main__":
    main()

