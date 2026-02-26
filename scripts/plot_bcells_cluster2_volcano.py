#!/usr/bin/env python3
"""
Make a volcano plot + Top-20 DE table for:
  results/paper2_bcells_scvi/tables/bcells_scarredCluster2_vs_others_wilcoxon.tsv

Outputs:
  - results/paper2_bcells_scvi/figures/volcano_bcells_cluster2_vs_others.png
  - results/paper2_bcells_scvi/figures/volcano_bcells_cluster2_vs_others.pdf
  - results/paper2_bcells_scvi/tables/bcells_scarredCluster2_vs_others_top20.tsv

This script is defensive about column names (Scanpy variants, your own TSV variants).
"""

import argparse
import math
import os
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


HIGHLIGHT_GENES_DEFAULT = [
    "CD74", "HLA-DRA", "HLA-DRB1", "HLA-DPA1", "HLA-DPB1", "HLA-DQB1", "HLA-DQA1",
    "IFI27", "IFIT1", "IFIT3", "ISG15", "MX1", "OAS1", "OASL",
    "NFKBIA", "TNFAIP3", "REL", "RELN", "IRF7", "STAT1",
]


def pick_col(df: pd.DataFrame, candidates):
    """Return the first matching column in candidates."""
    lower_map = {c.lower(): c for c in df.columns}
    for cand in candidates:
        if cand in df.columns:
            return cand
        lc = cand.lower()
        if lc in lower_map:
            return lower_map[lc]
    return None


def guess_gene_col(df: pd.DataFrame):
    return pick_col(df, ["gene", "genes", "symbol", "gene_symbol", "names", "name"])


def guess_logfc_col(df: pd.DataFrame):
    return pick_col(df, [
        "log2fc", "log2_fc", "logFC", "log_fc", "logfoldchange", "logfoldchanges",
        "avg_log2FC", "avg_logfc", "lfc",
    ])


def guess_padj_col(df: pd.DataFrame):
    return pick_col(df, [
        "p_adj", "padj", "pvals_adj", "pval_adj", "adj_pval", "fdr", "qval", "q_value",
    ])


def guess_pval_col(df: pd.DataFrame):
    return pick_col(df, ["pval", "p_value", "pvalue", "pvals"])


def safe_neglog10(x: np.ndarray):
    x = np.asarray(x, dtype=float)
    # avoid -inf
    x = np.where(np.isfinite(x) & (x > 0), x, np.nan)
    return -np.log10(x)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--de-tsv", required=True, help="DE TSV (wilcoxon) file")
    ap.add_argument("--out-dir", default="results/paper2_bcells_scvi/figures", help="Output figures directory")
    ap.add_argument("--out-top20", default="results/paper2_bcells_scvi/tables/bcells_scarredCluster2_vs_others_top20.tsv",
                    help="Output Top-20 TSV path")
    ap.add_argument("--padj-thresh", type=float, default=0.05, help="Significance threshold for adjusted p-value")
    ap.add_argument("--logfc-thresh", type=float, default=0.25, help="Effect-size threshold (abs logFC)")
    ap.add_argument("--highlight", nargs="*", default=HIGHLIGHT_GENES_DEFAULT, help="Genes to highlight")
    ap.add_argument("--title", default="B cells: cluster 2 (scarred) vs others (Wilcoxon)", help="Plot title")
    args = ap.parse_args()

    de_path = Path(args.de_tsv)
    if not de_path.exists():
        raise SystemExit(f"Missing --de-tsv: {de_path}")

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    out_png = out_dir / "volcano_bcells_cluster2_vs_others.png"
    out_pdf = out_dir / "volcano_bcells_cluster2_vs_others.pdf"
    out_top20 = Path(args.out_top20)
    out_top20.parent.mkdir(parents=True, exist_ok=True)

    df = pd.read_csv(de_path, sep="\t")

    gene_col = guess_gene_col(df)
    logfc_col = guess_logfc_col(df)
    padj_col = guess_padj_col(df)
    pval_col = guess_pval_col(df)

    if gene_col is None:
        raise SystemExit(f"Could not find a gene column. Columns: {df.columns.tolist()}")
    if logfc_col is None:
        raise SystemExit(f"Could not find a logFC column. Columns: {df.columns.tolist()}")

    # prefer padj; fall back to pval if needed
    p_col = padj_col or pval_col
    if p_col is None:
        raise SystemExit(f"Could not find a p-value / adjusted p-value column. Columns: {df.columns.tolist()}")

    # normalize gene symbols to string
    df[gene_col] = df[gene_col].astype(str)

    # numeric conversions
    df[logfc_col] = pd.to_numeric(df[logfc_col], errors="coerce")
    df[p_col] = pd.to_numeric(df[p_col], errors="coerce")

    # drop unusable rows
    df = df.dropna(subset=[logfc_col, p_col, gene_col]).copy()
    df = df[df[p_col] > 0].copy()

    df["neglog10p"] = safe_neglog10(df[p_col].values)

    # mark significance
    df["sig"] = (df[p_col] <= args.padj_thresh) & (df[logfc_col].abs() >= args.logfc_thresh)

    # Top-20 DE table (by smallest p then largest abs(logFC))
    top20 = (
        df.sort_values([p_col, df[logfc_col].abs().name], ascending=[True, False])
          .head(20)
          .copy()
    )

    # keep a clean subset of columns for the supplement table
    keep_cols = [gene_col, logfc_col, p_col]
    # add common useful columns if present
    for extra in ["score", "scores", "statistic", "stats", "pval", "p_value", "pvals", "pvals_adj"]:
        c = pick_col(df, [extra])
        if c and c not in keep_cols:
            keep_cols.append(c)

    top20_out = top20[keep_cols].rename(columns={
        gene_col: "gene",
        logfc_col: "logFC",
        p_col: "padj_or_pval",
    })
    top20_out.to_csv(out_top20, sep="\t", index=False)
    print(f"[OK] Wrote Top-20: {out_top20}")

    # --- Volcano plot ---
    x = df[logfc_col].values
    y = df["neglog10p"].values

    # base points
    plt.figure(figsize=(9, 7))
    plt.scatter(x[~df["sig"].values], y[~df["sig"].values], s=8, alpha=0.35)
    plt.scatter(x[df["sig"].values],  y[df["sig"].values],  s=10, alpha=0.75)

    # thresholds
    y_thresh = -math.log10(args.padj_thresh) if args.padj_thresh > 0 else None
    if y_thresh is not None and np.isfinite(y_thresh):
        plt.axhline(y_thresh, linewidth=1)
    plt.axvline(+args.logfc_thresh, linewidth=1)
    plt.axvline(-args.logfc_thresh, linewidth=1)

    # highlight genes
    highlights = set([g.upper() for g in args.highlight])
    hmask = df[gene_col].str.upper().isin(highlights)
    if hmask.any():
        plt.scatter(df.loc[hmask, logfc_col], df.loc[hmask, "neglog10p"], s=35, alpha=0.95)
        # label highlights (small offset to reduce overlap)
        for _, r in df.loc[hmask, [gene_col, logfc_col, "neglog10p"]].iterrows():
            plt.text(r[logfc_col] + 0.02, r["neglog10p"] + 0.05, r[gene_col], fontsize=9)

    plt.title(args.title)
    plt.xlabel("log fold-change (cluster 2 vs others)")
    plt.ylabel(f"-log10({p_col})")
    plt.tight_layout()

    plt.savefig(out_png, dpi=300)
    plt.savefig(out_pdf)
    plt.close()
    print(f"[OK] Wrote volcano: {out_png}")
    print(f"[OK] Wrote volcano: {out_pdf}")

    # quick console summary
    n = len(df)
    nsig = int(df["sig"].sum())
    print(f"[INFO] Rows={n}  Significant (padj<={args.padj_thresh} & abs(logFC)>={args.logfc_thresh}) = {nsig}")
    missing_highlights = sorted([g for g in args.highlight if g.upper() not in set(df[gene_col].str.upper())])
    if missing_highlights:
        print("[INFO] Highlight genes not found in DE table:", ", ".join(missing_highlights))


if __name__ == "__main__":
    main()

