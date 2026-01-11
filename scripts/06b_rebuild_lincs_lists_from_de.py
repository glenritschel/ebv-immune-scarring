#!/usr/bin/env python3
"""
Rebuild LINCS up/down gene lists from existing DE TSV tables.

This is a "repair" step when *.LINCS_*.up/down.txt were written empty.

Heuristics:
- gene column: first found among ["gene", "symbol", "names", "feature", "index"]
- logFC column: first found among ["logfoldchanges","logFC","log2FC","avg_log2FC"]
- adj p column: first found among ["pvals_adj","padj","qval","fdr","adj_pval"]

Outputs:
- results/tables/<DE basename>.LINCS.up.txt
- results/tables/<DE basename>.LINCS.down.txt
(one gene symbol per line)
"""

from __future__ import annotations
import argparse
from pathlib import Path
from typing import Optional, List

import pandas as pd


def pick_col(df: pd.DataFrame, candidates: List[str]) -> Optional[str]:
    cols = {c.lower(): c for c in df.columns}
    for c in candidates:
        if c.lower() in cols:
            return cols[c.lower()]
    return None


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--de", required=True, help="Path to DE TSV")
    ap.add_argument("--out-dir", default="results/tables", help="Output directory (default results/tables)")
    ap.add_argument("--padj-max", type=float, default=0.05, help="Adj p-value threshold (default 0.05)")
    ap.add_argument("--min-logfc", type=float, default=0.25, help="Min abs logFC (default 0.25)")
    ap.add_argument("--max-genes", type=int, default=150, help="Max genes per list (default 150)")
    args = ap.parse_args()

    de_path = Path(args.de)
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    df = pd.read_csv(de_path, sep="\t")
    if len(df) == 0:
        raise SystemExit(f"{de_path}: DE table has 0 rows")

    gene_col = pick_col(df, ["gene", "symbol", "names", "feature", "index"])
    lfc_col = pick_col(df, ["logfoldchanges", "logfc", "log2fc", "avg_log2fc"])
    padj_col = pick_col(df, ["pvals_adj", "padj", "qval", "fdr", "adj_pval"])

    if not gene_col or not lfc_col:
        raise SystemExit(f"{de_path}: missing gene/logFC columns. cols={list(df.columns)}")

    df = df.copy()
    df[gene_col] = df[gene_col].astype(str)

    df[lfc_col] = pd.to_numeric(df[lfc_col], errors="coerce")
    if padj_col:
        df[padj_col] = pd.to_numeric(df[padj_col], errors="coerce")

    # Filter
    filt = df.dropna(subset=[lfc_col]).copy()
    if padj_col:
        filt = filt[(filt[padj_col].isna()) | (filt[padj_col] <= args.padj_max)]
    filt = filt[filt[lfc_col].abs() >= args.min_logfc]

    if len(filt) == 0:
        raise SystemExit(
            f"{de_path}: no genes pass filters (padj<={args.padj_max}, abs(logFC)>={args.min_logfc}). "
            "Lower thresholds or verify DE."
        )

    up = filt.sort_values(lfc_col, ascending=False).head(args.max_genes)[gene_col].tolist()
    down = filt.sort_values(lfc_col, ascending=True).head(args.max_genes)[gene_col].tolist()

    base = de_path.name.replace(".tsv", "")
    up_path = out_dir / f"{base}.LINCS.up.txt"
    down_path = out_dir / f"{base}.LINCS.down.txt"

    up_path.write_text("\n".join(up) + "\n", encoding="utf-8")
    down_path.write_text("\n".join(down) + "\n", encoding="utf-8")

    print(f"[OK] Wrote UP:   {up_path} (n={len(up)})")
    print(f"[OK] Wrote DOWN: {down_path} (n={len(down)})")


if __name__ == "__main__":
    main()

