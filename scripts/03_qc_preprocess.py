# scripts/03_qc_preprocess.py
from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd
import scanpy as sc
import yaml

from src.qc import (
    add_qc_metrics,
    conservative_filter,
    normalize_log1p,
    params_from_config,
    qc_summary_row,
)


def main() -> int:
    ap = argparse.ArgumentParser(description="Conservative QC + preprocessing for GEO-derived .raw.h5ad datasets.")
    ap.add_argument("--config", default="configs/qc.yaml", help="QC config YAML")
    ap.add_argument("--gse", action="append", help="Dataset accession (repeatable). If omitted, run both.")
    ap.add_argument("--processed-root", default="data/processed", help="Directory containing *.raw.h5ad")
    ap.add_argument("--out-root", default="data/processed", help="Where to write *.qc.h5ad")
    ap.add_argument("--results-root", default="results/tables", help="Where to write qc_summary.tsv")
    args = ap.parse_args()

    cfg = yaml.safe_load(Path(args.config).read_text())
    processed_root = Path(args.processed_root)
    out_root = Path(args.out_root)
    results_root = Path(args.results_root)

    out_root.mkdir(parents=True, exist_ok=True)
    results_root.mkdir(parents=True, exist_ok=True)

    # default datasets if not provided
    datasets = args.gse or ["GSE158275", "GSE267750"]

    summary_rows = []

    for gse in datasets:
        gse = gse.strip().upper()
        in_path = processed_root / f"{gse}.raw.h5ad"
        if not in_path.exists():
            raise FileNotFoundError(f"Missing input: {in_path}")

        print(f"\n=== QC {gse} ===")
        params = params_from_config(cfg, gse)
        print(f"Params: min_genes={params.min_genes} min_counts={params.min_counts} "
              f"max_pct_mito={params.max_pct_mito} min_cells_per_gene={params.min_cells_per_gene}")

        adata = sc.read_h5ad(in_path)

        # QC metrics
        add_qc_metrics(adata, params.mito_prefixes)

        # keep a copy for summary
        adata_before = adata.copy()

        # filter (conservative)
        cells_removed, genes_removed = conservative_filter(
            adata,
            min_genes=params.min_genes,
            min_counts=params.min_counts,
            max_pct_mito=params.max_pct_mito,
            min_cells_per_gene=params.min_cells_per_gene,
        )
        print(f"Filtered: removed_cells={cells_removed} removed_genes={genes_removed}")
        print(f"Remaining: cells={adata.n_obs} genes={adata.n_vars}")

        # normalize + log
        normalize_log1p(adata, target_sum=params.target_sum, do_log1p=params.log1p)

        # summary row
        summary_rows.append(qc_summary_row(gse, adata_before, adata, params))

        out_path = out_root / f"{gse}.qc.h5ad"
        print(f"Writing {out_path} ...")
        adata.write_h5ad(out_path)

    # write summary table
    df = pd.DataFrame(summary_rows).sort_values("dataset")
    summary_path = results_root / "qc_summary.tsv"
    df.to_csv(summary_path, sep="\t", index=False)
    print(f"\nWrote QC summary: {summary_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

