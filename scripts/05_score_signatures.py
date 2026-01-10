# scripts/05_score_signatures.py
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

import warnings
warnings.filterwarnings("ignore", message="pkg_resources is deprecated as an API.*")

import scanpy as sc
import yaml


def score_signature(adata: sc.AnnData, name: str, genes: list[str]) -> tuple[int, int]:
    present = [g for g in genes if g in adata.var_names]
    missing = [g for g in genes if g not in adata.var_names]

    if len(present) < 3:
        # too few genes to be meaningful; still create column with NaNs
        adata.obs[f"score_{name}"] = np.nan
        return len(present), len(missing)

    sc.tl.score_genes(
        adata,
        gene_list=present,
        score_name=f"score_{name}",
        use_raw=False,  # we store normalized/log in X; counts in layers["counts"]
    )
    return len(present), len(missing)


def per_cluster_summary(adata: sc.AnnData, cluster_key: str, score_cols: list[str]) -> pd.DataFrame:
    rows = []
    clusters = adata.obs[cluster_key].astype(str)
    for cl in sorted(clusters.unique(), key=lambda x: int(x) if x.isdigit() else x):
        idx = clusters == cl
        row = {
            "cluster": cl,
            "n_cells": int(idx.sum()),
        }
        for c in score_cols:
            v = adata.obs.loc[idx, c].to_numpy()
            row[f"{c}_mean"] = float(np.nanmean(v))
            row[f"{c}_median"] = float(np.nanmedian(v))
        rows.append(row)
    return pd.DataFrame(rows).sort_values("n_cells", ascending=False)


def main() -> int:
    ap = argparse.ArgumentParser(description="Score gene signatures on scVI outputs and export per-cluster summaries.")
    ap.add_argument("--config", default="configs/signatures.yaml", help="Signature YAML")
    ap.add_argument("--gse", action="append", help="Dataset accession (repeatable). If omitted, run both.")
    ap.add_argument("--in-root", default="data/processed", help="Where *.scvi.h5ad live")
    ap.add_argument("--out-root", default="data/processed", help="Where to write *.scvi_scored.h5ad")
    ap.add_argument("--results-root", default="results/tables", help="Where to write score summary tables")
    ap.add_argument("--cluster-key", default="leiden_scvi", help="Cluster column in adata.obs")
    args = ap.parse_args()

    cfg = yaml.safe_load(Path(args.config).read_text())
    sigs: dict = cfg.get("signatures", {})
    if not sigs:
        raise ValueError("No signatures found under `signatures:` in config.")

    datasets = args.gse or ["GSE158275", "GSE267750"]
    in_root = Path(args.in_root)
    out_root = Path(args.out_root)
    results_root = Path(args.results_root)
    out_root.mkdir(parents=True, exist_ok=True)
    results_root.mkdir(parents=True, exist_ok=True)

    for gse in datasets:
        gse = gse.strip().upper()
        in_path = in_root / f"{gse}.scvi.h5ad"
        if not in_path.exists():
            raise FileNotFoundError(f"Missing input: {in_path}")

        print(f"\n=== Scoring {gse} ===")
        adata = sc.read_h5ad(in_path)

        if args.cluster_key not in adata.obs:
            raise ValueError(f"{gse}: cluster key {args.cluster_key!r} not found in obs")

        # score each signature
        score_cols = []
        report = []
        for name, genes in sigs.items():
            present, missing = score_signature(adata, name, list(genes))
            col = f"score_{name}"
            score_cols.append(col)
            report.append({"signature": name, "present": present, "missing": missing})
            print(f"  {name}: present={present} missing={missing}")

        # write signature coverage report
        cov_df = pd.DataFrame(report).sort_values("signature")
        cov_path = results_root / f"{gse}.signature_coverage.tsv"
        cov_df.to_csv(cov_path, sep="\t", index=False)
        print(f"Wrote: {cov_path}")

        # per-cluster summary
        summ = per_cluster_summary(adata, args.cluster_key, score_cols)
        summ_path = results_root / f"{gse}.cluster_signature_scores.tsv"
        summ.to_csv(summ_path, sep="\t", index=False)
        print(f"Wrote: {summ_path}")

        # save scored h5ad
        out_path = out_root / f"{gse}.scvi_scored.h5ad"
        print(f"Writing: {out_path}")
        adata.write_h5ad(out_path)

    print("\nDone.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

