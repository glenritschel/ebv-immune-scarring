#!/usr/bin/env python3
import argparse
import pandas as pd
import scanpy as sc
from scipy.stats import mannwhitneyu


def main():
    ap = argparse.ArgumentParser(description="Per-GSM fraction of scarred B-cell cluster (cluster 2) and case/control stats.")
    ap.add_argument("--h5ad", required=True, help="Input B-cell scVI h5ad (e.g., data/processed/GSE195452.bcells_scvi.h5ad)")
    ap.add_argument("--meta", required=True, help="GSM metadata TSV (must include gsm, condition, patient_id, selection_marker, tissue)")
    ap.add_argument("--cluster-key", default="leiden_scvi_bcell", help="Obs column containing B-cell clusters")
    ap.add_argument("--scar-cluster", default="2", help="Cluster label treated as scarred-like (default: '2')")
    ap.add_argument("--out", required=True, help="Output TSV path")
    args = ap.parse_args()

    adata = sc.read_h5ad(args.h5ad)
    meta = pd.read_csv(args.meta, sep="\t")

    need_cols = ["gsm", "condition", "patient_id", "selection_marker", "tissue"]
    missing_cols = [c for c in need_cols if c not in meta.columns]
    if missing_cols:
        raise SystemExit(f"[ERROR] metadata TSV missing columns: {missing_cols}")

    adata.obs = adata.obs.merge(meta[need_cols], on="gsm", how="left")

    if adata.obs["condition"].isna().any():
        missing = adata.obs.loc[adata.obs["condition"].isna(), "gsm"].nunique()
        print(f"[WARN] {missing} GSMs in adata not found in metadata table (condition missing).")

    if args.cluster_key not in adata.obs.columns:
        raise SystemExit(f"[ERROR] cluster key not found in adata.obs: {args.cluster_key}")

    adata.obs[args.cluster_key] = adata.obs[args.cluster_key].astype(str)
    scar_cluster = str(args.scar_cluster)

    adata.obs["is_scar_cluster"] = (adata.obs[args.cluster_key] == scar_cluster).astype(int)

    g = (
        adata.obs.groupby("gsm")
        .agg(
            n_cells=("is_scar_cluster", "size"),
            n_scar=("is_scar_cluster", "sum"),
            frac_scar=("is_scar_cluster", "mean"),
            condition=("condition", "first"),
            patient_id=("patient_id", "first"),
            selection_marker=("selection_marker", "first"),
            tissue=("tissue", "first"),
        )
        .reset_index()
    )

    g.to_csv(args.out, sep="\t", index=False)
    print("[OK] Wrote:", args.out)
    print(g["condition"].value_counts(dropna=False))

    # Case/control stats on per-GSM fraction
    g2 = g[g["condition"].isin(["case", "control"])].copy()
    case = g2.loc[g2["condition"] == "case", "frac_scar"].values
    ctrl = g2.loc[g2["condition"] == "control", "frac_scar"].values

    if len(case) and len(ctrl):
        stat = mannwhitneyu(case, ctrl, alternative="two-sided")
        print(f"[STAT] Mann–Whitney U: U={stat.statistic:.3g} p={stat.pvalue:.3g}")
        print(f"[STAT] Median frac: case={pd.Series(case).median():.6f}  control={pd.Series(ctrl).median():.6f}")
    else:
        print(f"[STAT] Not enough labeled GSMs for case/control comparison: n_case={len(case)} n_ctrl={len(ctrl)}")


if __name__ == "__main__":
    main()

