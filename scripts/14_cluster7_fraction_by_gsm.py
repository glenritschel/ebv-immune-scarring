#!/usr/bin/env python3
"""
Compute cluster7 fraction per GSM and compare case vs control using Mann–Whitney U.

Outputs:
- tables/cluster7_fraction_by_gsm.tsv
Prints MWU results for:
- All GSMs
- Skin only
- Skin + CD45+
- Skin + CD90+
"""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd
import scanpy as sc
from scipy.stats import mannwhitneyu


def mw(df: pd.DataFrame, label: str) -> None:
    df = df[df["condition"].isin(["case", "control"])].copy()
    case = df.loc[df["condition"] == "case", "frac_cluster7"].values
    ctrl = df.loc[df["condition"] == "control", "frac_cluster7"].values
    if len(case) < 5 or len(ctrl) < 5:
        print(label, "not enough samples:", len(case), len(ctrl))
        return
    stat = mannwhitneyu(case, ctrl, alternative="two-sided")
    print(
        f"{label}: n_case={len(case)} n_ctrl={len(ctrl)} "
        f" U={stat.statistic:.3g} p={stat.pvalue:.3g} "
        f" mean/median case={case.mean():.4f}/{pd.Series(case).median():.4f} "
        f" mean/median ctrl={ctrl.mean():.4f}/{pd.Series(ctrl).median():.4f}"
    )


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--in-h5ad", required=True)
    ap.add_argument("--meta-tsv", required=True)
    ap.add_argument("--out-tsv", required=True)
    ap.add_argument("--groupby", default="leiden_scvi")
    ap.add_argument("--cluster", default="7")
    args = ap.parse_args()

    adata = sc.read_h5ad(args.in_h5ad)
    meta = pd.read_csv(args.meta_tsv, sep="\t")

    if "gsm" not in adata.obs.columns:
        raise SystemExit("ERROR: adata.obs['gsm'] not found")
    if args.groupby not in adata.obs.columns:
        raise SystemExit(f"ERROR: adata.obs['{args.groupby}'] not found")

    adata.obs[args.groupby] = adata.obs[args.groupby].astype(str)

    # merge metadata onto obs
    adata.obs = adata.obs.merge(
        meta[["gsm", "patient_id", "condition", "selection_marker", "tissue"]],
        on="gsm",
        how="left",
    )

    missing = adata.obs.loc[adata.obs["condition"].isna(), "gsm"].nunique()
    if missing > 0:
        print(f"WARNING: {missing} GSMs in adata not found in metadata table.")

    adata.obs["is_cluster7"] = (adata.obs[args.groupby] == str(args.cluster)).astype(int)

    g = (
        adata.obs.groupby("gsm")
        .agg(
            n_cells=("is_cluster7", "size"),
            n_cluster7=("is_cluster7", "sum"),
            frac_cluster7=("is_cluster7", "mean"),
            condition=("condition", "first"),
            patient_id=("patient_id", "first"),
            selection_marker=("selection_marker", "first"),
            tissue=("tissue", "first"),
        )
        .reset_index()
    )

    out = Path(args.out_tsv)
    out.parent.mkdir(parents=True, exist_ok=True)
    g.to_csv(out, sep="\t", index=False)
    print("Wrote:", out)
    print(g["condition"].value_counts(dropna=False))

    # Stats
    mw(g, "All GSMs")
    if g["tissue"].notna().any():
        mw(g[g["tissue"].str.lower().eq("skin")], "Skin only")
        mw(g[(g["tissue"].str.lower().eq("skin")) & (g["selection_marker"].str.upper().eq("CD45+"))], "Skin + CD45+")
        mw(g[(g["tissue"].str.lower().eq("skin")) & (g["selection_marker"].str.upper().eq("CD90+"))], "Skin + CD90+")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())

