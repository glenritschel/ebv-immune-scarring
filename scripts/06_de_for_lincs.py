# scripts/06_de_for_lincs.py
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

import warnings
warnings.filterwarnings("ignore", message="pkg_resources is deprecated as an API.*")

import scanpy as sc


def run_de(
    adata: sc.AnnData,
    group_key: str,
    group: str,
    reference: str,
    method: str = "wilcoxon",
) -> pd.DataFrame:
    """
    Run DE for group vs reference using scanpy rank_genes_groups.
    Returns a tidy dataframe with gene, logfc, pvals, pvals_adj, scores.
    """
    sc.tl.rank_genes_groups(
        adata,
        groupby=group_key,
        groups=[group],
        reference=reference,
        method=method,
        pts=True,
        use_raw=False,
    )
    res = adata.uns["rank_genes_groups"]
    g = group

    df = pd.DataFrame({
        "gene": res["names"][g],
        "logfoldchanges": res["logfoldchanges"][g],
        "scores": res["scores"][g],
        "pvals": res["pvals"][g],
        "pvals_adj": res["pvals_adj"][g],
    })

    # Percent expressing
    if "pts" in res:
        df["pct_in_group"] = res["pts"][g]

    # For reference="rest", scanpy uses pts_rest. For a specific reference group, pts_rest may not exist.
    if "pts_rest" in res:
        df["pct_in_rest"] = res["pts_rest"][g]

    return df


def write_lincs_lists(df: pd.DataFrame, out_prefix: Path, top_n: int = 250) -> None:
    d = df.dropna(subset=["gene", "logfoldchanges"]).copy()

    # Optional: drop some common noise (keep or remove as you like)
    drop_prefixes = ("MT-", "RPL", "RPS")
    d = d[~d["gene"].astype(str).str.startswith(drop_prefixes)]

    up = d.sort_values("logfoldchanges", ascending=False).head(top_n)["gene"].astype(str).tolist()
    down = d.sort_values("logfoldchanges", ascending=True).head(top_n)["gene"].astype(str).tolist()

    up_path = Path(str(out_prefix) + ".up.txt")
    down_path = Path(str(out_prefix) + ".down.txt")
    up_path.write_text("\n".join(up) + "\n")
    down_path.write_text("\n".join(down) + "\n")


def main() -> int:
    ap = argparse.ArgumentParser(description="Run DE contrasts for LINCS input.")
    ap.add_argument("--in-root", default="data/processed", help="Where *.scvi_scored.h5ad live")
    ap.add_argument("--results-root", default="results/tables", help="Where to write DE outputs")
    ap.add_argument("--cluster-key", default="leiden_scvi", help="Cluster key in obs")
    args = ap.parse_args()

    in_root = Path(args.in_root)
    results_root = Path(args.results_root)
    results_root.mkdir(parents=True, exist_ok=True)

    contrasts = [
        # Clean B-cell internal contrasts
        ("GSE267750", "4", "18"),   # B activation / scarring
        ("GSE158275", "13", "2"),   # EBV high vs low LCL state
    ]

    for gse, group, ref in contrasts:
        in_path = in_root / f"{gse}.scvi_scored.h5ad"
        if not in_path.exists():
            raise FileNotFoundError(f"Missing input: {in_path}")

        print(f"\n=== DE {gse}: {group} vs {('rest' if ref is None else ref)} ===")
        adata = sc.read_h5ad(in_path)

        # Ensure group labels are strings
        adata.obs[args.cluster_key] = adata.obs[args.cluster_key].astype(str)

        reference = "rest" if ref is None else str(ref)
        df = run_de(adata, group_key=args.cluster_key, group=str(group), reference=reference)

        out_tsv = results_root / f"{gse}.DE_{args.cluster_key}_{group}_vs_{'rest' if ref is None else ref}.tsv"
        df.to_csv(out_tsv, sep="\t", index=False)
        print(f"Wrote: {out_tsv}")

        # LINCS lists
        out_prefix = results_root / f"{gse}.LINCS_{args.cluster_key}_{group}_vs_{'rest' if ref is None else ref}"
        drop_prefixes = ("MT-", "RPL", "RPS", "ACT", "EEF", "EIF")
        df = df[~df["gene"].astype(str).str.startswith(drop_prefixes)]

        if "pct_in_group" in df.columns:
            df = df[df["pct_in_group"] >= 0.1]

        write_lincs_lists(df, out_prefix, top_n=250)
        print(f"Wrote: {out_prefix}.up.txt and .down.txt")

    print("\nDone.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

