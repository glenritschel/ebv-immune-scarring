#!/usr/bin/env python3
"""
Cluster 7 vs Others differential expression (Wilcoxon) + LINCS UP/DN gene lists.

- Uses adata.obs[groupby] (default: leiden_scvi)
- Creates adata.layers['log1p'] from adata.layers['counts'] if missing
- Writes:
  - tables/cluster7_vs_others_wilcoxon.tsv
  - lincs/<tag>_UP150.txt
  - lincs/<tag>_DN150.txt
"""

from __future__ import annotations

import argparse
import os
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc


def ensure_log1p_layer(adata: sc.AnnData) -> str:
    if "log1p" in adata.layers:
        return "log1p"
    if "counts" not in adata.layers:
        # Fall back to X if no counts layer
        return None

    X = adata.layers["counts"]
    adata.layers["log1p"] = X.copy()
    adata.layers["log1p"] = adata.layers["log1p"].astype(np.float32)
    if hasattr(adata.layers["log1p"], "data"):  # sparse
        adata.layers["log1p"].data = np.log1p(adata.layers["log1p"].data)
    else:
        adata.layers["log1p"] = np.log1p(adata.layers["log1p"])
    return "log1p"


def write_lincs_lists(
    de: pd.DataFrame,
    out_up: Path,
    out_dn: Path,
    n: int = 150,
) -> None:
    # Sort by logFC descending for UP; ascending for DN
    de_up = de.sort_values(["logfoldchanges", "scores"], ascending=[False, False]).copy()
    de_dn = de.sort_values(["logfoldchanges", "scores"], ascending=[True, False]).copy()

    up = []
    for g in de_up["gene"].tolist():
        if g not in up:
            up.append(g)
        if len(up) >= n:
            break

    dn = []
    for g in de_dn["gene"].tolist():
        if g not in dn:
            dn.append(g)
        if len(dn) >= n:
            break

    # Ensure no overlap (rare, but enforce deterministically)
    up_set = set(up)
    dn = [g for g in dn if g not in up_set]
    dn = dn[:n]

    out_up.parent.mkdir(parents=True, exist_ok=True)
    out_dn.parent.mkdir(parents=True, exist_ok=True)

    out_up.write_text("\n".join(up) + "\n", encoding="utf-8")
    out_dn.write_text("\n".join(dn) + "\n", encoding="utf-8")

    print(f"UP n= {len(up)} unique= {len(set(up))}")
    print(f"DN n= {len(dn)} unique= {len(set(dn))}")
    print(f"UP∩DN overlap: {len(set(up).intersection(set(dn)))}")


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--in-h5ad", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--groupby", default="leiden_scvi")
    ap.add_argument("--cluster", default="7")
    ap.add_argument("--tag", default="GSE195452_cluster7_vs_others")
    ap.add_argument("--topn", type=int, default=150)
    args = ap.parse_args()

    in_h5ad = args.in_h5ad
    outdir = Path(args.outdir)
    lincs_dir = outdir / "lincs"
    tab_dir = outdir / "tables"

    outdir.mkdir(parents=True, exist_ok=True)
    lincs_dir.mkdir(parents=True, exist_ok=True)
    tab_dir.mkdir(parents=True, exist_ok=True)

    print(f"Loading: {in_h5ad}")
    adata = sc.read_h5ad(in_h5ad)

    if args.groupby not in adata.obs.columns:
        raise SystemExit(f"ERROR: adata.obs['{args.groupby}'] not found")

    adata.obs[args.groupby] = adata.obs[args.groupby].astype(str)

    layer = ensure_log1p_layer(adata)
    if layer is None:
        print("WARN: No log1p or counts layer found; DE will use adata.X")
    else:
        print(f"Using layer for DE: {layer}")

    # Create a binary group: scarred (cluster 7) vs other
    adata.obs["scarred_like"] = np.where(adata.obs[args.groupby] == str(args.cluster), "scarred", "other")
    print("Group sizes:")
    print(adata.obs["scarred_like"].value_counts())

    print(f"Running rank_genes_groups (wilcoxon): scarred={args.cluster} vs others; layer={layer}")
    sc.tl.rank_genes_groups(
        adata,
        groupby="scarred_like",
        groups=["scarred"],
        reference="other",
        method="wilcoxon",
        layer=layer,
        pts=True,
        use_raw=False,
    )

    # Extract DE table
    rg = adata.uns["rank_genes_groups"]
    genes = pd.DataFrame({
        "gene": rg["names"]["scarred"],
        "scores": rg["scores"]["scarred"],
        "logfoldchanges": rg["logfoldchanges"]["scarred"],
        "pvals": rg["pvals"]["scarred"],
        "pvals_adj": rg["pvals_adj"]["scarred"],
    })

    out_de = tab_dir / "cluster7_vs_others_wilcoxon.tsv"
    genes.to_csv(out_de, sep="\t", index=False)
    print(f"Wrote DE table: {out_de} (rows={genes.shape[0]})")

    out_up = lincs_dir / f"{args.tag}_UP{args.topn}.txt"
    out_dn = lincs_dir / f"{args.tag}_DN{args.topn}.txt"
    write_lincs_lists(genes, out_up, out_dn, n=args.topn)
    print(f"Wrote LINCS UP{args.topn}: {out_up}")
    print(f"Wrote LINCS DN{args.topn}: {out_dn}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())

