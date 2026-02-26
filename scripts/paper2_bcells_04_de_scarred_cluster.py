#!/usr/bin/env python3
import os
import numpy as np
import pandas as pd
import scanpy as sc

IN_H5AD = "data/processed/GSE195452.bcells_scvi.h5ad"
OUTDIR  = "results/paper2_bcells_scvi"
TABDIR  = os.path.join(OUTDIR, "tables")

# Option A: set explicitly (must match leiden_scvi_bcell labels exactly)
SCARRED_CLUSTER = None  # e.g. "3"

# Option B: infer from this table if present (recommended)
CLUSTER_SCORE_TSV = os.path.join(TABDIR, "bcells_cluster_scarring_scores.tsv")

def ensure_log1p_layer(adata):
    if "log1p" in adata.layers:
        return "log1p"
    if "counts" in adata.layers:
        X = adata.layers["counts"]
        adata.layers["log1p"] = X.copy()
        adata.layers["log1p"] = adata.layers["log1p"].astype(np.float32)
        if hasattr(adata.layers["log1p"], "data"):
            adata.layers["log1p"].data = np.log1p(adata.layers["log1p"].data)
        else:
            adata.layers["log1p"] = np.log1p(adata.layers["log1p"])
        return "log1p"
    return None

def pick_scarred_cluster_from_table():
    if not os.path.exists(CLUSTER_SCORE_TSV):
        return None
    df = pd.read_csv(CLUSTER_SCORE_TSV, sep="\t")
    # expects columns: leiden_scvi_bcell, score_scarring
    # but your script wrote "leiden_scvi_bcell" as "leiden_scvi_bcell" *or* as "leiden_scvi_bcell"?
    # Handle both common cases:
    if "leiden_scvi_bcell" in df.columns:
        cluster_col = "leiden_scvi_bcell"
    elif "leiden_scvi_bcell" not in df.columns and "leiden_scvi_bcell" in df.columns:
        cluster_col = "leiden_scvi_bcell"
    else:
        # fall back: first column
        cluster_col = df.columns[0]
    if "score_scarring" not in df.columns:
        return None

    top = df.sort_values("score_scarring", ascending=False).iloc[0][cluster_col]
    # normalize: 2.0 -> "2", "2" -> "2"
    try:
        return str(int(float(top)))
    except Exception:
        return str(top).strip()

    return str(top)

def main():
    os.makedirs(TABDIR, exist_ok=True)

    adata = sc.read_h5ad(IN_H5AD)

    if "leiden_scvi_bcell" not in adata.obs.columns:
        raise SystemExit(
            "Missing adata.obs['leiden_scvi_bcell']. "
            "Did you run paper2_bcells_02_scvi_cluster.py and save to the same IN_H5AD?"
        )

    adata.obs["leiden_scvi_bcell"] = adata.obs["leiden_scvi_bcell"].astype(str)
    clusters = sorted(adata.obs["leiden_scvi_bcell"].unique().tolist(), key=lambda x: (len(x), x))
    print("Available B-cell clusters (leiden_scvi_bcell):", clusters)

    scar = SCARRED_CLUSTER
    if scar is None:
        scar = pick_scarred_cluster_from_table()

    if scar is None:
        raise SystemExit(
            f"SCARRED_CLUSTER not set and {CLUSTER_SCORE_TSV} not found/usable. "
            "Set SCARRED_CLUSTER (string) to one of the available cluster IDs printed above."
        )

    scar = str(scar)
    if scar not in set(clusters):
        raise SystemExit(
            f"SCARRED_CLUSTER='{scar}' not present in leiden_scvi_bcell. "
            f"Available: {clusters}"
        )

    adata.obs["scarred_like"] = np.where(
        adata.obs["leiden_scvi_bcell"] == scar, "scarred", "other"
    )
    print("Group sizes:\n", adata.obs["scarred_like"].value_counts())

    if (adata.obs["scarred_like"] == "scarred").sum() < 10:
        raise SystemExit(
            f"Too few scarred cells in cluster {scar} for stable DE "
            f"(<10). Consider combining top 2–3 scarring clusters."
        )

    layer = ensure_log1p_layer(adata)

    sc.tl.rank_genes_groups(
        adata,
        groupby="scarred_like",
        groups=["scarred"],
        reference="other",
        method="wilcoxon",
        layer=layer,
        pts=True,
    )
    df = sc.get.rank_genes_groups_df(adata, group="scarred")
    out = os.path.join(TABDIR, f"bcells_scarredCluster{scar}_vs_others_wilcoxon.tsv")
    df.to_csv(out, sep="\t", index=False)
    print("Wrote:", out)

if __name__ == "__main__":
    main()

