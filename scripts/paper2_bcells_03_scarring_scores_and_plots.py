#!/usr/bin/env python3
import os
import numpy as np
import scanpy as sc
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

IN_H5AD = "data/processed/GSE195452.bcells_scvi.h5ad"
OUTDIR  = "results/paper2_bcells_scvi"
FIGDIR  = os.path.join(OUTDIR, "figures")
TABDIR  = os.path.join(OUTDIR, "tables")

SCARRING_GENES = [
    "ISG15","MX1","OAS1","IFIT1","IFIT3","IFI44L","IFI27",
    "HLA-DRA","HLA-DRB1","HLA-DPA1","HLA-DPB1","HLA-DQA1","HLA-DQB1",
    "CD74","CD83","NFKBIA","RELA","NFKB1","NFKB2","REL",
]

def ensure_log1p_layer(adata):
    if "log1p" in adata.layers:
        return
    if "counts" in adata.layers:
        X = adata.layers["counts"]
        adata.layers["log1p"] = X.copy()
        adata.layers["log1p"] = adata.layers["log1p"].astype(np.float32)
        if hasattr(adata.layers["log1p"], "data"):
            adata.layers["log1p"].data = np.log1p(adata.layers["log1p"].data)
        else:
            adata.layers["log1p"] = np.log1p(adata.layers["log1p"])
    else:
        adata.layers["log1p"] = adata.X.copy()

def mean_score(adata, genes, layer="log1p"):
    genes = [g for g in genes if g in adata.var_names]
    if not genes:
        return np.zeros(adata.n_obs, dtype=float)
    X = adata[:, genes].layers[layer]
    if hasattr(X, "toarray"):
        X = X.toarray()
    X = np.asarray(X)
    return X.mean(axis=1)

def main():
    os.makedirs(FIGDIR, exist_ok=True)
    os.makedirs(TABDIR, exist_ok=True)

    adata = sc.read_h5ad(IN_H5AD)
    ensure_log1p_layer(adata)

    # simple module scores
    isg = ["ISG15","MX1","OAS1","IFIT1","IFIT3","IFI44L","IFI27"]
    hla = ["HLA-DRA","HLA-DRB1","HLA-DPA1","HLA-DPB1","HLA-DQA1","HLA-DQB1","CD74"]
    nfk = ["NFKBIA","RELA","NFKB1","NFKB2","REL","CD83"]

    adata.obs["score_ISG"] = mean_score(adata, isg)
    adata.obs["score_HLA"] = mean_score(adata, hla)
    adata.obs["score_NFkB"] = mean_score(adata, nfk)
    adata.obs["score_scarring"] = mean_score(adata, SCARRING_GENES)

    # cluster-level summaries
    g = (
        adata.obs.groupby("leiden_scvi_bcell")[["score_ISG","score_HLA","score_NFkB","score_scarring"]]
        .mean()
        .sort_values("score_scarring", ascending=False)
        .reset_index()
    )
    g["leiden_scvi_bcell"] = g["leiden_scvi_bcell"].astype(str).str.replace(r"\.0$", "", regex=True)

    out_tsv = os.path.join(TABDIR, "bcells_cluster_scarring_scores.tsv")
    g.to_csv(out_tsv, sep="\t", index=False)

    # plots
    sc.pl.umap(adata, color=["leiden_scvi_bcell"], show=False)
    plt.savefig(os.path.join(FIGDIR, "umap_bcells_leiden_scvi.png"), dpi=200); plt.close()

    sc.pl.umap(adata, color=["score_ISG","score_HLA","score_NFkB","score_scarring"], show=False)
    plt.savefig(os.path.join(FIGDIR, "umap_bcells_scarring_scores.png"), dpi=200); plt.close()

    # heatmap-like matrixplot across clusters
    present = [g for g in SCARRING_GENES if g in adata.var_names]
    sc.pl.matrixplot(
        adata,
        var_names=present,
        groupby="leiden_scvi_bcell",
        layer="log1p",
        standard_scale="var",
        dendrogram=False,
        show=False
    )
    plt.gcf().set_size_inches(12, 6)
    plt.tight_layout()
    plt.savefig(os.path.join(FIGDIR, "matrixplot_bcells_scarring_genes.png"), dpi=200)
    plt.close()

    adata.write(IN_H5AD)  # save scores back into the same file
    print("Wrote:", out_tsv)
    print("Figures:", FIGDIR)
    print("\nTop clusters by scarring score:\n", g.head(10))

if __name__ == "__main__":
    main()

