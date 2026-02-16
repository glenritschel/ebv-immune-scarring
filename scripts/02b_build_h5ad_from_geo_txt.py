#!/usr/bin/env python

import gzip
import os
from pathlib import Path

import pandas as pd
import scipy.sparse as sp
import anndata as ad


def load_gsm_txt(path):
    print(f"Loading {path.name}")

    df = pd.read_csv(
        path,
        sep="\t",
        index_col=0,
        compression="gzip"
    )

    # genes x cells → cells x genes
    X = sp.csr_matrix(df.values.T)

    adata = ad.AnnData(
        X=X,
        obs=pd.DataFrame(index=df.columns),
        var=pd.DataFrame(index=df.index)
    )

    return adata


def main():
    raw_dir = Path("data/raw/GSE195452/unpacked")
    out_dir = Path("data/processed")
    out_dir.mkdir(parents=True, exist_ok=True)

    gsm_files = sorted(raw_dir.glob("GSM*.txt.gz"))
    if not gsm_files:
        raise RuntimeError("No GSM txt files found.")

    adatas = []

    for f in gsm_files:
        adata = load_gsm_txt(f)

        gsm_id = f.name.split("_")[0]
        adata.obs["gsm"] = gsm_id

        adatas.append(adata)

    print("Concatenating...")
    adata_all = ad.concat(adatas, join="outer", merge="same")

    out_path = out_dir / "GSE195452.raw.h5ad"
    print(f"Writing {out_path}")
    adata_all.write_h5ad(out_path)

    print("Done.")


if __name__ == "__main__":
    main()

