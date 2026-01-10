#!/usr/bin/env python3
from __future__ import annotations

import argparse
import gzip
import os
import re
import shutil
from pathlib import Path
from typing import Dict, List, Tuple

import scanpy as sc


GSM_RE = re.compile(r"^(GSM\d+)_")


def find_gsm_triples(unpacked_dir: Path) -> Dict[str, Dict[str, Path]]:
    """
    Finds per-GSM 10x MTX components in a 'flat' GEO unpacked directory.
    Accepts either genes.tsv(.gz) or features.tsv(.gz).
    Returns mapping: GSM -> {"barcodes": Path, "genes_or_features": Path, "matrix": Path}
    """
    if not unpacked_dir.exists():
        raise FileNotFoundError(f"Unpacked dir not found: {unpacked_dir}")

    files = [p for p in unpacked_dir.iterdir() if p.is_file()]
    gsm_map: Dict[str, Dict[str, Path]] = {}

    for p in files:
        m = GSM_RE.match(p.name)
        if not m:
            continue
        gsm = m.group(1)
        entry = gsm_map.setdefault(gsm, {})

        lname = p.name.lower()
        if "barcodes.tsv" in lname:
            entry["barcodes"] = p
        elif "genes.tsv" in lname or "features.tsv" in lname:
            entry["genes_or_features"] = p
        elif "matrix.mtx" in lname:
            entry["matrix"] = p

    # validate
    missing: List[str] = []
    for gsm, parts in sorted(gsm_map.items()):
        for k in ("barcodes", "genes_or_features", "matrix"):
            if k not in parts:
                missing.append(f"{gsm}:{k}")
    if missing:
        raise RuntimeError(
            "Missing expected 10x files for some GSMs. "
            "Expected barcodes.tsv(.gz), genes.tsv(.gz) or features.tsv(.gz), matrix.mtx(.gz).\n"
            f"Missing entries: {missing[:25]}" + (" ..." if len(missing) > 25 else "")
        )

    return gsm_map


def safe_symlink_or_copy(src: Path, dst: Path) -> None:
    dst.parent.mkdir(parents=True, exist_ok=True)
    if dst.exists():
        dst.unlink()
    try:
        os.symlink(src.resolve(), dst)
    except Exception:
        shutil.copy2(src, dst)


def build_sample_dir(tmp_root: Path, gsm: str, parts: Dict[str, Path]) -> Path:
    """
    Create standardized directory layout expected by scanpy.read_10x_mtx:
      barcodes.tsv.gz
      genes.tsv.gz (or features.tsv.gz)
      matrix.mtx.gz
    """
    sample_dir = tmp_root / gsm
    if sample_dir.exists():
        shutil.rmtree(sample_dir)
    sample_dir.mkdir(parents=True, exist_ok=True)

    # Standard names
    safe_symlink_or_copy(parts["barcodes"], sample_dir / "barcodes.tsv.gz")

    gpath = parts["genes_or_features"]

    # Scanpy >= 1.9 expects features.tsv(.gz). Many GEO exports use genes.tsv(.gz).
    # To be robust, always stage a features.tsv.gz file.
    safe_symlink_or_copy(gpath, sample_dir / "features.tsv.gz")

    # Also stage genes.tsv.gz for clarity/backward compatibility (optional)
    if "genes.tsv" in gpath.name.lower():
        safe_symlink_or_copy(gpath, sample_dir / "genes.tsv.gz")

    safe_symlink_or_copy(parts["matrix"], sample_dir / "matrix.mtx.gz")
    return sample_dir


def quick_peek_nrows(path: Path, n: int = 3) -> List[str]:
    """
    Quick peek into (gzipped) genes/features file to infer delimiter and column count.
    GEO often uses tab-delimited genes.tsv.gz with 1–3 columns.
    """
    lines: List[str] = []
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", encoding="utf-8", errors="replace") as f:
        for _ in range(n):
            line = f.readline()
            if not line:
                break
            lines.append(line.rstrip("\n"))
    return lines


def load_one_gsm(sample_dir: Path) -> sc.AnnData:
    """
    Robust loader for GEO 'flat 10x' exports where:
      - genes.tsv.gz may be 2-column (id, symbol)
      - features.tsv.gz may be 3-column (id, symbol, feature_type)
    """
    import pandas as pd
    import numpy as np
    from scipy.io import mmread

    # Files (we stage features.tsv.gz always; sometimes genes.tsv.gz too)
    mtx = sample_dir / "matrix.mtx.gz"
    barcodes = sample_dir / "barcodes.tsv.gz"
    features = sample_dir / "features.tsv.gz"

    if not (mtx.exists() and barcodes.exists() and features.exists()):
        raise FileNotFoundError(f"Missing one of matrix/barcodes/features in {sample_dir}")

    # Read matrix
    X = mmread(str(mtx)).tocsr().T

    # Read barcodes
    obs = pd.read_csv(barcodes, header=None, sep="\t")
    obs_names = obs[0].astype(str).tolist()

    # Read features/genes (2-col or 3-col)
    var = pd.read_csv(features, header=None, sep="\t")
    if var.shape[1] == 1:
        # extremely rare; treat as symbols
        gene_symbols = var[0].astype(str)
        gene_ids = gene_symbols
        feature_types = None
    elif var.shape[1] == 2:
        gene_ids = var[0].astype(str)
        gene_symbols = var[1].astype(str)
        feature_types = None
    else:
        gene_ids = var[0].astype(str)
        gene_symbols = var[1].astype(str)
        feature_types = var[2].astype(str)

    # Build AnnData
    adata = sc.AnnData(X=X)
    adata.obs_names = obs_names
    adata.var_names = gene_symbols.values
    adata.var["gene_ids"] = gene_ids.values
    if feature_types is not None:
        adata.var["feature_types"] = feature_types.values

    # Make unique var names (scanpy normally does this)
    adata.var_names_make_unique()

    return adata


def main() -> int:
    ap = argparse.ArgumentParser(
        description="Build a single .h5ad per GEO series from flattened 10x MTX triples in data/raw/<GSE>/unpacked/"
    )
    ap.add_argument("--gse", required=True, help="GEO Series accession, e.g. GSE158275")
    ap.add_argument(
        "--raw-root",
        default="data/raw",
        help="Root directory containing data/raw/<GSE>/unpacked/",
    )
    ap.add_argument(
        "--processed-root",
        default="data/processed",
        help="Output directory for processed .h5ad",
    )
    ap.add_argument(
        "--tmp-root",
        default="data/tmp/10x_stage",
        help="Temp dir for staging per-GSM 10x directories",
    )
    ap.add_argument(
        "--max-samples",
        type=int,
        default=0,
        help="For quick tests: limit number of GSM samples loaded (0 = all)",
    )
    args = ap.parse_args()

    gse = args.gse.strip().upper()
    unpacked_dir = Path(args.raw_root) / gse / "unpacked"
    out_dir = Path(args.processed_root)
    tmp_root = Path(args.tmp_root) / gse

    out_dir.mkdir(parents=True, exist_ok=True)
    tmp_root.mkdir(parents=True, exist_ok=True)

    gsm_map = find_gsm_triples(unpacked_dir)
    gsms = sorted(gsm_map.keys())
    if args.max_samples and args.max_samples > 0:
        gsms = gsms[: args.max_samples]

    print(f"Found {len(gsm_map)} GSM samples in {unpacked_dir}")
    if args.max_samples:
        print(f"Limiting to {len(gsms)} samples due to --max-samples")

    adatas = []
    for i, gsm in enumerate(gsms, start=1):
        parts = gsm_map[gsm]
        sample_dir = build_sample_dir(tmp_root, gsm, parts)

        # Small debug peek (optional but helpful)
        gf = parts["genes_or_features"]
        peek = quick_peek_nrows(gf, n=1)
        if peek:
            print(f"[{i}/{len(gsms)}] {gsm} genes/features peek: {peek[0][:120]}")

        print(f"[{i}/{len(gsms)}] Loading {gsm} ...")
        adata = load_one_gsm(sample_dir)
        adata.obs["gsm"] = gsm
        adata.obs_names = [f"{gsm}_{bc}" for bc in adata.obs_names]
        adata.obs["dataset"] = gse
        adatas.append(adata)

    print("Concatenating samples ...")
    # outer join on genes to be safe across samples
    adata_all = sc.concat(adatas, label="gsm", keys=gsms, join="outer", merge="same")
    adata_all.obs_names_make_unique()

    out_path = out_dir / f"{gse}.raw.h5ad"
    print(f"Writing {out_path} ...")
    adata_all.write_h5ad(out_path)

    print("Done.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

