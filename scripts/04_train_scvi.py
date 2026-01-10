# scripts/04_train_scvi.py
from __future__ import annotations

import argparse
from pathlib import Path

import warnings
warnings.filterwarnings("ignore", message="pkg_resources is deprecated as an API.*")


import scanpy as sc
import scvi


def train_scvi_for_dataset(
    gse: str,
    processed_root: Path,
    model_root: Path,
    out_root: Path,
    n_latent: int,
    max_epochs: int,
    batch_key: str,
    n_neighbors: int,
    leiden_resolution: float,
) -> Path:
    gse = gse.strip().upper()
    in_path = processed_root / f"{gse}.qc.h5ad"
    if not in_path.exists():
        raise FileNotFoundError(f"Missing input: {in_path}")

    print(f"\n=== scVI {gse} ===")
    print(f"Loading: {in_path}")
    adata = sc.read_h5ad(in_path)

    # Ensure counts layer exists
    if "counts" not in adata.layers:
        raise ValueError(f"{gse}: adata.layers['counts'] missing. QC step should have created it.")

    # Ensure batch key exists
    if batch_key not in adata.obs:
        raise ValueError(f"{gse}: batch_key='{batch_key}' not found in adata.obs columns: {list(adata.obs.columns)[:20]}...")

    # Setup for scVI
    scvi.model.SCVI.setup_anndata(adata, layer="counts", batch_key=batch_key)

    # Train model
    model = scvi.model.SCVI(adata, n_latent=n_latent)
    print(f"Training scVI: n_latent={n_latent} max_epochs={max_epochs} batch_key={batch_key}")
    model.train(max_epochs=max_epochs)

    # Save model
    out_model_dir = model_root / gse / "scvi"
    out_model_dir.mkdir(parents=True, exist_ok=True)
    print(f"Saving model to: {out_model_dir}")
    model.save(str(out_model_dir), overwrite=True)

    # Get latent representation
    print("Computing latent embedding ...")
    adata.obsm["X_scVI"] = model.get_latent_representation()

    # Neighbors + Leiden in latent space
    print("Computing neighbors + Leiden on X_scVI ...")
    sc.pp.neighbors(adata, use_rep="X_scVI", n_neighbors=n_neighbors)

    sc.tl.leiden(
    adata,
    resolution=leiden_resolution,
    key_added="leiden_scvi",
    flavor="igraph",
    n_iterations=2,
    directed=False,
)


    # Save output AnnData
    out_path = out_root / f"{gse}.scvi.h5ad"
    print(f"Writing: {out_path}")
    adata.write_h5ad(out_path)

    return out_path


def main() -> int:
    ap = argparse.ArgumentParser(description="Train scVI per dataset and write .scvi.h5ad outputs.")
    ap.add_argument("--gse", action="append", help="Dataset accession (repeatable). If omitted, run both.")
    ap.add_argument("--processed-root", default="data/processed", help="Where *.qc.h5ad files live")
    ap.add_argument("--out-root", default="data/processed", help="Where to write *.scvi.h5ad")
    ap.add_argument("--model-root", default="results/models", help="Where to save trained scVI models")

    ap.add_argument("--batch-key", default="gsm", help="Column in adata.obs to use as batch key")
    ap.add_argument("--n-latent", type=int, default=20, help="Latent dimension")
    ap.add_argument("--max-epochs", type=int, default=100, help="Training epochs")

    ap.add_argument("--n-neighbors", type=int, default=15, help="Neighbors for graph construction")
    ap.add_argument("--leiden-resolution", type=float, default=0.8, help="Leiden resolution for clustering in latent space")

    args = ap.parse_args()

    datasets = args.gse or ["GSE158275", "GSE267750"]

    processed_root = Path(args.processed_root)
    out_root = Path(args.out_root)
    model_root = Path(args.model_root)

    out_root.mkdir(parents=True, exist_ok=True)
    model_root.mkdir(parents=True, exist_ok=True)

    for gse in datasets:
        train_scvi_for_dataset(
            gse=gse,
            processed_root=processed_root,
            model_root=model_root,
            out_root=out_root,
            n_latent=args.n_latent,
            max_epochs=args.max_epochs,
            batch_key=args.batch_key,
            n_neighbors=args.n_neighbors,
            leiden_resolution=args.leiden_resolution,
        )

    print("\nDone.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

