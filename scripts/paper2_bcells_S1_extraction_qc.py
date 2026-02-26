#!/usr/bin/env python3
"""
Supplementary Figure S1 – B-cell extraction QC

Inputs:
  data/processed/GSE195452.bcells.h5ad

Outputs:
  results/paper2_bcells_scvi/figures/S1_bcells_marker_dotplot.png
  results/paper2_bcells_scvi/figures/S1_bcells_marker_dotplot.pdf
  results/paper2_bcells_scvi/tables/S1_bcells_cellcount_summary.tsv

What it does:
  1) Writes a simple cell-count summary (total cells, genes; plus optional counts by common metadata keys if present).
  2) Produces a marker-expression dotplot for (CD19, MS4A1, CD79A).
     - If cell_type-like columns exist, it groups by that.
     - Else if leiden_scvi_bcell exists, it groups by that.
     - Else it falls back to a single group ("all_cells") so the plot still renders.
"""

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt


DEFAULT_MARKERS = ["CD19", "MS4A1", "CD79A"]


def first_existing(obs, candidates):
    for c in candidates:
        if c in obs.columns:
            return c
    return None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--h5ad", default="data/processed/GSE195452.bcells.h5ad")
    ap.add_argument("--out-fig-dir", default="results/paper2_bcells_scvi/figures")
    ap.add_argument("--out-table-dir", default="results/paper2_bcells_scvi/tables")
    ap.add_argument("--markers", nargs="*", default=DEFAULT_MARKERS)
    ap.add_argument("--layer", default=None, help="Optional layer to use for expression (e.g., 'log1p'). Default: adata.X")
    args = ap.parse_args()

    in_path = Path(args.h5ad)
    if not in_path.exists():
        raise SystemExit(f"Missing input h5ad: {in_path}")

    out_fig_dir = Path(args.out_fig_dir); out_fig_dir.mkdir(parents=True, exist_ok=True)
    out_tbl_dir = Path(args.out_table_dir); out_tbl_dir.mkdir(parents=True, exist_ok=True)

    out_png = out_fig_dir / "S1_bcells_marker_dotplot.png"
    out_pdf = out_fig_dir / "S1_bcells_marker_dotplot.pdf"
    out_counts = out_tbl_dir / "S1_bcells_cellcount_summary.tsv"

    adata = sc.read_h5ad(in_path)

    # -----------------------------
    # 1) Cell-count summary table
    # -----------------------------
    rows = []

    rows.append({"metric": "n_cells", "value": int(adata.n_obs)})
    rows.append({"metric": "n_genes", "value": int(adata.n_vars)})

    # Optional: if common keys exist, add per-category counts
    common_keys = [
        "condition", "disease", "case_control",
        "patient_id", "donor", "sample_id", "gsm",
        "tissue", "selection_marker",
        "cell_type", "celltype", "cell_type_coarse", "cell_type_fine",
        "leiden", "leiden_scvi", "leiden_scvi_bcell",
    ]
    for k in common_keys:
        if k in adata.obs.columns:
            vc = adata.obs[k].astype(str).value_counts(dropna=False)
            for cat, n in vc.items():
                rows.append({"metric": f"count_by:{k}", "value": int(n), "category": cat})

    df_counts = pd.DataFrame(rows)
    # normalize columns
    if "category" not in df_counts.columns:
        df_counts["category"] = ""
    df_counts = df_counts[["metric", "category", "value"]]
    df_counts.to_csv(out_counts, sep="\t", index=False)
    print(f"[OK] Wrote: {out_counts}")

    # -----------------------------
    # 2) Marker expression dotplot
    # -----------------------------
    markers = [m for m in args.markers if m in adata.var_names]
    missing = [m for m in args.markers if m not in adata.var_names]
    if missing:
        print(f"[WARN] Markers not found in var_names (skipping): {missing}")
    if not markers:
        raise SystemExit(f"No requested markers found in var_names. Requested={args.markers}")

    # Choose a grouping column for the dotplot.
    groupby = first_existing(
        adata.obs,
        [
            # best: curated cell type labels
            "cell_type", "celltype", "cell_type_coarse", "cell_type_fine",
            # otherwise: bcell clustering
            "leiden_scvi_bcell",
            # otherwise: generic clusters
            "leiden_scvi", "leiden",
            # otherwise: sample-level grouping (still useful)
            "condition", "tissue", "selection_marker",
        ],
    )

    # If no grouping key exists, create a single group so Scanpy can plot.
    created_tmp_group = False
    if groupby is None:
        adata.obs["__all_cells__"] = "all_cells"
        groupby = "__all_cells__"
        created_tmp_group = True

    # Use an optional layer if provided and exists
    use_layer = None
    if args.layer:
        if args.layer in adata.layers:
            use_layer = args.layer
        else:
            print(f"[WARN] Requested layer '{args.layer}' not present. Using adata.X instead.")

    # Render and save: Scanpy uses current matplotlib figure.
    sc.set_figure_params(dpi=120, frameon=False)
    sc.pl.dotplot(
        adata,
        var_names=markers,
        groupby=groupby,
        layer=use_layer,
        standard_scale="var",   # helps comparability across markers
        dot_min=0.02,
        figsize=(8, 4),
        show=False,
        return_fig=True,
    ).savefig(out_png, dpi=300)

    # Also save PDF (re-render to ensure vector output is clean)
    fig = sc.pl.dotplot(
        adata,
        var_names=markers,
        groupby=groupby,
        layer=use_layer,
        standard_scale="var",
        dot_min=0.02,
        figsize=(8, 4),
        show=False,
        return_fig=True,
    )
    fig.savefig(out_pdf)
    plt.close("all")

    print(f"[OK] Wrote: {out_png}")
    print(f"[OK] Wrote: {out_pdf}")
    print(f"[INFO] groupby={groupby}  n_cells={adata.n_obs}  n_genes={adata.n_vars}")

    # cleanup
    if created_tmp_group:
        del adata.obs["__all_cells__"]


if __name__ == "__main__":
    main()

