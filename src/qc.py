# src/qc.py
from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, List, Tuple

import numpy as np
import pandas as pd
import scanpy as sc


@dataclass(frozen=True)
class QcParams:
    organism: str
    mito_prefixes: List[str]
    min_genes: int
    min_counts: int
    max_pct_mito: float
    min_cells_per_gene: int
    target_sum: float
    log1p: bool


def params_from_config(cfg: Dict, dataset: str) -> QcParams:
    defaults = cfg.get("defaults", {})
    overrides = (cfg.get("overrides", {}) or {}).get(dataset, {}) or {}

    def get(key, default=None):
        return overrides.get(key, defaults.get(key, default))

    return QcParams(
        organism=str(get("organism", "human")),
        mito_prefixes=list(get("mito_prefixes", ["MT-"])),
        min_genes=int(get("min_genes", 200)),
        min_counts=int(get("min_counts", 500)),
        max_pct_mito=float(get("max_pct_mito", 25.0)),
        min_cells_per_gene=int(get("min_cells_per_gene", 3)),
        target_sum=float(get("target_sum", 10000)),
        log1p=bool(get("log1p", True)),
    )


def add_qc_metrics(adata: sc.AnnData, mito_prefixes: List[str]) -> None:
    # mark mito genes
    varnames = adata.var_names.astype(str)
    mito = np.zeros(adata.n_vars, dtype=bool)
    upper = varnames.str.upper()
    for pref in mito_prefixes:
        mito |= upper.str.startswith(pref.upper())
    adata.var["mt"] = mito

    # compute qc metrics
    sc.pp.calculate_qc_metrics(
        adata,
        qc_vars=["mt"],
        percent_top=None,
        log1p=False,
        inplace=True,
    )


def conservative_filter(
    adata: sc.AnnData,
    min_genes: int,
    min_counts: int,
    max_pct_mito: float,
    min_cells_per_gene: int,
) -> Tuple[int, int]:
    """
    Returns: (n_cells_removed, n_genes_removed)
    """
    n_cells_before = adata.n_obs
    n_genes_before = adata.n_vars

    # cells
    keep = (
        (adata.obs["n_genes_by_counts"] >= min_genes)
        & (adata.obs["total_counts"] >= min_counts)
        & (adata.obs["pct_counts_mt"] <= max_pct_mito)
    )
    adata._inplace_subset_obs(keep)

    # genes
    sc.pp.filter_genes(adata, min_cells=min_cells_per_gene)

    n_cells_removed = n_cells_before - adata.n_obs
    n_genes_removed = n_genes_before - adata.n_vars
    return n_cells_removed, n_genes_removed


def normalize_log1p(adata: sc.AnnData, target_sum: float, do_log1p: bool) -> None:
    """
    Store raw counts in layers["counts"], normalize to target_sum, then log1p.
    """
    # Preserve raw counts
    if "counts" not in adata.layers:
        adata.layers["counts"] = adata.X.copy()

    sc.pp.normalize_total(adata, target_sum=target_sum)
    if do_log1p:
        sc.pp.log1p(adata)


def qc_summary_row(dataset: str, adata_before: sc.AnnData, adata_after: sc.AnnData, params: QcParams) -> Dict:
    def safe_median(col: str, ad: sc.AnnData):
        return float(np.median(ad.obs[col].to_numpy())) if col in ad.obs else float("nan")

    return {
        "dataset": dataset,
        "n_cells_before": int(adata_before.n_obs),
        "n_genes_before": int(adata_before.n_vars),
        "n_cells_after": int(adata_after.n_obs),
        "n_genes_after": int(adata_after.n_vars),
        "median_total_counts_before": safe_median("total_counts", adata_before),
        "median_genes_by_counts_before": safe_median("n_genes_by_counts", adata_before),
        "median_pct_mito_before": safe_median("pct_counts_mt", adata_before),
        "median_total_counts_after": safe_median("total_counts", adata_after),
        "median_genes_by_counts_after": safe_median("n_genes_by_counts", adata_after),
        "median_pct_mito_after": safe_median("pct_counts_mt", adata_after),
        "min_genes": params.min_genes,
        "min_counts": params.min_counts,
        "max_pct_mito": params.max_pct_mito,
        "min_cells_per_gene": params.min_cells_per_gene,
        "target_sum": params.target_sum,
        "log1p": params.log1p,
    }

