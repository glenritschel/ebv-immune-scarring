#!/usr/bin/env bash
set -euo pipefail

IN_H5AD="data/processed/GSE195452.scvi_scored.h5ad"
OUTDIR="results/paper2_cluster7_vs_others"
TAG="GSE195452_cluster7_vs_others"
CLUSTER="7"

LINCS_UP="${OUTDIR}/lincs/${TAG}_UP150.txt"
LINCS_DN="${OUTDIR}/lincs/${TAG}_DN150.txt"

UMAP_H5AD="data/processed/GSE195452.scvi_scored_umap.h5ad"
FIGDIR="${OUTDIR}/figures"
TABDIR="${OUTDIR}/tables"
META_TSV="${TABDIR}/gsm_metadata.tsv"
FRAC_TSV="${TABDIR}/cluster7_fraction_by_gsm.tsv"
L1000_DIR="${OUTDIR}/l1000cds2"

mkdir -p "${OUTDIR}" "${FIGDIR}" "${TABDIR}" "${L1000_DIR}"

echo "[1/6] DE + LINCS lists"
python scripts/06c_de_cluster7_vs_others.py \
  --in-h5ad "${IN_H5AD}" \
  --outdir "${OUTDIR}" \
  --groupby "leiden_scvi" \
  --cluster "${CLUSTER}" \
  --tag "${TAG}" \
  --topn 150

echo "[2/6] L1000CDS2 query (reverse mode)"
python scripts/08e_run_l1000cds2_for_cluster7.py \
  --up "${LINCS_UP}" \
  --down "${LINCS_DN}" \
  --out-dir "${L1000_DIR}" \
  --tag "${TAG}" \
  --db-version latest

echo "[3/6] UMAP + standard plots"
python scripts/09b_make_umap_and_plots.py \
  --in-h5ad "${IN_H5AD}" \
  --out-h5ad "${UMAP_H5AD}" \
  --figdir "${FIGDIR}" \
  --cluster "${CLUSTER}" \
  --groupby "leiden_scvi" \
  --use-rep "X_scVI"

echo "[4/6] Immune scarring matrixplot"
python scripts/11b_plot_immune_scarring_matrixplot.py \
  --in-h5ad "${UMAP_H5AD}" \
  --out-png "${FIGDIR}/immune_scarring_matrixplot_cluster_leiden_scvi.png" \
  --groupby "leiden_scvi" \
  --layer "log1p"

echo "[5/6] GSM metadata from GEO"
python scripts/13_gsm_metadata_from_geo.py \
  --gse "GSE195452" \
  --destdir "data/raw/GSE195452_meta" \
  --out-tsv "${META_TSV}"

echo "[6/6] Cluster7 fraction per GSM + case/control stats"
python scripts/14_cluster7_fraction_by_gsm.py \
  --in-h5ad "${UMAP_H5AD}" \
  --meta-tsv "${META_TSV}" \
  --out-tsv "${FRAC_TSV}" \
  --groupby "leiden_scvi" \
  --cluster "${CLUSTER}"

echo "DONE"
echo "Key outputs:"
echo "  DE:         ${OUTDIR}/tables/cluster7_vs_others_wilcoxon.tsv"
echo "  LINCS UP/DN: ${LINCS_UP} / ${LINCS_DN}"
echo "  L1000CDS2:   ${L1000_DIR}/${TAG}.l1000cds2.tsv"
echo "  Figures:     ${FIGDIR}/umap_*.png, ${FIGDIR}/immune_scarring_matrixplot_cluster_leiden_scvi.png"
echo "  GSM meta:    ${META_TSV}"
echo "  GSM frac:    ${FRAC_TSV}"

