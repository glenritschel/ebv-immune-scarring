#!/usr/bin/env bash
set -euo pipefail

# --- Inputs / outputs you may tweak ---
BCELL_H5AD="data/processed/GSE195452.bcells_scvi.h5ad"
META_TSV="results/paper2_cluster7_vs_others/tables/gsm_metadata.tsv"

OUT_ROOT="results/paper2_bcells_scvi"
DE_TSV="${OUT_ROOT}/tables/bcells_scarredCluster2_vs_others_wilcoxon.tsv"

LINCS_DIR="${OUT_ROOT}/lincs"
L1000_DIR="${OUT_ROOT}/l1000cds2"
L1000_TAG="GSE195452_bcells_cluster2"

UP_TXT="${LINCS_DIR}/bcells_scarredCluster2_vs_others_wilcoxon.LINCS.up.txt"
DN_TXT="${LINCS_DIR}/bcells_scarredCluster2_vs_others_wilcoxon.LINCS.down.txt"

L1000_TSV="${L1000_DIR}/${L1000_TAG}.tsv"
L1000_JSON="${L1000_DIR}/${L1000_TAG}.json"
OVERRIDES_TSV="${L1000_DIR}/drug_overrides.tsv"

ANNOT_TSV="${L1000_DIR}/l1000_hits_annotated.tsv"
ENRICH_TSV="${L1000_DIR}/l1000_class_enrichment.tsv"
FIG_DIR="${L1000_DIR}/figures"

FRAC_TSV="${OUT_ROOT}/tables/bcells_cluster2_fraction_by_gsm.tsv"

mkdir -p "${OUT_ROOT}/tables" "${LINCS_DIR}" "${L1000_DIR}" "${FIG_DIR}"

echo "[STEP] 1) Extract B cells"
python scripts/paper2_bcells_01_extract.py

echo "[STEP] 2) scVI train/integrate + cluster"
python scripts/paper2_bcells_02_scvi_cluster.py

echo "[STEP] 3) Score scarring signatures + plots"
python scripts/paper2_bcells_03_scarring_scores_and_plots.py

echo "[STEP] 4) DE: scarred cluster (2) vs others"
python scripts/paper2_bcells_04_de_scarred_cluster.py

echo "[STEP] 5) Build LINCS UP/DOWN lists (top 150)"
python scripts/06b_rebuild_lincs_lists_from_de.py \
  --de "${DE_TSV}" \
  --out-dir "${LINCS_DIR}" \
  --max-genes 150

echo "[STEP] 6) Query L1000CDS2"
python scripts/08_query_l1000cds2.py \
  --up "${UP_TXT}" \
  --down "${DN_TXT}" \
  --tag "${L1000_TAG}" \
  --out-json "${L1000_JSON}" \
  --out-tsv "${L1000_TSV}"

echo "[STEP] 7) Annotate hits (uses overrides if present)"
if [[ -f "${OVERRIDES_TSV}" ]]; then
  python scripts/08c_annotate_l1000_hits.py \
    --inputs "${L1000_TSV}" \
    --labels "${L1000_TAG}" \
    --override-tsv "${OVERRIDES_TSV}" \
    --out-dir "${L1000_DIR}"
else
  python scripts/08c_annotate_l1000_hits.py \
    --inputs "${L1000_TSV}" \
    --labels "${L1000_TAG}" \
    --out-dir "${L1000_DIR}"
fi

echo "[STEP] 8) Class enrichment"
python scripts/08d_l1000_class_enrichment.py \
  --hits "${ANNOT_TSV}" \
  --out-dir "${L1000_DIR}" \
  --min-class-count 2 \
  --topk-per-class 5 \
  --exclude-unmapped

echo "[STEP] 9) Plot class figures (score_mean)"
python scripts/12_plot_l1000_class_figures.py \
  --enrichment "${ENRICH_TSV}" \
  --out-dir "${FIG_DIR}" \
  --metric score_mean \
  --topn 15 \
  --min-signatures 1

echo "[STEP] 10) Plot class figures (score_max)"
python scripts/12_plot_l1000_class_figures.py \
  --enrichment "${ENRICH_TSV}" \
  --out-dir "${FIG_DIR}" \
  --metric score_max \
  --topn 15 \
  --min-signatures 1

echo "[STEP] 11) Per-GSM fraction of scarred B-cell cluster (cluster 2)"
python scripts/paper2_bcells_05_cluster2_fraction_by_gsm.py \
  --h5ad "${BCELL_H5AD}" \
  --meta "${META_TSV}" \
  --cluster-key leiden_scvi_bcell \
  --scar-cluster 2 \
  --out "${FRAC_TSV}"

echo "[DONE] Outputs:"
echo "  - DE:       ${DE_TSV}"
echo "  - LINCS:    ${UP_TXT}"
echo "  - LINCS:    ${DN_TXT}"
echo "  - L1000:    ${L1000_TSV}"
echo "  - Annot:    ${ANNOT_TSV}"
echo "  - Enrich:   ${ENRICH_TSV}"
echo "  - Figures:  ${FIG_DIR}/"
echo "  - Per-GSM:  ${FRAC_TSV}"

