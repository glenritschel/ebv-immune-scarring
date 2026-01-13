#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="./scripts"
CONFIG_DIR="./configs"
RESULTS_DIR="./results"
FIGURES_DIR="./figures"
LOG_DIR="./logs"

mkdir -p "$LOG_DIR" "$RESULTS_DIR" "$FIGURES_DIR"

RUN_ID="$(date +%Y%m%d_%H%M%S)"
LOG_FILE="$LOG_DIR/run_all_${RUN_ID}.log"
exec > >(tee -a "$LOG_FILE") 2>&1

echo "=== EBV Immune Scarring — Full Reproduction Run ==="
echo "Run ID: $RUN_ID"
echo "Date:   $(date -Is)"
echo

run () {
  echo
  echo ">>> $*"
  python "$@"
}

# -----------------------------
# DATA INGESTION
# -----------------------------
run scripts/01_fetch_data.py
run scripts/02_build_h5ad_from_geo_mtx.py
run scripts/03_qc_preprocess.py

# -----------------------------
# scVI MODELING
# -----------------------------
run scripts/04_train_scvi.py

# -----------------------------
# SIGNATURE ANALYSIS
# -----------------------------
run scripts/05_score_signatures.py

# -----------------------------
# DIFFERENTIAL EXPRESSION
# -----------------------------
run scripts/06_de_for_lincs.py
run scripts/06b_rebuild_lincs_lists_from_de.py

# -----------------------------
# LINCS / L1000 PREP
# -----------------------------
run scripts/07_harmonize_gene_symbols.py
run scripts/07b_make_clue_manifest.py

# -----------------------------
# LINCS / L1000 QUERY
# -----------------------------
# (May be skipped if cached)
run scripts/08b_clue_submit_and_download.py || true
run scripts/08_query_l1000cds2.py
run scripts/08_process_clue_results.py

# -----------------------------
# L1000 ANNOTATION (FROZEN)
# -----------------------------
ANNOT="results/l1000cds2/l1000_hits_annotated.tsv"
if [[ -f "$ANNOT" ]]; then
  echo "Skipping annotation — frozen artifact exists."
else
  run scripts/08c_annotate_l1000_hits.py
fi

run scripts/08d_l1000_class_enrichment.py

# -----------------------------
# FIGURE GENERATION
# -----------------------------
run scripts/09_plot_scvi_umap.py
run scripts/10_plot_signature_overlays.py
run scripts/11_plot_cluster_signature_heatmap.py
run scripts/12_plot_l1000_class_figures.py

echo
echo "=== PIPELINE COMPLETE ==="
echo "Log saved to: $LOG_FILE"

