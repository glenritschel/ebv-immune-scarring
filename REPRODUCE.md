# Reproducibility Guide — EBV Immune Scarring

This repository reproduces the analyses and figures for the manuscript:

**“Transcriptomic Definition of EBV Immune Scarring Reveals Drug-Reversible Signatures Across Immune Compartments”**

Target journal: *Communications Biology* (Nature Portfolio)

---

## 1) What this repo produces

Running the pipeline end-to-end will generate:

- Processed single-cell objects (`.h5ad`) with scVI latent representations and clustering
- Differential expression outputs and signature scoring summaries
- LINCS/L1000 connectivity query outputs and class-level summaries
- Vector-rendered figure PDFs (Figures 1–4)

**Note:** Public scRNA-seq datasets are re-downloaded and re-processed. Results should be qualitatively consistent with the manuscript, but exact numerical reproducibility may vary slightly depending on hardware and nondeterminism in GPU deep learning frameworks.

---

## 2) Inputs (public datasets)

This project uses the following public datasets:

- **GSE267750** — Primary human B cells
- **GSE158275** — EBV-transformed lymphoblastoid cell lines (LCLs)

Datasets are treated as **independent validation contexts** and are **not pooled**.

---

## 3) Environment setup

### Option A — Conda (recommended)
If you provide an `environment.yml`:

```bash
conda env create -f environment.yml
conda activate ebv-immune-scarring

## Pipeline execution order

The full analysis is executed via:

```bash
./run_all.sh

