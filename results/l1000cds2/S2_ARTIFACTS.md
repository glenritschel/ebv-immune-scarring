cat > results/l1000cds2/S2_ARTIFACTS.md <<'EOF'
# Sprint 2 Frozen Artifacts (Do Not Regenerate)

These files are frozen outputs from Sprint 2 for the EBV Immune Scarring project.
Do not edit manually. If regeneration is required, create a new artifact version
and record new checksums.

## Frozen artifact
- results/l1000cds2/l1000_hits_annotated.tsv

## How to verify
From repo root:

sha256sum -c results/l1000cds2/S2_ARTIFACTS.sha256.tsv

## Provenance
- Produced by: scripts/08c_annotate_l1000_hits.py
- Inputs:
  - results/l1000cds2/GSE267750_4_vs_18.reverse.top50.tsv
  - results/l1000cds2/GSE158275_13_vs_2.reverse.top50.tsv
  - configs/drug_overrides.tsv
- Output sanitized for TSV safety (tabs removed from annotation fields).
EOF

