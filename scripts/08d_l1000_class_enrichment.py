#!/usr/bin/env python3
"""
08d_l1000_class_enrichment.py

Summarize annotated L1000CDS2 reverse hits into:
- per-signature class enrichment-like summary
- cross-signature class overlap
- top compounds per class

Inputs:
  results/l1000cds2/l1000_hits_annotated.tsv  (from 08c)

Example:
  python scripts/08d_l1000_class_enrichment.py \
    --hits results/l1000cds2/l1000_hits_annotated.tsv \
    --out-dir results/l1000cds2 \
    --min-class-count 2 \
    --topk-per-class 5
"""

from __future__ import annotations

import argparse
from pathlib import Path
import pandas as pd


def _safe_float(x, default=0.0):
    try:
        return float(x)
    except Exception:
        return default


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--hits", required=True, help="Annotated hits TSV from 08c")
    ap.add_argument("--out-dir", required=True, help="Output directory")
    ap.add_argument("--min-class-count", type=int, default=2,
                    help="Minimum hits in a class within a signature to include in enrichment table")
    ap.add_argument("--topk-per-class", type=int, default=5,
                    help="How many top compounds to keep per (signature, class)")
    ap.add_argument("--exclude-unmapped", action="store_true", help="Exclude UNMAPPED class from summaries")
    args = ap.parse_args()

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    df = pd.read_csv(args.hits, sep="\t", dtype=str).fillna("")
    required = {"signature_name", "compound_raw", "class", "score"}
    missing = required - set(df.columns)
    if missing:
        raise SystemExit(f"[FATAL] Missing columns in hits file: {sorted(missing)}")

    # types
    df["score"] = df["score"].map(_safe_float)
    df["class"] = df["class"].astype(str).str.strip()
    df.loc[df["class"] == "", "class"] = "UNMAPPED"

    if args.exclude_unmapped:
        df = df[df["class"] != "UNMAPPED"].copy()

    # -------------------------
    # 1) Per-signature class summary ("enrichment-like")
    # -------------------------
    # n_hits: number of rows in that class
    # n_unique_compounds: unique compound_raw in that class
    # frac_hits: n_hits / total_hits_in_signature
    # frac_unique: n_unique_compounds / total_unique_compounds_in_signature
    # score stats
    sig_tot = df.groupby("signature_name").size().rename("sig_total_hits").reset_index()
    sig_tot_u = df.groupby("signature_name")["compound_raw"].nunique().rename("sig_total_unique").reset_index()

    grp = df.groupby(["signature_name", "class"], dropna=False)
    summ = grp.agg(
        n_hits=("compound_raw", "size"),
        n_unique_compounds=("compound_raw", "nunique"),
        score_mean=("score", "mean"),
        score_median=("score", "median"),
        score_max=("score", "max"),
    ).reset_index()

    summ = summ.merge(sig_tot, on="signature_name", how="left")
    summ = summ.merge(sig_tot_u, on="signature_name", how="left")

    summ["frac_hits"] = summ["n_hits"] / summ["sig_total_hits"].clip(lower=1)
    summ["frac_unique_compounds"] = summ["n_unique_compounds"] / summ["sig_total_unique"].clip(lower=1)

    # filter small classes
    summ_f = summ[summ["n_hits"] >= args.min_class_count].copy()

    # rank within signature by score_mean then n_hits
    summ_f["rank_within_signature"] = (
        summ_f.sort_values(["signature_name", "score_mean", "n_hits"], ascending=[True, False, False])
              .groupby("signature_name")
              .cumcount() + 1
    )

    out1 = out_dir / "l1000_class_enrichment.tsv"
    summ_f.sort_values(["signature_name", "rank_within_signature"]).to_csv(out1, sep="\t", index=False)
    print(f"[OK] Wrote: {out1}")

    # -------------------------
    # 2) Cross-signature class overlap
    # -------------------------
    # For each class: in how many signatures it appears; total hits; mean score over all rows
    cross = df.groupby(["class"]).agg(
        n_signatures=("signature_name", "nunique"),
        total_hits=("compound_raw", "size"),
        total_unique_compounds=("compound_raw", "nunique"),
        score_mean=("score", "mean"),
        score_median=("score", "median"),
        score_max=("score", "max"),
    ).reset_index()

    # Add which signatures contain it (useful for narrative)
    sigs = (df.groupby("class")["signature_name"]
              .apply(lambda s: ",".join(sorted(set(s))))
              .rename("signatures")
              .reset_index())
    cross = cross.merge(sigs, on="class", how="left")
    cross = cross.sort_values(["n_signatures", "score_mean", "total_hits"], ascending=[False, False, False])

    out2 = out_dir / "l1000_class_cross_signature.tsv"
    cross.to_csv(out2, sep="\t", index=False)
    print(f"[OK] Wrote: {out2}")

    # -------------------------
    # 3) Top compounds per class per signature
    # -------------------------
    # keep topK compounds (by best score) per (signature, class)
    # de-duplicate compound_raw within that group, keep best score row.
    df2 = df.sort_values(["signature_name", "class", "score"], ascending=[True, True, False]).copy()
    df2 = df2.drop_duplicates(subset=["signature_name", "class", "compound_raw"], keep="first")

    top = (df2.groupby(["signature_name", "class"], dropna=False)
             .head(args.topk_per_class)
             .copy())

    # include helpful fields if present
    keep_cols = ["signature_name", "class", "compound_raw", "score", "moa", "target", "source"]
    for extra in ["pert_id", "cell_id", "pert_dose", "pert_dose_unit", "pert_time", "pert_time_unit", "sig_id", "pubchem_id", "drugchem_id"]:
        if extra in top.columns:
            keep_cols.append(extra)

    top_out = top[keep_cols].copy()
    out3 = out_dir / "l1000_top_compounds_by_class.tsv"
    top_out.to_csv(out3, sep="\t", index=False)
    print(f"[OK] Wrote: {out3}")

    # Convenience: also write an "unmapped candidates" list to help you extend overrides
    unmapped = df[df["class"] == "UNMAPPED"].copy()
    if len(unmapped):
        cand = (unmapped.groupby("compound_raw")
                       .size()
                       .rename("n_hits")
                       .reset_index()
                       .sort_values(["n_hits", "compound_raw"], ascending=[False, True]))

        for cand in ["moa", "target", "source"]:
            if cand in top_out.columns:
                 top_out[cand] = top_out[cand].astype(str).str.replace("\t", " ", regex=False)

        out4 = out_dir / "l1000_unmapped_compounds.tsv"
        cand.to_csv(out4, sep="\t", index=False)
        print(f"[OK] Wrote: {out4} (candidates to add to overrides)")


if __name__ == "__main__":
    main()

