#!/usr/bin/env python3
"""
Process CLUE / LINCS query results (exported CSV/TSV) into reproducible tables.

Assumptions (flexible):
- You export results from CLUE for each query signature into a table.
- Columns vary by export type; this script tries to infer common ones.

Outputs:
- results/lincs/lincs_hits_full.tsv
- results/lincs/lincs_hits_filtered.tsv
- results/lincs/lincs_target_enrichment.tsv
- results/lincs/lincs_class_enrichment.tsv
- results/lincs/lincs_overlap_lcl_vs_primary.tsv

You can run this multiple times as you add more exports; it will concatenate.
"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import pandas as pd


def sniff_delimiter(path: Path) -> str:
    sample = path.read_text(encoding="utf-8")[:4096]
    try:
        dialect = csv.Sniffer().sniff(sample, delimiters=",\t")
        return dialect.delimiter
    except Exception:
        # Default to tab (CLUE often exports TSV)
        return "\t"


def load_table(path: Path) -> pd.DataFrame:
    delim = sniff_delimiter(path)
    df = pd.read_csv(path, sep=delim)
    df.columns = [c.strip() for c in df.columns]
    return df


def first_existing(df: pd.DataFrame, candidates: List[str]) -> Optional[str]:
    cols = {c.lower(): c for c in df.columns}
    for cand in candidates:
        if cand.lower() in cols:
            return cols[cand.lower()]
    return None


def canonicalize(df: pd.DataFrame, source_file: str, signature: str) -> pd.DataFrame:
    """
    Try to map various CLUE column names to a canonical schema.
    """
    # Common columns
    col_name = first_existing(df, ["pert_iname", "compound", "drug", "name"])
    col_moa = first_existing(df, ["moa", "mechanism_of_action", "mechanism"])
    col_target = first_existing(df, ["target", "targets", "primary_target"])
    col_score = first_existing(df, ["score", "wtcs", "cs", "connectivity_score"])
    col_fdr = first_existing(df, ["fdr", "q_value", "qvalue", "adj_p", "padj"])
    col_cell = first_existing(df, ["cell_id", "cell_line", "cell", "cellline"])
    col_time = first_existing(df, ["pert_time", "time", "treatment_time"])
    col_dose = first_existing(df, ["pert_dose", "dose"])

    # Require a name + score
    if col_name is None or col_score is None:
        raise ValueError(
            f"{source_file}: could not find required columns (name + score). "
            f"Have: {list(df.columns)}"
        )

    out = pd.DataFrame()
    out["signature"] = signature
    out["source_file"] = source_file
    out["compound"] = df[col_name].astype(str)

    out["score"] = pd.to_numeric(df[col_score], errors="coerce")
    out["fdr"] = pd.to_numeric(df[col_fdr], errors="coerce") if col_fdr else pd.NA

    out["moa"] = df[col_moa].astype(str) if col_moa else pd.NA
    out["target"] = df[col_target].astype(str) if col_target else pd.NA
    out["cell_line"] = df[col_cell].astype(str) if col_cell else pd.NA
    out["time"] = df[col_time].astype(str) if col_time else pd.NA
    out["dose"] = df[col_dose].astype(str) if col_dose else pd.NA

    # Drop rows without score
    out = out.dropna(subset=["score"])
    return out


def enrich_counts(df: pd.DataFrame, col: str) -> pd.DataFrame:
    tmp = df.dropna(subset=[col]).copy()
    tmp[col] = tmp[col].astype(str).str.strip()
    tmp = tmp[tmp[col] != ""]
    counts = tmp.groupby([col]).size().reset_index(name="n")
    counts = counts.sort_values("n", ascending=False)
    return counts


def overlap_table(df: pd.DataFrame, sig_a: str, sig_b: str, key: str = "compound") -> pd.DataFrame:
    a = set(df.loc[df["signature"] == sig_a, key].astype(str))
    b = set(df.loc[df["signature"] == sig_b, key].astype(str))
    inter = sorted(a.intersection(b))
    out = pd.DataFrame({key: inter})
    out["in_both"] = True
    out["in_" + sig_a] = out[key].isin(a)
    out["in_" + sig_b] = out[key].isin(b)
    return out


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument(
        "--clue-dir",
        default="results/clue_exports",
        help="Directory containing CLUE export files (.csv/.tsv). Default: results/clue_exports",
    )
    ap.add_argument(
        "--out-dir",
        default="results/lincs",
        help="Output directory. Default: results/lincs",
    )
    ap.add_argument(
        "--signature-from-filename",
        action="store_true",
        help="Derive signature name from filename prefix (recommended).",
    )
    ap.add_argument(
        "--sig-a",
        default="GSE267750_primary_Bcell",
        help="Signature A name for overlap table",
    )
    ap.add_argument(
        "--sig-b",
        default="GSE158275_LCL",
        help="Signature B name for overlap table",
    )
    ap.add_argument("--fdr-max", type=float, default=0.05, help="Filter hits by FDR <= this (default 0.05)")
    ap.add_argument("--score-max", type=float, default=-90.0, help="Filter hits by score <= this (default -90)")

    args = ap.parse_args()
    clue_dir = Path(args.clue_dir)
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    files = sorted([p for p in clue_dir.glob("*") if p.suffix.lower() in (".csv", ".tsv", ".txt")])
    if not files:
        raise SystemExit(f"No CLUE export files found in {clue_dir}")

    all_rows: List[pd.DataFrame] = []

    for fp in files:
        df = load_table(fp)
        signature = fp.stem.split("__")[0] if args.signature_from_filename else "UNKNOWN_SIGNATURE"
        canon = canonicalize(df, source_file=fp.name, signature=signature)
        all_rows.append(canon)
        print(f"[OK] Loaded {fp.name}: {len(canon):,} rows (signature={signature})")

    full = pd.concat(all_rows, ignore_index=True)
    full_path = out_dir / "lincs_hits_full.tsv"
    full.to_csv(full_path, sep="\t", index=False)

    # Filter for reversal hits (negative / strong)
    filt = full.copy()
    if "fdr" in filt.columns:
        filt = filt[(filt["fdr"].isna()) | (filt["fdr"] <= args.fdr_max)]
    filt = filt[filt["score"] <= args.score_max]

    filt_path = out_dir / "lincs_hits_filtered.tsv"
    filt.to_csv(filt_path, sep="\t", index=False)

    # Enrich MOA + target (simple frequency; conservative)
    moa_enr = enrich_counts(filt, "moa")
    target_enr = enrich_counts(filt, "target")

    moa_path = out_dir / "lincs_class_enrichment.tsv"
    target_path = out_dir / "lincs_target_enrichment.tsv"
    moa_enr.to_csv(moa_path, sep="\t", index=False)
    target_enr.to_csv(target_path, sep="\t", index=False)

    # Overlap
    overlap = overlap_table(filt, args.sig_a, args.sig_b, key="compound")
    overlap_path = out_dir / "lincs_overlap_lcl_vs_primary.tsv"
    overlap.to_csv(overlap_path, sep="\t", index=False)

    print("[DONE]")
    print(f"  full:     {full_path}")
    print(f"  filtered: {filt_path}")
    print(f"  moa:      {moa_path}")
    print(f"  target:   {target_path}")
    print(f"  overlap:  {overlap_path}")


if __name__ == "__main__":
    main()

