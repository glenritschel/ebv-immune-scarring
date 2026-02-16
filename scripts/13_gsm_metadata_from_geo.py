#!/usr/bin/env python3
"""
Download/parse GSE family SOFT via GEOparse, create a compact GSM metadata table.

Heuristic condition labeling:
- if patient_id contains 'Ctrl' (case-insensitive) OR equals 'Control' => control
- else if patient_id is non-empty => case
- else => unknown

Writes TSV with columns:
gsm, title, source_name_ch1, tissue, selection_marker, patient_id, condition
"""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd
import GEOparse


def first_or_empty(v):
    if isinstance(v, list) and len(v) > 0:
        return v[0]
    if v is None:
        return ""
    return str(v)


def extract_characteristics(gsm) -> dict:
    # GSM metadata stores "characteristics_ch1" as list of "key: value" strings
    ch = gsm.metadata.get("characteristics_ch1", []) or []
    out = {}
    for item in ch:
        if ":" in item:
            k, val = item.split(":", 1)
            out[k.strip().lower()] = val.strip()
    return out


def label_condition(patient_id: str) -> str:
    pid = (patient_id or "").strip()
    if pid == "":
        return "unknown"
    low = pid.lower()
    if "ctrl" in low or pid.lower() == "control":
        return "control"
    return "case"


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--gse", default="GSE195452")
    ap.add_argument("--destdir", default="data/raw/GSE195452_meta")
    ap.add_argument("--out-tsv", required=True)
    args = ap.parse_args()

    gse = GEOparse.get_GEO(args.gse, destdir=args.destdir)

    rows = []
    for gsm_id, gsm in gse.gsms.items():
        meta = gsm.metadata
        ch = extract_characteristics(gsm)

        title = first_or_empty(meta.get("title"))
        source = first_or_empty(meta.get("source_name_ch1"))
        tissue = ch.get("tissue", source)  # fall back to source_name
        selection_marker = ch.get("selection marker", "")
        patient_id = ch.get("patient id", "")

        rows.append({
            "gsm": gsm_id,
            "title": title,
            "source_name_ch1": source,
            "tissue": tissue,
            "selection_marker": selection_marker,
            "patient_id": patient_id,
            "condition": label_condition(patient_id),
        })

    df = pd.DataFrame(rows).sort_values("gsm")
    out = Path(args.out_tsv)
    out.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(out, sep="\t", index=False)
    print("Wrote:", out)
    print(df["condition"].value_counts(dropna=False))
    print("Example rows:")
    print(df.head(10).to_string(index=False))

    return 0


if __name__ == "__main__":
    raise SystemExit(main())

