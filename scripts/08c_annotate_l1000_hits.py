#!/usr/bin/env python3
"""
Annotate L1000CDS2 top hits with drug classes / MOA.

Inputs:
- One or more L1000CDS2 TSV exports (the "topMeta" parsed TSVs from 08_query_l1000cds2.py)
  Expected columns include:
    rank, score, pert_desc, pert_id, cell_id, pert_dose, pert_dose_unit, pert_time, pert_time_unit, sig_id, pubchem_id

Outputs:
- l1000_hits_annotated.tsv: row-level annotated hits (compound -> class/moa)
- l1000_class_summary.tsv: per-signature class summary (counts + score stats)
- l1000_overlap_by_compound.tsv: overlap of compounds across signatures

Annotation strategy:
1) Optional user override mapping file (TSV) takes precedence.
2) Curated exact-name map for common immunology / signaling drugs.
3) Keyword/rule-based mapping (e.g., "trichostatin" -> HDAC inhibitor).
4) Otherwise class="UNMAPPED".

Override file format (TSV, case-insensitive match on compound_normalized):
  compound    class    moa    target
Example:
  dorsomorphin    TGFb_BMP_pathway    BMP/AMPK inhibitor    ALK2/3/6, AMPK

Usage:
  python scripts/08c_annotate_l1000_hits.py \
    --inputs results/l1000cds2/GSE267750_4_vs_18.reverse.top50.tsv results/l1000cds2/GSE158275_13_vs_2.reverse.top50.tsv \
    --labels GSE267750_primary_Bcell_4_vs_18 GSE158275_LCL_13_vs_2 \
    --out-dir results/l1000cds2

If you used different filenames, just pass them in.
"""

from __future__ import annotations

import argparse
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import pandas as pd


# ---------------------------
# Normalization helpers
# ---------------------------

_WS_RE = re.compile(r"\s+")
_PUNCT_RE = re.compile(r"[^a-z0-9\-\+ ]+")


def normalize_name(s: str) -> str:
    """
    Normalize compound names to improve matching:
    - lower-case
    - remove most punctuation
    - collapse whitespace
    """
    import unicodedata

    if s is None:
        return ""
    # normalize unicode (fix NBSP, odd hyphens, etc.)
    x = unicodedata.normalize("NFKC", str(s))
    x = x.replace("\u00a0", " ")  # NBSP -> space
    x = x.strip().lower()

    x = x.replace("_", " ")

    x = _PUNCT_RE.sub("", x)

    # normalize common separators
    x = x.replace("–", "-").replace("—", "-")
    x = x.replace("  ", " ")

    x = x.strip('"').strip("'")

    # strip common salt / formulation suffixes
    for suffix in [
        " dihydrochloride", " hydrochloride", " hydrobromide", " sulfate",
        " sodium salt", " potassium salt", " phosphate", " acetate",
        " tartrate", " mesylate", " fumarate", " maleate"
    ]:
        if x.endswith(suffix):
            x = x[: -len(suffix)].strip()


    x = _WS_RE.sub(" ", x).strip()
    return x


# ---------------------------
# Curated mappings + rules
# ---------------------------

# Exact-name (normalized) mapping for common compounds
# You can extend this list over time; overrides file is recommended for personalization.
CURATED: Dict[str, Dict[str, str]] = {
    # JAK/STAT
    "ruxolitinib": {"class": "JAK_inhibitor", "moa": "JAK1/2 inhibitor", "target": "JAK1, JAK2"},
    "tofacitinib": {"class": "JAK_inhibitor", "moa": "JAK inhibitor", "target": "JAK1/3 (and others)"},
    "baricitinib": {"class": "JAK_inhibitor", "moa": "JAK1/2 inhibitor", "target": "JAK1, JAK2"},
    "fedratinib": {"class": "JAK_inhibitor", "moa": "JAK2 inhibitor", "target": "JAK2"},
    "cucurbitacin i": {"class": "JAK_STAT_modulator", "moa": "STAT3 pathway inhibitor", "target": "STAT3 (pathway)"},
    "cucurbitacin": {"class": "JAK_STAT_modulator", "moa": "STAT pathway inhibitor", "target": "STAT (pathway)"},

    # NF-kB / inflammatory
    "bms 345541": {"class": "NFkB_inhibitor", "moa": "IKK inhibitor", "target": "IKK"},
    "parthenolide": {"class": "NFkB_inhibitor", "moa": "NF-kB pathway inhibitor", "target": "NF-kB (pathway)"},

    # HDAC / epigenetic
    "trichostatin a": {"class": "HDAC_inhibitor", "moa": "HDAC inhibitor", "target": "HDAC"},
    "vorinostat": {"class": "HDAC_inhibitor", "moa": "HDAC inhibitor", "target": "HDAC"},
    "panobinostat": {"class": "HDAC_inhibitor", "moa": "HDAC inhibitor", "target": "HDAC"},
    "romidepsin": {"class": "HDAC_inhibitor", "moa": "HDAC inhibitor", "target": "HDAC"},
    "entinostat": {"class": "HDAC_inhibitor", "moa": "HDAC inhibitor", "target": "HDAC"},
    "valproic acid": {"class": "HDAC_inhibitor", "moa": "HDAC inhibitor (weak)", "target": "HDAC (weak)"},

    # PI3K/mTOR
    "rapamycin": {"class": "mTOR_inhibitor", "moa": "mTOR inhibitor", "target": "mTOR"},
    "sirolimus": {"class": "mTOR_inhibitor", "moa": "mTOR inhibitor", "target": "mTOR"},
    "everolimus": {"class": "mTOR_inhibitor", "moa": "mTOR inhibitor", "target": "mTOR"},
    "temsirolimus": {"class": "mTOR_inhibitor", "moa": "mTOR inhibitor", "target": "mTOR"},
    "wortmannin": {"class": "PI3K_inhibitor", "moa": "PI3K inhibitor", "target": "PI3K"},
    "ly294002": {"class": "PI3K_inhibitor", "moa": "PI3K inhibitor", "target": "PI3K"},

    # BTK
    "ibrutinib": {"class": "BTK_inhibitor", "moa": "BTK inhibitor", "target": "BTK"},
    "acalabrutinib": {"class": "BTK_inhibitor", "moa": "BTK inhibitor", "target": "BTK"},

    # TGFβ/BMP/SMAD-related
    "dorsomorphin": {"class": "TGFb_BMP_pathway", "moa": "BMP/AMPK inhibitor", "target": "BMP (ALK2/3/6), AMPK"},
    "ldn193189": {"class": "TGFb_BMP_pathway", "moa": "BMP type I receptor inhibitor", "target": "ALK2/3"},
    "sb431542": {"class": "TGFb_BMP_pathway", "moa": "TGFβR1/ALK5 inhibitor", "target": "ALK5"},
    "teniposide": {"class": "Topoisomerase_inhibitor", "moa": "Topoisomerase II inhibitor", "target": "TOP2"},
    "mitoxantrone": {"class": "Topoisomerase_inhibitor", "moa": "Topoisomerase II poison", "target": "TOP2"},
    "bi-2536": {"class": "PLK_inhibitor", "moa": "PLK1 inhibitor", "target": "PLK1"},
    "at-7519": {"class": "CDK_inhibitor", "moa": "CDK inhibitor", "target": "CDK1/2/5/9 (varies)"},
    "tg101348": {"class": "JAK_inhibitor", "moa": "JAK2 inhibitor", "target": "JAK2"},
    "bms-387032": {"class": "JAK_inhibitor", "moa": "JAK inhibitor", "target": "JAK (family)"},
    "gsk-461364": {"class": "PLK_inhibitor", "moa": "PLK1 inhibitor", "target": "PLK1"},
    "cgp-60474": {"class": "CDK_inhibitor", "moa": "CDK inhibitor", "target": "CDK (family)"},
    "pf-431396": {"class": "FAK_SRC_inhibitor", "moa": "FAK inhibitor", "target": "PTK2 (FAK)"},
    "gsk-3-inhibitor-ii": {"class": "GSK3_inhibitor", "moa": "GSK3 inhibitor", "target": "GSK3A/B"},
    "l-sulforophane": {"class": "NRF2_activator", "moa": "NRF2 pathway activator", "target": "NRF2 (pathway)"},

}

# Keyword/rule-based mapping
RULES: List[Tuple[re.Pattern, Dict[str, str]]] = [
    # HDAC
    (re.compile(r"\btrichostatin\b"), {"class": "HDAC_inhibitor", "moa": "HDAC inhibitor", "target": "HDAC"}),
    (re.compile(r"\bhdac\b"), {"class": "HDAC_inhibitor", "moa": "HDAC inhibitor", "target": "HDAC"}),
    (re.compile(r"\b(saha|vorinostat|panobinostat|romidepsin|entinostat)\b"),
     {"class": "HDAC_inhibitor", "moa": "HDAC inhibitor", "target": "HDAC"}),

    # Bromodomain / BET
    (re.compile(r"\b(jq1|i-bet|ibet|bet)\b"), {"class": "BET_inhibitor", "moa": "BET bromodomain inhibitor", "target": "BRD2/3/4"}),

    # Kinase families
    (re.compile(r"\bjak\b"), {"class": "JAK_inhibitor", "moa": "JAK inhibitor", "target": "JAK (family)"}),
    (re.compile(r"\bstat\b"), {"class": "JAK_STAT_modulator", "moa": "STAT pathway modulator", "target": "STAT (pathway)"}),
    (re.compile(r"\bpi3k\b"), {"class": "PI3K_inhibitor", "moa": "PI3K inhibitor", "target": "PI3K"}),
    (re.compile(r"\bmtor\b"), {"class": "mTOR_inhibitor", "moa": "mTOR inhibitor", "target": "mTOR"}),
    (re.compile(r"\bbtk\b"), {"class": "BTK_inhibitor", "moa": "BTK inhibitor", "target": "BTK"}),

    # TGF/BMP
    (re.compile(r"\btgf\b|\balk5\b|\bsmad\b"), {"class": "TGFb_BMP_pathway", "moa": "TGFβ pathway modulator", "target": "TGFβ/SMAD (pathway)"}),
    (re.compile(r"\bbmp\b|\balk2\b|\balk3\b|\balk6\b"), {"class": "TGFb_BMP_pathway", "moa": "BMP pathway modulator", "target": "BMP (pathway)"}),

    # Generic LINCS BRD compounds
    (re.compile(r"\bbrd-[a-z0-9]+\b"), {"class": "LINCS_BRD_unannotated", "moa": "LINCS compound (BRD)", "target": ""}),
    (re.compile(r"\bbi-?2536\b"), {"class": "PLK_inhibitor", "moa": "PLK1 inhibitor", "target": "PLK1"}),
    (re.compile(r"\bat-?7519\b"), {"class": "CDK_inhibitor", "moa": "CDK inhibitor", "target": "CDK (family)"}),
    (re.compile(r"\btg101348\b|\btg-?101348\b"), {"class": "JAK_inhibitor", "moa": "JAK2 inhibitor", "target": "JAK2"}),

]


#vvvv
#def load_overrides(path: Optional[Path]) -> Dict[str, Dict[str, str]]:
#    if not path:
#        return {}
#    if not path.exists():
#        raise SystemExit(f"Override file not found: {path}")
#
#    df = pd.read_csv(path, sep="\t")
#    needed = {"compound", "class"}
#    if not needed.issubset(df.columns):
#        raise SystemExit(f"Override TSV must contain columns: {sorted(needed)} (found {list(df.columns)})")
#
#    out: Dict[str, Dict[str, str]] = {}
#    for _, r in df.iterrows():
#        key = normalize_name(r["compound"])
#        out[key] = {
#            "class": str(r.get("class", "") or ""),
#            "moa": str(r.get("moa", "") or ""),
#            "target": str(r.get("target", "") or ""),
#        }
#    return out
#^^^^

def load_overrides(path: Path) -> dict:
    """
    Load manual compound overrides from a TSV with columns:
    compound, class, moa, target

    Returns dict keyed by normalize_name(compound).
    """
    if path is None:
        return {}

    df = pd.read_csv(path, sep="\t", dtype=str, comment="#").fillna("")

    required = {"compound", "class", "moa", "target"}
    missing = required - set(df.columns)
    if missing:
        raise SystemExit(f"[FATAL] override TSV missing columns: {sorted(missing)}. Found: {df.columns.tolist()}")

    overrides = {}
    for _, r in df.iterrows():
        compound = str(r["compound"]).strip()
        if not compound:
            continue

        key = normalize_name(compound)
        overrides[key] = {
            "class": str(r["class"]).strip() or "UNMAPPED",
            "moa": str(r["moa"]).strip(),
            "target": str(r["target"]).strip(),
            "source": "override",
        }

    return overrides


def annotate_one(compound: str, overrides: Dict[str, Dict[str, str]]) -> Dict[str, str]:
    key = normalize_name(compound)


    if key in overrides:
        return {**overrides[key], "source": "override"}

    if key in CURATED:
        return {**CURATED[key], "source": "curated"}

    for pat, ann in RULES:
        if pat.search(key):
            return {**ann, "source": f"rule:{pat.pattern}"}

    return {"class": "UNMAPPED", "moa": "", "target": "", "source": "unmapped"}


def load_hits(tsv_path: Path, label: str) -> pd.DataFrame:
    df = pd.read_csv(tsv_path, sep="\t")
    # Minimal expected columns
    for c in ["rank", "score", "pert_desc"]:
        if c not in df.columns:
            raise SystemExit(f"{tsv_path}: missing required column '{c}'. Found: {list(df.columns)}")

    df = df.copy()
    df["signature_name"] = label
    df["input_file"] = str(tsv_path)
    df["compound_raw"] = df["pert_desc"].astype(str).str.strip()
    df["compound_normalized"] = df["compound_raw"].map(normalize_name)
    return df


def make_class_summary(df: pd.DataFrame) -> pd.DataFrame:
    # score: higher=stronger reversal in L1000CDS2
    g = df.groupby(["signature_name", "class"], dropna=False)
    out = g.agg(
        n_hits=("compound_raw", "count"),
        n_unique_compounds=("compound_normalized", "nunique"),
        score_mean=("score", "mean"),
        score_median=("score", "median"),
        score_max=("score", "max"),
    ).reset_index()

    # Sort to show most “supported” classes per signature
    out = out.sort_values(["signature_name", "n_hits", "score_mean"], ascending=[True, False, False])
    return out


def make_overlap(df: pd.DataFrame) -> pd.DataFrame:
    # Which compounds/classes appear across signatures?
    # Pivot by compound_normalized
    base = df[["signature_name", "compound_normalized", "compound_raw", "class", "score"]].copy()
    # keep a representative raw name per compound per signature (first)
    base = base.sort_values(["signature_name", "score"], ascending=[True, False])

    # Aggregate per signature + compound
    per_sig = base.groupby(["signature_name", "compound_normalized"], as_index=False).agg(
        compound_raw=("compound_raw", "first"),
        class_=("class", "first"),
        best_score=("score", "max"),
    )

    # Pivot to wide: columns per signature best_score
    wide = per_sig.pivot(index="compound_normalized", columns="signature_name", values="best_score")
    wide = wide.reset_index()

    # Attach a representative name/class from the full set
    rep = base.groupby("compound_normalized", as_index=False).agg(
        representative_name=("compound_raw", "first"),
        representative_class=("class", "first"),
    )

    out = rep.merge(wide, on="compound_normalized", how="left")
    # count how many signatures this compound appears in
    sig_cols = [c for c in out.columns if c not in {"compound_normalized", "representative_name", "representative_class"}]
    out["n_signatures_present"] = out[sig_cols].notna().sum(axis=1)

    out = out.sort_values(["n_signatures_present"], ascending=False)
    return out


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--inputs", nargs="+", required=True, help="One or more L1000CDS2 topMeta TSV files")
    ap.add_argument(
        "--labels",
        nargs="+",
        required=True,
        help="One label per input (e.g., GSE267750_primary_Bcell_4_vs_18)",
    )
    ap.add_argument("--override-tsv", default=None, help="Optional override TSV (compound->class/moa/target)")
    ap.add_argument("--out-dir", default="results/l1000cds2", help="Output directory")
    args = ap.parse_args()

    if len(args.inputs) != len(args.labels):
        raise SystemExit("--inputs and --labels must have the same number of items")

    override_path = Path(args.override_tsv) if args.override_tsv else None
    overrides = load_overrides(override_path)

    dfs: List[pd.DataFrame] = []
    for p, label in zip(args.inputs, args.labels):
        dfs.append(load_hits(Path(p), label))

    df = pd.concat(dfs, ignore_index=True)

    # Drop obvious junk / sentinel rows
    df["compound_raw"] = df["pert_desc"].astype(str).str.strip()

    # common bad values in L1000CDS2 exports (seen in the wild)
    BAD_COMPOUNDS = {"-666", "-1", "nan", "none", "null", ""}
    df = df[~df["compound_raw"].str.lower().isin(BAD_COMPOUNDS)].copy()

    # Annotate (ensure row alignment)
    df = df.reset_index(drop=True)

    anns = df["compound_raw"].map(lambda x: annotate_one(x, overrides))
    ann_df = pd.DataFrame(list(anns)).reset_index(drop=True)

    # Assign by position (avoids index alignment issues)
    for c in ["class", "moa", "target", "source"]:
        df[c] = ann_df[c].values

    df["class"] = df["class"].fillna("").astype(str).str.strip()
    df.loc[df["class"] == "", "class"] = "UNMAPPED"

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    hits_path = out_dir / "l1000_hits_annotated.tsv"
    class_path = out_dir / "l1000_class_summary.tsv"
    overlap_path = out_dir / "l1000_overlap_by_compound.tsv"

    # Stable, paper-friendly column order
    cols = [
        "signature_name",
        "rank",
        "score",
        "compound_raw",
        "class",
        "moa",
        "target",
        "source",
        "pert_id",
        "cell_id",
        "pert_dose",
        "pert_dose_unit",
        "pert_time",
        "pert_time_unit",
        "sig_id",
        "pubchem_id",
        "input_file",
    ]
    cols = [c for c in cols if c in df.columns] + [c for c in df.columns if c not in cols]
    
    # Sanitize text fields to avoid embedded tabs breaking TSV output
    for cand in ["class", "moa", "target", "source"]:
        if cand in df.columns:
            df[cand] = (
                 df[cand]
                 .astype(str)
                 .str.replace("\t", " ", regex=False)
                .str.strip()
             )


    df.to_csv(hits_path, sep="\t", index=False, columns=cols)

    class_summary = make_class_summary(df)
    class_summary.to_csv(class_path, sep="\t", index=False)

    overlap = make_overlap(df)
    overlap.to_csv(overlap_path, sep="\t", index=False)

    # Console summary
    print(f"[OK] Wrote: {hits_path}")
    print(f"[OK] Wrote: {class_path}")
    print(f"[OK] Wrote: {overlap_path}")
    unmapped = (df["class"] == "UNMAPPED").sum()
    print(f"[INFO] Total hits: {len(df)} | UNMAPPED: {unmapped} ({unmapped/len(df):.1%})")
    top_unmapped = (
        df.loc[df["class"] == "UNMAPPED", "compound_raw"]
        .value_counts()
        .head(15)
        .to_string()
    )
    if unmapped:
        print("[INFO] Top unmapped (by frequency):")
        print(top_unmapped)

    # Write an override template for unmapped compounds
    template_path = out_dir / "drug_overrides_template.tsv"
    unmapped_names = (
        df.loc[df["class"] == "UNMAPPED", "compound_raw"]
          .drop_duplicates()
          .sort_values()
    )
    tmpl = pd.DataFrame({"compound": unmapped_names, "class": "", "moa": "", "target": ""})
    tmpl.to_csv(template_path, sep="\t", index=False)
    print(f"[OK] Wrote override template: {template_path}")


if __name__ == "__main__":
    main()

