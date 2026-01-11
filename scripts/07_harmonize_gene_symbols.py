#!/usr/bin/env python3
"""
Harmonize gene symbols for LINCS/CLUE submission.

- Reads *.LINCS_*.up.txt and *.LINCS_*.down.txt (one gene symbol per line)
- Normalizes common scRNA-seq symbol artifacts (case, whitespace, version suffixes)
- Optionally applies an alias->approved mapping file if provided
- Writes harmonized lists + a mapping report

Design note:
  This script does not call external services (offline/reproducible).
  If you later add a curated mapping table, results become stronger without changing code.
"""

from __future__ import annotations

import argparse
import csv
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Tuple


VERSION_SUFFIX_RE = re.compile(r"\.\d+$")
ENSEMBL_LIKE_RE = re.compile(r"^ENS[A-Z]{0,3}G\d{6,}$", re.IGNORECASE)


@dataclass(frozen=True)
class MapStats:
    total_in: int
    total_out: int
    dropped_blank: int
    dropped_ensembl: int
    mapped_alias: int
    unchanged: int
    duplicates_removed: int


def read_gene_list(path: Path) -> List[str]:
    """
    Robustly read gene lists that may include extra columns:
      GENE<TAB>score
      GENE score
      "GENE",score
    We keep only the first token (tab/space/comma separated).
    """
    genes: List[str] = []
    with path.open("r", encoding="utf-8") as f:
        for line in f:
            s = line.strip()
            if not s:
                genes.append("")
                continue
            parts = re.split(r"[\t, ]+", s)
            genes.append(parts[0] if parts else "")
    return genes



def normalize_symbol(sym: str) -> str:
    s = sym.strip()
    s = VERSION_SUFFIX_RE.sub("", s)  # e.g. "STAT1.1" -> "STAT1"
    s = s.replace('"', "").replace("'", "")
    s = s.replace("\ufeff", "")  # BOM
    # Common scRNA artifacts: "HLA-DRA " or " hla-dra"
    s = s.upper()
    # Sometimes gene lists include "GENE:..." or "gene=..."
    for prefix in ("GENE:", "GENE=", "SYMBOL:", "SYMBOL="):
        if s.startswith(prefix):
            s = s[len(prefix):].strip()
    return s


def load_alias_map(path: Path) -> Dict[str, str]:
    """
    Expects TSV/CSV with columns:
      alias, approved

    Accepts header variants; ignores extra columns.
    """
    mapping: Dict[str, str] = {}
    with path.open("r", encoding="utf-8") as f:
        sample = f.read(4096)
        f.seek(0)
        dialect = csv.Sniffer().sniff(sample, delimiters=",\t")
        reader = csv.DictReader(f, dialect=dialect)
        # Find likely column names
        cols = {c.lower(): c for c in reader.fieldnames or []}
        alias_col = cols.get("alias") or cols.get("prev_symbol") or cols.get("previous_symbol") or cols.get("synonym")
        appr_col = cols.get("approved") or cols.get("approved_symbol") or cols.get("symbol")
        if not alias_col or not appr_col:
            raise ValueError(
                f"Alias map {path} must have columns like alias + approved (got: {reader.fieldnames})"
            )
        for row in reader:
            alias = normalize_symbol(row[alias_col])
            appr = normalize_symbol(row[appr_col])
            if alias and appr:
                mapping[alias] = appr
    return mapping


def harmonize(
    genes_in: Iterable[str],
    alias_map: Optional[Dict[str, str]] = None,
) -> Tuple[List[str], List[Tuple[str, str, str]] , MapStats]:
    """
    Returns:
      genes_out: unique, ordered
      report_rows: (original, normalized, final)
      stats
    """
    alias_map = alias_map or {}
    report_rows: List[Tuple[str, str, str]] = []

    out: List[str] = []
    seen = set()

    dropped_blank = 0
    dropped_ensembl = 0
    mapped_alias = 0
    unchanged = 0
    dup_removed = 0

    total_in = 0
    for g in genes_in:
        total_in += 1
        orig = g
        norm = normalize_symbol(g)

        if not norm:
            dropped_blank += 1
            report_rows.append((orig, norm, ""))
            continue

        # If some lists accidentally contain Ensembl IDs, drop (CLUE expects symbols).
        if ENSEMBL_LIKE_RE.match(norm):
            dropped_ensembl += 1
            report_rows.append((orig, norm, ""))
            continue

        final = alias_map.get(norm, norm)
        if final != norm:
            mapped_alias += 1
        else:
            unchanged += 1

        if final in seen:
            dup_removed += 1
            report_rows.append((orig, norm, final))
            continue

        seen.add(final)
        out.append(final)
        report_rows.append((orig, norm, final))

    stats = MapStats(
        total_in=total_in,
        total_out=len(out),
        dropped_blank=dropped_blank,
        dropped_ensembl=dropped_ensembl,
        mapped_alias=mapped_alias,
        unchanged=unchanged,
        duplicates_removed=dup_removed,
    )
    return out, report_rows, stats


def write_list(path: Path, genes: List[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as f:
        for g in genes:
            f.write(f"{g}\n")


def write_report(path: Path, rows: List[Tuple[str, str, str]], stats: MapStats, src: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["source_file", "original", "normalized", "final"])
        for orig, norm, final in rows:
            w.writerow([src, orig, norm, final])

    # Append a summary alongside (same base name)
    summary_path = path.with_suffix(".summary.tsv")
    with summary_path.open("w", encoding="utf-8", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["metric", "value"])
        for k, v in stats.__dict__.items():
            w.writerow([k, v])


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument(
        "--in-dir",
        default="results",
        help="Root to search for *.LINCS_*.up.txt/*.down.txt (default: results)",
    )
    ap.add_argument(
        "--out-dir",
        default="results/lincs",
        help="Output directory (default: results/lincs)",
    )
    ap.add_argument(
        "--alias-map",
        default="",
        help="Optional TSV/CSV alias map with columns alias,approved (e.g., configs/hgnc_aliases.tsv)",
    )
    ap.add_argument(
        "--min-genes",
        type=int,
        default=25,
        help="Warn if harmonized list is smaller than this (default: 25)",
    )

    args = ap.parse_args()
    in_dir = Path(args.in_dir)
    out_dir = Path(args.out_dir)
    alias_map_path = Path(args.alias_map) if args.alias_map else None

    alias_map: Optional[Dict[str, str]] = None
    if alias_map_path:
        alias_map = load_alias_map(alias_map_path)

    lincs_files = (
        sorted(in_dir.rglob("*.LINCS_*.up.txt")) +
        sorted(in_dir.rglob("*.LINCS_*.down.txt")) +
        sorted(in_dir.rglob("*.LINCS.up.txt")) +
        sorted(in_dir.rglob("*.LINCS.down.txt"))
    )

    if not lincs_files:
        raise SystemExit(f"No LINCS gene-list files found under: {in_dir} (expected *.LINCS_*.up.txt/*.down.txt)")

    # Consolidated report across all files
    consolidated_rows: List[Tuple[str, str, str, str]] = []  # (src, orig, norm, final)
    consolidated_summary: List[Tuple[str, int]] = []

    for fp in lincs_files:
        genes_in = read_gene_list(fp)
        genes_out, report_rows, stats = harmonize(genes_in, alias_map=alias_map)

        out_name = fp.name.replace(".txt", ".harmonized.txt")
        out_path = out_dir / out_name
        write_list(out_path, genes_out)

        report_path = out_dir / (fp.name.replace(".txt", ".mapping_report.tsv"))
        write_report(report_path, report_rows, stats, src=fp.name)

        if len(genes_out) < args.min_genes:
            print(f"[WARN] {fp.name}: harmonized list has only {len(genes_out)} genes (<{args.min_genes})")

        consolidated_summary.append((fp.name, stats.total_out))
        for orig, norm, final in report_rows:
            consolidated_rows.append((fp.name, orig, norm, final))

        print(f"[OK] {fp.name} -> {out_path} (out={stats.total_out}, alias_mapped={stats.mapped_alias})")

    # Write consolidated mapping report
    out_dir.mkdir(parents=True, exist_ok=True)
    consolidated_report = out_dir / "gene_mapping_report.tsv"
    with consolidated_report.open("w", encoding="utf-8", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["source_file", "original", "normalized", "final"])
        for src, orig, norm, final in consolidated_rows:
            w.writerow([src, orig, norm, final])

    print(f"[DONE] Consolidated report: {consolidated_report}")


if __name__ == "__main__":
    main()

