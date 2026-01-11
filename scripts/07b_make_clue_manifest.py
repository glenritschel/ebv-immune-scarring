#!/usr/bin/env python3
"""
Create a reproducible CLUE submission manifest from harmonized LINCS gene lists.

Looks in an input directory (default: results/lincs) for:
  - *.LINCS.up.harmonized.txt / *.LINCS.down.harmonized.txt
  - *.LINCS_FIXED.up.harmonized.txt / *.LINCS_FIXED.down.harmonized.txt (if you used that naming)
  - *.LINCS_*.up.harmonized.txt / *.LINCS_*.down.harmonized.txt (older style)

Writes:
  - results/lincs/CLUE_submission_manifest.tsv

The manifest includes:
  signature_name, dataset, contrast, direction_files, n_up, n_down, notes
"""

from __future__ import annotations

import argparse
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Tuple

# Parse dataset/contrast from common filenames
# Examples:
#   GSE267750.DE_leiden_scvi_4_vs_18.LINCS.up.harmonized.txt
#   GSE158275.DE_leiden_scvi_13_vs_2.LINCS.down.harmonized.txt
DATASET_RE = re.compile(r"^(GSE\d+)\.")
CONTRAST_RE = re.compile(r"(?:leiden_scvi_)?(\d+_vs_\d+)")
DE_CONTRAST_RE = re.compile(r"DE_leiden_scvi_(\d+_vs_\d+)")


@dataclass
class Entry:
    signature_name: str
    dataset: str
    contrast: str
    up_file: Path
    down_file: Path
    n_up: int
    n_down: int
    notes: str


def count_lines(p: Path) -> int:
    with p.open("r", encoding="utf-8") as f:
        return sum(1 for _ in f)


def infer_dataset(name: str) -> str:
    m = DATASET_RE.search(name)
    return m.group(1) if m else "UNKNOWN"


def infer_contrast(name: str) -> str:
    m = DE_CONTRAST_RE.search(name)
    if m:
        return m.group(1)
    m = CONTRAST_RE.search(name)
    return m.group(1) if m else "UNKNOWN"


def infer_signature(dataset: str, contrast: str) -> str:
    # You can customize this naming scheme to match your manuscript.
    # Keeps it stable and human-readable.
    if dataset == "GSE267750":
        return f"{dataset}_primary_Bcell_{contrast}"
    if dataset == "GSE158275":
        return f"{dataset}_LCL_{contrast}"
    return f"{dataset}_{contrast}"


def pair_files(files: List[Path]) -> List[Tuple[Path, Path]]:
    """
    Pair up/down files based on shared stem before '.up.' / '.down.'
    """
    up_map: Dict[str, Path] = {}
    down_map: Dict[str, Path] = {}

    for p in files:
        n = p.name
        if ".up." in n:
            key = n.split(".up.")[0]
            up_map[key] = p
        elif ".down." in n:
            key = n.split(".down.")[0]
            down_map[key] = p

    pairs: List[Tuple[Path, Path]] = []
    for key, up in sorted(up_map.items()):
        down = down_map.get(key)
        if down:
            pairs.append((up, down))
    return pairs


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--in-dir", default="results/lincs", help="Directory containing harmonized lists")
    ap.add_argument(
        "--out",
        default="results/lincs/CLUE_submission_manifest.tsv",
        help="Output manifest TSV path",
    )
    ap.add_argument(
        "--min-genes",
        type=int,
        default=25,
        help="Warn if up/down has fewer than this many genes",
    )
    args = ap.parse_args()

    in_dir = Path(args.in_dir)
    out_path = Path(args.out)

    if not in_dir.exists():
        raise SystemExit(f"Input dir not found: {in_dir}")

    # Harmonized files
    files = sorted(in_dir.glob("*.harmonized.txt"))
    if not files:
        raise SystemExit(f"No harmonized files found in {in_dir} (expected *.harmonized.txt)")

    pairs = pair_files(files)
    if not pairs:
        raise SystemExit(f"Found harmonized files but could not pair up/down lists in {in_dir}")

    entries: List[Entry] = []

    for up, down in pairs:
        dataset = infer_dataset(up.name)
        contrast = infer_contrast(up.name)

        n_up = count_lines(up)
        n_down = count_lines(down)

        signature_name = infer_signature(dataset, contrast)

        notes = ""
        if n_up < args.min_genes or n_down < args.min_genes:
            notes = f"WARN: small list (n_up={n_up}, n_down={n_down})"

        entries.append(
            Entry(
                signature_name=signature_name,
                dataset=dataset,
                contrast=contrast,
                up_file=up,
                down_file=down,
                n_up=n_up,
                n_down=n_down,
                notes=notes,
            )
        )

    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w", encoding="utf-8") as f:
        f.write("\t".join(["signature_name", "dataset", "contrast", "up_file", "down_file", "n_up", "n_down", "notes"]) + "\n")
        for e in entries:
            f.write(
                "\t".join(
                    [
                        e.signature_name,
                        e.dataset,
                        e.contrast,
                        str(e.up_file),
                        str(e.down_file),
                        str(e.n_up),
                        str(e.n_down),
                        e.notes,
                    ]
                )
                + "\n"
            )

    print(f"[OK] Wrote: {out_path}")
    print(f"[OK] Signatures: {len(entries)}")
    for e in entries:
        print(f" - {e.signature_name}: up={e.n_up}, down={e.n_down}")


if __name__ == "__main__":
    main()

