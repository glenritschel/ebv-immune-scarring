#!/usr/bin/env python3
"""
Query L1000CDS2 (Ma'ayan Lab) with up/down gene sets to identify mimickers/reversers.

API docs + payload format:
  POST https://maayanlab.cloud/L1000CDS2/query
  - data.upGenes / data.dnGenes: arrays of gene symbols
  - config.searchMethod: "geneSet"
  - config.aggravate: False => reverse mode, True => aggravate mode
  - config['db-version']: e.g. "latest"
Docs: https://maayanlab.cloud/L1000CDS2/help/  (API section)
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import List, Dict, Any, Optional

import requests


API_URL = "https://maayanlab.cloud/L1000CDS2/query"


def read_genes(path: Path) -> List[str]:
    genes: List[str] = []
    with path.open("r", encoding="utf-8") as f:
        for line in f:
            g = line.strip()
            if not g:
                continue
            # keep only first token if accidental whitespace/columns exist
            genes.append(g.split()[0])
    return genes


def upper_genes(genes: List[str]) -> List[str]:
    return [g.upper() for g in genes]


def write_tsv_topmeta(topmeta: List[Dict[str, Any]], out_tsv: Path) -> None:
    # Flatten to stable columns (some may be missing depending on hit)
    cols = [
        "rank",
        "score",
        "pert_desc",
        "pert_id",
        "cell_id",
        "pert_dose",
        "pert_dose_unit",
        "pert_time",
        "pert_time_unit",
        "sig_id",
        "pubchem_id",
        "drugchem_id",
    ]
    out_tsv.parent.mkdir(parents=True, exist_ok=True)
    with out_tsv.open("w", encoding="utf-8") as f:
        f.write("\t".join(cols) + "\n")
        for i, row in enumerate(topmeta, start=1):
            flat = {**row}
            flat["rank"] = i
            f.write("\t".join(str(flat.get(c, "")) for c in cols) + "\n")


def main() -> int:
    ap = argparse.ArgumentParser(description="Query L1000CDS2 with up/down gene sets")
    ap.add_argument("--up", required=True, help="Path to UP genes file (one gene per line)")
    ap.add_argument("--down", required=True, help="Path to DOWN genes file (one gene per line)")
    ap.add_argument("--tag", default="EBV_Immune_Scarring", help="Tag stored in metadata")
    ap.add_argument("--cell", default=None, help="Optional cell line metadata (e.g., MCF7)")
    ap.add_argument("--db-version", default="latest",
                    help="Database version (default: latest). Docs mention e.g. 'cpcd-gse70138-v1.0' or 'cpcd-v1.0'")
    ap.add_argument("--aggravate", action="store_true",
                    help="If set: aggravate mode (mimic/aggravate). Default is reverse mode (find reversers).")
    ap.add_argument("--combination", action="store_true",
                    help="If set: also ask for drug combinations.")
    ap.add_argument("--share", action="store_true",
                    help="If set: allow service to generate a shareId (public sharing). Default off.")
    ap.add_argument("--timeout", type=int, default=120, help="HTTP timeout seconds")
    ap.add_argument("--out-json", required=True, help="Where to write raw JSON response")
    ap.add_argument("--out-tsv", required=True, help="Where to write parsed topMeta TSV")
    args = ap.parse_args()

    up_path = Path(args.up)
    dn_path = Path(args.down)
    out_json = Path(args.out_json)
    out_tsv = Path(args.out_tsv)

    up = upper_genes(read_genes(up_path))
    dn = upper_genes(read_genes(dn_path))

    if len(up) < 5 or len(dn) < 5:
        print(f"[WARN] Very small gene sets: up={len(up)} down={len(dn)}", file=sys.stderr)

    data: Dict[str, Any] = {"upGenes": up, "dnGenes": dn}
    config: Dict[str, Any] = {
        "aggravate": bool(args.aggravate),
        "searchMethod": "geneSet",
        "share": bool(args.share),
        "combination": bool(args.combination),
        "db-version": args.db_version,
    }
    meta = [{"key": "Tag", "value": args.tag}]
    if args.cell:
        meta.append({"key": "Cell", "value": args.cell})

    payload: Dict[str, Any] = {"data": data, "config": config, "meta": meta}
    headers = {"content-type": "application/json"}

    out_json.parent.mkdir(parents=True, exist_ok=True)
    out_tsv.parent.mkdir(parents=True, exist_ok=True)

    r = requests.post(API_URL, data=json.dumps(payload), headers=headers, timeout=args.timeout)
    r.raise_for_status()
    res = r.json()

    with out_json.open("w", encoding="utf-8") as f:
        json.dump(res, f, indent=2, sort_keys=True)

    topmeta = res.get("topMeta", [])
    if not isinstance(topmeta, list) or len(topmeta) == 0:
        print("[WARN] Response contained no topMeta results.", file=sys.stderr)
        # still exit 0 because JSON is saved; you can inspect errors
        return 0

    write_tsv_topmeta(topmeta, out_tsv)

    mode = "AGGRAVATE" if args.aggravate else "REVERSE"
    print(f"[OK] mode={mode} up={len(up)} down={len(dn)}")
    print(f"[OK] wrote JSON: {out_json}")
    print(f"[OK] wrote TSV:  {out_tsv} (n={len(topmeta)})")
    if "shareId" in res:
        print(f"[OK] shareId: {res.get('shareId')}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())

