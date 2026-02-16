#!/usr/bin/env python3
"""
Convenience wrapper to run scripts/08_query_l1000cds2.py with the required out paths.
"""

from __future__ import annotations

import argparse
import subprocess
from pathlib import Path


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--up", required=True)
    ap.add_argument("--down", required=True)
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--tag", required=True)
    ap.add_argument("--db-version", default="latest")
    ap.add_argument("--aggravate", action="store_true")
    ap.add_argument("--timeout", type=int, default=120)
    args = ap.parse_args()

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    out_json = out_dir / f"{args.tag}.l1000cds2.json"
    out_tsv = out_dir / f"{args.tag}.l1000cds2.tsv"

    cmd = [
        "python", "scripts/08_query_l1000cds2.py",
        "--up", args.up,
        "--down", args.down,
        "--tag", args.tag,
        "--db-version", args.db_version,
        "--timeout", str(args.timeout),
        "--out-json", str(out_json),
        "--out-tsv", str(out_tsv),
    ]
    if args.aggravate:
        cmd.append("--aggravate")

    print("Running:", " ".join(cmd))
    subprocess.check_call(cmd)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

