# scripts/01_fetch_data.py
from __future__ import annotations

import argparse
from pathlib import Path

import yaml

from src.geo_download import download_suppl_files, write_manifest


def main() -> int:
    ap = argparse.ArgumentParser(description="Download GEO supplementary files for configured datasets.")
    ap.add_argument("--config", default="configs/datasets.yaml", help="Path to datasets YAML")
    ap.add_argument("--dataset", action="append", help="Dataset key to download (repeatable). If omitted, downloads all.")
    args = ap.parse_args()

    cfg_path = Path(args.config)
    if not cfg_path.exists():
        raise FileNotFoundError(f"Config not found: {cfg_path}")

    cfg = yaml.safe_load(cfg_path.read_text())
    datasets = cfg.get("datasets", {})
    if not datasets:
        raise ValueError("No datasets found under `datasets:` in config.")

    selected_keys = args.dataset or list(datasets.keys())

    for key in selected_keys:
        if key not in datasets:
            raise KeyError(f"Dataset {key!r} not found in config. Available: {sorted(datasets.keys())}")

        d = datasets[key]
        gse = d["geo_accession"]
        raw_dir = Path(d["paths"]["raw_dir"])
        inc = d.get("download", {}).get("include_patterns", []) or []
        exc = d.get("download", {}).get("exclude_patterns", []) or []

        print(f"\n=== {key} ({gse}) ===")
        print(f"Downloading to: {raw_dir}")
        print(f"include_patterns={inc} exclude_patterns={exc}")

        downloaded, skipped = download_suppl_files(
            gse=gse,
            out_dir=raw_dir,
            include_patterns=inc,
            exclude_patterns=exc,
        )
        manifest_path = write_manifest(gse, raw_dir, downloaded, skipped, inc, exc)

        print(f"Downloaded: {len(downloaded)} files")
        if skipped:
            print(f"Skipped (already present): {len(skipped)} files")
        print(f"Manifest: {manifest_path}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())

