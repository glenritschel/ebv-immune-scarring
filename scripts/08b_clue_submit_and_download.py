#!/usr/bin/env python3
"""
Submit CLUE (Connectivity Map) L1000 gene-set queries via the CLUE Query API and download results.

- Reads a manifest TSV (default: results/lincs/CLUE_submission_manifest.tsv) produced by 07b.
- For each row, submits an L1000 gene-set query job to: https://api.clue.io/api/jobs  (tool_id sig_fastgutc_tool by default)
- Polls job status using: https://api.clue.io/api/jobs/findByJobId/<job_id>  (as documented)
- Downloads the result bundle from download_url (usually an S3 URL) and extracts it.

Auth:
- Requires a CLUE API key in env var: CLUE_USER_KEY
  (You can find it in CLUE -> Account Settings -> API Key)

Note on gene identifiers:
- CLUE docs recommend Entrez gene IDs for queries. This script assumes your harmonized lists are acceptable as-is.
  If you need strict Entrez conversion, we can add a symbol->entrez conversion step using /api/genes. :contentReference[oaicite:1]{index=1}
"""

from __future__ import annotations

import argparse
import csv
import os
import sys
import time
import tarfile
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import requests


API_JOBS = "https://api.clue.io/api/jobs"
API_POLL = "https://api.clue.io/api/jobs/findByJobId/{job_id}"  # documented path :contentReference[oaicite:2]{index=2}


@dataclass
class ManifestRow:
    signature_name: str
    dataset: str
    contrast: str
    up_file: Path
    down_file: Path
    n_up: int
    n_down: int
    notes: str


def read_gene_list(path: Path) -> List[str]:
    genes: List[str] = []
    with path.open("r", encoding="utf-8") as f:
        for line in f:
            g = line.strip()
            if not g:
                continue
            # If user accidentally has extra columns, take first token
            g = g.split()[0]
            genes.append(g)
    # de-dup while preserving order
    seen = set()
    out = []
    for g in genes:
        if g not in seen:
            seen.add(g)
            out.append(g)
    return out


def make_stringified_gmt(tag: str, genes: List[str]) -> str:
    # CLUE examples show "TAG\t\t<ID1>\t<ID2>..." :contentReference[oaicite:3]{index=3}
    # We'll keep the double tab after the set name.
    return f"{tag}\t\t" + "\t".join(genes)


def normalize_url(url: str) -> str:
    # CLUE often returns URLs beginning with //s3.amazonaws.com/... :contentReference[oaicite:4]{index=4}
    if url.startswith("//"):
        return "https:" + url
    if url.startswith("http://"):
        return url.replace("http://", "https://", 1)
    return url


def load_manifest(path: Path) -> List[ManifestRow]:
    rows: List[ManifestRow] = []
    with path.open("r", encoding="utf-8") as f:
        reader = csv.DictReader(f, delimiter="\t")
        required = {"signature_name", "dataset", "contrast", "up_file", "down_file", "n_up", "n_down", "notes"}
        missing = required - set(reader.fieldnames or [])
        if missing:
            raise SystemExit(f"Manifest missing columns: {sorted(missing)}")

        for r in reader:
            rows.append(
                ManifestRow(
                    signature_name=r["signature_name"],
                    dataset=r["dataset"],
                    contrast=r["contrast"],
                    up_file=Path(r["up_file"]),
                    down_file=Path(r["down_file"]),
                    n_up=int(r["n_up"]),
                    n_down=int(r["n_down"]),
                    notes=r.get("notes", "") or "",
                )
            )
    return rows


def clue_submit_job(
    session: requests.Session,
    user_key: str,
    name: str,
    up_gmt: str,
    down_gmt: str,
    tool_id: str,
    dataset: str,
    data_type: str,
    ignore_warnings: bool,
) -> str:
    headers = {
        "Accept": "application/json",
        "Content-Type": "application/json",
        "user_key": user_key,  # documented header :contentReference[oaicite:5]{index=5}
    }
    payload = {
        "tool_id": tool_id,                 # "sig_gutc_tool" or "sig_fastgutc_tool" :contentReference[oaicite:6]{index=6}
        "name": name,
        "uptag-cmapfile": up_gmt,
        "dntag-cmapfile": down_gmt,
        "data_type": data_type,             # default "L1000" :contentReference[oaicite:7]{index=7}
        "dataset": dataset,                 # "Touchstone" :contentReference[oaicite:8]{index=8}
        "ignoreWarnings": ignore_warnings,
    }
    resp = session.post(API_JOBS, headers=headers, json=payload, timeout=120)
    if resp.status_code >= 400:
        raise RuntimeError(f"Submit failed ({resp.status_code}): {resp.text[:800]}")
    j = resp.json()

    # Expected: j["result"]["job_id"] :contentReference[oaicite:9]{index=9}
    job_id = None
    if isinstance(j, dict):
        job_id = (j.get("result") or {}).get("job_id")
    if not job_id:
        raise RuntimeError(f"Submit response missing job_id: {j}")
    return str(job_id)


def clue_poll_job(
    session: requests.Session,
    user_key: str,
    job_id: str,
    poll_seconds: int,
    timeout_seconds: int,
) -> Dict:
    headers = {
        "Accept": "application/json",
        "user_key": user_key,  # documented header :contentReference[oaicite:10]{index=10}
    }
    t0 = time.time()
    while True:
        resp = session.get(API_POLL.format(job_id=job_id), headers=headers, timeout=60)
        if resp.status_code >= 400:
            raise RuntimeError(f"Poll failed ({resp.status_code}): {resp.text[:800]}")
        j = resp.json()

        status = (j.get("status") or j.get("result", {}).get("status") or "").lower()
        download_status = (j.get("download_status") or j.get("result", {}).get("download_status") or "").lower()

        if download_status == "completed" or status == "completed":
            return j

        if time.time() - t0 > timeout_seconds:
            raise TimeoutError(f"Job {job_id} did not complete within {timeout_seconds}s. Last response: {j}")

        time.sleep(poll_seconds)


def download_file(session: requests.Session, url: str, out_path: Path) -> None:
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with session.get(url, stream=True, timeout=300) as r:
        r.raise_for_status()
        with out_path.open("wb") as f:
            for chunk in r.iter_content(chunk_size=1024 * 1024):
                if chunk:
                    f.write(chunk)


def extract_tar_gz(tar_path: Path, out_dir: Path) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)
    with tarfile.open(tar_path, "r:gz") as tf:
        tf.extractall(out_dir)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--manifest", default="results/lincs/CLUE_submission_manifest.tsv")
    ap.add_argument("--out-dir", default="results/clue_downloads")
    ap.add_argument("--tool-id", default="sig_fastgutc_tool")  # faster tool
    ap.add_argument("--dataset", default="Touchstone")
    ap.add_argument("--data-type", default="L1000")
    ap.add_argument("--ignore-warnings", action="store_true", default=True)
    ap.add_argument("--poll-seconds", type=int, default=15)
    ap.add_argument("--timeout-seconds", type=int, default=60 * 30)  # 30 min
    ap.add_argument("--limit", type=int, default=0, help="If >0, only process first N signatures")
    args = ap.parse_args()

    user_key = os.environ.get("CLUE_USER_KEY")
    if not user_key:
        raise SystemExit("Missing env var CLUE_USER_KEY (your CLUE API key).")

    manifest_path = Path(args.manifest)
    rows = load_manifest(manifest_path)
    if args.limit and args.limit > 0:
        rows = rows[: args.limit]

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    session = requests.Session()

    print(f"[INFO] Loaded {len(rows)} signatures from {manifest_path}")
    print(f"[INFO] Submitting to {API_JOBS} using tool_id={args.tool_id}, dataset={args.dataset}, data_type={args.data_type}")

    for row in rows:
        print(f"\n=== {row.signature_name} ===")

        up_genes = read_gene_list(row.up_file)
        down_genes = read_gene_list(row.down_file)

        up_gmt = make_stringified_gmt("UP", up_genes)
        down_gmt = make_stringified_gmt("DOWN", down_genes)

        job_id = clue_submit_job(
            session=session,
            user_key=user_key,
            name=row.signature_name,
            up_gmt=up_gmt,
            down_gmt=down_gmt,
            tool_id=args.tool_id,
            dataset=args.dataset,
            data_type=args.data_type,
            ignore_warnings=args.ignore_warnings,
        )
        print(f"[OK] Submitted job_id={job_id}")

        j = clue_poll_job(
            session=session,
            user_key=user_key,
            job_id=job_id,
            poll_seconds=args.poll_seconds,
            timeout_seconds=args.timeout_seconds,
        )

        # download_url may be at top-level or in result :contentReference[oaicite:11]{index=11}
        download_url = j.get("download_url") or (j.get("result") or {}).get("download_url")
        if not download_url:
            raise RuntimeError(f"Job completed but no download_url found. Response keys: {list(j.keys())}")

        download_url = normalize_url(str(download_url))
        sig_dir = out_dir / row.signature_name
        tar_path = sig_dir / f"{row.signature_name}.tar.gz"

        print(f"[INFO] Downloading: {download_url}")
        download_file(session, download_url, tar_path)
        print(f"[OK] Downloaded: {tar_path}")

        extract_dir = sig_dir / "extracted"
        extract_tar_gz(tar_path, extract_dir)
        print(f"[OK] Extracted -> {extract_dir}")

        # write a tiny provenance file
        (sig_dir / "job_id.txt").write_text(job_id + "\n", encoding="utf-8")
        (sig_dir / "download_url.txt").write_text(download_url + "\n", encoding="utf-8")

    print("\n[DONE] All signatures processed.")


if __name__ == "__main__":
    main()

