# src/geo_download.py
from __future__ import annotations

import ftplib
import json
import os
import re
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, List, Optional, Tuple


GEO_FTP_HOST = "ftp.ncbi.nlm.nih.gov"


@dataclass(frozen=True)
class GeoSupplFile:
    filename: str
    size_bytes: Optional[int] = None


def _ensure_dir(path: Path) -> None:
    path.mkdir(parents=True, exist_ok=True)


def geo_series_ftp_dir(gse: str) -> str:
    """
    GEO FTP structure:
      /geo/series/GSE158nnn/GSE158275/suppl/
    """
    m = re.fullmatch(r"GSE(\d+)", gse.strip().upper())
    if not m:
        raise ValueError(f"Invalid GEO Series accession: {gse!r} (expected like 'GSE158275')")
    n = int(m.group(1))
    # group into nn nnn (e.g. 158275 -> 158nnn)
    prefix = f"GSE{n // 1000:03d}nnn"
    return f"/geo/series/{prefix}/{gse.strip().upper()}/suppl"


def list_suppl_files(gse: str, timeout: int = 60) -> List[GeoSupplFile]:
    """
    List supplementary files via FTP LIST. We try to parse file sizes when possible.
    """
    remote_dir = geo_series_ftp_dir(gse)
    files: List[GeoSupplFile] = []

    with ftplib.FTP(GEO_FTP_HOST, timeout=timeout) as ftp:
        ftp.login()  # anonymous
        ftp.cwd(remote_dir)

        # Use MLSD if supported (gives sizes reliably). Fall back to LIST.
        try:
            for name, facts in ftp.mlsd():
                if facts.get("type") != "file":
                    continue
                size = facts.get("size")
                files.append(GeoSupplFile(filename=name, size_bytes=int(size) if size else None))
            return sorted(files, key=lambda f: f.filename)
        except Exception:
            pass

        lines: List[str] = []
        ftp.retrlines("LIST", lines.append)

    # Parse LIST (best effort)
    for line in lines:
        # typical format: "-rw-r--r-- 1 ftp ftp 12345 Jan 01 00:00 filename"
        parts = line.split()
        if len(parts) < 9:
            continue
        name = parts[-1]
        try:
            size = int(parts[4])
        except Exception:
            size = None
        files.append(GeoSupplFile(filename=name, size_bytes=size))

    return sorted(files, key=lambda f: f.filename)


def _matches_patterns(name: str, include: Iterable[str], exclude: Iterable[str]) -> bool:
    lname = name.lower()
    inc = [p.lower() for p in include if p]
    exc = [p.lower() for p in exclude if p]
    if inc and not any(p in lname for p in inc):
        return False
    if exc and any(p in lname for p in exc):
        return False
    return True


def download_suppl_files(
    gse: str,
    out_dir: Path,
    include_patterns: Optional[List[str]] = None,
    exclude_patterns: Optional[List[str]] = None,
    timeout: int = 60,
    retries: int = 3,
    sleep_seconds: float = 1.0,
) -> Tuple[List[Path], List[str]]:
    """
    Download supplementary files for a GEO series accession into out_dir.
    Skips files that already exist with matching size (if size known).
    Returns: (downloaded_paths, skipped_filenames)
    """
    include_patterns = include_patterns or []
    exclude_patterns = exclude_patterns or []
    _ensure_dir(out_dir)

    remote_dir = geo_series_ftp_dir(gse)
    entries = list_suppl_files(gse, timeout=timeout)

    selected = [e for e in entries if _matches_patterns(e.filename, include_patterns, exclude_patterns)]
    if not selected:
        raise RuntimeError(
            f"No supplementary files matched filters for {gse}. "
            f"include={include_patterns} exclude={exclude_patterns}"
        )

    downloaded: List[Path] = []
    skipped: List[str] = []

    with ftplib.FTP(GEO_FTP_HOST, timeout=timeout) as ftp:
        ftp.login()
        ftp.cwd(remote_dir)

        for e in selected:
            dest = out_dir / e.filename

            # Skip if exists and size matches (when size known)
            if dest.exists() and e.size_bytes is not None and dest.stat().st_size == e.size_bytes:
                skipped.append(e.filename)
                continue

            # Download with retry
            ok = False
            last_err = None
            for attempt in range(1, retries + 1):
                try:
                    # Download to temp then rename atomically
                    tmp = dest.with_suffix(dest.suffix + ".part")
                    if tmp.exists():
                        tmp.unlink()

                    with open(tmp, "wb") as f:
                        ftp.retrbinary(f"RETR {e.filename}", f.write, blocksize=1024 * 1024)

                    # Optional size check
                    if e.size_bytes is not None and tmp.stat().st_size != e.size_bytes:
                        raise IOError(
                            f"Size mismatch for {e.filename}: expected {e.size_bytes}, got {tmp.stat().st_size}"
                        )

                    os.replace(tmp, dest)
                    downloaded.append(dest)
                    ok = True
                    break
                except Exception as ex:
                    last_err = ex
                    time.sleep(sleep_seconds)

            if not ok:
                raise RuntimeError(f"Failed downloading {e.filename} after {retries} attempts: {last_err}")

    return downloaded, skipped


def write_manifest(
    gse: str,
    out_dir: Path,
    downloaded: List[Path],
    skipped: List[str],
    include_patterns: List[str],
    exclude_patterns: List[str],
) -> Path:
    manifest = {
        "geo_accession": gse,
        "timestamp_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "out_dir": str(out_dir),
        "include_patterns": include_patterns,
        "exclude_patterns": exclude_patterns,
        "downloaded_files": [p.name for p in downloaded],
        "skipped_files": skipped,
    }
    path = out_dir / "manifest.json"
    path.write_text(json.dumps(manifest, indent=2))
    return path

