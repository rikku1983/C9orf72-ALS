#!/usr/bin/env python3
"""
Download GSE288365 — 10x Visium spatial transcriptomics of human cervical
spinal cord (9 control + 4 sporadic ALS + 5 C9orf72 ALS; paired RNA + protein).

Populates data/raw/GSE288365/ with the series supplementary files from GEO FTP.
Data is git-ignored; this script is the reproducible recovery path.

Usage:
    python code/download/download_GSE288365_visium.py            # all files
    python code/download/download_GSE288365_visium.py --matrices  # skip images
"""
import argparse
import re
import sys
import urllib.request
from pathlib import Path

GSE = "GSE288365"
BASE = f"https://ftp.ncbi.nlm.nih.gov/geo/series/{GSE[:-3]}nnn/{GSE}/suppl/"
PROJECT_ROOT = Path(__file__).resolve().parents[2]
OUTDIR = PROJECT_ROOT / "data" / "raw" / GSE

# file suffixes needed for a standard Scanpy/Squidpy Visium load
MATRIX_SUFFIXES = (
    "matrix.mtx.gz", "barcodes.tsv.gz", "features.tsv.gz",
    "tissue_positions.csv.gz", "scalefactors_json.json.gz",
    "tissue_hires_image.png.gz", "tissue_lowres_image.png.gz",
)


def list_remote_files():
    html = urllib.request.urlopen(BASE, timeout=60).read().decode(errors="ignore")
    return sorted(set(re.findall(r'href="([^"?/][^"]+)"', html)))


def download(fname, dest):
    url = BASE + fname
    tmp = dest.with_suffix(dest.suffix + ".part")
    with urllib.request.urlopen(url, timeout=300) as r, open(tmp, "wb") as f:
        total = int(r.headers.get("Content-Length", 0))
        got = 0
        while True:
            chunk = r.read(1 << 20)
            if not chunk:
                break
            f.write(chunk)
            got += len(chunk)
    tmp.rename(dest)
    return got


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--matrices", action="store_true",
                    help="download only count matrices + spatial coords (skip large images)")
    args = ap.parse_args()

    OUTDIR.mkdir(parents=True, exist_ok=True)
    files = list_remote_files()
    if args.matrices:
        # only the flat 10x files needed for a Scanpy/Squidpy load;
        # SeuratObject.rds (200-500 MB each) are NOT used by the python pipeline
        files = [f for f in files if f.endswith(MATRIX_SUFFIXES)]
    print(f"{GSE}: {len(files)} files to fetch into {OUTDIR}")

    for i, fname in enumerate(files, 1):
        dest = OUTDIR / fname
        if dest.exists():
            print(f"[{i}/{len(files)}] skip (exists) {fname}")
            continue
        try:
            n = download(fname, dest)
            print(f"[{i}/{len(files)}] {fname}  ({n/1e6:.1f} MB)")
        except Exception as e:
            print(f"[{i}/{len(files)}] FAILED {fname}: {e}", file=sys.stderr)


if __name__ == "__main__":
    main()
