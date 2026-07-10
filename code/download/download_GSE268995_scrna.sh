#!/usr/bin/env bash
# Download GSE268995 (PBMC scRNA-seq, 40 samples) count triplets only (skip 38GB RDS tar).
set -uo pipefail
REPO=/home/sunli/C9orf72-ALS
OUT=$REPO/data/raw/GSE268995
MAP=$REPO/data/metadata/GSE268995_sample_condition_map.csv
mkdir -p "$OUT"
LOG=$OUT/download.log; : > "$LOG"
tail -n +2 "$MAP" | while IFS=, read -r sample gsm disease condition; do
  pfx="${gsm}_${sample}"
  # GSM nnn dir: first 7 chars + nnn
  stem=$(echo "$gsm" | sed 's/...$//')   # GSM8304 -> GSM8304? no; need GSM8304nnn
  gnnn="${gsm:0:7}nnn"                    # GSM8304338 -> GSM8304nnn
  base="https://ftp.ncbi.nlm.nih.gov/geo/samples/${gnnn}/${gsm}/suppl"
  for suf in matrix.mtx.gz barcodes.tsv.gz features.tsv.gz; do
    url="${base}/${pfx}${suf}"
    dst="$OUT/${pfx}${suf}"
    if [[ -s "$dst" ]]; then echo "skip $dst" >>"$LOG"; continue; fi
    curl -sfL "$url" -o "$dst" && echo "OK  $pfx$suf $(stat -c%s "$dst")" >>"$LOG" || echo "FAIL $url" >>"$LOG"
  done
  echo "done $sample ($condition)" >>"$LOG"
done
echo "ALL DONE" >>"$LOG"
