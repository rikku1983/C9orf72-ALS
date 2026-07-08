#!/bin/bash
# Sequentially download each per-sample author Seurat RDS, extract meta.data, delete RDS.
# Keeps peak disk ~1 GB. Metadata CSVs (~few MB each) accumulate in seurat_extract/.
set -u
REPO=/home/sunli/C9orf72-ALS
BASE="https://ftp.ncbi.nlm.nih.gov/geo/series/GSE288nnn/GSE288365/suppl"
RDSDIR="$REPO/data/raw/GSE288365_rds"
OUT="$REPO/data/processed/seurat_extract"
mkdir -p "$RDSDIR" "$OUT"
cd "$RDSDIR"
n=0; ok=0
while read fn; do
  n=$((n+1))
  tag=$(echo "$fn" | sed -E 's/^GSE288365_//; s/_SeuratObject.rds.gz$//')
  metaout="$OUT/ALL_${tag}_meta.csv"
  if [ -f "$metaout" ]; then echo "[$n] SKIP (meta exists) $tag"; ok=$((ok+1)); continue; fi
  echo "[$n] fetching $fn"
  curl -s -o "$fn" "$BASE/$fn"
  if [ ! -s "$fn" ]; then echo "[$n] DOWNLOAD FAILED $fn"; continue; fi
  gunzip -f "$fn"                       # -> .rds (removes .gz)
  rdsf="${fn%.gz}"
  Rscript "$REPO/code/analysis/09b_extract_meta_only.R" "$rdsf" "$metaout" 2>&1 | tail -1
  if [ -f "$metaout" ]; then ok=$((ok+1)); fi
  rm -f "$rdsf" "$fn"                   # reclaim disk
  echo "[$n] done $tag  (disk: $(du -sh $RDSDIR 2>/dev/null | cut -f1))"
done < /home/sunli/C9orf72-ALS/data/metadata/rds_files.txt
echo "ALLDONE extracted=$ok / $n"
