#!/usr/bin/env Rscript
# Extract author annotations from GSE268995 per-sample Seurat RDS objects.
# For each of 40 samples: download RDS -> pull meta.data + harmony.umap + predicted_ADT -> delete RDS.
# Outputs (merged across all samples):
#   results/tables/GSE268995_author_meta.csv.gz    (cell_id + meta cols + umap coords)
#   results/tables/GSE268995_predicted_ADT.csv.gz  (cell_id + 228 surface proteins)
suppressMessages({library(SeuratObject); library(Matrix); library(data.table)})
options(timeout=900)  # large RDS files; default 60s times out mid-download

REPO <- "/home/sunli/C9orf72-ALS"
RDSDIR <- file.path(REPO, "data/raw/GSE268995_rds")
CACHE <- file.path(RDSDIR, "cache"); dir.create(CACHE, showWarnings=FALSE, recursive=TRUE)
dir.create(RDSDIR, showWarnings=FALSE, recursive=TRUE)
MAP <- fread(file.path(REPO, "data/metadata/GSE268995_sample_condition_map.csv"))
LOG <- file.path(RDSDIR, "extract.log")
logmsg <- function(...) { cat(sprintf(...), "\n", file=LOG, append=TRUE); cat(sprintf(...), "\n") }

# robust download with retries
dl_retry <- function(url, dest, tries=4) {
  for (k in seq_len(tries)) {
    ok <- tryCatch(download.file(url, dest, quiet=TRUE, mode="wb")==0, error=function(e) FALSE)
    if (ok && file.exists(dest) && file.size(dest) > 1e6) return(TRUE)
    Sys.sleep(5*k)
  }
  FALSE
}

# meta columns we keep (annotation + repertoire summary; drop verbose per-chain gene cols to keep file small)
KEEP <- c("orig.ident","sample_id","nCount_RNA","nFeature_RNA","percent.mt",
          "seurat_clusters","harmony_clusters","DF.classifications",
          "predicted.celltype.l1","predicted.celltype.l1.score",
          "predicted.celltype.l2","predicted.celltype.l2.score",
          "diagnosis","diagnosis_general","age","sex","race",
          "tcr_clonotype_id","tcr_frequency","tra_cdr3s","trb_cdr3s",
          "inkt_evidence","mait_evidence",
          "bcr_clonotype_id","bcr_frequency","igh_cdr3s","igl_cdr3s")

for (i in seq_len(nrow(MAP))) {
  s <- MAP$sample[i]; gsm <- MAP$gsm[i]; cond <- MAP$condition[i]
  meta_cache <- file.path(CACHE, sprintf("%s_meta.rds", s))
  adt_cache  <- file.path(CACHE, sprintf("%s_adt.rds", s))
  if (file.exists(meta_cache) && file.exists(adt_cache)) { logmsg("SKIP %s (cached)", s); next }
  gnnn <- paste0(substr(gsm,1,7),"nnn")
  url <- sprintf("https://ftp.ncbi.nlm.nih.gov/geo/samples/%s/%s/suppl/%s_%s.rds.gz", gnnn, gsm, gsm, s)
  gz <- file.path(RDSDIR, sprintf("%s_%s.rds.gz", gsm, s))
  rds <- sub("\\.gz$","",gz)
  if (!dl_retry(url, gz)) { logmsg("FAIL download %s after retries", s); next }
  # .rds.gz is double-gzipped; readRDS auto-decompresses only for a FILENAME, not a connection.
  system2("gunzip", c("-kf", gz))   # -> rds (keeps .gz)
  obj <- readRDS(rds)
  md <- obj@meta.data
  bc <- sub(sprintf("^%s_", s), "", rownames(md))   # RDS uses SAMPLE_BARCODE -> strip prefix
  cell_id <- paste0(s, "::", bc)                     # match our AnnData SAMPLE::BARCODE
  cols_present <- intersect(KEEP, colnames(md))
  m <- data.table(cell_id=cell_id, sample=s, condition=cond, md[, cols_present, drop=FALSE])
  if ("harmony.umap" %in% names(obj@reductions)) {
    um <- obj@reductions$harmony.umap@cell.embeddings
    m[, `:=`(umap1=um[,1], umap2=um[,2])]
  }
  saveRDS(m, meta_cache)
  has_adt <- "predicted_ADT" %in% names(obj@assays)
  if (has_adt) {
    adt <- as.matrix(obj@assays$predicted_ADT@data)  # proteins x cells (direct slot; avoids Seurat load)
    saveRDS(data.table(cell_id=cell_id, t(adt)), adt_cache)
  } else saveRDS(data.table(cell_id=cell_id), adt_cache)
  logmsg("OK %s (%s): %d cells, %d meta cols, ADT=%s", s, cond, length(cell_id),
         length(cols_present), ifelse(has_adt,"yes","no"))
  rm(obj, md); if (exists("adt")) rm(adt); gc()
  file.remove(rds); if (file.exists(gz)) file.remove(gz)  # conserve disk (re-downloadable from GEO)
}

# merge all cached per-sample outputs
meta_files <- file.path(CACHE, sprintf("%s_meta.rds", MAP$sample))
meta_files <- meta_files[file.exists(meta_files)]
META <- rbindlist(lapply(meta_files, readRDS), fill=TRUE, use.names=TRUE)
fwrite(META, file.path(REPO,"results/tables/GSE268995_author_meta.csv.gz"))
logmsg("wrote author_meta: %d cells x %d cols (%d samples)", nrow(META), ncol(META), length(meta_files))

adt_files <- file.path(CACHE, sprintf("%s_adt.rds", MAP$sample))
adt_files <- adt_files[file.exists(adt_files)]
ADT <- rbindlist(lapply(adt_files, readRDS), fill=TRUE, use.names=TRUE)
fwrite(ADT, file.path(REPO,"results/tables/GSE268995_predicted_ADT.csv.gz"))
logmsg("wrote predicted_ADT: %d cells x %d proteins", nrow(ADT), ncol(ADT)-1)
logmsg("ALL DONE")
