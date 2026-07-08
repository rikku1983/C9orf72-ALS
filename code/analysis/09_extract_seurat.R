#!/usr/bin/env Rscript
# 09_extract_seurat.R — extract meta.data + reductions + key-gene expression
# from author Seurat .rds objects to CSV, for validation against our AnnData.
#
# Usage: Rscript 09_extract_seurat.R <path_to.rds.gz> <out_prefix>
# Writes: <out_prefix>_meta.csv, <out_prefix>_genes.csv, <out_prefix>_struct.txt

suppressMessages({library(Matrix)})
args <- commandArgs(trailingOnly=TRUE)
rds <- args[1]; outp <- args[2]
KEY_GENES <- c("GPNMB","CHIT1","C1QA","C1QB","C1QC","CHAT","GFAP","NEFH",
               "AIF1","CD68","MSR1","TREM2","TYROBP","APOE","P2RY12","MBP")

cat("Reading", rds, "\n")
obj <- readRDS(rds)   # full Seurat S4 object; access via slots to avoid needing Seurat pkg

struct <- file(paste0(outp,"_struct.txt"), open="wt")
w <- function(...) { cat(..., "\n", file=struct, append=TRUE); cat(..., "\n") }
w("class:", class(obj)[1])
w("assay slots:", paste(names(obj@assays), collapse=", "))
w("active.assay:", obj@active.assay)
w("reduction slots:", paste(names(obj@reductions), collapse=", "))
w("meta.data cols:", paste(colnames(obj@meta.data), collapse=", "))

# --- meta.data ---
md <- obj@meta.data
md$barcode <- rownames(md)
write.csv(md, paste0(outp,"_meta.csv"), row.names=FALSE)
w("meta.data written, n cells:", nrow(md))

# --- reductions (embeddings) via slot ---
for (rd in names(obj@reductions)) {
  emb <- obj@reductions[[rd]]@cell.embeddings
  if (is.null(emb) || nrow(emb)==0) next
  df <- as.data.frame(emb[, 1:min(3,ncol(emb)), drop=FALSE])
  df$barcode <- rownames(df)
  write.csv(df, paste0(outp,"_",rd,".csv"), row.names=FALSE)
  w("reduction", rd, "written dims:", paste(dim(emb),collapse="x"))
}

# --- key-gene expression via assay slot ---
da <- obj@active.assay
assay <- obj@assays[[da]]
# Assay (v3/v4): slots counts/data ; Assay5 (v5): layers list
getmat <- function(a) {
  if (.hasSlot(a, "data") && length(a@data) > 0) return(a@data)          # v3/v4 normalized
  if (.hasSlot(a, "layers")) {                                            # v5
    ln <- names(a@layers)
    pick <- if ("data" %in% ln) "data" else ln[1]
    return(a@layers[[pick]])
  }
  if (.hasSlot(a, "counts")) return(a@counts)
  NULL
}
m <- getmat(assay)
# rownames may live on the assay for v5 layers
rn <- tryCatch(rownames(assay), error=function(e) NULL)
if (is.null(rn) && !is.null(rownames(m))) rn <- rownames(m)
if (!is.null(m) && !is.null(rn)) {
  rownames(m) <- rn
  present <- KEY_GENES[KEY_GENES %in% rn]
  w("key genes present:", paste(present, collapse=", "))
  if (length(present) > 0) {
    sub <- as.matrix(m[present, , drop=FALSE])
    cn <- colnames(m); if (is.null(cn)) cn <- rownames(obj@meta.data)
    colnames(sub) <- cn
    gdf <- as.data.frame(t(sub)); gdf$barcode <- rownames(gdf)
    write.csv(gdf, paste0(outp,"_genes.csv"), row.names=FALSE)
    w("gene matrix written:", paste(dim(sub),collapse="x"))
  }
} else {
  w("WARN: could not extract expression matrix")
}
close(struct)
cat("DONE", outp, "\n")
