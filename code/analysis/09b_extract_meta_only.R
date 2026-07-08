#!/usr/bin/env Rscript
# Extract meta.data ONLY from an author Seurat RDS (fast, light).
# Usage: Rscript 09b_extract_meta_only.R <rds> <out_meta.csv>
suppressMessages({library(Matrix)})
args <- commandArgs(trailingOnly=TRUE)
rds <- args[1]; outcsv <- args[2]
obj <- readRDS(rds)              # full Seurat S4; access slots directly
md <- obj@meta.data
md$barcode <- rownames(md)
write.csv(md, outcsv, row.names=FALSE)
cat("wrote", outcsv, nrow(md), "spots x", ncol(md), "cols\n")
