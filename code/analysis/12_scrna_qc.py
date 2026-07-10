#!/usr/bin/env python3
"""GSE268995 PBMC scRNA-seq QC (memory-streaming): per-sample QC + concat_on_disk."""
import os,gc
for v in ["OMP_NUM_THREADS","OPENBLAS_NUM_THREADS","MKL_NUM_THREADS","NUMBA_NUM_THREADS"]:
    os.environ[v]="1"
os.environ["KMP_AFFINITY"]="disabled"; os.environ["OMP_PROC_BIND"]="FALSE"
import numpy as np, pandas as pd, scipy.io
import scanpy as sc, anndata as ad
from anndata.experimental import concat_on_disk

REPO="/home/sunli/C9orf72-ALS"
RAW=f"{REPO}/data/raw/GSE268995"
TMP=f"{REPO}/data/processed/tmp_scrna"; os.makedirs(TMP,exist_ok=True)
MAP=pd.read_csv(f"{REPO}/data/metadata/GSE268995_sample_condition_map.csv")

obs_all=[]; tmp_paths=[]
for _,r in MAP.iterrows():
    s=r['sample']; pfx=f"{r.gsm}_{s}"
    M=scipy.io.mmread(f"{RAW}/{pfx}matrix.mtx.gz").T.tocsr().astype(np.float32)
    bc=pd.read_csv(f"{RAW}/{pfx}barcodes.tsv.gz",header=None)[0].values
    ft=pd.read_csv(f"{RAW}/{pfx}features.tsv.gz",header=None,sep="\t")
    a=ad.AnnData(X=M)
    a.obs_names=[f"{s}::{b}" for b in bc]
    a.var_names=ft[1].values; a.var["ensembl"]=ft[0].values
    a.var_names_make_unique()
    a.obs["sample"]=s; a.obs["gsm"]=r.gsm; a.obs["condition"]=r.condition; a.obs["disease"]=r.disease
    a.var["mt"]=a.var_names.str.startswith("MT-")
    a.var["ribo"]=a.var_names.str.startswith(("RPS","RPL"))
    sc.pp.calculate_qc_metrics(a,qc_vars=["mt","ribo"],inplace=True,percent_top=None)
    obs_all.append(a.obs[["sample","gsm","condition","disease","n_genes_by_counts","total_counts","pct_counts_mt","pct_counts_ribo"]].copy())
    p=f"{TMP}/{s}.h5ad"; a.write(p); tmp_paths.append(p)
    print(f"{s} {r.condition}: {a.n_obs} cells, medGenes={a.obs.n_genes_by_counts.median():.0f} medMT={a.obs.pct_counts_mt.median():.1f}%",flush=True)
    del a,M; gc.collect()

OBS=pd.concat(obs_all); OBS.to_parquet(f"{REPO}/data/processed/GSE268995_percell_qc.parquet")
print(f"\nTOTAL: {len(OBS)} cells x 36601 genes")

# per-sample summary
g=OBS.groupby("sample")
summ=pd.DataFrame({"condition":g["condition"].first(),"n_cells":g.size(),
    "median_genes":g["n_genes_by_counts"].median(),"median_counts":g["total_counts"].median(),
    "median_pct_mt":g["pct_counts_mt"].median(),"median_pct_ribo":g["pct_counts_ribo"].median()}).round(1)
summ.to_csv(f"{REPO}/results/tables/GSE268995_qc_per_sample.csv")
print("\n=== PER-SAMPLE QC ==="); print(summ.to_string())
print("\n=== PER-CONDITION ===")
print(OBS.groupby("condition").agg(n_samples=("sample","nunique"),n_cells=("sample","size"),
    med_genes=("n_genes_by_counts","median"),med_mt=("pct_counts_mt","median")).round(1).to_string())

# merge on disk (memory-safe)
concat_on_disk(tmp_paths, f"{REPO}/data/processed/GSE268995_raw_qc.h5ad", join="outer")
print("\nconcat_on_disk done -> GSE268995_raw_qc.h5ad")
