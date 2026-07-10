#!/usr/bin/env python
"""Pseudobulk aggregation: sum raw counts per (sample x L1 celltype).
Streams per-sample from backed h5ad (memory-safe). Output: one CSV per celltype
(genes x samples) + a sample metadata table. Enables within-celltype DE by condition."""
import os
for v in ("OMP_NUM_THREADS","OPENBLAS_NUM_THREADS","MKL_NUM_THREADS","NUMBA_NUM_THREADS"):
    os.environ[v]="1"
os.environ["KMP_AFFINITY"]="disabled"; os.environ["OMP_PROC_BIND"]="FALSE"
import numpy as np, pandas as pd, scanpy as sc, scipy.sparse as sp, gc, json

ROOT="/home/sunli/C9orf72-ALS"
H5=f"{ROOT}/data/processed/GSE268995_raw_qc.h5ad"
OBS=f"{ROOT}/results/tables/GSE268995_annotated_obs.parquet"
OUTDIR=f"{ROOT}/results/tables/pseudobulk_GSE268995"; os.makedirs(OUTDIR,exist_ok=True)

obs=pd.read_parquet(OBS).set_index("cell_id")
obs=obs[obs["author_retained"]]
L1="predicted.celltype.l1"
celltypes=sorted(obs[L1].dropna().unique())
print("celltypes:",celltypes,flush=True)

A=sc.read_h5ad(H5,backed="r")
genes=A.var_names.to_numpy()
samples=sorted(obs["sample"].unique())

# accumulator: {celltype: DataFrame(genes x samples)}
acc={ct:np.zeros((len(genes),len(samples)),dtype=np.float64) for ct in celltypes}
ncells={ct:{s:0 for s in samples} for ct in celltypes}
name_to_pos={n:i for i,n in enumerate(A.obs_names.to_numpy())}

for si,s in enumerate(samples):
    cells_s=obs.index[obs["sample"]==s]
    pos=np.array([name_to_pos[c] for c in cells_s if c in name_to_pos])
    if len(pos)==0: continue
    order=np.argsort(pos); pos_sorted=pos[order]
    X=A[pos_sorted].to_memory().X  # cells x genes sparse
    X=sp.csr_matrix(X)
    cells_sorted=[cells_s[i] for i in order]
    ct_s=obs.loc[cells_sorted,L1].to_numpy()
    for ct in celltypes:
        mask=ct_s==ct
        if mask.sum()==0: continue
        acc[ct][:,si]=np.asarray(X[mask].sum(axis=0)).ravel()
        ncells[ct][s]=int(mask.sum())
    del X; gc.collect()
    print(f"  {s}: {len(pos)} cells",flush=True)
A.file.close()

for ct in celltypes:
    dfm=pd.DataFrame(acc[ct],index=genes,columns=samples).astype(int)
    safe=ct.replace(" ","_").replace("/","_")
    dfm.to_csv(f"{OUTDIR}/pb_{safe}.csv")
nc=pd.DataFrame(ncells).T  # celltype x sample
nc.to_csv(f"{OUTDIR}/ncells_per_celltype_sample.csv")
# sample meta
smeta=obs.groupby("sample").agg(condition=("condition","first"),
        sample_label_conflict=("sample_label_conflict","first")).reset_index()
smeta.to_csv(f"{OUTDIR}/sample_meta.csv",index=False)
print("DONE. wrote",len(celltypes),"pseudobulk matrices to",OUTDIR,flush=True)
