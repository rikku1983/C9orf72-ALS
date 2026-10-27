#!/usr/bin/env python
"""
Stage 19: REDO pseudobulk DE (C9-ALS vs sALS) using the AUTHORS' manual histological
annotation (auth_manual) instead of our marker-score regions.

WHY: stage 18 showed our marker-score `region` column disagrees badly with the authors'
histology -- only 15% of our 'ventral_horn_MN' spots lie in the authors' Anterior_horns
(33% are white matter, 21% meninges). Every 'ventral horn' claim from stages 05-08
therefore carries a wrong anatomical label and must be recomputed.

Aggregates RAW counts to one profile per donor x author-region, then DESeq2 across donors.
"""
import os
for v in ("OMP_NUM_THREADS","OPENBLAS_NUM_THREADS","MKL_NUM_THREADS","NUMBA_NUM_THREADS"): os.environ[v]="1"
os.environ["KMP_AFFINITY"]="disabled"
import numpy as np, pandas as pd, h5py, scipy.sparse as sp
ROOT="/home/sunli/C9orf72-ALS"; TAB=f"{ROOT}/results/tables"
MIN_SPOTS=50; MIN_DONORS=3
REGIONS=["Anterior_horns","White_matter","Posterior_white_matter","Grey_matter"]

def dec(g):
    if isinstance(g,h5py.Group):
        c=g["categories"][:].astype(str); k=g["codes"][:]; return np.where(k>=0,c[np.clip(k,0,None)],None)
    x=g[:]; return x.astype(str) if x.dtype.kind in "SO" else x

# obs from annotated (has auth_manual), raw counts from qc h5ad
with h5py.File(f"{ROOT}/data/processed/GSE288365_annotated_authmeta.h5ad","r") as h:
    o=h["obs"]; spots=o[o.attrs["_index"]][:].astype(str)
    meta=pd.DataFrame({"spot":spots,"donor_id":dec(o["donor_id"]),"condition":dec(o["condition"]),
                       "anat":dec(o["auth_manual"])})
with h5py.File(f"{ROOT}/data/processed/GSE288365_qc.h5ad","r") as h:
    qspots=h["obs"][h["obs"].attrs["_index"]][:].astype(str)
    genes=h["var"][h["var"].attrs["_index"]][:].astype(str)
    pos=pd.Series(np.arange(len(qspots)),index=qspots).reindex(meta.spot).values
    assert not np.isnan(pos).any(), "spot mismatch between annotated and qc h5ad"
    meta["row"]=pos.astype(int)
    X=h["X"]; indptr=X["indptr"][:]
    groups=[(d,r,c) for (d,r,c),g in meta.groupby(["donor_id","anat","condition"],observed=True) if len(g)>=MIN_SPOTS and r in REGIONS]
    gidx={k:i for i,k in enumerate(groups)}
    PB=np.zeros((len(groups),len(genes)),dtype=np.float64)
    key=list(zip(meta.donor_id,meta.anat,meta.condition))
    rows_by_g={}
    for i,k in enumerate(key):
        if k in gidx: rows_by_g.setdefault(gidx[k],[]).append(meta.row.values[i])
    for gi,rws in rows_by_g.items():
        rws=np.sort(np.array(rws))
        for s in range(0,len(rws),4000):
            blk=rws[s:s+4000]
            for rr in blk:
                a0,a1=indptr[rr],indptr[rr+1]
                np.add.at(PB[gi],X["indices"][a0:a1],X["data"][a0:a1])
pbdf=pd.DataFrame(PB,columns=genes)
gmeta=pd.DataFrame(groups,columns=["donor_id","anat","condition"])
gmeta["n_spots"]=[len(rows_by_g[i]) for i in range(len(groups))]
pbdf.index=gmeta.donor_id+"|"+gmeta.anat
gmeta.to_csv(f"{TAB}/pseudobulk_authoranat_groups.csv",index=False)
print(gmeta.groupby(["anat","condition"]).size().unstack(fill_value=0).to_string())

from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats
for reg in REGIONS:
    m=gmeta[gmeta.anat==reg].copy()
    sub=m[m.condition.isin(["C9-ALS","sALS"])]
    if sub.groupby("condition").size().min()<MIN_DONORS:
        print(f"[skip] {reg}: n per group {dict(sub.groupby('condition').size())}"); continue
    cnt=pbdf.loc[(m.donor_id+"|"+m.anat).values][ (sub.index-m.index[0]).map(lambda x:True) if False else slice(None)]
    cnt=pbdf.loc[(sub.donor_id+"|"+sub.anat).values].round().astype(int)
    keepg=((cnt>=10).sum(0)>=3)
    cnt=cnt.loc[:,keepg]
    md=sub[["condition"]].copy(); md.index=cnt.index
    dds=DeseqDataSet(counts=cnt,metadata=md,design_factors="condition",ref_level=["condition","sALS"],quiet=True)
    dds.deseq2()
    res=DeseqStats(dds,contrast=["condition","C9-ALS","sALS"],quiet=True); res.summary()
    out=res.results_df.sort_values("padj")
    out.to_csv(f"{TAB}/pseudobulk_authoranat_C9vsSALS_{reg}.csv")
    nsig=int((out.padj<0.05).sum())
    print(f"\n=== {reg}: n={dict(sub.groupby('condition').size())} genes={cnt.shape[1]} FDR<0.05: {nsig}")
    print(out.head(8)[["log2FoldChange","pvalue","padj"]].round(4).to_string())
