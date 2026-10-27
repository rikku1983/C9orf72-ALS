#!/usr/bin/env python
"""Stage 19b: GSEA (Hallmark) on C9-ALS vs sALS rankings computed per AUTHOR histological
region (stage 19). Replaces the stage-07 'ventral horn' GSEA, whose region label was wrong.

RESULT: TNF-alpha/NF-kB and Inflammatory Response survive in the TRUE anterior horn and in
other grey matter. Unfolded Protein Response and mTORC1 do NOT survive -- they were artifacts
of our mislabelled marker-score region (which was 33% white matter, 21% meninges).
"""
import os
for v in ("OMP_NUM_THREADS","OPENBLAS_NUM_THREADS","MKL_NUM_THREADS","NUMBA_NUM_THREADS"): os.environ[v]="1"
os.environ["KMP_AFFINITY"]="disabled"
import pandas as pd, gseapy as gp
ROOT="/home/sunli/C9orf72-ALS"; TAB=f"{ROOT}/results/tables"
REGIONS=["Anterior_horns","White_matter","Posterior_white_matter","Grey_matter"]
allr=[]
for reg in REGIONS:
    d=pd.read_csv(f"{TAB}/pseudobulk_authoranat_C9vsSALS_{reg}.csv",index_col=0).dropna(subset=["stat"])
    r=gp.prerank(rnk=d["stat"].sort_values(ascending=False),gene_sets="MSigDB_Hallmark_2020",
                 permutation_num=1000,seed=0,threads=1,min_size=10,max_size=500,outdir=None,verbose=False).res2d
    r["region"]=reg
    r=r.rename(columns={"Term":"term","NES":"nes","FDR q-val":"fdr","NOM p-val":"pval"})
    r["nes"]=pd.to_numeric(r["nes"]); r["fdr"]=pd.to_numeric(r["fdr"])
    allr.append(r[["region","term","nes","pval","fdr"]])
    print(f"{reg}: Hallmark FDR<0.05 = {(r.fdr<0.05).sum()}")
pd.concat(allr).to_csv(f"{TAB}/gsea_authoranat_C9vsSALS_Hallmark.csv",index=False)
