#!/usr/bin/env python
"""
07_pathway_enrichment.py — GSEA prerank on pseudobulk DE results.

For each contrast x scope, rank genes by (sign of log2FC * -log10(pvalue)) and
run gseapy.prerank against Hallmark, GO-BP, and Reactome (Enrichr libraries).

Input : results/tables/pseudobulk_<contrast>_<scope>.csv
Output: results/tables/gsea_<contrast>_<scope>_<lib>.csv
        results/figures/gsea_summary.png
"""
import os
os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("KMP_AFFINITY", "disabled")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
import warnings; warnings.filterwarnings("ignore")
import numpy as np, pandas as pd, gseapy as gp
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

REPO="/home/sunli/C9orf72-ALS"; TAB=f"{REPO}/results/tables"; FIG=f"{REPO}/results/figures"
LIBS={"Hallmark":"MSigDB_Hallmark_2020",
      "GOBP":"GO_Biological_Process_2023",
      "Reactome":"Reactome_2022"}
CONTRASTS=["C9ALS_vs_control","sALS_vs_control","C9ALS_vs_sALS"]
SCOPES=["global","ventral_horn_MN"]

def rank_df(path):
    d=pd.read_csv(path,index_col=0)
    d=d.dropna(subset=["log2FoldChange","pvalue"])
    d=d[~d.index.duplicated(keep="first")]
    d["score"]=np.sign(d["log2FoldChange"])*-np.log10(d["pvalue"].clip(lower=1e-300))
    r=d["score"].sort_values(ascending=False)
    return r

def main():
    summary=[]
    for contrast in CONTRASTS:
        for scope in SCOPES:
            path=f"{TAB}/pseudobulk_{contrast}_{scope}.csv"
            if not os.path.exists(path):
                print("skip missing",path); continue
            rnk=rank_df(path)
            for libname,lib in LIBS.items():
                try:
                    pre=gp.prerank(rnk=rnk.reset_index(),gene_sets=lib,
                                   min_size=5,max_size=500,permutation_num=1000,
                                   seed=0,no_plot=True,outdir=None,threads=4)
                    res=pre.res2d.copy()
                    res.to_csv(f"{TAB}/gsea_{contrast}_{scope}_{libname}.csv",index=False)
                    res["contrast"]=contrast;res["scope"]=scope;res["lib"]=libname
                    summary.append(res)
                    sig=res[res["FDR q-val"].astype(float)<0.25]
                    print(f"{contrast} {scope} {libname}: {len(sig)} sets FDR<0.25")
                except Exception as e:
                    print("ERR",contrast,scope,libname,type(e).__name__,str(e)[:120])
    if summary:
        allres=pd.concat(summary,ignore_index=True)
        allres.to_csv(f"{TAB}/gsea_all_summary.csv",index=False)
        print("saved gsea_all_summary.csv rows=",len(allres))

if __name__=="__main__":
    main()
