#!/usr/bin/env python
"""
Stage 17 (spatial arm): Does the myeloid/complement/DAM/NF-kB program spatially
co-localize with TDP-43 pathology, and does the coupling differ C9-ALS vs sALS?

Uses author spot annotations in GSE288365_annotated_authmeta.h5ad:
  auth_annotation_general: {Not_selected, TDP43_distant, TDP43_adjacent, Ring}
  auth_TDP43_mean:         continuous TDP-43 IF intensity

CONCLUSION: The dramatic spot-level divergence (sALS recruits programs to
pathology, C9 does not) is DRIVEN BY PSEUDOREPLICATION. Donor-level aggregation
(only 3 donors/genotype have both pathology-bearing and unaffected spots) shows
heterogeneous, sign-inconsistent deltas -> NOT a reproducible genotype effect.
Reported as an honest negative.
"""
import os
for v in ("OMP_NUM_THREADS","OPENBLAS_NUM_THREADS","MKL_NUM_THREADS","NUMBA_NUM_THREADS"):
    os.environ[v]="1"
os.environ["KMP_AFFINITY"]="disabled"; os.environ["OMP_PROC_BIND"]="FALSE"
import numpy as np, pandas as pd, scanpy as sc
from scipy import stats
ROOT="/home/sunli/C9orf72-ALS"

a=sc.read_h5ad(f"{ROOT}/data/processed/GSE288365_annotated_authmeta.h5ad")
programs={
 "myeloid":["AIF1","CD68","ITGAX","TYROBP","FCER1G","CST3","C1QA","C1QB","C1QC","LAPTM5","CTSS","LY86"],
 "complement":["C1QA","C1QB","C1QC","C3","C1R","C1S","CFB","C2","C4A","C4B"],
 "DAM":["GPNMB","APOE","TREM2","TYROBP","CTSB","CTSD","LGALS3","CD9","CD63","FTL","APOC1","LPL","SPP1","MSR1"],
 "NFkB_TNF":["NFKB1","NFKB2","RELB","TNF","TNFAIP3","NFKBIA","CXCL8","CCL2","IL1B","TNFAIP2","BIRC3","TRAF1","CD44","ICAM1"],
}
for name,genes in programs.items():
    sc.tl.score_genes(a,[g for g in genes if g in a.var_names],score_name=f"sc_{name}",ctrl_size=50)
scores=[f"sc_{k}" for k in programs]

o=a.obs.copy()
o["is_ALS"]=o["condition"].isin(["C9-ALS","sALS"])
grade_order=["Not_selected","TDP43_distant","TDP43_adjacent","Ring"]
o["tdp_grade"]=pd.Categorical(o["auth_annotation_general"],categories=grade_order,ordered=True)
als=o[o["is_ALS"] & o["tdp_grade"].notna()].copy()
als["path_bearing"]=als["tdp_grade"].isin(["Ring","TDP43_adjacent","TDP43_distant"])

# spot-level (pseudoreplicated -- reported for transparency, NOT trusted)
rows=[]
for cond in ["C9-ALS","sALS"]:
    sub=als[als["condition"]==cond]
    for s in scores:
        pb=sub.loc[sub["path_bearing"],s].values; ns=sub.loc[~sub["path_bearing"],s].values
        u,p=stats.mannwhitneyu(pb,ns,alternative="two-sided")
        rows.append({"genotype":cond,"program":s.replace("sc_",""),"path_bearing":pb.mean(),
                     "not_selected":ns.mean(),"diff":pb.mean()-ns.mean(),"p":p,"n_path":len(pb)})
res=pd.DataFrame(rows)
res.to_csv(f"{ROOT}/results/tables/GSE288365_TDP43_coupling_spotlevel.csv",index=False)

# donor-level (the honest test)
rows=[]
for cond in ["C9-ALS","sALS"]:
    sub=als[als["condition"]==cond]
    for donor,dg in sub.groupby("donor_id",observed=True):
        pb=dg[dg["path_bearing"]]; ns=dg[~dg["path_bearing"]]
        if len(pb)<10 or len(ns)<10: continue
        for s in scores:
            rows.append({"genotype":cond,"donor":donor,"program":s.replace("sc_",""),
                         "delta":pb[s].mean()-ns[s].mean(),"n_path":len(pb)})
dd=pd.DataFrame(rows)
dd.to_csv(f"{ROOT}/results/tables/GSE288365_TDP43_coupling_donorlevel.csv",index=False)
als.groupby("tdp_grade",observed=True)[scores].mean().to_csv(f"{ROOT}/results/tables/GSE288365_program_by_TDPgrade.csv")
print("stage 17 done; donor-level mean deltas:")
print(dd.pivot_table(index="program",columns="genotype",values="delta").round(4))
