#!/usr/bin/env python
"""Point 1: broaden the CNS myeloid-program test in peripheral blood.
(a) Surface proteins (imputed CITE-seq / predicted_ADT, 228 markers) within monocytes,
    C9 vs control and C9 vs pooled sALS -> CD68 protein trends up in C9 (nominal p=0.033 vs sALS).
(b) Per-gene DAM breakdown in monocyte pseudobulk -> tetraspanin/integrin arm (CD9, ITGAX, CD63)
    trends up; CNS lipid core (GPNMB/APOE/TREM2/SPP1) not expressed peripherally (filtered out).
(c) DAM module also tested in DC compartment -> no signal.
Same author-RDS barcode-shift remap applied to ADT cell_ids. Nothing survives FDR (4 C9 donors)."""
import os
for v in ("OMP_NUM_THREADS","OPENBLAS_NUM_THREADS","MKL_NUM_THREADS","NUMBA_NUM_THREADS"):
    os.environ[v]="1"
os.environ["KMP_AFFINITY"]="disabled"; os.environ["OMP_PROC_BIND"]="FALSE"
import numpy as np, pandas as pd
from scipy import stats

ROOT="/home/sunli/C9orf72-ALS"
PB=f"{ROOT}/results/tables/pseudobulk_GSE268995"
SHIFT={"A1":"C1","B1":"D1","C1":"E1","D1":"F1","E1":"A1","F1":"B1"}

# --- (a) ADT surface proteins ---
adt=pd.read_csv(f"{ROOT}/results/tables/GSE268995_predicted_ADT.csv.gz")
adt["barcode"]=adt["cell_id"].str.split("::").str[1]
adt["our_sample"]=adt["cell_id"].str.split("::").str[0].map(lambda s:SHIFT.get(s,s))
adt["phys_cell"]=adt["our_sample"]+"::"+adt["barcode"]
obs=pd.read_parquet(f"{ROOT}/results/tables/GSE268995_annotated_obs.parquet").set_index("cell_id")
obs=obs[obs["author_retained"]]
adt=adt.set_index("phys_cell")
adt["L1"]=obs["predicted.celltype.l1"].reindex(adt.index).values
adt["sample"]=obs["sample"].reindex(adt.index).values
protein_cols=[c for c in adt.columns if c not in ("cell_id","barcode","our_sample","L1","sample")]
smeta=pd.read_csv(f"{PB}/sample_meta.csv").set_index("sample")

mono_adt=adt[adt["L1"]=="Mono"]
persamp=mono_adt.groupby("sample")[protein_cols].mean().join(smeta[["condition"]])
markers=[m for m in ["CD68","CD11c","HLA-DR","CD14","CD16","CD163","CD64","CD86","CD11b-1","CD11b-2","CD206","CD169","CD115"] if m in protein_cols]
def test(tbl,cols,gA,gB):
    r=[]
    for m in cols:
        a=tbl.loc[tbl.condition.isin(gA),m].values; b=tbl.loc[tbl.condition.isin(gB),m].values
        _,p=stats.mannwhitneyu(a,b,alternative="two-sided")
        r.append({"protein":m,"mean_C9":a.mean(),"mean_other":b.mean(),"diff":a.mean()-b.mean(),"p":p})
    r=pd.DataFrame(r).sort_values("p"); r["p_fdr"]=stats.false_discovery_control(r["p"]); return r
test(persamp,markers,["C9-ALS"],["control"]).to_csv(f"{ROOT}/results/tables/GSE268995_mono_ADT_C9vsControl.csv",index=False)
test(persamp,markers,["C9-ALS"],["sALS_slow","sALS_fast"]).to_csv(f"{ROOT}/results/tables/GSE268995_mono_ADT_C9vsSALS.csv",index=False)

# --- (b) per-gene DAM in monocyte pseudobulk ---
dam=["GPNMB","APOE","TREM2","TYROBP","CD68","ITGAX","CTSB","CTSD","LGALS3","SPP1","FTL","APOC1","LPL","CD9","MSR1","CD63"]
def cpm(ct):
    m=pd.read_csv(f"{PB}/pb_{ct}.csv",index_col=0)[smeta.index]; mf=m[(m>=10).sum(axis=1)>=8]
    return np.log1p(mf.div(mf.sum(axis=0),axis=1)*1e6)
cm=cpm("Mono"); rows=[]
for g in dam:
    if g not in cm.index: rows.append({"gene":g,"in_mono":False}); continue
    v=cm.loc[g]
    a=v[smeta[smeta.condition=="C9-ALS"].index].values;c=v[smeta[smeta.condition=="control"].index].values
    s=v[smeta[smeta.condition.isin(["sALS_slow","sALS_fast"])].index].values
    _,pc=stats.mannwhitneyu(a,c);_,psl=stats.mannwhitneyu(a,s)
    rows.append({"gene":g,"in_mono":True,"C9":a.mean(),"ctrl":c.mean(),"sALS":s.mean(),
                 "d_vs_ctrl":a.mean()-c.mean(),"p_vs_ctrl":pc,"d_vs_sALS":a.mean()-s.mean(),"p_vs_sALS":psl})
pd.DataFrame(rows).to_csv(f"{ROOT}/results/tables/GSE268995_mono_DAMgenes_percond.csv",index=False)
print("point-1 tables written")
