#!/usr/bin/env python
"""
Stage 22 (Phase II): spatial context for the mediator modules -- motor-neuron proximity,
MAP2 neuronal immunofluorescence, and TDP-43 pathology.

DESIGN NOTE. Every C9-vs-sALS test in this project compares 4 vs 5 donors and is therefore
underpowered. Motor-neuron distance and MAP2 are WITHIN-DONOR gradients: each donor contributes
its own slope / correlation, which are then tested across all 18 donors. Far better powered.

CONFOUND CONTROL. 91% of near-motor-neuron spots lie in the anterior horn, so the distance
gradient is computed WITHIN THE ANTERIOR HORN ONLY -- otherwise it re-measures the grey/white axis.

Requires module assignments from stage 21 (mediator_modules.csv).
Outputs: results/figures/mediator_MNdistance_MAP2.png
         results/tables/mediator_{MNdistance_gradient,MAP2_correlation,TDP43_donorlevel}.csv
"""
import os
for v in ("OMP_NUM_THREADS","OPENBLAS_NUM_THREADS","MKL_NUM_THREADS","NUMBA_NUM_THREADS"): os.environ[v]="1"
os.environ["KMP_AFFINITY"]="disabled"
import numpy as np, pandas as pd, h5py, scipy.sparse as sp
from scipy import stats
ROOT="/home/sunli/C9orf72-ALS"; PROC=f"{ROOT}/data/processed"; TAB=f"{ROOT}/results/tables"
FA=f"{PROC}/GSE288365_annotated_authmeta.h5ad"; FL=f"{PROC}/GSE288365_annotated.h5ad"
GRPS=["group 1","group 2","group 3","group 4"]

def dec(g):
    if isinstance(g,h5py.Group):
        c=g["categories"][:].astype(str); k=g["codes"][:]; return np.where(k>=0,c[np.clip(k,0,None)],None)
    x=g[:]; return x.astype(str) if x.dtype.kind in "SO" else x

mod=pd.read_csv(f"{TAB}/mediator_modules.csv").set_index("gene")
genes_by_mod={c:mod[mod.cluster==c].index.tolist() for c in sorted(mod.cluster.unique())}
allg=sorted(mod.index)

with h5py.File(FA,"r") as h:
    o=h["obs"]
    obs=pd.DataFrame({"donor_id":dec(o["donor_id"]),"condition":dec(o["condition"]),"anat":dec(o["auth_manual"]),
                      "mn":dec(o["auth_MN_dist_group"]),"tdp":dec(o["auth_annotation_general"])})
    obs["MAP2"]=o["auth_MAP2_norm"][:]
with h5py.File(FL,"r") as h:
    gn=h["var"][h["var"].attrs["_index"]][:].astype(str); X=h["X"]; ip=X["indptr"][:]
    cols=np.array([np.where(gn==g)[0][0] for g in allg]); n=len(ip)-1
    E=np.zeros((n,len(cols)),np.float32)
    for s in range(0,n,20000):
        e=min(s+20000,n); a0,a1=ip[s],ip[e]
        m=sp.csr_matrix((X["data"][a0:a1],X["indices"][a0:a1],ip[s:e+1]-a0),shape=(e-s,len(gn)))
        E[s:e]=m[:,cols].toarray()
Edf=pd.DataFrame(E,columns=allg); z=(Edf-Edf.mean())/Edf.std()
MODS=[f"M{c}" for c in genes_by_mod]
for c,gl in genes_by_mod.items(): obs[f"M{c}"]=z[gl].mean(1).values

# (1) MN-distance gradient, anterior horn only
ah=obs[(obs.anat=="Anterior_horns")&(obs.mn.isin(GRPS))].copy(); ah["gi"]=ah.mn.str[-1].astype(int)
cnt=ah.groupby(["donor_id","gi"]).size().unstack(fill_value=0)
ah=ah[ah.donor_id.isin(cnt.index[(cnt>=15).sum(1)>=3])]
dm_=ah.groupby(["donor_id","condition","gi"],observed=True)[MODS].mean().reset_index()
nn=ah.groupby(["donor_id","gi"]).size().rename("n").reset_index()
dm_=dm_.merge(nn,on=["donor_id","gi"]); dm_=dm_[dm_.n>=15]
rows=[]; slopes=[]
for m in MODS:
    sl=[{"module":m,"donor_id":d,"condition":g.condition.iloc[0],"slope":np.polyfit(g.gi,g[m],1)[0]}
        for d,g in dm_.groupby("donor_id") if g.gi.nunique()>=3]
    sl=pd.DataFrame(sl); slopes.append(sl)
    c9=sl[sl.condition=="C9-ALS"].slope; sa=sl[sl.condition=="sALS"].slope; ct=sl[sl.condition=="control"].slope
    rows.append({"module":m,"mean_slope":sl.slope.mean(),"n_donors":len(sl),"n_negative":int((sl.slope<0).sum()),
                 "p_slope_vs0":stats.wilcoxon(sl.slope).pvalue,
                 "ALSvsCtrl_slope_p":stats.mannwhitneyu(pd.concat([c9,sa]),ct).pvalue,
                 "C9vssALS_slope_p":stats.mannwhitneyu(c9,sa).pvalue if min(len(c9),len(sa))>=3 else np.nan})
pd.DataFrame(rows).to_csv(f"{TAB}/mediator_MNdistance_gradient.csv",index=False)
pd.concat(slopes).to_csv(f"{TAB}/mediator_MNdistance_slopes_perdonor.csv",index=False)

# (2) MAP2 within-donor correlation, anterior horn
sub=obs[(obs.anat=="Anterior_horns")&obs.MAP2.notna()]
rows=[]
for m in MODS:
    cs=np.array([stats.spearmanr(g[m],g.MAP2)[0] for d,g in sub.groupby("donor_id") if len(g)>=50])
    rows.append({"module":m,"mean_rho_vs_MAP2":cs.mean(),"n_donors":len(cs),
                 "n_positive":int((cs>0).sum()),"p_vs0":stats.wilcoxon(cs).pvalue})
pd.DataFrame(rows).to_csv(f"{TAB}/mediator_MAP2_correlation.csv",index=False)

# (3) TDP-43 donor-level
obs["pathb"]=obs.tdp.isin(["Ring","TDP43_adjacent","TDP43_distant"])
rows=[]
for m in MODS:
    d_=[{"delta":g[g.pathb][m].mean()-g[(~g.pathb)&g.tdp.notna()][m].mean()}
        for (d,c),g in obs[obs.condition!="control"].groupby(["donor_id","condition"],observed=True)
        if (g.pathb.sum()>=10) and ((~g.pathb)&g.tdp.notna()).sum()>=10]
    d_=pd.DataFrame(d_)
    rows.append({"module":m,"n_donors":len(d_),"mean_delta":d_.delta.mean(),
                 "n_positive":int((d_.delta>0).sum()),
                 "p":stats.wilcoxon(d_.delta).pvalue if len(d_)>=5 else np.nan})
pd.DataFrame(rows).to_csv(f"{TAB}/mediator_TDP43_donorlevel.csv",index=False)
print("stage 22 done")
