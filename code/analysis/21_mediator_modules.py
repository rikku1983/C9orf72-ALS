#!/usr/bin/env python
"""
Stage 21 (Phase II): map the 118 bulk causal-mediator genes onto spinal-cord cell types
and author histological regions.

The bulk data were not available; this characterises the gene list in GSE288365 rather
than re-testing the mediation model.

Pipeline: panel coverage -> pseudobulk co-expression modules -> Cell2location cell-type
attribution -> regional profile -> spatial maps -> condition contrasts -> PBMC detectability
-> annotated table.

NOTE on module detection: spot-level correlation is useless here (zero-inflation collapses
everything into one cluster of 99/104). Modules are computed on donor x region pseudobulk
log-CPM with the region main effect removed.

Outputs: results/figures/mediator_*.png, results/tables/mediator_*.csv
"""
import os
for v in ("OMP_NUM_THREADS","OPENBLAS_NUM_THREADS","MKL_NUM_THREADS","NUMBA_NUM_THREADS"): os.environ[v]="1"
os.environ["KMP_AFFINITY"]="disabled"
import numpy as np, pandas as pd, h5py, scipy.sparse as sp
from scipy import stats
from scipy.cluster.hierarchy import linkage, fcluster, leaves_list
from scipy.spatial.distance import squareform
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

ROOT="/home/sunli/C9orf72-ALS"; PROC=f"{ROOT}/data/processed"; TAB=f"{ROOT}/results/tables"; FIG=f"{ROOT}/results/figures"
FA=f"{PROC}/GSE288365_annotated_authmeta.h5ad"; FL=f"{PROC}/GSE288365_annotated.h5ad"
REG4=["Anterior_horns","Grey_matter","White_matter","Posterior_white_matter"]
RL={"Anterior_horns":"Anterior horn","Grey_matter":"Other grey","White_matter":"Anterolat. WM","Posterior_white_matter":"Posterior WM"}
CT=['Astrocytes','Endothelial','Ependymal.Cells','Lymphocytes','Macrophages','Meninges','Microglia',
    'Neurons','OPC','Oligodendrocytes','Pericytes','ProliferatingMicroglia','Schwann']
PAL={"control":"#4477AA","sALS":"#DDAA33","C9-ALS":"#CC6677"}; ORDER=["control","sALS","C9-ALS"]
MODNAME={1:"M1 Oligodendrocyte/myelin",2:"M2 Astrocyte",3:"M3 Neuronal",4:"M4 Vascular",5:"M5 Myeloid/DAM"}
MODCOL={1:"#88CCEE",2:"#44AA99",3:"#DDCC77",4:"#AA4499",5:"#CC6677"}
ACOL={"Anterior_horns":"#332288","Grey_matter":"#44AA99","White_matter":"#882255","Posterior_white_matter":"#88CCEE","Meninges":"#BBBBBB"}

MED="""NFIX MAPK8IP2 GABARAPL2 CAPG CYBA SOAT1 PYGM BRINP1 ADD1 SMOX FXYD5 PAPOLA TREM2 DSP CIRBP HCK
PPP1R16B PLLP PYCARD CORO2B SIGLEC8 SLC38A1 SOD2 B3GAT2 NCF2 CTSD LDB3 PLAU NCKAP1L RPS10 MYRF RPL23
WDFY4 CD68 RPL27 LLGL1 LGALS3 CHIT1 SLC38A2 RPS15A SLC37A2 CD63 ITM2C GMPR CHRNA1 ARRB2 PIK3R5 EMP3
CD53 PYCR2 LIX1 PSD2 PLA2G7 DAAM2 TMEM47 NADK2 GJA1 OTOA SLC24A2 C1QC ITGB2 LSS NBEAL2 RPL26 PAQR4
LAPTM5 SDC3 IGSF8 CTSS S100A11 NFASC EDIL3 SYK HEPACAM NDRG2 CMTM5 CKB APBB1 STX3 DPEP2 KLK6 CD300A
LGALS9 P2RY1 GPR37L1 DCLK2 GPR37 C3AR1 IL17D CARNS1 ANO6 DNHD1 ALOX15B RPL24P4 FAM181B ANXA2 RPS17
PRKX APOBR OPALIN SLC5A3 IFI30 GPX1 CEBPA TRIL TNFRSF14 AC006011.2 RPL7A ADORA3 SPI1 SLC24A4 APOC1
GRN TYROBP CNTN1 CHCHD2 GPNMB IQGAP2""".split()

def dec(g):
    if isinstance(g,h5py.Group):
        c=g["categories"][:].astype(str); k=g["codes"][:]; return np.where(k>=0,c[np.clip(k,0,None)],None)
    x=g[:]; return x.astype(str) if x.dtype.kind in "SO" else x

with h5py.File(FA,"r") as h:
    o=h["obs"]
    obs=pd.DataFrame({"donor_id":dec(o["donor_id"]),"condition":dec(o["condition"]),
                      "anat":dec(o["auth_manual"]),"section":dec(o["section"])})
    for c in CT: obs[f"ct_{c}"]=o[f"auth_norm_ct_{c}"][:]
    xy=h["obsm/spatial"][:]
obs["x"],obs["y"]=xy[:,0],xy[:,1]
CTcols=[f"ct_{c}" for c in CT]

with h5py.File(FL,"r") as h:
    genes=h["var"][h["var"].attrs["_index"]][:].astype(str)
    found=[g for g in MED if g in genes]; missing=[g for g in MED if g not in genes]
    X=h["X"]; ip=X["indptr"][:]; cols=np.array([np.where(genes==g)[0][0] for g in found])
    n=len(ip)-1; E=np.zeros((n,len(cols)),np.float32)
    for s in range(0,n,20000):
        e=min(s+20000,n); a0,a1=ip[s],ip[e]
        m=sp.csr_matrix((X["data"][a0:a1],X["indices"][a0:a1],ip[s:e+1]-a0),shape=(e-s,len(genes)))
        E[s:e]=m[:,cols].toarray()
Edf=pd.DataFrame(E,columns=found)

base=pd.DataFrame({"gene":found,"pct_spots_expressing":(Edf>0).mean(0).values*100,
                   "mean_lognorm":Edf.mean(0).values})
base["detected"]=base.pct_spots_expressing>=1.0
base["artifact_flag"]=np.where(base.gene.isin(["PYGM","CHRNA1","CARNS1","LDB3"]),"skeletal_muscle","")
miss=pd.DataFrame({"gene":missing}); miss["detected"]=False
miss["artifact_flag"]=np.where(miss.gene.str.startswith(("RPL","RPS")),"ribosomal_not_in_probe_panel","not_in_probe_panel")
base_full=pd.concat([base,miss],ignore_index=True); base_full["in_spatial_panel"]=~base_full.gene.isin(missing)
base_full.to_csv(f"{TAB}/mediator_genes_base.csv",index=False)
use=base.loc[base.detected,"gene"].tolist()

# ---- modules on donor x region pseudobulk ----
pb=pd.read_csv(f"{TAB}/pseudobulk_authoranat_counts.csv",index_col=0)
gm=pd.read_csv(f"{TAB}/pseudobulk_authoranat_groups.csv")
cpm=np.log1p(pb.div(pb.sum(1),axis=0)*1e6)
use2=[g for g in use if g in cpm.columns and cpm[g].std()>0]
Mr=cpm[use2].copy(); Mr["anat"]=gm.anat.values
Mc=Mr.groupby("anat").transform(lambda v:v-v.mean())
R=Mc.corr()
D=1-R.values; np.fill_diagonal(D,0); D=(D+D.T)/2
Z=linkage(squareform(D,checks=False),method="average")
mod=pd.DataFrame({"gene":use2,"cluster":fcluster(Z,5,criterion="maxclust")}).set_index("gene")
mod["module"]=mod.cluster.map(MODNAME)

# ---- cell-type attribution on pseudobulk ----
pbkey=gm.donor_id+"|"+gm.anat
ctpb=obs[obs.anat.isin(gm.anat.unique())].groupby(["donor_id","anat"],observed=True)[CTcols].mean()
ctpb.index=[f"{a}|{b}" for a,b in ctpb.index]; ctpb=ctpb.reindex(pbkey.values)
ctc=ctpb.copy(); ctc["anat"]=gm.anat.values; ctc=ctc.groupby("anat").transform(lambda v:v-v.mean())
CORR2=pd.DataFrame({c.replace("ct_",""):Mc.corrwith(ctc[c]) for c in CTcols})
mod["best_celltype_pb"]=CORR2.idxmax(1); mod["best_rho_pb"]=CORR2.max(1)
mod.reset_index().to_csv(f"{TAB}/mediator_modules.csv",index=False)
CORR2.to_csv(f"{TAB}/mediator_celltype_correlation.csv")

# ---- module scores, regional profile, contrasts ----
zc=(Edf[use2]-Edf[use2].mean())/Edf[use2].std()
MODS=[f"M{c}" for c in sorted(mod.cluster.unique())]
for c in sorted(mod.cluster.unique()): obs[f"M{c}"]=zc[mod[mod.cluster==c].index.tolist()].mean(1).values
sub=obs[obs.anat.isin(REG4)]
dr=sub.groupby(["donor_id","condition","anat"],observed=True)[MODS].mean().reset_index()
dr["n"]=sub.groupby(["donor_id","condition","anat"],observed=True).size().values
dr=dr[dr.n>=50]; dr.to_csv(f"{TAB}/mediator_module_donor_region.csv",index=False)
rows=[]; als=dr[dr.condition!="control"]
for m in MODS:
    pk=dr.groupby("anat")[m].mean().reindex(REG4).idxmax()
    p=als.pivot_table(index="donor_id",columns="anat",values=m)
    for r2 in REG4:
        if r2==pk: continue
        d_=(p[pk]-p[r2]).dropna()
        rows.append({"module":m,"peak_region":pk,"vs_region":r2,"mean_diff":d_.mean(),
                     "n_pos":int((d_>0).sum()),"n":len(d_),"p":stats.wilcoxon(d_).pvalue})
reg=pd.DataFrame(rows); reg.to_csv(f"{TAB}/mediator_region_profile.csv",index=False)
rank={"control":0,"sALS":1,"C9-ALS":2}; rows=[]
for m in MODS:
    for rg in REG4:
        q=dr[dr.anat==rg]; g=lambda c:q.loc[q.condition==c,m]
        alls=pd.concat([g("sALS"),g("C9-ALS")]); rho,pt=stats.spearmanr(q.condition.map(rank),q[m])
        rows.append({"module":m,"region":rg,"ALSvsCtrl_diff":alls.mean()-g("control").mean(),
                     "ALSvsCtrl_p":stats.mannwhitneyu(alls,g("control")).pvalue,
                     "C9vssALS_diff":g("C9-ALS").mean()-g("sALS").mean(),
                     "C9vssALS_p":stats.mannwhitneyu(g("C9-ALS"),g("sALS")).pvalue,
                     "trend_rho":rho,"trend_p":pt,"n_donors":len(q)})
con=pd.DataFrame(rows); con.to_csv(f"{TAB}/mediator_module_contrasts.csv",index=False)

# ---- PBMC detectability ----
blood=None
PBD=f"{TAB}/pseudobulk_GSE268995"
if os.path.isdir(PBD):
    mono=pd.read_csv(f"{PBD}/pb_Mono.csv",index_col=0)
    cpm_m=np.log1p(mono.div(mono.sum(0),axis=1)*1e6)
    blood=pd.DataFrame({"gene":cpm_m.index,"blood_mean_logCPM":cpm_m.mean(1).values,
                        "blood_n_samples_detected":(mono>0).sum(1).values})
    blood=blood[blood.gene.isin(MED)]
    blood["blood_detectable"]=(blood.blood_mean_logCPM>1.0)&(blood.blood_n_samples_detected>=20)
    blood.to_csv(f"{TAB}/mediator_blood_detectability.csv",index=False)

# ---- annotated table ----
pk_true=reg.drop_duplicates("module").set_index("module")["peak_region"]
mstat=con.loc[con.groupby("module")["ALSvsCtrl_p"].idxmin()].set_index("module")
ann=base_full.merge(mod.reset_index()[["gene","module","best_celltype_pb","best_rho_pb"]],on="gene",how="left")
ann["module_peak_region"]=ann.module.map({MODNAME[int(m[1])]:pk_true[m] for m in MODS})
for col,src in [("module_ALSvsCtrl_p","ALSvsCtrl_p"),("module_C9vssALS_p","C9vssALS_p"),
                ("module_trend_rho","trend_rho"),("module_trend_p","trend_p")]:
    ann[col]=ann.module.map({MODNAME[int(m[1])]:mstat.loc[m,src] for m in MODS})
if blood is not None: ann=ann.merge(blood[["gene","blood_mean_logCPM","blood_detectable"]],on="gene",how="left")
ann["status"]=np.where(~ann.in_spatial_panel,"not in spatial probe panel",
              np.where(~ann.detected,"below detection in spatial",""))
ann.sort_values(["module","gene"],na_position="last").to_csv(f"{TAB}/mediator_genes_annotated.csv",index=False)
print("modules:"); print(mod.module.value_counts().sort_index().to_string())
print("\nALS-vs-control significant modules:",sorted(con[con.ALSvsCtrl_p<0.05].module.unique()))
print("C9-vs-sALS significant modules:",sorted(con[con.C9vssALS_p<0.05].module.unique()) or "none")
print("stage 21 done (figures are produced by the companion notebook cells; tables written)")
