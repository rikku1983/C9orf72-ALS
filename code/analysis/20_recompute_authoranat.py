#!/usr/bin/env python
"""
Stage 20: FULL recompute of every region-dependent analysis using the AUTHORS' manual
histology (auth_manual) instead of our marker-score regions.

Supersedes the region-stratified parts of stages 05, 06, 07, 08, 10, 11.
Reason: our marker-score 'ventral_horn_MN' was only 15% concordant with the authors'
Anterior_horns (33% white matter, 21% meninges) -- see region_annotation_disagreement.csv.

Outputs (all suffixed _authoranat):
  pseudobulk_authoranat_counts.csv            donor x region raw-count pseudobulk
  de_authoranat_<contrast>_<region>.csv       DESeq2, 3 contrasts x 4 regions
  gsea_authoranat_<contrast>_<region>.csv     Hallmark GSEA per contrast x region
  composition_authoranat_<contrast>.csv       C2L cell-type proportions per region
  protein_authoranat_<contrast>.csv           isotype-subtracted CLR protein per region
  programscores_authoranat_donor_region.csv   DAM/complement/homeostatic/NFkB per donor x region
  trend_authoranat.csv                        ordinal control<sALS<C9 trend per region
"""
import os
for v in ("OMP_NUM_THREADS","OPENBLAS_NUM_THREADS","MKL_NUM_THREADS","NUMBA_NUM_THREADS"): os.environ[v]="1"
os.environ["KMP_AFFINITY"]="disabled"
import numpy as np, pandas as pd, h5py, scipy.sparse as sp
from scipy import stats
ROOT="/home/sunli/C9orf72-ALS"; TAB=f"{ROOT}/results/tables"; PROC=f"{ROOT}/data/processed"
REGIONS=["Anterior_horns","Grey_matter","White_matter","Posterior_white_matter"]
CONTRASTS=[("C9-ALS","sALS"),("C9-ALS","control"),("sALS","control")]
MIN_SPOTS=50; MIN_DONORS=3
PROGRAMS={
 "DAM":["GPNMB","CHIT1","TREM2","TYROBP","APOE","ITGAX","CST7","LPL","CD9","CD63","SPP1","CTSB","CTSD","LGALS3","FABP5","MSR1","CD68"],
 "complement":["C1QA","C1QB","C1QC","C3"],
 "homeostatic":["P2RY12","TMEM119","CX3CR1","CSF1R","SALL1"],
 "NFkB_C9":["BIRC3","SOCS2","CALHM6"]}
CT=['Astrocytes','Endothelial','Ependymal.Cells','Lymphocytes','Macrophages','Meninges','Microglia',
    'Neurons','OPC','Oligodendrocytes','Pericytes','ProliferatingMicroglia','Schwann']
ISO=["IgG2a","IgG1","IgG2b","Rat-IgG2a"]

def dec(g):
    if isinstance(g,h5py.Group):
        c=g["categories"][:].astype(str); k=g["codes"][:]; return np.where(k>=0,c[np.clip(k,0,None)],None)
    x=g[:]; return x.astype(str) if x.dtype.kind in "SO" else x

# ---------- obs + author metadata ----------
with h5py.File(f"{PROC}/GSE288365_annotated_authmeta.h5ad","r") as h:
    o=h["obs"]; spots=o[o.attrs["_index"]][:].astype(str)
    meta=pd.DataFrame({"spot":spots,"donor_id":dec(o["donor_id"]),"condition":dec(o["condition"]),
                       "anat":dec(o["auth_manual"])})
    for c in CT: meta[f"ct_{c}"]=o[f"auth_norm_ct_{c}"][:]
    prot=h["obsm/protein"][:]; pnames=[x.decode() if isinstance(x,bytes) else str(x) for x in h["uns/protein_names"][:]]
print("spots:",len(meta),"| protein features:",len(pnames))

# ---------- protein: CLR + isotype subtraction ----------
P=pd.DataFrame(prot,columns=pnames)
clr=np.log1p(P.values); clr=clr-clr.mean(1,keepdims=True)
C=pd.DataFrame(clr,columns=pnames)
iso=[c for c in pnames if any(i.lower() in c.lower() for i in ISO)]
bg=C.loc[:,iso].mean(1).values if iso else 0.0
# PTPRC appears twice (two CD45 clones) -> average duplicate columns
MYE=sorted({c for c in pnames if c.upper() in {"CD68","CD163","ITGAX","ITGAM","CD14","HLA-DRA","PTPRC"}})
for m in MYE:
    col=C.loc[:,m]
    v=col.mean(1).values if getattr(col,"ndim",1)>1 else col.values
    meta[f"prot_{m}"]=v-bg
print("isotypes:",iso,"| myeloid protein:",MYE)

# ---------- program scores (lognorm X) ----------
with h5py.File(f"{PROC}/GSE288365_annotated.h5ad","r") as h:
    genes=h["var"][h["var"].attrs["_index"]][:].astype(str); X=h["X"]; ip=X["indptr"][:]
    allg=sorted({g for v in PROGRAMS.values() for g in v})
    have=[g for g in allg if g in genes]; cols=np.array([np.where(genes==g)[0][0] for g in have])
    n=len(ip)-1; E=np.zeros((n,len(cols)),np.float32)
    for s in range(0,n,20000):
        e=min(s+20000,n); a0,a1=ip[s],ip[e]
        m=sp.csr_matrix((X["data"][a0:a1],X["indices"][a0:a1],ip[s:e+1]-a0),shape=(e-s,len(genes)))
        E[s:e]=m[:,cols].toarray()
Edf=pd.DataFrame(E,columns=have)
z=(Edf-Edf.mean())/Edf.std()
for nm,gl in PROGRAMS.items():
    g=[x for x in gl if x in have]; meta[f"prog_{nm}"]=z[g].mean(1).values
for g in ["GPNMB","CHIT1","APOE","C1QA","C1QC","TREM2"]:
    if g in have: meta[f"gene_{g}"]=Edf[g].values

# ---------- donor x region summaries ----------
sub=meta[meta.anat.isin(REGIONS)]
valcols=[c for c in meta.columns if c.startswith(("ct_","prot_","prog_","gene_"))]
dr=sub.groupby(["donor_id","condition","anat"],observed=True)[valcols].mean().reset_index()
dr["n_spots"]=sub.groupby(["donor_id","condition","anat"],observed=True).size().values
dr=dr[dr.n_spots>=MIN_SPOTS]
dr.to_csv(f"{TAB}/programscores_authoranat_donor_region.csv",index=False)

def mwu(a,b):
    a=np.asarray(a,float); b=np.asarray(b,float)
    if len(a)<2 or len(b)<2: return np.nan
    return stats.mannwhitneyu(a,b).pvalue

# composition + protein + program contrasts per region
rows=[]
for reg in REGIONS:
    q=dr[dr.anat==reg]
    for A,B in CONTRASTS+[("ALS","control")]:
        gA=q[q.condition.isin(["C9-ALS","sALS"])] if A=="ALS" else q[q.condition==A]
        gB=q[q.condition==B]
        for c in valcols:
            p=mwu(gA[c],gB[c])
            rows.append({"region":reg,"contrast":f"{A}_vs_{B}","feature":c,
                         "mean_A":gA[c].mean(),"mean_B":gB[c].mean(),"diff":gA[c].mean()-gB[c].mean(),
                         "p":p,"nA":len(gA),"nB":len(gB)})
con=pd.DataFrame(rows)
con[con.feature.str.startswith("ct_")].to_csv(f"{TAB}/composition_authoranat_allcontrasts.csv",index=False)
con[con.feature.str.startswith("prot_")].to_csv(f"{TAB}/protein_authoranat_allcontrasts.csv",index=False)
con[con.feature.str.startswith(("prog_","gene_"))].to_csv(f"{TAB}/programs_authoranat_allcontrasts.csv",index=False)

# ordinal trend
rank={"control":0,"sALS":1,"C9-ALS":2}
tr=[]
for reg in REGIONS:
    q=dr[dr.anat==reg].copy(); q["ord"]=q.condition.map(rank)
    for c in valcols:
        rho,p=stats.spearmanr(q["ord"],q[c]); tr.append({"region":reg,"feature":c,"rho":rho,"p":p,"n":len(q)})
pd.DataFrame(tr).to_csv(f"{TAB}/trend_authoranat.csv",index=False)
print("composition/protein/program/trend tables written")

# ---------- pseudobulk raw counts per donor x region ----------
with h5py.File(f"{PROC}/GSE288365_qc.h5ad","r") as h:
    qs=h["obs"][h["obs"].attrs["_index"]][:].astype(str)
    qgenes=h["var"][h["var"].attrs["_index"]][:].astype(str)
    pos=pd.Series(np.arange(len(qs)),index=qs).reindex(meta.spot).values.astype(int)
    X=h["X"]; ip=X["indptr"][:]
    keys=[(d,r,c) for (d,r,c),g in sub.groupby(["donor_id","anat","condition"],observed=True) if len(g)>=MIN_SPOTS]
    gi={k:i for i,k in enumerate(keys)}
    PB=np.zeros((len(keys),len(qgenes)))
    kk=list(zip(meta.donor_id,meta.anat,meta.condition))
    bucket={}
    for i,k in enumerate(kk):
        if k in gi: bucket.setdefault(gi[k],[]).append(pos[i])
    for g_,rws in bucket.items():
        for rr in np.sort(np.array(rws)):
            a0,a1=ip[rr],ip[rr+1]; np.add.at(PB[g_],X["indices"][a0:a1],X["data"][a0:a1])
pb=pd.DataFrame(PB,columns=qgenes)
gm=pd.DataFrame(keys,columns=["donor_id","anat","condition"]); gm["n_spots"]=[len(bucket[i]) for i in range(len(keys))]
pb.index=gm.donor_id+"|"+gm.anat
pb.to_csv(f"{TAB}/pseudobulk_authoranat_counts.csv"); gm.to_csv(f"{TAB}/pseudobulk_authoranat_groups.csv",index=False)
print("pseudobulk:",pb.shape); print(gm.groupby(["anat","condition"]).size().unstack(fill_value=0).to_string())

# ---------- DESeq2 + GSEA, 3 contrasts x 4 regions ----------
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats
import gseapy as gp
summ=[]; gsea_all=[]
for reg in REGIONS:
    m=gm[gm.anat==reg]
    for A,B in CONTRASTS:
        s=m[m.condition.isin([A,B])]
        if s.groupby("condition").size().min()<MIN_DONORS: continue
        cnt=pb.loc[(s.donor_id+"|"+s.anat).values].round().astype(int)
        cnt=cnt.loc[:,(cnt>=10).sum(0)>=3]
        md=s[["condition"]].copy(); md.index=cnt.index
        dds=DeseqDataSet(counts=cnt,metadata=md,design="~condition",quiet=True); dds.deseq2()
        res=DeseqStats(dds,contrast=["condition",A,B],quiet=True); res.summary()
        out=res.results_df.sort_values("padj")
        tag=f"{A.replace('-','')}_vs_{B}"
        out.to_csv(f"{TAB}/de_authoranat_{tag}_{reg}.csv")
        nsig=int((out.padj<0.05).sum())
        summ.append({"region":reg,"contrast":f"{A} vs {B}","n_genes":cnt.shape[1],"n_FDR05":nsig,
                     "nA":int((s.condition==A).sum()),"nB":int((s.condition==B).sum()),
                     "top_gene":out.index[0],"top_padj":out.padj.iloc[0]})
        r=gp.prerank(rnk=out.dropna(subset=["stat"])["stat"].sort_values(ascending=False),
                     gene_sets="MSigDB_Hallmark_2020",permutation_num=1000,seed=0,threads=1,
                     min_size=10,max_size=500,outdir=None,verbose=False).res2d
        r=r.rename(columns={"Term":"term","NES":"nes","FDR q-val":"fdr","NOM p-val":"pval"})
        r["nes"]=pd.to_numeric(r["nes"]); r["fdr"]=pd.to_numeric(r["fdr"])
        r["region"]=reg; r["contrast"]=f"{A} vs {B}"
        gsea_all.append(r[["region","contrast","term","nes","pval","fdr"]])
        print(f"{reg:24s} {A} vs {B}: genes={cnt.shape[1]} FDR05={nsig} | top Hallmark: "
              f"{r.loc[r.nes.abs().idxmax(),'term']} NES={r.loc[r.nes.abs().idxmax(),'nes']:+.2f} FDR={r.loc[r.nes.abs().idxmax(),'fdr']:.3f}")
pd.DataFrame(summ).to_csv(f"{TAB}/de_authoranat_summary.csv",index=False)
pd.concat(gsea_all).to_csv(f"{TAB}/gsea_authoranat_all.csv",index=False)
print("\nDONE stage 20")
