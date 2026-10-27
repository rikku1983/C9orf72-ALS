#!/usr/bin/env python
"""
Stage 18 (spatial arm): spatial distribution of the DAM program, GPNMB and CHIT1 in
control / sALS / C9-ALS spinal cord, using the AUTHORS' manual histological annotation
(auth_manual) instead of our marker-score regions (which mislabel white-matter rim spots).

DAM core = mean z-score of the 15-gene DAM set with GPNMB and CHIT1 REMOVED (non-circular).
All condition / region statistics are at DONOR level (spot = pseudoreplicate).
Representative section per condition = section whose mean DAM core is closest to the
condition median of section means (fixed rule, not hand-picked).

Outputs: results/figures/spatial_DAM_GPNMB_CHIT1_{maps,quant}.png
         results/tables/spatial_DAMcore_GPNMB_CHIT1_donor_region.csv
         results/tables/spatial_DAMcore_GPNMB_CHIT1_stats.csv
         results/tables/spatial_GPNMB_CHIT1_DAMcore_within_section_rho.csv
"""
import os
for v in ("OMP_NUM_THREADS","OPENBLAS_NUM_THREADS","MKL_NUM_THREADS","NUMBA_NUM_THREADS"): os.environ[v]="1"
os.environ["KMP_AFFINITY"]="disabled"; os.environ["OMP_PROC_BIND"]="FALSE"
import numpy as np, pandas as pd, h5py, scipy.sparse as sp
from scipy import stats
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm
from matplotlib.patches import Patch
from matplotlib.lines import Line2D

ROOT="/home/sunli/C9orf72-ALS"; FIG=f"{ROOT}/results/figures"; TAB=f"{ROOT}/results/tables"
F=f"{ROOT}/data/processed/GSE288365_annotated.h5ad"
FA=f"{ROOT}/data/processed/GSE288365_annotated_authmeta.h5ad"
DAM17=['GPNMB','CHIT1','TREM2','TYROBP','APOE','ITGAX','CST7','LPL','CD9','CD63','SPP1','CTSB','CTSD','LGALS3','FABP5','MSR1','CD68']
DAMCORE=[g for g in DAM17 if g not in ("GPNMB","CHIT1")]
ORDER=["control","sALS","C9-ALS"]; FEAT=["DAM_core","GPNMB","CHIT1"]
PAL={"control":"#4477AA","sALS":"#DDAA33","C9-ALS":"#CC6677"}
ANAT=["Anterior_horns","Grey_matter","White_matter","Posterior_white_matter","Meninges"]
REG4=ANAT[:4]
ALAB={"Anterior_horns":"Anterior horn","Grey_matter":"Other grey matter","White_matter":"Anterolateral WM",
      "Posterior_white_matter":"Posterior WM (dorsal columns)","Meninges":"Meninges"}
ACOL={"Anterior_horns":"#332288","Grey_matter":"#44AA99","White_matter":"#882255","Posterior_white_matter":"#88CCEE","Meninges":"#BBBBBB"}

def dec(g):
    if isinstance(g,h5py.Group):
        cats=g["categories"][:].astype(str); codes=g["codes"][:]
        return np.where(codes>=0,cats[np.clip(codes,0,None)],None)
    x=g[:]; return x.astype(str) if x.dtype.kind in "SO" else x

# ---- extract genes per spot (chunked; low memory) ----
with h5py.File(F,"r") as h:
    genes=h["var"][h["var"].attrs["_index"]][:].astype(str); obs=h["obs"]
    df=pd.DataFrame({"spot":obs[obs.attrs["_index"]][:].astype(str)})
    for c in ["section","condition","donor_id"]: df[c]=dec(obs[c])
    xy=h["obsm/spatial"][:]; X=h["X"]; cols=np.array([np.where(genes==g)[0][0] for g in DAM17])
    indptr=X["indptr"][:]; n=len(indptr)-1; out=np.zeros((n,len(cols)),np.float32)
    for s in range(0,n,20000):
        e=min(s+20000,n); a0,a1=indptr[s],indptr[e]
        m=sp.csr_matrix((X["data"][a0:a1],X["indices"][a0:a1],indptr[s:e+1]-a0),shape=(e-s,len(genes)))
        out[s:e]=m[:,cols].toarray()
df=pd.concat([df,pd.DataFrame(out,columns=DAM17)],axis=1); df["x"]=xy[:,0]; df["y"]=xy[:,1]
z=(df[DAMCORE]-df[DAMCORE].mean())/df[DAMCORE].std(); df["DAM_core"]=z.mean(axis=1)
with h5py.File(FA,"r") as h:
    obs=h["obs"]; names=obs[obs.attrs["_index"]][:].astype(str)
    assert (names==df.spot.values).all()
    df["anat"]=dec(obs["auth_manual"]); df["macro"]=obs["auth_norm_ct_Macrophages"][:]; df["micro"]=obs["auth_norm_ct_Microglia"][:]

# ---- donor x anatomy means + stats ----
d2=df[df.anat.isin(ANAT)]; FEAT2=FEAT+["macro","micro"]
dra=d2.groupby(["donor_id","condition","anat"],observed=True)[FEAT2].mean().reset_index()
dra.to_csv(f"{TAB}/spatial_DAMcore_GPNMB_CHIT1_donor_region.csv",index=False)
rows=[]
for f in FEAT2:
    p=dra[dra.condition!="control"].pivot_table(index="donor_id",columns="anat",values=f)
    for r in ["Anterior_horns","Grey_matter","Posterior_white_matter","Meninges"]:
        d=(p["White_matter"]-p[r]).dropna()
        rows.append({"test":"ALS donors: anterolateral WM minus region","feature":f,"region":r,
                     "mean_diff":d.mean(),"n_pos":int((d>0).sum()),"n":len(d),"p":stats.wilcoxon(d).pvalue})
    for r in REG4:
        q=dra[dra.anat==r]; g=lambda c:q.loc[q.condition==c,f]
        rows.append({"test":"ALS vs control (MWU donors)","feature":f,"region":r,"mean_diff":pd.concat([g('sALS'),g('C9-ALS')]).mean()-g('control').mean(),
                     "n":len(q),"p":stats.mannwhitneyu(pd.concat([g('sALS'),g('C9-ALS')]),g('control')).pvalue})
        rows.append({"test":"C9 vs sALS (MWU donors)","feature":f,"region":r,"mean_diff":g('C9-ALS').mean()-g('sALS').mean(),
                     "n":len(g('C9-ALS'))+len(g('sALS')),"p":stats.mannwhitneyu(g('C9-ALS'),g('sALS')).pvalue})
st=pd.DataFrame(rows)
ws=pd.DataFrame([{"section":s,"condition":c,"gene":f,"rho":stats.spearmanr(g[f],g["DAM_core"])[0],"n_spots":len(g)}
                 for (s,c),g in df.groupby(["section","condition"],observed=True) for f in ["GPNMB","CHIT1"]])
ws.to_csv(f"{TAB}/spatial_GPNMB_CHIT1_DAMcore_within_section_rho.csv",index=False)
for gname in ["GPNMB","CHIT1"]:
    w=ws[ws.gene==gname]
    st.loc[len(st)]={"test":"within-section rho vs DAM core: ALS vs control sections (MWU)","feature":gname,"region":"whole section",
                     "mean_diff":w[w.condition!="control"].rho.median()-w[w.condition=="control"].rho.median(),"n":len(w),
                     "p":stats.mannwhitneyu(w[w.condition!="control"].rho,w[w.condition=="control"].rho).pvalue}
q=dra[dra.anat=="White_matter"]
for lab,qq in [("all donors",q),("ALS donors",q[q.condition!="control"])]:
    for gname in ["GPNMB","CHIT1"]:
        rho,p=stats.spearmanr(qq.DAM_core,qq[gname])
        st.loc[len(st)]={"test":f"donor-level Spearman vs DAM core, anterolateral WM ({lab})","feature":gname,"region":"White_matter","mean_diff":rho,"n":len(qq),"p":p}
st.to_csv(f"{TAB}/spatial_DAMcore_GPNMB_CHIT1_stats.csv",index=False)

# ---- representative sections ----
sec=df.groupby(["section","condition"],observed=True)["DAM_core"].mean().reset_index()
rep={c:sec[sec.condition==c].iloc[(sec[sec.condition==c].DAM_core-sec[sec.condition==c].DAM_core.median()).abs().argmin()]["section"] for c in ORDER}
donor_of=df.drop_duplicates("section").set_index("section")["donor_id"].to_dict()
print("representative:",{c:donor_of[s] for c,s in rep.items()})

# ---- Figure 1: maps ----
try: apply_figure_style(sizes=(8,7,6))
except NameError: plt.rcParams.update({"font.size":8,"axes.titlesize":8,"axes.labelsize":8,"xtick.labelsize":6,"ytick.labelsize":6,"legend.fontsize":7,"axes.spines.top":False,"axes.spines.right":False})
sub={c:df[df.section==rep[c]] for c in ORDER}; allrep=pd.concat(sub.values())
lim={"DAM_core":np.percentile(np.abs(allrep.DAM_core),99),"GPNMB":np.percentile(allrep.GPNMB,99.5),"CHIT1":np.percentile(allrep.CHIT1[allrep.CHIT1>0],90)}
fig,axes=plt.subplots(4,3,figsize=(7.2,8.6)); fig.subplots_adjust(left=0.1,right=0.88,top=0.94,bottom=0.07,hspace=0.08,wspace=0.04)
mp={}; s=4.0
for j,c in enumerate(ORDER):
    d=sub[c]
    for i,f in enumerate(["anat","DAM_core","GPNMB","CHIT1"]):
        ax=axes[i,j]
        if f=="anat": ax.scatter(d.x,d.y,c=d.anat.map(ACOL).fillna("#F2F2F2"),s=s,linewidths=0,rasterized=True)
        elif f=="DAM_core":
            o=d.sort_values(f,key=np.abs); mp[f]=ax.scatter(o.x,o.y,c=o[f],cmap="RdBu_r",norm=TwoSlopeNorm(0,-lim[f],lim[f]),s=s,linewidths=0,rasterized=True)
        else:
            zz=d[d[f]==0]; nz=d[d[f]>0].sort_values(f)
            ax.scatter(zz.x,zz.y,c="#E6E6E6",s=s,linewidths=0,rasterized=True)
            mp[f]=ax.scatter(nz.x,nz.y,c=nz[f],cmap="magma_r",vmin=0,vmax=lim[f],s=s*1.3,linewidths=0,rasterized=True)
        ax.set_aspect("equal"); ax.invert_yaxis(); ax.set_xticks([]); ax.set_yticks([])
        for sp_ in ax.spines.values(): sp_.set_visible(False)
    axes[0,j].set_title(f"{c}\n{donor_of[rep[c]]}",color=PAL[c],loc="center")
for i,lab in enumerate(["Histological region","DAM core score","GPNMB (log-norm)","CHIT1 (log-norm)"]): axes[i,0].set_ylabel(lab)
for i,f in enumerate(["DAM_core","GPNMB","CHIT1"],start=1):
    b=axes[i,2].get_position(); cb=fig.colorbar(mp[f],cax=fig.add_axes([0.9,b.y0+0.02,0.012,b.height-0.04])); cb.ax.tick_params(labelsize=6)
fig.legend(handles=[Patch(color=ACOL[r],label=ALAB[r]) for r in ANAT],loc="lower center",bbox_to_anchor=(0.49,0.005),ncol=3,frameon=False,fontsize=7)
fig.savefig(f"{FIG}/spatial_DAM_GPNMB_CHIT1_maps.png",dpi=300); plt.close(fig)

# ---- Figure 2: donor-level quantification ----
fig2,ax2=plt.subplots(2,3,figsize=(7.2,5.0)); fig2.subplots_adjust(left=0.09,right=0.98,top=0.93,bottom=0.1,hspace=0.75,wspace=0.42)
off={"control":-0.25,"sALS":0,"C9-ALS":0.25}; rng=np.random.default_rng(0)
YL={"DAM_core":"DAM core score","GPNMB":"GPNMB (log-norm)","CHIT1":"CHIT1 (log-norm)","macro":"Macrophages (% of spot cells)"}
TT=["Highest in anterolateral WM","GPNMB: spares dorsal columns","CHIT1: same pattern","Macrophages track the signal"]
for k,f in enumerate(["DAM_core","GPNMB","CHIT1","macro"]):
    ax=ax2.flat[k]
    for i,rg in enumerate(REG4):
        for c in ORDER:
            v=dra[(dra.anat==rg)&(dra.condition==c)][f].values; x=i+off[c]
            ax.scatter(x+rng.uniform(-0.05,0.05,len(v)),v,s=10,color=PAL[c],alpha=0.85,linewidths=0,zorder=3)
            ax.hlines(np.median(v),x-0.1,x+0.1,color="black",lw=1.2,zorder=4)
    ax.set_xticks(range(4)); ax.set_xticklabels(["Anterior horn","Other grey","Anterolat. WM","Posterior WM"],rotation=35,ha="right",rotation_mode="anchor",fontsize=6)
    ax.set_ylabel(YL[f]); ax.set_title(TT[k],loc="left"); ax.margins(0.05)
ax2.flat[0].legend(handles=[Line2D([0],[0],marker="o",lw=0,color=PAL[c],label=c,ms=4) for c in ORDER],frameon=False,fontsize=6,loc="upper left",handletextpad=0.2)
ax=ax2.flat[4]
for i,c in enumerate(ORDER):
    v=ws[(ws.gene=="GPNMB")&(ws.condition==c)].rho.values
    ax.scatter(i+rng.uniform(-0.08,0.08,len(v)),v,s=14,color=PAL[c],linewidths=0,zorder=3); ax.hlines(np.median(v),i-0.18,i+0.18,color="black",lw=1.2)
ax.axhline(0,color="grey",lw=0.6,ls=":"); ax.set_xticks(range(3)); ax.set_xticklabels(ORDER,fontsize=6)
ax.set_ylabel("Spearman ρ, GPNMB vs\nDAM core (per section)"); ax.set_title("Coupling stronger in ALS",loc="left"); ax.margins(0.08)
ax=ax2.flat[5]
for c in ORDER:
    qq=q[q.condition==c]; ax.scatter(qq.DAM_core,qq.GPNMB,s=18,color=PAL[c],linewidths=0,zorder=3)
r_all=stats.spearmanr(q.DAM_core,q.GPNMB)[0]; r_als=stats.spearmanr(q[q.condition!="control"].DAM_core,q[q.condition!="control"].GPNMB)[0]
ax.text(0.03,0.97,f"all ρ={r_all:.2f}\nALS only ρ={r_als:.2f}",transform=ax.transAxes,va="top",fontsize=7)
ax.set_xlabel("DAM core score"); ax.set_ylabel("GPNMB (log-norm)"); ax.set_title("Donor level, anterolateral WM",loc="left"); ax.margins(0.08)
for a_,L in zip(ax2.flat,"abcdef"): a_.text(-0.22,1.08,L,transform=a_.transAxes,fontweight="bold",fontsize=9,va="bottom")
fig2.savefig(f"{FIG}/spatial_DAM_GPNMB_CHIT1_quant.png",dpi=300); plt.close(fig2)
print("stage 18 done")
