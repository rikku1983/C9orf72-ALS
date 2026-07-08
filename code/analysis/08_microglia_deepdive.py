#!/usr/bin/env python
"""
08_microglia_deepdive.py — GPNMB / CHIT1 and the microglia program in C9-ALS.

Analyses:
  (A) DAM (disease-associated microglia) vs homeostatic microglia module scores
  (B) GPNMB / CHIT1 expression by condition x region (donor-level means)
  (C) Co-expression of GPNMB/CHIT1 with the DAM program (spot-level correlation)
  (D) Microglia state shift: DAM/homeostatic ratio by condition
  (E) Spatial localization: which anatomical region carries the signal

Input : data/processed/GSE288365_annotated.h5ad  (X = log-normalized)
Output: results/figures/microglia_deepdive.png
        results/tables/microglia_gene_by_condition_region.csv
        results/tables/microglia_dam_correlations.csv
"""
import os
os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("KMP_AFFINITY", "disabled")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
import warnings; warnings.filterwarnings("ignore")
import numpy as np, pandas as pd, scanpy as sc
from scipy.stats import mannwhitneyu, spearmanr
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

REPO="/home/sunli/C9orf72-ALS"; PROC=f"{REPO}/data/processed"
FIG=f"{REPO}/results/figures"; TAB=f"{REPO}/results/tables"
ORDER=["control","sALS","C9-ALS"]; COL={"control":"#4C72B0","sALS":"#DD8452","C9-ALS":"#C44E52"}

# DAM / lipid-droplet-associated / phagocytic microglia program (human ALS + neurodegen literature)
DAM=["GPNMB","CHIT1","TREM2","TYROBP","APOE","ITGAX","CST7","LPL","CD9","CD63",
     "SPP1","CTSB","CTSD","LGALS3","FABP5","MSR1","CD68"]
HOMEO=["P2RY12","TMEM119","CX3CR1","CSF1R","SALL1"]
FOCUS=["GPNMB","CHIT1"]

def main():
    a=sc.read_h5ad(f"{PROC}/GSE288365_annotated.h5ad")
    dam=[g for g in DAM if g in a.var_names]; homeo=[g for g in HOMEO if g in a.var_names]
    sc.tl.score_genes(a,dam,score_name="DAM_score",use_raw=False)
    sc.tl.score_genes(a,homeo,score_name="homeo_score",use_raw=False)

    # ---- (B) donor-level gene expression by condition x region ----
    genes=FOCUS+["DAM_score","homeo_score"]
    expr=pd.DataFrame(a[:,FOCUS].X.toarray(),columns=FOCUS,index=a.obs_names)
    expr["DAM_score"]=a.obs["DAM_score"].values; expr["homeo_score"]=a.obs["homeo_score"].values
    expr["donor_id"]=a.obs["donor_id"].values; expr["condition"]=a.obs["condition"].values
    expr["region"]=a.obs["region"].values
    donor_glob=expr.groupby(["donor_id","condition"],observed=True)[genes].mean().reset_index()
    donor_vh=expr[expr.region=="ventral_horn_MN"].groupby(["donor_id","condition"],observed=True)[genes].mean().reset_index()
    rec=[]
    for scope,dat in [("global",donor_glob),("ventral_horn_MN",donor_vh)]:
        for g in genes:
            row={"gene":g,"scope":scope}
            for cond in ORDER:
                row[f"{cond}_mean"]=dat.loc[dat.condition==cond,g].mean()
            try: row["p_C9_vs_ctrl"]=mannwhitneyu(dat.loc[dat.condition=='C9-ALS',g],dat.loc[dat.condition=='control',g]).pvalue
            except Exception: row["p_C9_vs_ctrl"]=np.nan
            try: row["p_C9_vs_sALS"]=mannwhitneyu(dat.loc[dat.condition=='C9-ALS',g],dat.loc[dat.condition=='sALS',g]).pvalue
            except Exception: row["p_C9_vs_sALS"]=np.nan
            rec.append(row)
    pd.DataFrame(rec).to_csv(f"{TAB}/microglia_gene_by_condition_region.csv",index=False)

    # ---- (C) spot-level co-expression of GPNMB/CHIT1 with DAM markers (in ALS spots) ----
    als=a[a.obs.condition.isin(["C9-ALS","sALS"])]
    Xd=pd.DataFrame(als[:,dam].X.toarray(),columns=dam)
    corr=[]
    for f in FOCUS:
        for g in dam:
            if g==f: continue
            rho,p=spearmanr(Xd[f],Xd[g])
            corr.append({"focus":f,"partner":g,"spearman_rho":rho,"p":p})
    cdf=pd.DataFrame(corr); cdf.to_csv(f"{TAB}/microglia_dam_correlations.csv",index=False)

    # ---- (E) region localization (donor-level, C9-ALS only) ----
    reg_order=["ventral_horn_MN","dorsal_horn","gray_neuropil","white_matter"]
    c9=expr[expr.condition=="C9-ALS"]
    reg_means=c9.groupby("region",observed=True)[FOCUS+["DAM_score"]].mean().reindex(reg_order)

    # =================== FIGURE ===================
    fig=plt.figure(figsize=(18,9))
    gs=fig.add_gridspec(2,4,hspace=.42,wspace=.36)

    # Row 1: GPNMB, CHIT1, DAM_score, homeo_score by condition (ventral horn donor means)
    for i,g in enumerate(["GPNMB","CHIT1","DAM_score","homeo_score"]):
        ax=fig.add_subplot(gs[0,i])
        for k,cond in enumerate(ORDER):
            vals=donor_vh.loc[donor_vh.condition==cond,g].values
            ax.scatter(np.full(len(vals),k)+np.random.uniform(-.08,.08,len(vals)),vals,
                       color=COL[cond],s=55,edgecolor="k",lw=.5,zorder=3)
            ax.hlines(np.mean(vals),k-.25,k+.25,color="k",lw=2,zorder=4)
        ax.set_xticks(range(3));ax.set_xticklabels(ORDER,fontsize=8,rotation=15)
        ax.set_title(f"{g}  (ventral horn)",fontsize=10)
        try:
            p=mannwhitneyu(donor_vh.loc[donor_vh.condition=='C9-ALS',g],donor_vh.loc[donor_vh.condition=='control',g]).pvalue
            ax.text(.97,.03,f"C9vsCtrl p={p:.3f}",transform=ax.transAxes,ha="right",va="bottom",fontsize=8)
        except Exception: pass

    # Row 2, panel 0-1: co-expression bars
    ax=fig.add_subplot(gs[1,0:2])
    cd=cdf.pivot(index="partner",columns="focus",values="spearman_rho").drop(index=[g for g in FOCUS if g in cdf.partner.unique()],errors="ignore")
    cd=cd.sort_values("GPNMB",ascending=True)
    y=np.arange(len(cd)); w=.4
    ax.barh(y-w/2,cd["GPNMB"],w,label="GPNMB",color="#C44E52",edgecolor="k",lw=.4)
    ax.barh(y+w/2,cd["CHIT1"],w,label="CHIT1",color="#8172B3",edgecolor="k",lw=.4)
    ax.set_yticks(y);ax.set_yticklabels(cd.index,fontsize=8)
    ax.set_xlabel("Spearman rho (ALS spots)",fontsize=9)
    ax.set_title("Co-expression of GPNMB/CHIT1 with DAM program",fontsize=10)
    ax.legend(fontsize=8);ax.axvline(0,color="k",lw=.6)

    # Row 2, panel 2: DAM vs homeostatic ratio by condition
    ax=fig.add_subplot(gs[1,2])
    donor_vh["log_ratio"]=np.log2((donor_vh["DAM_score"]-donor_vh["DAM_score"].min()+.01)/
                                  (donor_vh["homeo_score"]-donor_vh["homeo_score"].min()+.01))
    for k,cond in enumerate(ORDER):
        vals=donor_vh.loc[donor_vh.condition==cond,"log_ratio"].values
        ax.scatter(np.full(len(vals),k)+np.random.uniform(-.08,.08,len(vals)),vals,
                   color=COL[cond],s=55,edgecolor="k",lw=.5,zorder=3)
        ax.hlines(np.mean(vals),k-.25,k+.25,color="k",lw=2,zorder=4)
    ax.set_xticks(range(3));ax.set_xticklabels(ORDER,fontsize=8,rotation=15)
    ax.set_title("DAM / homeostatic balance\n(ventral horn)",fontsize=10)
    ax.set_ylabel("log2 DAM:homeo (shifted)",fontsize=8)

    # Row 2, panel 3: region localization heatmap (C9-ALS)
    ax=fig.add_subplot(gs[1,3])
    hm=reg_means[FOCUS].T
    im=ax.imshow(hm.values,cmap="Reds",aspect="auto")
    ax.set_xticks(range(len(reg_order)));ax.set_xticklabels([r.replace('_','\n') for r in reg_order],fontsize=7)
    ax.set_yticks(range(len(FOCUS)));ax.set_yticklabels(FOCUS,fontsize=9)
    ax.set_title("Regional expression\n(C9-ALS, mean log-norm)",fontsize=10)
    for ii in range(hm.shape[0]):
        for jj in range(hm.shape[1]):
            ax.text(jj,ii,f"{hm.values[ii,jj]:.2f}",ha="center",va="center",fontsize=7,
                    color="white" if hm.values[ii,jj]>hm.values.max()*.6 else "black")
    fig.colorbar(im,ax=ax,fraction=.046,pad=.04)

    fig.suptitle("GPNMB / CHIT1 and the microglia program in C9-ALS spinal cord",fontsize=13,y=.98)
    fig.savefig(f"{FIG}/microglia_deepdive.png",dpi=150,bbox_inches="tight");plt.close(fig)

    # console summary
    print("=== GPNMB/CHIT1/DAM/homeo by condition (ventral horn donor means) ===")
    print(pd.DataFrame(rec).query("scope=='ventral_horn_MN'")[["gene","control_mean","sALS_mean","C9-ALS_mean","p_C9_vs_ctrl","p_C9_vs_sALS"]].to_string(index=False))
    print("\n=== top DAM co-expression partners (mean of GPNMB,CHIT1 rho) ===")
    top=cdf.groupby("partner")["spearman_rho"].mean().sort_values(ascending=False).head(10)
    print(top.to_string())
    print("\n=== regional means (C9-ALS) ===")
    print(reg_means.to_string())

if __name__=="__main__":
    main()
