#!/usr/bin/env python
"""
06_scores_by_condition.py — Donor-level cell-type/pathway signature scores by
condition and region. Aggregates spot module scores to donor means (avoids
pseudoreplication) and plots control / sALS / C9-ALS.

Input : data/processed/GSE288365_annotated.h5ad
Output: results/figures/scores_by_condition.png
        results/tables/scores_donor_region.csv
"""
import os
os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("KMP_AFFINITY", "disabled")
import warnings; warnings.filterwarnings("ignore")
import numpy as np, pandas as pd, scanpy as sc
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

REPO="/home/sunli/C9orf72-ALS"; PROC=f"{REPO}/data/processed"
FIG=f"{REPO}/results/figures"; TAB=f"{REPO}/results/tables"
ORDER=["control","sALS","C9-ALS"]; COL={"control":"#4C72B0","sALS":"#DD8452","C9-ALS":"#C44E52"}

def main():
    a=sc.read_h5ad(f"{PROC}/GSE288365_annotated.h5ad")
    score_cols=[c for c in a.obs.columns if c.startswith("ct_")]
    df=a.obs[["donor_id","condition","region"]+score_cols].copy()
    # donor-level means, whole-section and ventral-horn-only
    glob=df.groupby(["donor_id","condition"],observed=True)[score_cols].mean().reset_index()
    vh=df[df["region"]=="ventral_horn_MN"].groupby(["donor_id","condition"],observed=True)[score_cols].mean().reset_index()
    glob.to_csv(f"{TAB}/scores_donor_region.csv",index=False)

    key=["ct_motor_neuron","ct_complement","ct_microglia_mye","ct_astrocyte"]
    fig,axes=plt.subplots(2,4,figsize=(17,8))
    from scipy.stats import mannwhitneyu
    for j,(data,scope) in enumerate([(glob,"whole section"),(vh,"ventral horn")]):
        for i,sc_col in enumerate(key):
            ax=axes[j,i]
            for k,cond in enumerate(ORDER):
                vals=data.loc[data["condition"]==cond,sc_col].values
                ax.scatter(np.full(len(vals),k)+np.random.uniform(-.08,.08,len(vals)),
                           vals,color=COL[cond],s=55,edgecolor="k",lw=.5,zorder=3)
                ax.hlines(np.mean(vals),k-.25,k+.25,color="k",lw=2,zorder=4)
            ax.set_xticks(range(3));ax.set_xticklabels(ORDER,fontsize=9)
            ax.set_title(f"{sc_col.replace('ct_','')}\n({scope})",fontsize=10)
            # C9 vs control p
            try:
                c9=data.loc[data['condition']=='C9-ALS',sc_col];ct=data.loc[data['condition']=='control',sc_col]
                p=mannwhitneyu(c9,ct).pvalue
                ax.text(.98,.02,f"C9vsCtrl p={p:.3f}",transform=ax.transAxes,ha="right",va="bottom",fontsize=8)
            except Exception: pass
    fig.suptitle("Cell-type signature scores by condition (donor-level means)",fontsize=13)
    fig.tight_layout();fig.savefig(f"{FIG}/scores_by_condition.png",dpi=150);plt.close(fig)
    print("done.")

if __name__=="__main__":
    main()
