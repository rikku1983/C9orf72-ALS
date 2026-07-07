#!/usr/bin/env python
"""
03_annotate_regions.py — Annotate spinal-cord anatomical regions on Visium
spots via marker-gene scoring, and score cell-type signatures per spot.

Input : data/processed/GSE288365_clustered.h5ad
Output: data/processed/GSE288365_annotated.h5ad
        results/figures/region_scores_umap.png
        results/figures/spatial_regions_<section>.png (a few example sections)
        results/tables/region_by_cluster.csv

Spinal cord regions (Visium spots are ~55um = multi-cell, so these are
dominant-signal region labels, not pure cell types):
- Gray matter — ventral horn (motor neurons): CHAT, ISL1, MNX1, PRPH, SLC5A7
- Gray matter — dorsal horn (sensory):        PAX2, TLX3, SLC17A6, POU4F1
- White matter (myelinated tracts):           MBP, PLP1, MOBP, MOG, MAG
- Gray matter neuropil (general neuronal):     RBFOX3, SNAP25, SYT1, MAP2

Cell-type signatures scored per spot (module scores):
- Motor neuron:   CHAT, ISL1, MNX1, PRPH, SLC5A7, NEFH, NEFL
- Astrocyte:      GFAP, AQP4, SLC1A2, SLC1A3, ALDH1L1
- Microglia/mye:  AIF1, CSF1R, C1QA, C1QB, C1QC, TYROBP, TREM2, P2RY12, CX3CR1
- Oligodendro:    MBP, PLP1, MOBP, MOG, MAG, CLDN11
- OPC:            PDGFRA, CSPG4, OLIG1, OLIG2
- Endothelial:    CLDN5, PECAM1, FLT1, VWF
- Complement:     C1QA, C1QB, C1QC, C3, C4A, C4B, CR1, ITGAM, ITGAX  (progression-linked)
"""
import os
os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("KMP_AFFINITY", "disabled")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
import warnings; warnings.filterwarnings("ignore")
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

REPO = "/home/sunli/C9orf72-ALS"
PROC = f"{REPO}/data/processed"
FIG  = f"{REPO}/results/figures"
TAB  = f"{REPO}/results/tables"

REGION_MARKERS = {
    "ventral_horn_MN": ["CHAT", "ISL1", "MNX1", "PRPH", "SLC5A7"],
    "dorsal_horn":     ["PAX2", "TLX3", "SLC17A6", "POU4F1"],
    "white_matter":    ["MBP", "PLP1", "MOBP", "MOG", "MAG"],
    "gray_neuropil":   ["RBFOX3", "SNAP25", "SYT1", "MAP2"],
}
CELLTYPE_SIGS = {
    "motor_neuron":  ["CHAT", "ISL1", "MNX1", "PRPH", "SLC5A7", "NEFH", "NEFL"],
    "astrocyte":     ["GFAP", "AQP4", "SLC1A2", "SLC1A3", "ALDH1L1"],
    "microglia_mye": ["AIF1", "CSF1R", "C1QA", "C1QB", "C1QC", "TYROBP",
                      "TREM2", "P2RY12", "CX3CR1"],
    "oligodendro":   ["MBP", "PLP1", "MOBP", "MOG", "MAG", "CLDN11"],
    "OPC":           ["PDGFRA", "CSPG4", "OLIG1", "OLIG2"],
    "endothelial":   ["CLDN5", "PECAM1", "FLT1", "VWF"],
    "complement":    ["C1QA", "C1QB", "C1QC", "C3", "C4A", "C4B",
                      "CR1", "ITGAM", "ITGAX"],
}

def score_sets(a, sets, prefix=""):
    present = {}
    for name, genes in sets.items():
        g = [x for x in genes if x in a.var_names]
        if not g:
            print(f"  [warn] {name}: no markers present"); continue
        sc.tl.score_genes(a, g, score_name=f"{prefix}{name}", use_raw=False)
        present[name] = g
    return present

def main():
    a = sc.read_h5ad(f"{PROC}/GSE288365_clustered.h5ad")
    print(f"[load] {a.shape}, {a.obs['leiden'].nunique()} clusters")

    reg = score_sets(a, REGION_MARKERS, prefix="reg_")
    cts = score_sets(a, CELLTYPE_SIGS, prefix="ct_")

    # assign each spot to the region with the highest z-scored region score
    reg_cols = [f"reg_{k}" for k in reg]
    Z = a.obs[reg_cols].apply(lambda c: (c - c.mean()) / (c.std() + 1e-9))
    a.obs["region"] = (Z.idxmax(axis=1).str.replace("reg_", "", regex=False)
                       .astype("category"))
    print("[region] spot counts:\n", a.obs["region"].value_counts())

    # region vs cluster crosstab
    ct = pd.crosstab(a.obs["leiden"], a.obs["region"])
    ct.to_csv(f"{TAB}/region_by_cluster.csv")

    a.write(f"{PROC}/GSE288365_annotated.h5ad")

    # UMAP of region + key cell-type scores
    panels = ["region", "ct_motor_neuron", "ct_microglia_mye",
              "ct_complement", "ct_astrocyte", "reg_white_matter"]
    fig, axes = plt.subplots(2, 3, figsize=(16, 9))
    for ax, p in zip(axes.ravel(), panels):
        sc.pl.umap(a, color=p, ax=ax, show=False, size=3,
                   cmap=None if p == "region" else "viridis")
    fig.tight_layout(); fig.savefig(f"{FIG}/region_scores_umap.png", dpi=140)
    plt.close(fig)
    print("done.")

if __name__ == "__main__":
    main()
