#!/usr/bin/env python
"""
02_integrate_cluster.py — Normalize, batch-integrate (Harmony), and cluster
the merged GSE288365 Visium spots.

Input : data/processed/GSE288365_qc.h5ad
Output: data/processed/GSE288365_clustered.h5ad
        results/figures/umap_condition.png, umap_leiden.png
        results/tables/cluster_composition.csv

Design notes:
- log1p-normalize (not SCT): low RNA depth + protein-focused assay.
- HVGs on RNA only (2000, seurat_v3 on raw counts).
- Harmony batch key = section (20 sections); protects against per-slide
  technical variation while preserving biology across donors.
- Leiden at resolution 1.0 for a first pass (regions annotated in stage 03).
"""
import os, gc
os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("KMP_AFFINITY", "disabled")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
import warnings; warnings.filterwarnings("ignore")
import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

REPO = "/home/sunli/C9orf72-ALS"
PROC = f"{REPO}/data/processed"
FIG  = f"{REPO}/results/figures"
TAB  = f"{REPO}/results/tables"
sc.settings.n_jobs = 4
RES = 1.0
N_HVG = 2000
N_PCS = 30

def step1_lognorm():
    """Load QC counts, compute HVGs, normalize+log1p, persist. Idempotent."""
    if os.path.exists(f"{PROC}/GSE288365_lognorm.h5ad"):
        print("[step1] lognorm exists, skip")
        return
    a = sc.read_h5ad(f"{PROC}/GSE288365_qc.h5ad")
    print(f"[step1] load {a.shape}")
    sc.pp.highly_variable_genes(a, flavor="seurat_v3", n_top_genes=N_HVG,
                                batch_key="section")
    sc.pp.normalize_total(a, target_sum=1e4)
    sc.pp.log1p(a)
    a.write(f"{PROC}/GSE288365_lognorm.h5ad")
    print("[step1] wrote lognorm")

def step2_cluster():
    """Cluster on HVG subset only; save tiny embeddings (no expression matrix).

    Memory-lean (7 GB box): load full sparse lognorm (~1.2 GB), slice HVGs,
    free the full object, then run PCA DIRECTLY on the sparse log-norm matrix.
    We deliberately skip sc.pp.scale — zero-centering there densifies the
    132k x 2000 matrix (~2 GB) and PCA/Harmony on top exceeds RAM. PCA with
    zero_center=True centers implicitly on sparse input without densifying.
    """
    a = sc.read_h5ad(f"{PROC}/GSE288365_lognorm.h5ad")
    hvg = a.var_names[a.var["highly_variable"]].tolist()
    ah = a[:, hvg].copy()
    del a; gc.collect()
    # ensure float32 sparse
    ah.X = ah.X.astype(np.float32)
    print(f"[step2] HVG subset {ah.shape}")
    sc.tl.pca(ah, n_comps=N_PCS, svd_solver="arpack", zero_center=True)
    # Call harmonypy directly: this version returns Z_corr as (cells x PCs);
    # scanpy's harmony_integrate wrapper assumes (PCs x cells) and transposes,
    # corrupting the shape. Assign shape-robustly.
    import harmonypy
    ho = harmonypy.run_harmony(ah.obsm["X_pca"], ah.obs, ["section"],
                               max_iter_harmony=20)
    Z = np.asarray(ho.Z_corr)
    if Z.shape[0] != ah.n_obs:      # got (PCs x cells) -> transpose
        Z = Z.T
    assert Z.shape == (ah.n_obs, N_PCS), f"harmony shape {Z.shape}"
    ah.obsm["X_pca_harmony"] = Z.astype(np.float32)
    del ho, Z; gc.collect()
    sc.pp.neighbors(ah, use_rep="X_pca_harmony", n_neighbors=15, n_pcs=N_PCS)
    sc.tl.leiden(ah, resolution=RES, key_added="leiden", flavor="igraph",
                 n_iterations=2, directed=False)
    sc.tl.umap(ah, min_dist=0.3)
    print(f"[step2] {ah.obs['leiden'].nunique()} leiden clusters")
    # persist ONLY embeddings + labels (tiny)
    np.savez(f"{PROC}/_embeddings.npz",
             obs_names=ah.obs_names.values.astype(str),
             leiden=ah.obs["leiden"].values.astype(str),
             X_pca=ah.obsm["X_pca"], X_pca_harmony=ah.obsm["X_pca_harmony"],
             X_umap=ah.obsm["X_umap"])
    del ah; gc.collect()

def step3_merge_and_plot():
    """Attach labels+embeddings to full lognorm, save clustered, make figures."""
    a = sc.read_h5ad(f"{PROC}/GSE288365_lognorm.h5ad")
    emb = np.load(f"{PROC}/_embeddings.npz", allow_pickle=True)
    assert list(emb["obs_names"]) == list(a.obs_names), "obs order mismatch"
    a.obs["leiden"] = pd.Categorical(emb["leiden"])
    a.obsm["X_pca"] = emb["X_pca"]
    a.obsm["X_pca_harmony"] = emb["X_pca_harmony"]
    a.obsm["X_umap"] = emb["X_umap"]
    a.write(f"{PROC}/GSE288365_clustered.h5ad")
    comp = (pd.crosstab(a.obs["leiden"], a.obs["condition"])
            .reindex(columns=["control", "sALS", "C9-ALS"]))
    comp.to_csv(f"{TAB}/cluster_composition.csv")
    for color, fn in [("condition", "umap_condition.png"),
                      ("leiden", "umap_leiden.png")]:
        fig, ax = plt.subplots(figsize=(7, 6))
        sc.pl.umap(a, color=color, ax=ax, show=False, size=3,
                   legend_loc="on data" if color == "leiden" else "right margin")
        fig.tight_layout(); fig.savefig(f"{FIG}/{fn}", dpi=150)
        plt.close(fig)
    print("[step3] done.")

STEPS = {"1": step1_lognorm, "2": step2_cluster, "3": step3_merge_and_plot}

if __name__ == "__main__":
    import sys
    sel = sys.argv[1] if len(sys.argv) > 1 else "all"
    if sel == "all":
        step1_lognorm(); step2_cluster(); step3_merge_and_plot()
    else:
        STEPS[sel]()
