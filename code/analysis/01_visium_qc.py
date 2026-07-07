#!/usr/bin/env python
"""
01_visium_qc.py — Load all GSE288365 Visium spinal-cord sections, attach
condition metadata, run per-spot QC, and merge into one AnnData.

Input : data/raw/GSE288365/  (flat 10x files, one set per section)
        data/metadata/GSE288365_sample_condition_map.csv
Output: data/processed/GSE288365_qc.h5ad   (merged, QC-filtered, raw counts in .X)
        results/figures/qc_*.png
        results/tables/qc_per_section.csv

Run:  python code/analysis/01_visium_qc.py
"""
import os
# sandbox: disable OpenMP thread-affinity (pthread_setaffinity_np not permitted)
os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("KMP_AFFINITY", "disabled")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
import json, gzip, warnings
warnings.filterwarnings("ignore")
import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
from scipy.io import mmread
from scipy.sparse import csr_matrix
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

REPO = "/home/sunli/C9orf72-ALS"
RAW  = f"{REPO}/data/raw/GSE288365"
META = f"{REPO}/data/metadata/GSE288365_sample_condition_map.csv"
PROC = f"{REPO}/data/processed"
FIG  = f"{REPO}/results/figures"
TAB  = f"{REPO}/results/tables"
for d in (PROC, FIG, TAB):
    os.makedirs(d, exist_ok=True)

# ---- QC thresholds (RNA modality only; low-depth postmortem cord Visium) ----
# This assay is Visium GEX + 35-plex Antibody Capture (protein). Protein counts
# dominate (~99.7% of total), so ALL QC is computed on the RNA sub-matrix.
# RNA median ~850 counts / ~630 genes per spot -> lenient thresholds.
MIN_COUNTS = 200      # remove empty / off-tissue-edge / debris spots
MIN_GENES  = 100
MAX_PCT_MT = 30.0     # postmortem CNS runs high; 30% removes lysed/debris spots
MIN_CELLS_PER_GENE = 3

def load_section(prefix):
    """Build an AnnData for one section from flat 10x files.

    The matrix carries two modalities: 'Gene Expression' (RNA, ~18k genes) and
    'Antibody Capture' (35 protein markers). We keep RNA in .X and stash the
    35 protein markers in .obsm['protein'] (+ names in .uns['protein_names']).
    """
    base = f"{RAW}/GSE288365_{prefix}"
    # matrix: features x spots -> transpose to spots x features
    Xall = mmread(f"{base}_matrix.mtx.gz").T.tocsr()
    barcodes = pd.read_csv(f"{base}_barcodes.tsv.gz", header=None)[0].values
    feats = pd.read_csv(f"{base}_features.tsv.gz", header=None, sep="\t")
    feats.columns = ["ensembl", "symbol", "ftype"][:feats.shape[1]]
    is_rna  = (feats["ftype"] == "Gene Expression").values
    is_prot = (feats["ftype"] == "Antibody Capture").values
    # --- RNA modality -> .X ---
    A = ad.AnnData(X=csr_matrix(Xall[:, is_rna]))
    A.obs_names = barcodes
    A.var_names = feats.loc[is_rna, "symbol"].astype(str).values
    A.var_names_make_unique()
    A.var["ensembl"] = feats.loc[is_rna, "ensembl"].values
    # --- protein modality -> .obsm['protein'] ---
    prot = np.asarray(Xall[:, is_prot].todense())
    A.obsm["protein"] = prot
    A.uns["protein_names"] = list(feats.loc[is_prot, "symbol"].astype(str).values)
    # spatial positions
    pos = pd.read_csv(f"{base}_tissue_positions.csv.gz")
    pos = pos.set_index("barcode")
    pos = pos.reindex(A.obs_names)
    A.obs["in_tissue"]  = pos["in_tissue"].values
    A.obs["array_row"]  = pos["array_row"].values
    A.obs["array_col"]  = pos["array_col"].values
    A.obsm["spatial"]   = pos[["pxl_col_in_fullres", "pxl_row_in_fullres"]].values.astype(float)
    # keep only in-tissue spots
    A = A[A.obs["in_tissue"] == 1].copy()
    # scalefactors (store for later plotting)
    with gzip.open(f"{base}_scalefactors_json.json.gz", "rt") as fh:
        A.uns["scalefactors"] = json.load(fh)
    return A

def main():
    meta = pd.read_csv(META)
    print(f"[meta] {len(meta)} sections")
    per_section = []
    adatas = []
    for _, row in meta.iterrows():
        sec = row["geo_section"]
        A = load_section(sec)
        # QC metrics
        A.var["mt"] = A.var_names.str.startswith("MT-")
        sc.pp.calculate_qc_metrics(A, qc_vars=["mt"], inplace=True, percent_top=None)
        n_before = A.n_obs
        keep = (
            (A.obs["total_counts"] >= MIN_COUNTS) &
            (A.obs["n_genes_by_counts"] >= MIN_GENES) &
            (A.obs["pct_counts_mt"] <= MAX_PCT_MT)
        )
        A = A[keep].copy()
        # annotate condition
        for c in ["donor_id", "condition", "subdiag", "age", "sex"]:
            A.obs[c] = row[c]
        A.obs["section"] = sec
        A.obs_names = [f"{sec}::{b}" for b in A.obs_names]
        adatas.append(A)
        per_section.append(dict(
            section=sec, donor=row["donor_id"], condition=row["condition"],
            n_spots_raw=n_before, n_spots_qc=A.n_obs,
            median_counts=float(np.median(A.obs["total_counts"])),
            median_genes=float(np.median(A.obs["n_genes_by_counts"])),
            median_pct_mt=float(np.median(A.obs["pct_counts_mt"])),
        ))
        print(f"[{sec}] {row['condition']:8s} spots {n_before}->{A.n_obs}  "
              f"med_counts={per_section[-1]['median_counts']:.0f}")

    # merge (outer join on genes; keep raw counts)
    merged = ad.concat(adatas, join="outer", label="section_key",
                       keys=[a.obs["section"][0] for a in adatas], index_unique=None, fill_value=0)
    merged.uns["protein_names"] = adatas[0].uns["protein_names"]
    sc.pp.filter_genes(merged, min_cells=MIN_CELLS_PER_GENE)
    merged.obs["condition"] = pd.Categorical(
        merged.obs["condition"], categories=["control", "sALS", "C9-ALS"], ordered=True)
    print(f"[merged] {merged.n_obs} spots x {merged.n_vars} genes")
    print(merged.obs["condition"].value_counts())

    merged.write(f"{PROC}/GSE288365_qc.h5ad")
    qc = pd.DataFrame(per_section)
    qc.to_csv(f"{TAB}/qc_per_section.csv", index=False)

    # ---- QC figure ----
    fig, axes = plt.subplots(1, 3, figsize=(13, 4))
    order = ["control", "sALS", "C9-ALS"]
    colors = {"control": "#4C72B0", "sALS": "#DD8452", "C9-ALS": "#C44E52"}
    for ax, (col, lab) in zip(axes, [("median_counts", "Median counts/spot"),
                                     ("median_genes", "Median genes/spot"),
                                     ("n_spots_qc", "QC-passing spots")]):
        for i, cond in enumerate(order):
            sub = qc[qc["condition"] == cond]
            ax.scatter(np.full(len(sub), i) + np.random.uniform(-.1, .1, len(sub)),
                       sub[col], color=colors[cond], s=40, alpha=.8, edgecolor="k", lw=.4)
        ax.set_xticks(range(3)); ax.set_xticklabels(order); ax.set_title(lab)
    fig.suptitle("GSE288365 Visium — per-section QC")
    fig.tight_layout()
    fig.savefig(f"{FIG}/qc_per_section.png", dpi=150)
    print("done.")

if __name__ == "__main__":
    main()
