#!/usr/bin/env python
"""
05_pseudobulk_de.py — Pseudobulk differential expression across conditions.

Why pseudobulk: 132k spots come from only 18 donors. Spot-level tests treat
spots as independent replicates (pseudoreplication) and grossly inflate
significance. We aggregate RAW counts to one profile per donor (and per
donor x region), then run DESeq2 across donors -> donor-level statistics.

Input : data/processed/GSE288365_annotated.h5ad  (needs .layers or raw counts)
        data/processed/GSE288365_qc.h5ad          (raw counts fallback)
Output: results/tables/pseudobulk_<contrast>_<scope>.csv
        data/processed/pseudobulk_counts.csv   (donor x gene matrix + meta)

Contrasts (each vs the reference):
- C9-ALS  vs control
- sALS    vs control
- C9-ALS  vs sALS      <-- C9-specific signal (primary MOA question)
Scopes: 'global' (all spots per donor) and 'ventral_horn_MN' (MN region only).
"""
import os
os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("KMP_AFFINITY", "disabled")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
import warnings; warnings.filterwarnings("ignore")
import numpy as np
import pandas as pd
import scanpy as sc
from scipy.sparse import issparse

REPO = "/home/sunli/C9orf72-ALS"
PROC = f"{REPO}/data/processed"
TAB  = f"{REPO}/results/tables"

MIN_SPOTS_PER_PSEUDOBULK = 50   # a donor x region group needs >=50 spots
MIN_DONORS_PER_GROUP     = 3    # need >=3 donors per condition to test

def get_raw_counts():
    """Load raw counts + region labels aligned by spot."""
    ann = sc.read_h5ad(f"{PROC}/GSE288365_annotated.h5ad")
    region = ann.obs[["donor_id", "condition", "region", "section"]].copy()
    del ann
    qc = sc.read_h5ad(f"{PROC}/GSE288365_qc.h5ad")  # raw counts in .X
    # align (qc has all QC spots; annotated is subset/same order)
    qc = qc[region.index].copy()
    X = qc.X
    genes = qc.var_names.to_numpy()
    return X, genes, region

def build_pseudobulk(X, genes, meta, by_region=False):
    """Sum raw counts within donor (optionally donor x region)."""
    rows, index = [], []
    if by_region:
        groups = meta.groupby(["donor_id", "region"], observed=True)
    else:
        groups = meta.groupby("donor_id", observed=True)
    for key, idx in groups.groups.items():
        pos = meta.index.get_indexer(idx)
        if len(pos) < MIN_SPOTS_PER_PSEUDOBULK:
            continue
        sub = X[pos]
        s = np.asarray(sub.sum(axis=0)).ravel() if issparse(sub) else sub.sum(0)
        rows.append(s)
        donor = key[0] if by_region else key
        region = key[1] if by_region else "global"
        cond = meta.loc[idx, "condition"].iloc[0]
        index.append((donor, region, cond, len(pos)))
    mat = pd.DataFrame(np.vstack(rows), columns=genes)
    info = pd.DataFrame(index, columns=["donor_id", "region", "condition", "n_spots"])
    return mat.reset_index(drop=True), info

def run_deseq(counts, info, cond_a, cond_b):
    """DESeq2: cond_a vs cond_b (positive LFC = up in cond_a)."""
    from pydeseq2.dds import DeseqDataSet
    from pydeseq2.ds import DeseqStats
    keep = info["condition"].isin([cond_a, cond_b]).values
    if keep.sum() == 0:
        return None
    c = counts.loc[keep].copy()
    m = info.loc[keep].copy().reset_index(drop=True)
    # need enough donors on each side
    vc = m["condition"].value_counts()
    if vc.get(cond_a, 0) < MIN_DONORS_PER_GROUP or vc.get(cond_b, 0) < MIN_DONORS_PER_GROUP:
        return None
    # filter low-count genes (>=10 counts total)
    c = c.loc[:, c.sum(0) >= 10]
    c.index = [f"s{i}" for i in range(len(c))]
    m.index = c.index
    m["condition"] = pd.Categorical(m["condition"], categories=[cond_b, cond_a])
    dds = DeseqDataSet(counts=c.astype(int), metadata=m, design_factors="condition",
                       ref_level=["condition", cond_b], quiet=True)
    dds.deseq2()
    st = DeseqStats(dds, contrast=["condition", cond_a, cond_b], quiet=True)
    st.summary()
    res = st.results_df.sort_values("padj")
    return res

CONTRASTS = [("C9-ALS", "control"), ("sALS", "control"), ("C9-ALS", "sALS")]

def main():
    X, genes, meta = get_raw_counts()
    print(f"[load] {X.shape[0]} spots, {len(genes)} genes, {meta['donor_id'].nunique()} donors")

    # GLOBAL pseudobulk (one profile per donor)
    gmat, ginfo = build_pseudobulk(X, genes, meta, by_region=False)
    gmat.to_csv(f"{PROC}/pseudobulk_counts.csv")
    ginfo.to_csv(f"{PROC}/pseudobulk_info.csv", index=False)
    print(f"[global] {len(ginfo)} donor pseudobulks")
    print(ginfo.groupby("condition", observed=True).size())

    # VENTRAL HORN pseudobulk (donor x ventral_horn_MN)
    rmat, rinfo = build_pseudobulk(X, genes, meta, by_region=True)
    vh = rinfo["region"] == "ventral_horn_MN"
    print(f"[ventral_horn] {vh.sum()} donor pseudobulks")

    for scope, (mat, info, mask) in {
        "global": (gmat, ginfo, np.ones(len(ginfo), bool)),
        "ventral_horn_MN": (rmat, rinfo, vh.values),
    }.items():
        m_mat = mat.loc[mask].reset_index(drop=True)
        m_info = info.loc[mask].reset_index(drop=True)
        for a, b in CONTRASTS:
            res = run_deseq(m_mat, m_info, a, b)
            if res is None:
                print(f"  [skip] {a} vs {b} ({scope}): too few donors")
                continue
            tag = f"{a}_vs_{b}".replace("-", "").replace(" ", "")
            fn = f"{TAB}/pseudobulk_{tag}_{scope}.csv"
            res.to_csv(fn)
            n_sig = (res["padj"] < 0.05).sum()
            print(f"  [{scope}] {a} vs {b}: {n_sig} genes padj<0.05  -> {os.path.basename(fn)}")
    print("done.")

if __name__ == "__main__":
    main()
