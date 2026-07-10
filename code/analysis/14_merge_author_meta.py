#!/usr/bin/env python
"""Merge author annotations onto our QC'd cells -> obs-level parquet (memory-safe).
The expression matrix stays in GSE268995_raw_qc.h5ad; this table joins to it by
cell_id (== our obs_names, SAMPLE::BARCODE using our matrix-sample labels).

Resolves the column-1 barcode-shift in the author RDS deposit:
  author-file A1->our-matrix C1, B1->D1, C1->E1, D1->F1, E1->A1, F1->B1 (two 3-cycles).
GEO series-matrix condition is AUTHORITATIVE. Author donor labels conflict for
5 samples (A1,B1,C1,D1,E1) -> flagged, not trusted. All 4 C9-ALS samples concordant.
"""
import os
for v in ("OMP_NUM_THREADS","OPENBLAS_NUM_THREADS","MKL_NUM_THREADS","NUMBA_NUM_THREADS"):
    os.environ[v]="1"
os.environ["KMP_AFFINITY"]="disabled"; os.environ["OMP_PROC_BIND"]="FALSE"
import numpy as np, pandas as pd, h5py

ROOT="/home/sunli/C9orf72-ALS"
META=f"{ROOT}/results/tables/GSE268995_author_meta.csv.gz"
H5=f"{ROOT}/data/processed/GSE268995_raw_qc.h5ad"
OUT=f"{ROOT}/results/tables/GSE268995_annotated_obs.parquet"

# --- our cells: read obs_names + our per-cell QC from h5ad without loading X ---
with h5py.File(H5,"r") as f:
    # obs index
    obs_grp=f["obs"]
    idx_key=obs_grp.attrs.get("_index","_index")
    our_names=obs_grp[idx_key][:].astype(str)
    def col(name):
        d=obs_grp[name]
        if isinstance(d,h5py.Group):  # categorical
            cats=d["categories"][:].astype(str); codes=d["codes"][:]
            return pd.Categorical.from_codes(codes,cats)
        return d[:]
    our=pd.DataFrame({"cell_id":our_names})
    for c in ["sample","gsm","condition","disease","n_genes_by_counts","total_counts","pct_counts_mt","pct_counts_ribo"]:
        try: our[c]=col(c)
        except Exception as e: print("skip",c,e)
print(f"our cells: {len(our)}",flush=True)

# --- author meta + barcode-shift remap ---
m=pd.read_csv(META,low_memory=False)
m["author_sample"]=m["cell_id"].str.split("::").str[0]
m["barcode"]=m["cell_id"].str.split("::").str[1]
SHIFT={"A1":"C1","B1":"D1","C1":"E1","D1":"F1","E1":"A1","F1":"B1"}
m["our_sample"]=m["author_sample"].map(lambda s: SHIFT.get(s,s))
m["phys_cell"]=m["our_sample"]+"::"+m["barcode"]
assert m["phys_cell"].is_unique
m=m.set_index("phys_cell")
print(f"author cells: {len(m)}",flush=True)

ann_cols=['predicted.celltype.l1','predicted.celltype.l1.score','predicted.celltype.l2',
          'predicted.celltype.l2.score','seurat_clusters','harmony_clusters','DF.classifications',
          'diagnosis','diagnosis_general','age','sex','race',
          'tcr_clonotype_id','tcr_frequency','mait_evidence','inkt_evidence',
          'bcr_clonotype_id','bcr_frequency','umap1','umap2']

# --- join author annotations onto our cells (inner = author-retained) ---
maligned=m.reindex(our["cell_id"].values)
for c in ann_cols:
    our[c]=maligned[c].values
diag_map={"als_slow":"sALS_slow","als_fast":"sALS_fast","healthy_control":"control","als_c9orf72":"C9-ALS"}
our["author_condition"]=maligned["diagnosis"].map(diag_map).values
our["author_retained"]=our["predicted.celltype.l1"].notna()
conflict_samples={"A1","B1","C1","D1","E1"}
our["sample_label_conflict"]=our["sample"].isin(conflict_samples)

our.to_parquet(OUT,index=False)
print("WROTE",OUT,flush=True)
print("author_retained:",int(our["author_retained"].sum()),"/",len(our),flush=True)
sub=our[our["author_retained"]]
print("\nL1 (author-retained):\n",sub['predicted.celltype.l1'].value_counts().to_string(),flush=True)
print("\ncondition (GEO authoritative, author-retained):\n",sub['condition'].value_counts().to_string(),flush=True)
print("\nconflict cells (author-retained):",int(sub['sample_label_conflict'].sum()),flush=True)
