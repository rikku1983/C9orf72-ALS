#!/usr/bin/env python
"""Stage 10: integrate author RDS metadata (deconvolution + TDP-43 pathology +
MN-distance + imaging) into our AnnData, then compositional & pathology analysis.

Prereqs:
  - data/processed/GSE288365_annotated.h5ad  (our object; stages 01-03)
  - data/processed/seurat_extract/ALL_*_meta.csv  (author metadata; via
    code/download/download_extract_all_rds.sh + 09b_extract_meta_only.R)

Outputs: results/figures/author_meta_integration.png + composition/tdp43 tables,
and data/processed/GSE288365_annotated_authmeta.h5ad (enriched obs).
"""
import os
os.environ.setdefault('OMP_NUM_THREADS','1'); os.environ.setdefault('KMP_AFFINITY','disabled')
os.environ.setdefault('OPENBLAS_NUM_THREADS','1'); os.environ.setdefault('MKL_NUM_THREADS','1')
import re, glob, numpy as np, pandas as pd, anndata as ad
from scipy.stats import mannwhitneyu

REPO="/home/sunli/C9orf72-ALS"
base=f"{REPO}/data/processed/seurat_extract"
A=ad.read_h5ad(f"{REPO}/data/processed/GSE288365_annotated.h5ad")

# ---- load & stack author metadata; key = capture-area '||' cell barcode ----
def parse_bc(bc): left,cell=bc.rsplit("__",1); return left,cell
auth=[]
for f in sorted(glob.glob(f"{base}/ALL_*_meta.csv")):
    df=pd.read_csv(f)
    df['left']=[parse_bc(b)[0] for b in df['barcode']]
    df['cell']=[parse_bc(b)[1] for b in df['barcode']]
    auth.append(df)
auth=pd.concat(auth,ignore_index=True)
auth['key']=auth['left']+"||"+auth['cell']
auth_idx=auth.set_index('key')

# ---- map our section_key -> author capture area (normalized string equality) ----
def norm(s): return re.sub(r'[^A-Za-z0-9]+','-',s).upper().strip('-')
lefts=sorted(auth['left'].unique())
mapping={s:[l for l in lefts if norm(l)==norm(s)][0] for s in A.obs['section_key'].unique()}

CT_RAW=['Astrocytes','Endothelial','Ependymal.Cells','Lymphocytes','Macrophages','Meninges',
        'Microglia','Neurons','OPC','Oligodendrocytes','Pericytes','ProliferatingMicroglia','Schwann','Total_Cells']
CT_NORM=[f'norm_ct_{c}' for c in CT_RAW if c!='Total_Cells']
IMG=['TDP43_mean','TDP43_stdev','tdp43_a2_d2','MAP2_norm','MAP2_raw_mean','MAP2_90_pctile','cleaned_IBA1_binary','IBA1_positive']
SPATIAL=['MN_dist_map','MN_adjacent','MN_dist_group','annotation_general','manual','manual_3']
WANT=[c for c in CT_RAW+CT_NORM+IMG+SPATIAL if c in auth.columns]

our_key=A.obs['section_key'].map(mapping).astype(str).values+"||"+np.array([b.split("::")[-1] for b in A.obs_names])
aligned=auth_idx.reindex(our_key)[WANT]; aligned.index=A.obs_names
for c in aligned.columns: A.obs['auth_'+c]=aligned[c].values
print("per-spot match rate:", round(A.obs['auth_Total_Cells'].notna().mean(),4))
A.write(f"{REPO}/data/processed/GSE288365_annotated_authmeta.h5ad")

# ---- compositional analysis (donor-level, ventral horn) ----
ctn=['auth_'+c for c in CT_NORM]
matched=A.obs[A.obs['auth_Total_Cells'].notna()]
vh=matched[matched['region']=='ventral_horn_MN']
donor_vh=vh.groupby('donor_id',observed=True).agg({**{c:'mean' for c in ctn},'condition':'first'})

def contrast(dd,g1,g2):
    out=[]
    for c in ctn:
        a=dd[dd['condition']==g1][c].dropna(); b=dd[dd['condition']==g2][c].dropna()
        if len(a)>=2 and len(b)>=2:
            out.append({'celltype':c.replace('auth_norm_ct_',''),'diff':a.mean()-b.mean(),
                        'p':mannwhitneyu(a,b)[1]})
    return pd.DataFrame(out).sort_values('p')
dv=donor_vh.copy(); dv['condition']=dv['condition'].map({'control':'Control','C9-ALS':'ALS','sALS':'ALS'})
contrast(dv,'ALS','Control').to_csv(f"{REPO}/results/tables/composition_ALSvsControl_VH.csv",index=False)
contrast(donor_vh,'C9-ALS','sALS').to_csv(f"{REPO}/results/tables/composition_C9vsSALS_VH.csv",index=False)

# ---- TDP-43 pathology gradient (all ALS) ----
order=['Not_selected','TDP43_distant','TDP43_adjacent','Ring']
GENES=[g for g in ['GPNMB','CHIT1','C1QA','C1QB','C1QC','CD68','TYROBP','APOE','MSR1','AIF1',
       'TREM2','BIRC3','SOCS2','CALHM6','GFAP','CHAT','NEFH','NEFL','P2RY12'] if g in A.var_names]
sub=A[A.obs['auth_annotation_general'].isin(order)]
Xg=sub[:,GENES].X; Xg=np.asarray(Xg.todense()) if hasattr(Xg,'todense') else np.asarray(Xg)
E=pd.DataFrame(Xg,columns=GENES); E['anno']=pd.Categorical(sub.obs['auth_annotation_general'].values,order,ordered=True)
E.groupby('anno',observed=True)[GENES].mean().T.to_csv(f"{REPO}/results/tables/tdp43_gradient_expression.csv")
print("stage 10 done")
