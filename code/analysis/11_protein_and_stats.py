#!/usr/bin/env python
"""Stage 11: (A) protein-layer validation of myeloid activation using the 35-plex
antibody-capture data with isotype-control background subtraction; (B) properly-
powered statistics for the disease-severity gradient — donor random-intercept
mixed models (composite program scores) + Jonckheere-Terpstra monotonic trend
(donor pseudobulk, robust to zero-inflation) for single genes.

Input:  data/processed/GSE288365_annotated_authmeta.h5ad  (stage 10)
Outputs: results/figures/protein_stats_validation.png
         results/tables/protein_{ALSvsControl,C9vsSALS}_VH.csv, protein_C9vsSALS_lesions.csv
         results/tables/mixedmodel_gradient_VH.csv, jonckheere_trend_VH.csv

Key caveats:
 - The 35-plex panel has NO GPNMB/CHIT1/TREM2/APOE, so protein cannot directly
   validate the RNA program switch; it validates general myeloid activation
   (CD68/ITGAX/CD163/ITGAM/CD14/HLA-DRA).
 - protein obsm lost column names on h5ad round-trip; names are in uns['protein_names'].
 - Gaussian LMM on sparse single genes -> singular (donor_var->0); use JT trend there.
 - C9-vs-sALS lesion contrast underpowered (2-3 donors/group); reported, not significant.
"""
import os
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMBA_NUM_THREADS']: os.environ.setdefault(k,'1')
os.environ.setdefault('KMP_AFFINITY','disabled'); os.environ.setdefault('OMP_PROC_BIND','FALSE')
import numpy as np, pandas as pd, anndata as ad, scanpy as sc
from scipy.stats import mannwhitneyu, norm
from itertools import combinations
import statsmodels.formula.api as smf
import warnings; warnings.filterwarnings('ignore')

REPO="/home/sunli/C9orf72-ALS"
A=ad.read_h5ad(f"{REPO}/data/processed/GSE288365_annotated_authmeta.h5ad")

# ---------- (A) protein layer ----------
pnames=list(A.uns['protein_names']); cols=[]; seen={}
for n in pnames:
    if n in seen: seen[n]+=1; cols.append(f"{n}.{seen[n]}")
    else: seen[n]=0; cols.append(n)
Praw=pd.DataFrame(np.asarray(A.obsm['protein'],float),columns=cols,index=A.obs_names)
iso=['mouse_IgG2a','mouse_IgG1k','mouse_IgG2bk','rat_IgG2a']
xl=np.log1p(Praw); Pclr=xl.sub(xl.mean(1),axis=0)          # CLR
Pbg=Pclr.sub(Pclr[iso].mean(1),axis=0)                      # isotype background subtraction
prot_markers=['CD68','CD163','ITGAX','ITGAM','CD14','HLA-DRA','PTPRC']  # FCGR3A drop (at bg)
obs=A.obs.copy()
for m in prot_markers: obs['prot_'+m]=Pbg[m].values
obs['condgrp']=obs['condition'].map({'control':'Control','C9-ALS':'ALS','sALS':'ALS'})
vh=obs['region']=='ventral_horn_MN'

def donor_contrast(mask,g1,g2,gcol='condition'):
    dd=obs[mask].groupby('donor_id',observed=True).agg({**{'prot_'+m:'mean' for m in prot_markers},gcol:'first'})
    r=[]
    for m in prot_markers:
        a=dd[dd[gcol]==g1]['prot_'+m].dropna(); b=dd[dd[gcol]==g2]['prot_'+m].dropna()
        if len(a)>=2 and len(b)>=2:
            r.append({'marker':m,'diff':a.mean()-b.mean(),'p':mannwhitneyu(a,b)[1],'n1':len(a),'n2':len(b)})
    return pd.DataFrame(r).sort_values('p')
donor_contrast(vh,'ALS','Control','condgrp').to_csv(f"{REPO}/results/tables/protein_ALSvsControl_VH.csv",index=False)
donor_contrast(vh,'C9-ALS','sALS').to_csv(f"{REPO}/results/tables/protein_C9vsSALS_VH.csv",index=False)
lesion=vh & obs['auth_annotation_general'].isin(['Ring','TDP43_adjacent'])
donor_contrast(lesion,'C9-ALS','sALS').to_csv(f"{REPO}/results/tables/protein_C9vsSALS_lesions.csv",index=False)

# ---------- (B) gradient statistics ----------
DAM=['GPNMB','CHIT1','TREM2','TYROBP','APOE','ITGAX','CST7','LPL','CD9','CD63','SPP1','CTSB','CTSD','LGALS3','FABP5','MSR1','CD68']
for nm,gl in [('DAM',DAM),('COMPLEMENT',['C1QA','C1QB','C1QC','C3']),('C9PROG',['BIRC3','SOCS2','CALHM6'])]:
    sc.tl.score_genes(A,[g for g in gl if g in A.var_names],score_name='score_'+nm,use_raw=False)

tgts=['score_DAM','score_COMPLEMENT','score_C9PROG','GPNMB','CHIT1','SOCS2','BIRC3','TREM2','APOE','C1QA','C1QC']
gi={g:A.var_names.get_loc(g) for g in tgts if g in A.var_names}
Xg=A[:,list(gi)].X; Xg=np.asarray(Xg.todense()) if hasattr(Xg,'todense') else np.asarray(Xg)
md=A.obs[['donor_id','condition','region']].copy().reset_index(drop=True)
for i,g in enumerate(gi): md[g]=Xg[:,i]
for s in ['score_DAM','score_COMPLEMENT','score_C9PROG']: md[s]=A.obs[s].values
vh_all=md[md['region']=='ventral_horn_MN'].copy()
vh_all['cond_ord']=vh_all['condition'].map({'control':0,'sALS':1,'C9-ALS':2}).astype(float)

# mixed model (donor random intercept), flag unreliable singular fits
rows=[]
for v in tgts:
    m=smf.mixedlm(f"{v} ~ cond_ord",vh_all,groups=vh_all['donor_id']).fit(method='lbfgs')
    dv=m.cov_re.iloc[0,0]; se=m.bse['cond_ord']
    rows.append({'target':v,'slope_per_step':m.params['cond_ord'],'se':se,'p_mixed':m.pvalues['cond_ord'],
                 'donor_var':dv,'reliable':(dv>1e-6) and np.isfinite(se) and se<10})
pd.DataFrame(rows).to_csv(f"{REPO}/results/tables/mixedmodel_gradient_VH.csv",index=False)

# Jonckheere-Terpstra on donor pseudobulk (robust for single genes)
def jonckheere(groups):
    J=sum((y>x)+0.5*(y==x) for i,j in combinations(range(len(groups)),2) for x in groups[i] for y in groups[j])
    ns=[len(g) for g in groups]; N=sum(ns)
    mu=(N**2-sum(n**2 for n in ns))/4
    var=(N**2*(2*N+3)-sum(n**2*(2*n+3) for n in ns))/72
    z=(J-mu)/np.sqrt(var); return z,2*(1-norm.cdf(abs(z)))
dpb=vh_all.groupby('donor_id',observed=True).agg({**{t:'mean' for t in tgts},'condition':'first'})
r=[]
for t in tgts:
    z,p=jonckheere([dpb[dpb['condition']==g][t].values for g in ['control','sALS','C9-ALS']])
    r.append({'target':t,'ctrl_mean':dpb[dpb['condition']=='control'][t].mean(),
              'sALS_mean':dpb[dpb['condition']=='sALS'][t].mean(),
              'C9_mean':dpb[dpb['condition']=='C9-ALS'][t].mean(),'JT_z':z,'p_trend':p})
pd.DataFrame(r).sort_values('p_trend').to_csv(f"{REPO}/results/tables/jonckheere_trend_VH.csv",index=False)
print("stage 11 done")
