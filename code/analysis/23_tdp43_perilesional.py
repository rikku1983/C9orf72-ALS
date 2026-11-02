#!/usr/bin/env python
"""
Stage 23 (spatial arm): perilesional structure around TDP-43 pathology.

MOTIVATION
Stages 17 and 22 reported TDP-43 coupling as an honest negative. Both tested a
MONOTONE relationship (grade ladder; proximal vs matched-distant), which cannot
detect a peak at an intermediate radius. This stage re-examines the question at
fine distance resolution inside the authors' annotated neighbourhood.

KEY FACT ESTABLISHED HERE (geometry, not assumption)
The authors' `auth_tdp43_a2_d2` levels are concentric dilation rings around
manually marked TDP-43+ spots, identical to the micron in every section:
    TDP43_positive   0 um    the spot itself (90 spots cohort-wide)
    TDP43_adjacent   112 um  ring 1
    Ring             200 um  ring 2      <- NOT a severity tier; a dilation ring
    TDP43_distant    330 um  ring 3
    Not_selected     ~2600 um median
`auth_annotation_general` merges positive+adjacent into one "TDP43_adjacent"
level (472 spots) = the paper's binary "TDP-43+" class.

RESULT
Within the anterior horn of the 4 lesion-richest sections, myeloid cell fraction
and the Myeloid/DAM program form a HALO: the TDP-43+ spot itself is neuron-
dominated (+6.4 pp neurons, -1.2 pp myeloid vs section mean) while myeloid
peaks at ring 3-4 (300-400 um, +2.0 pp) and troughs at rings 9-11 (-1.6 pp).
The gain is myeloid-SPECIFIC: at ring 3 neurons fall 2.0 pp and myeloid gains
2.0 pp while oligodendrocytes still LOSE 0.9 pp and astrocytes gain 0.6 pp.
Present in 3 of 4 sections; reversed in S4_GW35, whose lesions are 63% in the
POSTERIOR horn so its anterior-horn distance field is anchored on few centres.

CAVEATS
- Spot-level and pseudoreplicated: the evidence is section consistency, not p.
- Peak at 300-400 um is further than microglial process reach (~50-100 um), so
  this is a local reactive field, not direct perineuronal ensheathment.
- The rightmost distance bin (17+ rings) rises again for BOTH myeloid and
  oligodendrocytes: that is the anterior-horn periphery abutting white matter,
  an anatomical edge effect, not lesion-related. Do not read it as a 2nd peak.
- Ring labels are nearest-lesion, not unique-lesion; in dense lesion clusters
  the shells of neighbouring lesions merge.

NOT USED for the SfN talk (user decision, Nov 2026); kept for later follow-up.
"""
import os
for v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMBA_NUM_THREADS"):
    os.environ[v] = "1"
os.environ["KMP_AFFINITY"] = "disabled"; os.environ["OMP_PROC_BIND"] = "FALSE"
import numpy as np, pandas as pd, h5py, scipy.sparse as spx
from scipy import stats
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = "/home/sunli/C9orf72-ALS"
PROC = f"{ROOT}/data/processed"
TAB, FIG = f"{ROOT}/results/tables", f"{ROOT}/results/figures"
QC_H5 = f"{PROC}/GSE288365_qc.h5ad"                      # RAW integer UMI counts
AM_H5 = f"{PROC}/GSE288365_annotated_authmeta.h5ad"      # log-normalised .X + author obs
OBS_CSV = f"{TAB}/author_meta_merged_obs.csv"            # author obs, NaN-correct

# ---------- gene sets: bulk118_ward_k3, spatial-detectable members ----------
MYE = ["ADORA3","APOBR","APOC1","ARRB2","C1QC","C3AR1","CAPG","CD300A","CD53","CD68","CEBPA",
       "CHIT1","CTSD","CTSS","CYBA","FXYD5","GPNMB","GRN","HCK","IFI30","IQGAP2","ITGB2",
       "LAPTM5","LGALS3","LGALS9","NBEAL2","NCF2","NCKAP1L","PIK3R5","PLA2G7","PLAU","PYCARD",
       "SIGLEC8","SLC37A2","SOAT1","SOD2","SPI1","STX3","SYK","TNFRSF14","TREM2","TYROBP",
       "WDFY4","ANO6","ANXA2","CD63","DSP","EMP3","PAPOLA","S100A11","SLC38A2","SLC5A3"]
AST = ["B3GAT2","BRINP1","CHCHD2","CKB","CNTN1","DAAM2","DCLK2","FAM181B","GABARAPL2","GJA1",
       "GMPR","GPR37L1","HEPACAM","IL17D","ITM2C","LIX1","NADK2","NDRG2","NFIX","P2RY1","PRKX",
       "PSD2","PYCR2","PYGM","SDC3","SLC24A4","SLC38A1","SMOX","TMEM47","TRIL"]
OLI = ["ADD1","APBB1","CARNS1","CIRBP","CMTM5","CORO2B","DNHD1","EDIL3","GPR37","IGSF8","KLK6",
       "LDB3","LLGL1","LSS","MAPK8IP2","MYRF","NFASC","OPALIN","PAQR4","PLLP","PPP1R16B","SLC24A2"]
NEED = sorted(set(MYE + AST + OLI + ["GPNMB", "CHIT1"]))
FEATS = ["Myeloid/DAM", "Astrocyte", "Oligo/myelin", "GPNMB", "CHIT1"]
R4 = ["Anterior_horns", "Grey_matter", "White_matter", "Posterior_white_matter"]
TDPORD = ["TDP43_positive", "TDP43_adjacent", "Ring", "TDP43_distant", "Not_selected"]

# ---------- author obs (CSV holds NaN correctly; the h5ad categorical uses code -1,
#            which `categories[codes]` silently decodes as the LAST category) ----------
O = pd.read_csv(OBS_CSV, low_memory=False)

# ---------- expression: raw counts -> CP10K -> log1p -> z across ALL spots ----------
with h5py.File(QC_H5, "r") as h:
    gv = h["var"][h["var"].attrs["_index"]][:].astype(str)
    bq = h["obs"][h["obs"].attrs["_index"]][:].astype(str)
    idx = {g: int(np.where(gv == g)[0][0]) for g in NEED if g in set(gv)}
    X = h["X"]; indptr = X["indptr"][:]; n = len(indptr) - 1
    keys = list(idx); M = np.zeros((n, len(keys)), np.float32); tot = np.zeros(n)
    for s0 in range(0, n, 20000):
        e = min(s0 + 20000, n); a0, a1 = indptr[s0], indptr[e]
        m = spx.csr_matrix((X["data"][a0:a1], X["indices"][a0:a1], indptr[s0:e+1] - a0),
                           shape=(e - s0, len(gv)))
        tot[s0:e] = np.asarray(m.sum(1)).ravel()
        M[s0:e] = m[:, [idx[k] for k in keys]].toarray()
with h5py.File(AM_H5, "r") as h:
    assert np.array_equal(h["obs"][h["obs"].attrs["_index"]][:].astype(str), bq), "spot order differs"
    sp_xy = h["obsm"]["spatial"][:]
    def dec(k):
        g = h["obs"][k]
        return np.asarray(g["categories"][:].astype(str))[g["codes"][:]] if isinstance(g, h5py.Group) else g[:]
    arr, acc = dec("array_row").astype(int), dec("array_col").astype(int)
assert len(O) == n, "obs csv / h5ad length mismatch"

Z = pd.DataFrame((lambda L: (L - L.mean(0)) / L.std(0))(np.log1p(M / tot[:, None] * 1e4)), columns=keys)
A = pd.DataFrame({"section_key": O.section_key.values, "donor": O.donor_id.values,
                  "cond": O.condition.values, "anat": O.auth_manual.values,
                  "tdp": O.auth_tdp43_a2_d2.values, "ar": arr, "ac": acc,
                  "px": sp_xy[:, 0], "py": sp_xy[:, 1]})
for nm, gl in [("Myeloid/DAM", MYE), ("Astrocyte", AST), ("Oligo/myelin", OLI)]:
    A[nm] = Z[[g for g in gl if g in keys]].mean(1).values
A["GPNMB"] = Z["GPNMB"].values; A["CHIT1"] = Z["CHIT1"].values
A["neuron_frac"] = O["auth_norm_ct_Neurons"].values
A["myeloid_frac"] = (O["auth_norm_ct_Microglia"].values
                     + O["auth_norm_ct_ProliferatingMicroglia"].values
                     + O["auth_norm_ct_Macrophages"].values)
for ct in ["Oligodendrocytes", "Astrocytes"]:
    A[ct] = O["auth_norm_ct_" + ct].values

# ---------- verify the ring geometry (reproduces the dilation radii) ----------
geo = []
for sk, s in A[A.tdp.isin(TDPORD[:4])].groupby("section_key", observed=True):
    pos = A[(A.section_key == sk) & (A.tdp == "TDP43_positive")]
    if len(pos) == 0: continue
    P0 = np.c_[pos.ac * 50.0, pos.ar * 100.0]
    for lev in TDPORD[:4]:
        q = s[s.tdp == lev]
        if not len(q): continue
        Q = np.c_[q.ac * 50.0, q.ar * 100.0]
        d = np.sqrt(((Q[:, None, :] - P0[None, :, :]) ** 2).sum(-1)).min(1)
        geo.append(dict(section=sk, level=lev, n=len(q), median_um=np.median(d)))
pd.DataFrame(geo).to_csv(f"{TAB}/tdp43_grade_geometry.csv", index=False)

# ---------- the 4 lesion-richest sections, anterior horn only ----------
top4 = (A[A.tdp.isin(TDPORD[:4])].groupby("section_key", observed=True).size()
        .sort_values(ascending=False).head(4).index.tolist())
AH = A[(A.anat == "Anterior_horns") & A.section_key.isin(top4)].copy()
for sk, s in AH.groupby("section_key", observed=True):
    pos = A[(A.section_key == sk) & (A.tdp == "TDP43_positive")]
    P0 = np.c_[pos.ac * 50.0, pos.ar * 100.0]
    Q = np.c_[s.ac * 50.0, s.ar * 100.0]
    AH.loc[s.index, "dist"] = np.sqrt(((Q[:, None, :] - P0[None, :, :]) ** 2).sum(-1)).min(1)

EDG = [0, 60, 140, 240, 340, 450, 600, 800, 1100, 1600, 1e9]
LBL = ["0", "1", "2", "3", "4", "5-6", "7-8", "9-11", "12-16", "17+"]
AH["bin"] = pd.cut(AH.dist, EDG, labels=LBL, right=False)
CEN = ["Myeloid/DAM", "GPNMB", "myeloid_frac", "neuron_frac", "Oligodendrocytes", "Astrocytes"]
for c in CEN:                      # centre within section so section level cannot drive the profile
    AH[c + "_c"] = AH[c] - AH.groupby("section_key", observed=True)[c].transform("mean")

g = AH.groupby("bin", observed=True)
P = pd.DataFrame({"n": g.size()})
for c, lab in [("neuron_frac_c", "Neurons"), ("myeloid_frac_c", "Myeloid"),
               ("Oligodendrocytes_c", "Oligodendrocytes"), ("Astrocytes_c", "Astrocytes"),
               ("Myeloid/DAM_c", "Myeloid/DAM score"), ("GPNMB_c", "GPNMB")]:
    P[lab] = g[c].mean(); P[lab + "_sem"] = g[c].sem()
P = P.reindex(LBL).reset_index().rename(columns={"index": "bin"})
P.to_csv(f"{TAB}/tdp43_perilesional_profile.csv", index=False)

# per-section halo test: rings 2-4 (200-450 um) vs rings 7-11 (800-1600 um)
rows = []
for sk, s in AH.groupby("section_key", observed=True):
    near, mid = s[s.bin.isin(["2", "3", "4"])], s[s.bin.isin(["7-8", "9-11"])]
    rows.append(dict(section=sk, cond=s.cond.iloc[0], n_near=len(near), n_mid=len(mid),
                     myeloid_near=near.myeloid_frac.mean(), myeloid_mid=mid.myeloid_frac.mean(),
                     delta_pp=near.myeloid_frac.mean() - mid.myeloid_frac.mean(),
                     p_spotlevel=stats.mannwhitneyu(near.myeloid_frac, mid.myeloid_frac)[1]))
pd.DataFrame(rows).to_csv(f"{TAB}/tdp43_halo_persection.csv", index=False)

# categorical version in the authors' own levels, rings 0-3 only, all sections with >=10 such spots
AT = A[(A.anat == "Anterior_horns") & A.tdp.isin(TDPORD)].copy()
AT["code"] = AT.tdp.map({v: i for i, v in enumerate(TDPORD)})
R03 = AT[AT.tdp != "Not_selected"]
keep = R03.groupby("section_key", observed=True).size()
rows = []
for sk in keep[keep >= 40].index:
    s = R03[R03.section_key == sk]
    r = dict(section=sk, cond=s.cond.iloc[0], n=len(s),
             rho_neuron=stats.spearmanr(s.neuron_frac, s.code)[0],
             rho_myeloid_frac=stats.spearmanr(s.myeloid_frac, s.code)[0])
    for f in FEATS:
        r["rho_" + f.replace("/", "_")] = (stats.spearmanr(s[f], s.code)[0]
                                           if s[f].nunique() > 1 else np.nan)
    rows.append(r)
pd.DataFrame(rows).to_csv(f"{TAB}/tdp43_rings03_gradient.csv", index=False)

# ---------- figure ----------
CT4 = {"Neurons": "#332288", "Myeloid": "#CC6677",
       "Oligodendrocytes": "#88CCEE", "Astrocytes": "#44AA99"}
plt.rcParams.update({"font.size": 7, "axes.labelsize": 7, "xtick.labelsize": 6,
                     "ytick.labelsize": 6, "figure.dpi": 150})
fig, axes = plt.subplots(2, 1, figsize=(5.4, 5.2), sharex=True)
x = np.arange(len(P))
for lab, col in CT4.items():
    axes[0].errorbar(x, P[lab], yerr=P[lab + "_sem"], color=col, marker="o", ms=3.4, lw=1.2,
                     capsize=1.6, elinewidth=0.6, label=lab)
for ax in axes:
    ax.axhline(0, color="0.75", lw=0.6, zorder=0)
    ax.axvspan(1.5, 4.5, color="#CC6677", alpha=0.07, zorder=0)
    for sp_ in ("top", "right"): ax.spines[sp_].set_visible(False)
axes[0].set_ylabel("cell fraction, section-centred (pp)", fontsize=7.4)
axes[0].legend(frameon=False, fontsize=6.4, ncol=2, loc="upper right",
               handletextpad=0.3, columnspacing=1.0)
for lab, col, ls in [("Myeloid/DAM score", "#CC6677", "-"), ("GPNMB", "#882255", "--")]:
    axes[1].errorbar(x, P[lab], yerr=P[lab + "_sem"], color=col, marker="o", ms=3.4, lw=1.2,
                     capsize=1.6, elinewidth=0.6, ls=ls, label=lab)
axes[1].set_ylabel("expression, section-centred (z)", fontsize=7.4)
axes[1].legend(frameon=False, fontsize=6.4, loc="upper left", handletextpad=0.4)
axes[1].set_xticks(x); axes[1].set_xticklabels(P["bin"])
axes[1].set_xlabel("distance from nearest TDP-43+ spot (spot steps, 100 um each)")
fig.subplots_adjust(left=0.145, right=0.985, top=0.985, bottom=0.11, hspace=0.10)
fig.savefig(f"{FIG}/tdp43_perilesional_profile.png", dpi=300)

print("sections used:", [s[:34] for s in top4])
print(P[["bin", "n", "Neurons", "Myeloid", "Oligodendrocytes", "Astrocytes"]].round(2).to_string(index=False))
print("\nwrote tdp43_{grade_geometry,perilesional_profile,halo_persection,rings03_gradient}.csv"
      " and tdp43_perilesional_profile.png")
