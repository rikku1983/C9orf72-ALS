# C9orf72-ALS Project — State / Resume Guide

**Last updated:** end of session 11 (Phase III: C9-versus-sALS resolution, stage 24).
**Repo:** github.com/rikku1983/C9orf72-ALS (branch `main`). Local: `/home/sunli/C9orf72-ALS`.

Read this block first, then the session-8 guide below it for the older context.

---

## SESSION 11 — Phase III: resolving the C9-versus-sALS direction

**Code:** `code/analysis/24_c9_vs_sals_resolution.py` (self-contained, verified to
reproduce every number in REPORT.md §10). **Report:** new §10, plus a correction
appended to §4.6.

### The five things established

1. **Three-group ordering resolves the apparent TREM2/APOE conflict only partly.**
   GPNMB and CHIT1 run control < sALS < C9 in all four regions; TREM2, APOE and
   SPI1 run control < C9 < sALS. Both genotypes engage the myeloid programme —
   the dissociation is which arm. But the bulk causal model compares C9 against
   sALS **with no control arm**, so the conflict is not dissolved (see item 5).

2. **The Myeloid/DAM cluster is two opposing arms.** 15 genes higher in C9, 26
   higher in sALS, 11 mixed. The sALS arm is a coherent receptor cascade
   (TREM2–TYROBP–SYK–NCKAP1L–SPI1 5/5 down in C9; NADPH oxidase 2/2; inflammasome
   4/4; lineage TFs 2/2). The C9 arm is lipid storage and membrane remodelling
   (SOAT1, GPNMB, APOC1; ANO6, STX3, NBEAL2, ANXA2, S100A11, EMP3; SLC5A3,
   SLC38A2; SOD2). **Both arms are induced versus control in both genotypes** —
   each genotype favours one (crossover, descriptive p 1.4e-7).
   C9orf72, SMCR8 and WDR41 themselves show no reduction in C9 donors, so the
   lysosomal loss-of-function interpretation is a hypothesis, not a finding.

3. **CORRECTION — the mean-z score and pseudobulk DE weight genes differently.**
   They agree gene by gene (r +0.48 to +0.71) but can differ in sign at cluster
   level. The score is a mean dominated by a few genes: top 5 of 52 contribute
   87% of the difference in posterior WM, SOD2 alone 36%. Cause: a typical
   myeloid gene is detected in only 3.8–6.6% of spots. **Z-scoring equalises each
   gene's contribution to the score LEVEL, not to the between-group DIFFERENCE** —
   an earlier claim to the contrary was wrong. Figure 12 is unaffected (no
   C9-vs-sALS comparison in it is significant; smallest p 0.190).

4. **Do NOT regress out Cell2location composition.** Collinearity rho(myeloid %,
   DAM expression) = +0.94 in anterior/lateral WM (12% variance left); the
   covariate is a projection of the same matrix; and abundance is a mediator, not
   a confounder. State leakage into the deconvolution is real but not dominant
   (lineage markers +0.62 vs activation +0.65). IBA1 immunofluorescence cannot
   arbitrate — the categorical is saturated and the continuous intensity does not
   track C2L microglia (rho -0.01 within section). Report marginal AND
   abundance-conditional side by side instead.

5. **Aggregation geometry explains much of the myeloid cross-platform conflict.**
   Whole-section pseudobulk (bulk-like) moves Myeloid/DAM from -0.106 to **+0.012**;
   Astrocyte (-0.283 -> -0.208) and Oligodendrocyte/myelin (-0.223 -> -0.280) do
   not reconcile. Composition accounts for only 10-14%; the rest is log-of-a-sum
   weighting. Winner's curse is ruled out (it predicts shrinkage, not sign flip).

### NEXT STAGE — priority order

**On the user's side (highest value first):**
  1. Score the **two arms separately** in the 156-sample bulk cohort. If the
     C9-preferring arm's a-path is stronger than the sALS-preferring arm's, that
     replicates the dissociation across platforms without either dataset
     reversing. Gene lists: `results/tables/myeloid_two_arms_genes.csv`.
  2. Test the **denominator effect**: regress out a neuronal-content proxy
     (NEFL, SNAP25, SYT1, RBFOX3) and see whether the 118-gene a-paths attenuate.
  3. Re-estimate a-paths for TREM2, APOE, SPI1, TYROBP, SYK with **parents
     restricted to non-gene covariates**.
  4. Check whether MuSiC proportions were nodes or covariates in the network
     (potential circularity — derived from the same expression matrix).
  5. RIN/PMI regression on the ribosomal eigengene; rerun mediation with cluster
     eigengenes; report homeostatic:activation ratio alongside MuSiC proportions.

**On this side:**
  - Request from the Gate lab: **per-cell IBA1 area and roughness values with
    donor identity**, split by genotype. Expression-independent state measure,
    already acquired, never reported by subtype. The published IBA1/CD68/DAPI
    counts cannot serve (n = 2 vs 4, fields-of-view pseudoreplication, ALS-vs-control
    only, and CD68 is itself an activation marker).
  - Optional stage 25: limma-voom with duplicateCorrelation(block = donor) for
    pooled-region inference; per-region `~myeloid_fraction + condition` as a
    sensitivity analysis only (NOT as the primary model — see item 4).
  - Housekeeping: `719e52c` and `a5f2883` (2026-10-28) are labelled "Stage 24"
    and "Stage 23"; stage numbers 23 and 24 are each used twice. Renumber once.
  - Rerun stage 17 with mean-z instead of `sc.tl.score_genes` for consistency.
  - Rebuild stage-18 `spatial_DAM_GPNMB_CHIT1_quant.png` GPNMB panel on pseudobulk.
  - Compositional reanalysis of cell proportions (CLR / Dirichlet / scCODA).

### Talk-ready vs not

**Ready:** GPNMB and CHIT1 exceeding sALS in all four regions (agrees with bulk);
ALS-versus-control myeloid expansion (four independent lines); the CD68
pre-specified protein test (d +0.82, one-sided p 0.095).
**Not ready:** the TREM2/APOE direction; any single-number summary of the myeloid
cluster; the two-arm structure (manuscript, not talk).

### Technical notes that cost time this session
  - OpenMP guards MUST be set before importing numpy in a fresh kernel
    (OMP/OPENBLAS/MKL/NUMBA_NUM_THREADS=1, KMP_AFFINITY=disabled,
    OMP_PROC_BIND=FALSE) or the kernel dies with OMP error #179.
  - `auth_cleaned_IBA1_binary` is **not binary** — it is a continuous intensity.
    `auth_IBA1_positive` is saturated at 96–99% and carries no information.
  - Decode categoricals as `np.where(codes >= 0, categories[codes], None)`;
    `auth_manual` stores 3,121 missing values as code -1.

---

**Last updated:** end of session 8 (spatial TDP-43 coupling stage 17 + PBMC arm complete + literature review on cortex vs spinal cord).
**Repo:** github.com/rikku1983/C9orf72-ALS (branch `main`). Local: `/home/sunli/C9orf72-ALS`.
**Latest commit before this note:** `358056b` (Stage 17).

This file is the single entry point to resume. Read it top to bottom for full context.

---


## !! ANNOTATION CORRECTION (session 9) -- READ BEFORE USING STAGES 05-08 !!
Our marker-score `region` column is UNRELIABLE. Only 15% of spots we labelled
'ventral_horn_MN' lie in the authors' Anterior_horns (33% white matter, 21% meninges).
USE `auth_manual` (authors' manual histology) FOR ALL REGION WORK.
Recomputed on correct anatomy (stages 18/19/19b):
  - SURVIVES: TNF-alpha/NF-kB C9-vs-sALS (anterior horn NES +1.86 FDR 0.014;
    other grey matter NES +2.23 FDR <0.001). Inflammatory Response (ant. horn FDR 0.009).
    Hypoxia (grey matter FDR <0.001). IL-6/JAK/STAT3 (ant. horn FDR 0.016).
  - DOES NOT SURVIVE: Unfolded Protein Response + mTORC1 (were FDR 0.000 under the wrong
    label, now FDR 0.21 in true anterior horn) -> the 'proteostasis' story was an artifact.
  - Gene-level C9-vs-sALS: 0 genes FDR<0.05 in true anterior horns (CANX padj 0.0007 -> 0.164).
    39 genes in anterolateral WM, 63 in posterior WM, 2 in other grey matter.
  - DAM/GPNMB/CHIT1 peak = ANTEROLATERAL WHITE MATTER, not ventral horn; dorsal columns spared.
REPORT.md still carries the OLD 'ventral horn' labels and must be revised.

## 1. Research goal
Decode the mechanism of action (MOA) by which the **C9orf72 hexanucleotide repeat expansion**
drives ALS progression — identify genes / pathways / cell populations mediating C9-specific
effects. PRIMARY contrast = **C9-ALS vs sALS** (what makes C9 distinct, not just ALS vs control).

## 2. Datasets (two compartments, UNPAIRED cohorts — no donor matching between them)
- **GSE288365** — 10x Visium **spatial, spinal cord** (postmortem). 20 sections / 18 donors:
  9 control, 5 sALS, 4 C9-ALS. IF-guided (IBA1/MAP2/TDP-43). PRIMARY CNS-lesion dataset.
- **GSE268995** — scRNA-seq **PBMC / peripheral blood** (CITE-seq, 228 imputed surface proteins).
  40 samples: 18 control, 9 sALS-slow, 9 sALS-fast, 4 C9-ALS. Systemic-immune proxy, NOT CNS.
- **User's bulk spinal-cord RNA-seq DEGs** — STILL PENDING, not yet provided. Long-standing
  next step: cross-validate the C9-specific NF-kB/GPNMB program.

## 3. Analysis stages (all in code/analysis/, committed)
**Spatial arm (spinal cord):**
- 01 visium QC; 02 integrate+cluster; 03 region annotation (ventral horn/MN, dorsal horn, WM, gray)
- 05 pseudobulk DE; 06 scores by condition; 07 GSEA/pathway; 08 microglia deep-dive
- 09/09b extract Seurat author meta (R); 10 integrate author meta; 11 protein + stats
- 17 TDP-43 pathology <-> program coupling (this session)

**PBMC arm (blood):**
- 12 scRNA QC; 13 extract author meta (R); 14 merge author meta (barcode-shift resolved)
- 15 pseudobulk aggregate + monocyte DE; 16 myeloid module + surface-protein (ADT) test

## 4. KEY FINDINGS (the through-line)
Central result: a **C9-associated myeloid/DAM/complement/NF-kB program strongest in SPINAL CORD**,
with a **weak directionally-consistent peripheral echo in blood monocytes**.

Spatial (ROBUST when stratified by anatomy):
- DAM composite score rises monotonically control->sALS->C9 (slope +0.158/step, p<1e-4; APOE p=8e-4).
- GSEA C9-vs-sALS ventral horn: NF-kB / complement / myeloid themes (NES up to 2.30).
- ALS ventral-horn spots: GPNMB corr w/ complement (rho=+0.22) and myeloid (rho=+0.19).
- Protein layer ALS-vs-control VH: CD68 D+0.56 (p=0.013), ITGAX D+0.52 (p=0.027).

PBMC (WEAK, mostly NS — 4 C9 donors, underpowered):
- Composition FLAT across conditions (nothing survives FDR).
- Monocyte gene-level DE: nothing survives FDR.
- Peripheral echo: CD68 protein up in C9 monocytes (p=0.033 vs sALS) + CD9/ITGAX/CD63 mRNA up.
  BUT only the activation-marker arm; CNS DAM lipid core (GPNMB/APOE/TREM2/SPP1) NOT expressed
  in blood monocytes -> peripheral signal is partial, not the full CNS program.

Honest negatives:
- Stage 17: TDP-43 pathology does NOT couple to myeloid/complement programs reproducibly.
  A dramatic spot-level effect was PSEUDOREPLICATION; gone at donor level (3 donors/genotype).

## 5. LITERATURE CONTEXT (session 8 review — cortex vs spinal cord)
- Field trends toward **dying-forward / corticomotoneuronal** model (cortex-first, secondary
  spinal degeneration via glutamate excitotoxicity + trans-synaptic TDP-43 spread) — supported
  by TMS cortical hyperexcitability, TDP-43 restricted to corticofugal neurons, Brettschneider
  4-stage pathology staging, and a causal mouse model. NOT settled consensus (origin formally
  unresolved; axonal-spread evidence indirect; hybrid "dying-outward" model exists).
- Implication: cortex may capture UPSTREAM/earlier biology, spinal cord the effector stage.
  C9 carries extra cortical weight (ALS-FTD spectrum, DPR pathology heavier in cortex).
- IMPORTANT TENSION vs our finding: published single-microglia work (Nat Neurosci 2025,
  s41593-025-02075-1) reports C9 haploinsufficiency IMPAIRS microglial DAM/reactive transition
  vs sALS — opposite direction to our spot-level DAM rise. Our Visium can't resolve
  "more microglia" vs "more-activated microglia"; snRNA-seq would adjudicate.

## 6. CANDIDATE NEXT DATASETS (C9 CNS single-nucleus — scoped, not downloaded)
- Spinal cord + C9 (ideal match): microglia paper data is EGA CONTROLLED-ACCESS (EGAD00001009686);
  comparative 85-subject snRNAseq is bioRxiv 2025.08.26.672029 (preprint, likely unreleased).
- Cortex + C9 (openly analyzable, but NOT spinal cord): Li et al. Nat Commun s41467-023-41033-y
  (motor+frontal cortex C9-ALS/C9-FTD/control); Pineda et al. large cortex atlas (33 ALS incl 16 C9).
- ACTION IF PURSUED: verify GEO/Synapse/CELLxGENE download status; test our C9 programs
  (DAM/myeloid/complement/NF-kB) at true single-microglia resolution -> adjudicate the
  "more cells vs impaired activation" tension above.

## 7. PENDING TASKS (priority order)
1. Integrate user's bulk spinal-cord RNA-seq DEGs (awaiting file) — strongest cross-compartment test.
2. Optionally refresh REPORT.md to fold in PBMC arm + stage-17 negative (currently covers 01-11).
3. Optionally pursue a C9 CNS snRNA-seq dataset (cortex openly, spinal cord access-gated).

## 8. CRITICAL TECHNICAL NOTES (preserve — see also compaction archives)
- 7 GB RAM ceiling. DO NOT .to_memory()/.copy() the 280k-cell subset (OOM exit 137).
  For obs-level work read via h5py (decode categoricals); for matrix ops stream per-sample from backed h5ad.
- OpenMP guards in EVERY python entry BEFORE numpy/scanpy: OMP_NUM_THREADS=1 etc; KMP_AFFINITY=disabled.
  AVOID sklearn (OpenMP abort). scipy.stats for FDR/tests.
- env `c9als` (python): scanpy, anndata, h5py, pydeseq2, scipy, pandas, pyarrow, matplotlib,
  seaborn, decoupler, gseapy, statsmodels, squidpy, leidenalg, harmonypy. env `r`: R 4.5.3 + Seurat.
- Push requires token: git push "https://x-access-token:${GITHUB_TOKEN}@github.com/rikku1983/C9orf72-ALS.git" main
- GSE268995 barcode-shift: author RDS file-labels are cyclically shifted vs GSM matrices.
  SHIFT {A1:C1,B1:D1,C1:E1,D1:F1,E1:A1,F1:B1}. GEO series-matrix condition is AUTHORITATIVE.
  All 4 C9 samples (C2,D3,F3,H5) perfectly concordant. Details in code/analysis/14.
- Large intermediate tables (parquet, predicted_ADT.csv.gz ~500M, pseudobulk/) are GITIGNORED —
  regenerate from stage scripts 14-16; key summaries saved as artifacts.

## 9. KEY PATHS
- Spatial processed: data/processed/GSE288365_annotated.h5ad,
  data/processed/GSE288365_annotated_authmeta.h5ad (auth_ TDP-43 pathology cols).
- PBMC processed: data/processed/GSE268995_raw_qc.h5ad (backed, 296k x 36601).
- Metadata: data/metadata/GSE288365_sample_condition_map.csv, GSE268995_sample_condition_map.csv.
- REPORT.md = full methods+results writeup (currently stages 01-11).
