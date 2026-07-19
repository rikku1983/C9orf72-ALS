# C9orf72-ALS Project — State / Resume Guide

**Last updated:** end of session 8 (spatial TDP-43 coupling stage 17 + PBMC arm complete + literature review on cortex vs spinal cord).
**Repo:** github.com/rikku1983/C9orf72-ALS (branch `main`). Local: `/home/sunli/C9orf72-ALS`.
**Latest commit before this note:** `358056b` (Stage 17).

This file is the single entry point to resume. Read it top to bottom for full context.

---

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
