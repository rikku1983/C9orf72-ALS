# C9orf72-ALS Spatial Transcriptomics — Project State / Resume Guide

**Last updated:** end of session 2 (GPNMB/CHIT1 microglia deep-dive).
**Repo:** github.com/rikku1983/C9orf72-ALS (branch `main`). Local: `/home/sunli/C9orf72-ALS`.

This file is the single entry point to resume. Read it top to bottom and you have the full context.

---

## 1. Research goal

Decode the mechanism of action (MOA) by which the **C9orf72 hexanucleotide repeat
expansion** drives ALS progression — identify genes / pathways / cell populations that
mediate C9-specific effects in the **spinal cord** (cervical, postmortem). Plan: analyze
spatial transcriptomics first, then integrate with the user's existing **bulk RNA-seq**
spinal-cord DEGs (bulk list NOT yet provided).

**Primary contrast for MOA = C9-ALS vs sALS** — isolates what the repeat expansion adds
*beyond ALS itself*. C9-vs-control and sALS-vs-control capture the shared disease program.

**Framing caveat:** Visium here is an *end-stage* snapshot — it cannot measure progression
*rate* directly. The logic is to test whether C9 carries more of the molecular
"fast-progression signature" (complement/myeloid, per the companion paper) than sporadic ALS.

---

## 2. Primary dataset — GSE288365

- 10x **Visium spatial** + **35-plex Antibody Capture protein** (CITE-like), human cervical spinal cord.
- **20 sections / 18 donors** (2 donors have paired TOP+BOT sections).
- **9 control + 4 C9-ALS + 5 sALS** donors.
- Downloaded flat files only (466 MB); SeuratObject.rds files EXCLUDED (~8 GB) to save space.
- **Modality split:** the matrix mixes 18,085 Gene-Expression genes + 35 protein markers.
  Protein = 99.7% of raw counts, so it MUST be split out or QC is meaningless. Stage 01
  puts RNA in `.X` and the 35 protein markers in `obsm['protein']`.
- RNA depth is low (~850 median counts, ~630 genes/spot) — expected for postmortem cord +
  protein-focused assay. Workable but the reason we use donor-level pseudobulk, not spot-level tests.
- 35 protein markers include CD163, CD68, HLA-DRA, PTPRC, CR2 (immune/myeloid/complement)
  + 4 isotype controls (mouse_IgG2a/IgG1k/IgG2bk, rat_IgG2a) for background subtraction.
  **Protein layer NOT yet analyzed** (candidate next step).

**Companion paper** (Nat Neurosci, DOI 10.1038/s41593-026-02300-5): most unique
transcriptional alterations occur in C9-ALS; complement + lipid-programmed myeloid states
converge at motor-neuron loss / TDP-43 pathology; complement genes track fast progression.

**Other datasets (not used yet):** GSE268995 (scRNA-seq companion), GSE272626 (bulk 266
spinal cord validation), bioRxiv 2025.08.26.672029 (embargoed snRNA-seq).

---

## 3. Pipeline — all scripts in `code/analysis/`, run in conda env `c9als`

| Stage | Script | What it does | Key outputs |
|---|---|---|---|
| 01 | `01_visium_qc.py` | Load 20 sections, split RNA/protein, QC filter, merge | `GSE288365_qc.h5ad`, qc_per_section.{png,csv} |
| 02 | `02_integrate_cluster.py` | Log-norm, HVG, PCA, Harmony (batch=section), Leiden | `GSE288365_lognorm.h5ad`, `_clustered.h5ad`, umap_*.png |
| 03 | `03_annotate_regions.py` | Region labels + cell-type module scores | `GSE288365_annotated.h5ad`, region_by_cluster.csv |
| 05 | `05_pseudobulk_de.py` | Donor-level DESeq2, 3 contrasts × 2 scopes | pseudobulk_*.csv (6 files) |
| 06 | `06_scores_by_condition.py` | Cell-type scores by condition (donor means) | scores_by_condition.png, scores_donor_region.csv |
| 07 | `07_pathway_enrichment.py` | GSEA prerank (Hallmark/GO-BP/Reactome) × all contrasts | gsea_*.csv (18 files) + gsea_all_summary.csv |
| 08 | `08_microglia_deepdive.py` | GPNMB/CHIT1 + DAM vs homeostatic microglia program | microglia_deepdive.png, microglia_*.csv |

(There is no stage 04 — numbering skipped.)

**To fully reproduce from scratch:** re-run `code/download/download_GSE288365_visium.py`
to fetch raw flat files → `01` → `02` → `03` → then `05`,`06`,`07`,`08` in any order (all
read `GSE288365_annotated.h5ad` or the pseudobulk tables).

---

## 4. Key findings so far

### The shared ALS signature (present in BOTH C9 and sALS, vs control)
Concentrated in the **ventral horn** (motor-neuron zone):
- **Motor-neuron loss:** CHAT (−2.3), NEFH (−2.7), NEFL (−3.3), all padj ≤ 0.005.
- **Complement activation:** C1QA/B/C up (padj ≤ 0.01).
- **Myeloid/DAM activation:** MSR1, CD68, AIF1, TYROBP, GPNMB, CHIT1, CHI3L1/2 up.
- In the direct C9-vs-sALS contrast these are **n.s.** → shared ALS engine, not C9-specific.

### The C9-SPECIFIC signature (C9-ALS vs sALS)
Much smaller (128 DE genes global / 37 ventral horn; 74–89 GSEA sets vs >1000 for either-vs-control).
- **UP in C9:** TNF-α signaling via NF-κB (NES 2.07, q<0.001; genes BIRC3, SOCS2, CALHM6);
  HSF1-mediated heat-shock / proteostasis (q≈0.015); RIPK1-mediated necroptosis (q≈0.01).
- **DOWN in C9:** synaptic vesicle exocytosis / neurotransmitter secretion (matches CPLX3
  down); glial + oligodendrocyte differentiation.
- Gene-level C9-specific: **BIRC3, SOCS2, CALHM6 UP; ITGAX (CD11c) DOWN.**

### GPNMB / CHIT1 microglia deep-dive (stage 08)
- Both strongly UP in ALS vs control, steepest in ventral horn (CHIT1 log2FC +4.3, GPNMB +3.1).
- **Shared, not C9-specific:** C9-vs-sALS n.s. (padj 0.3–0.6); ordering control<sALS<C9 but
  within donor scatter (only 4 C9 / 5 sALS donors → low power for this contrast).
- Co-express with lipid/phagocytic DAM program: APOE, CTSB, CTSD, CD68, LGALS3, TYROBP, SPP1.
- DAM score rises in ALS while homeostatic (P2RY12/TMEM119/CX3CR1) does not → activation shift.
- **GPNMB spatially tracks degeneration:** in ALS ventral-horn spots, +correlated with
  complement (rho +0.22) and myeloid (+0.19), **−correlated with motor-neuron signature (−0.19)**.

### One-line synthesis
Neurodegeneration + complement + DAM myeloid activation = the common ALS engine (both C9 and
sporadic). What the C9 repeat adds on top = **NF-κB/TNF inflammation + HSF1 proteostasis stress
+ RIPK1 necroptosis, plus a sharper synaptic-transmission and glial-differentiation deficit.**

---

## 5. Critical technical solutions (MUST preserve — these cost real debugging time)

1. **Sandbox env guards in EVERY script**, at top before numpy/scanpy import (else OpenMP
   `pthread_setaffinity` abort, exit 134):
   ```python
   os.environ.setdefault('OMP_NUM_THREADS','1')
   os.environ.setdefault('KMP_AFFINITY','disabled')
   os.environ.setdefault('OPENBLAS_NUM_THREADS','1')
   ```
2. **Memory ceiling ~7 GB** for 132k spots × 18k genes. (a) never hold full matrix + dense HVG
   at once; (b) SKIP `sc.pp.scale` (densifies → OOM) — run `sc.tl.pca` directly on sparse
   log-norm with `zero_center=True`; (c) stage 02 is split into 3 idempotent steps persisting
   only tiny embeddings (`_embeddings.npz`) between clustering and merge.
3. **Harmony bug:** this harmonypy version returns Z_corr as (cells×PCs); scanpy's
   `harmony_integrate` wrapper assumes (PCs×cells) and transposes → corrupts to shape (30,).
   FIX: call `harmonypy.run_harmony()` directly, assign shape-robustly (transpose if
   `Z.shape[0] != n_obs`).
4. **Pseudobulk, not spot-level:** 132k spots from only 18 donors → spot-level Wilcoxon =
   pseudoreplication. Aggregate raw counts per donor (global) and donor×region, DESeq2 across
   donors. `MIN_SPOTS_PER_PSEUDOBULK=50`, `MIN_DONORS_PER_GROUP=3`.
5. **Git push (sandbox can't persist git remotes — .git/config read-only):** in `repl`,
   `tok = host.credentials.get('GitHub')['token']`, base64-encode `x-access-token:{tok}`,
   write to `handoff/gh.json`; in `bash` read it and
   `git -c http.extraheader="AUTHORIZATION: basic $B64" push https://github.com/rikku1983/C9orf72-ALS.git main`;
   then `rm` the token file. GitHub username `rikku1983`, name LI SUN.
6. **save_artifacts:** copy repo files into the workspace dir first, then save by basename.
7. **PubMed MCP:** tools are `search_articles` (arg `query`, `max_results`) and
   `get_article_metadata` (arg `pmids=[...]`, a LIST). Enrichr (`maayanlab.cloud`) required a
   `request_network_access` grant for GSEA gene-set libraries — already granted this session.

---

## 6. Data files (all git-ignored — large, on ephemeral workspace disk)

`data/processed/`: GSE288365_qc.h5ad (1.3G), _lognorm.h5ad (904M), _clustered.h5ad (936M),
**_annotated.h5ad (947M) ← the main analysis object** (has region labels + ct_ scores),
pseudobulk_counts.csv, pseudobulk_info.csv, _embeddings.npz.
`data/raw/GSE288365/` (466M, 140 flat files).

**IMPORTANT:** these are NOT in git and the workspace disk is swept after long idle gaps.
If they're gone on resume, regenerate: re-download (script in `code/download/`) then run
stages 01→02→03. All result CSVs and figures ARE committed to git, so findings are safe even
if the .h5ad objects are lost.

---

## 7. Pending / candidate next steps (pick up here)

1. **Protein-layer validation** — analyze `obsm['protein']` (35 markers): CD68, CD163,
   HLA-DRA, PTPRC, CR2 as direct myeloid/complement readouts, with isotype-control
   (mouse/rat IgG) background subtraction. Corroborate RNA DAM/complement scores at protein level.
2. **Bulk RNA-seq integration (the original plan)** — awaiting user's spinal-cord DEG CSV
   (gene + log2FC + p-value + contrast). Test overlap of bulk C9 signature with spatial DE;
   check whether bulk hits localize to specific regions/cell types.
3. **Mixed-effects re-test** — does the control<sALS<C9 DAM/complement gradient survive a
   proper donor-random-effect model? (current pseudobulk = fixed-effect DESeq2.)
4. **Cell-type deconvolution** — Visium spots are multi-cell; consider deconvolution
   (cell2location / RCTD) using the GSE268995 scRNA-seq companion as reference for finer
   cell-population attribution.
5. **Write-up** — assemble stages 01–08 into a methods/results narrative.
6. Optional: track bioRxiv 2025.08.26.672029 for the snRNA-seq release.

---

## 8. Environment

Conda env **`c9als`**: scanpy, squidpy, anndata, leidenalg, python-igraph, harmonypy,
scikit-misc, decoupler, gseapy (1.3.0), pydeseq2 (0.5.4), matplotlib, seaborn, pandas, numpy,
scipy, h5py. Spec in `env/environment.yml`. Always pass `environment="c9als"` to python/bash/r.

## 9. Saved artifacts (Claude Science project)

Figures: qc_per_section, umap_leiden, umap_condition, region_scores_umap,
scores_by_condition, gsea_c9_specific, microglia_deepdive.
Tables: pseudobulk_* (6), gsea_all_summary, scores_donor_region, microglia_* (2),
cluster_composition, region_by_cluster, dataset inventory, sample_condition_map.
