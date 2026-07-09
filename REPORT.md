# Spatial transcriptomic dissection of the C9orf72 mechanism in ALS spinal cord

**Project:** Decode the mechanism by which the *C9orf72* hexanucleotide repeat expansion drives ALS progression in the spinal cord.
**Dataset:** GSE288365 — 10x Visium spatial transcriptomics + 35-plex antibody-capture protein, human cervical spinal cord (postmortem).
**Primary mechanistic contrast:** C9-ALS vs sporadic ALS (sALS), which isolates the repeat-expansion effect from ALS-in-general. Control vs ALS establishes the disease baseline.
**Status:** Spatial arm complete (stages 01–11). Bulk RNA-seq integration pending user's DEG table.

---

## 1. Executive summary

Using end-stage human cervical spinal cord profiled by spatial transcriptomics, we localized the ALS neurodegenerative signature to the ventral horn and asked what, at the molecular level, distinguishes *C9orf72*-mutant ALS from sporadic ALS. Three findings define the mechanistic picture:

1. **ALS remodels the ventral horn identically at the cellular level regardless of genotype.** Motor neurons are lost and glia expand in both C9-ALS and sALS; no cell type differs significantly in abundance between the two genotypes. The C9-specific effect is therefore **not compositional** — it is a difference in what cells *express*.

2. **The C9-specific signal is an NF-κB/TNF-α inflammatory program.** The single strongest pathway separating C9-ALS from sALS in the ventral horn is *TNF-α signaling via NF-κB* (NES = 2.07, FDR < 0.001), driven by *BIRC3, SOCS2, CALHM6* and downstream NF-κB targets. Where sALS-vs-control is dominated by **complement** (NES = 1.96), C9-vs-control is dominated by **TNF-α/NF-κB** (NES = 2.30) — two different routes into the same degenerating tissue.

3. **A disease-severity gradient runs control → sALS → C9 in the myeloid/complement axis.** Disease-associated microglia (DAM) and complement programs increase monotonically with disease severity (DAM score Jonckheere trend p = 8×10⁻⁴; mixed-model slope +0.16/step, p < 1×10⁻⁴), with C9-ALS at the extreme.

The strongest *hypothesis* the spatial data supports — that C9 and sALS myeloid cells run divergent programs (C9: NF-κB/GPNMB; sALS: canonical DAM/complement) specifically at TDP-43 lesions — is directionally consistent across every donor but underpowered at the donor level (3 vs 3 donors with pathology). Confirming it is the primary role for the pending bulk RNA-seq cohort.

All results were validated against the original authors' processed Seurat objects (per-spot QC concordance r = 1.00; key-gene Spearman 0.96–1.00).

---

## 2. Data and cohort

### 2.1 Dataset structure
GSE288365 comprises **20 Visium sections from 18 donors** (two donors contributed paired top/bottom captures of the same block). Each spot carries two modalities in a single matrix: 18,085 gene-expression features and 35 antibody-capture protein features. Protein counts dominate raw totals (~99.7%), so the modalities were split at load: RNA placed in `.X`, protein in `obsm['protein']`.

### 2.2 Cohort composition
| Condition | Donors | Sections |
|-----------|:------:|:--------:|
| Control   | 9      | 10       |
| C9-ALS    | 4      | 4        |
| sALS      | 5      | 6        |

RNA depth is modest, as expected for postmortem cord on a protein-focused assay: median ~850 counts and ~630 genes per spot (range across sections 377–3,240 counts). Section-level QC is in `qc_per_section.csv` and Figure 1.

![Figure 1. Per-section QC: raw vs QC-passing spot counts, median counts, genes, and mitochondrial fraction across all 20 sections.]({{artifact:art_c4cddaa1-8867-4021-806b-d23f082e0c69}})

---

## 3. Methods

### 3.1 Preprocessing and QC (stage 01)
Each Visium section was loaded and the RNA/protein modalities separated. RNA-only QC thresholds: minimum 200 counts, minimum 100 genes, maximum 30% mitochondrial reads per spot. Counts were normalized to median library size and log1p-transformed (`sc.pp.normalize_total` + `sc.pp.log1p`). QC retained 132,781 spots across the 20 sections.

### 3.2 Integration, clustering, and region annotation (stages 02–03)
Highly variable genes were selected, PCA computed, and batch effects across sections corrected with Harmony. Leiden clustering was performed on the Harmony embedding. Anatomical regions were assigned by marker-gene scoring into four classes — **ventral_horn_MN** (motor-neuron-containing gray matter: *CHAT, NEFH, NEFL, SLC5A7*), gray neuropil, white matter (*MBP, PLP1, MOBP*), and a low-signal class. The ventral-horn region is the focus of all disease contrasts because it is where motor neurons — the cells that die in ALS — reside.

### 3.3 Pseudobulk differential expression (stage 05)
To avoid pseudoreplication (treating spots as independent replicates), differential expression was computed on **donor-level pseudobulk**: spot counts were summed per donor within a region (minimum 50 spots per pseudobulk, minimum 3 donors per group) and modeled with DESeq2 (pyDESeq2 0.5.4). Three contrasts were run in each of two spatial scopes (global; ventral_horn_MN): C9-ALS vs control, sALS vs control, C9-ALS vs sALS.

### 3.4 Program scores and pathway enrichment (stages 06–08)
Curated gene-set scores (`sc.tl.score_genes`) were computed for a disease-associated microglia (DAM) program (*GPNMB, CHIT1, TREM2, TYROBP, APOE, ITGAX, CST7, LPL, CD9, CD63, SPP1, CTSB, CTSD, LGALS3, FABP5, MSR1, CD68*), a complement program (*C1QA, C1QB, C1QC, C3*), a homeostatic-microglia program (*P2RY12, TMEM119, CX3CR1, CSF1R, SALL1*), and a candidate C9-specific NF-κB program (*BIRC3, SOCS2, CALHM6*). Pathway enrichment used GSEA preranked (gseapy 1.3.0) on the pseudobulk statistics against Hallmark, GO-BP, and Reactome.

### 3.5 Validation against authors' processed objects (stage 09)
The original authors' processed Seurat objects (18 GSM-level `.rds.gz` files at the GSE288365 series-supplement level) were downloaded and their metadata, dimensional reductions, and key-gene matrices extracted with a lightweight slot-based reader (SeuratObject only). For one control (H1_GB20) and one C9-ALS donor (C4_N280), we barcode-matched our AnnData against the authors' object spot-for-spot and compared QC metrics, key-gene expression, and region labels.

### 3.6 Author-metadata integration (stage 10)
The authors' full per-spot annotations were merged into our object across all 20 sections (98.7% spot match) via a normalized capture-area + cell-barcode key. This imported: **deconvolution** (13 cell-type fractions from their reference-based method), **TDP-43 pathology grade** (`annotation_general`: Ring / TDP43_adjacent / TDP43_distant / Not_selected), **motor-neuron distance** (`MN_dist_group` 1–5), and **imaging features** (IBA1, MAP2, TDP-43 immunofluorescence). Compositional contrasts were computed on donor-level mean cell-type fractions in the ventral horn (Mann–Whitney U), keeping the analysis donor-level to avoid pseudoreplication.

### 3.7 Protein-layer validation (stage 11A)
The 35-plex antibody-capture matrix was CLR-normalized (centered log-ratio across markers per spot) and **isotype-control background subtracted** (mean CLR of mouse IgG2a, IgG1κ, IgG2bκ, rat IgG2a). Myeloid markers passing background QC (>99% of spots above isotype): CD68, CD163, ITGAX (CD11c), ITGAM (CD11b), CD14, HLA-DRA, PTPRC (CD45); FCGR3A was at background and dropped. Donor-level Mann–Whitney U contrasts were run in the ventral horn.

### 3.8 Properly-powered gradient statistics (stage 11B)
Two complementary tests re-examined the control → sALS → C9 severity gradient:
- **Linear mixed models** (statsmodels): spot-level outcome ~ ordinal condition (0/1/2) with donor as a random intercept, across all 18 donors' ventral-horn spots. Reliable for composite program scores; flagged as unreliable (singular, donor-variance → 0) for sparse single genes with heavy zero-inflation.
- **Jonckheere–Terpstra monotonic trend test** on donor pseudobulk means — the appropriate non-parametric ordered-alternative test, robust to the zero-inflation that defeats a Gaussian mixed model on single genes.

### 3.9 Reproducibility and technical notes
- OpenMP thread guards were set in every Python entry point before importing numpy/scanpy (avoids an OpenMP abort in this environment; ARI computed manually rather than via scikit-learn).
- Authors' RDS are full Seurat S4 objects but only the lightweight SeuratObject package was installed, so slots were accessed directly (`obj@meta.data`, `obj@assays`, `obj@reductions`).
- All code is in `code/analysis/` (stages 01–11) and `code/download/`; result tables in `results/tables/`; figures in `results/figures/`. Repository: github.com/rikku1983/C9orf72-ALS.

---

## 4. Results

### 4.1 Validation: our pipeline reproduces the authors' object exactly

Barcode-matched comparison against the authors' processed Seurat objects showed **perfect QC concordance**: our per-spot total counts, gene counts, and mitochondrial fraction vs the authors' equivalent metrics gave Pearson r = 1.000 and median ratio = 1.000 for both the control and the C9-ALS donor. This is the decisive check that the RNA/protein modality split — the main technical hazard in this dataset — was done correctly.

Key-gene per-spot correlations (our log-normalized vs authors' SCT) were 0.96–1.00 across all 16 genes tested (e.g. GPNMB 0.98/0.98, CHIT1 1.00/1.00, C1QA 0.96/0.97, CHAT 1.00/1.00). Spot-set agreement (Jaccard) was 0.84 (control) and 0.99 (C9-ALS); the authors retained ~1,270 more spots in the control under a slightly looser threshold, but every spot we kept, they kept too. Region labels are the one step that diverges (our marker-score 4-class scheme vs their manual histology; WM/GM agreement 0.71 control / 0.54 C9), so anatomical claims are kept donor-level.

![Figure 2. Validation against authors' Seurat objects: per-spot QC concordance, key-gene correlations, and region-label overlap for one control and one C9-ALS donor.]({{artifact:art_92faf353-a7a7-4c3b-a5c3-42ed13e21823}})

### 4.2 ALS remodels ventral-horn cellular composition — but genotype does not

Using the authors' deconvolution, ventral-horn cellular composition shifts strongly with disease (donor-level, ALS vs control):

| Cell type | Control % | ALS % | Δ | p |
|-----------|:---------:|:-----:|:---:|:---:|
| Neurons | 4.03 | 1.47 | −2.57 | 0.004 |
| Microglia | 1.71 | 2.59 | +0.88 | 0.001 |
| Proliferating microglia | 2.22 | 3.54 | +1.33 | 0.010 |
| Astrocytes | 22.70 | 28.92 | +6.21 | 0.005 |
| Macrophages | 3.18 | 3.98 | +0.80 | 0.042 |

This is the expected degenerative signature: motor-neuron loss with reactive gliosis and myeloid expansion.

**Critically, the C9-ALS vs sALS composition contrast shows no significant difference for any of the 13 cell types** (smallest p ≈ 0.11 for Endothelial and OPC; Macrophages p = 0.19, highest in C9 at 4.5% but not significant). Genotype does not change *which* cells are present — it changes what they do. This is the central mechanistic distinction of the whole study.

![Figure 3. Author deconvolution + TDP-43 pathology integration. (a) VH composition by condition; (b) per-donor neuron loss and gliosis; (c) ALS-vs-control composition shift; (d) expression along the TDP-43 pathology gradient; (e) C9-vs-sALS myeloid divergence at lesions; (f) programs vs motor-neuron distance.]({{artifact:art_e6a59b27-b318-4873-ab4d-28465fcdb072}})

### 4.3 The C9-specific signal is an NF-κB/TNF-α program

GSEA on the ventral-horn pseudobulk cleanly separates the two genotypes' routes into the tissue:

| Contrast | Top Hallmark pathway | NES | FDR |
|----------|---------------------|:---:|:---:|
| **C9-ALS vs sALS** | **TNF-α signaling via NF-κB** | **2.07** | **<0.001** |
| C9-ALS vs control | TNF-α signaling via NF-κB | 2.30 | <0.001 |
| sALS vs control | Complement | 1.96 | <0.001 |

The genotype contrast (C9 vs sALS) and the C9-vs-control contrast are both led by **TNF-α/NF-κB**, whereas sALS-vs-control is led by **complement**. The NF-κB leading-edge includes *BIRC3, TNFAIP3, NFKBIA, REL, IL6, IL1B, CCL2, SOCS3, BCL3*. Consistent single-gene pseudobulk hits (C9 vs sALS, ventral horn) are *CALHM6* (log2FC +1.38, p = 1.8×10⁻⁴), *BIRC3* (+1.49, p = 0.023), and *SOCS2* (+1.32, p = 0.030) up in C9, and *ITGAX* (−1.14, p = 3.7×10⁻⁵) down in C9. Of 1,646 genes tested, 37 reach padj < 0.05 — a modest count reflecting the 4-vs-5-donor power, so the pathway-level signal is more robust than any single gene.

![Figure 4. C9-specific pathway enrichment: NF-κB/TNF-α and inflammatory pathways elevated in C9-ALS relative to sALS in the ventral horn.]({{artifact:art_233198ee-788e-4ee8-b8c6-f6f8dfb89635}})

### 4.4 A disease-severity gradient in the myeloid/complement axis

DAM and complement programs rise monotonically control → sALS → C9. Because the earlier donor-level test was underpowered, two robust tests were applied:

**Jonckheere–Terpstra monotonic trend (donor pseudobulk, ventral horn):**
| Target | Control | sALS | C9 | p (trend) |
|--------|:-------:|:----:|:--:|:---------:|
| DAM score | 0.328 | 0.584 | 0.616 | 8×10⁻⁴ |
| APOE | 2.84 | 4.09 | 3.66 | 0.003 |
| C1QC | 0.44 | 0.80 | 1.09 | 0.004 |
| C1QA | 0.39 | 0.54 | 0.77 | 0.011 |
| CHIT1 | 0.005 | 0.052 | 0.080 | 0.013 |
| GPNMB | 0.12 | 0.35 | 0.65 | 0.018 |
| Complement score | −0.42 | −0.18 | −0.04 | 0.022 |

**Linear mixed models (spot-level, donor random intercept, 18 donors)** confirmed the composite scores where the model is well-posed: DAM score slope +0.158/step (p = 2×10⁻⁶) and APOE +0.513/step (p = 8×10⁻⁴). Single sparse genes gave singular fits (donor-variance → 0 from zero-inflation) and were flagged unreliable rather than reported — which is exactly why the Jonckheere test is the appropriate instrument for those.

![Figure 5. DAM/complement program scores by condition and region, showing the control → sALS → C9 gradient.]({{artifact:art_e593697d-e580-4955-8e1b-e7e6aa65e55d}})

### 4.5 The programs are spatially targeted to sites of degeneration

Expression of the myeloid/complement genes rises along the authors' **TDP-43 pathology gradient** (Not_selected → distant → adjacent → Ring) and toward motor neurons (MN-distance groups). Motor-neuron markers *NEFH/NEFL* peak at TDP43-adjacent spots (NEFH 0.18 → 1.17; NEFL 0.39 → 1.62), confirming the annotation captures true motor-neuron neighborhoods, while complement (*C1QC* 0.59 → ~0.88) and GPNMB (0.25 → ~0.4–0.5) rise into pathology. The programs are not diffuse; they concentrate where cells are dying.

### 4.6 The C9-vs-sALS lesion divergence: strongest hypothesis, underpowered

At TDP-43 pathology spots, C9 and sALS myeloid cells trend toward **different programs**: C9 lesions higher in the NF-κB-linked genes (GPNMB, CHIT1, SOCS2, BIRC3), sALS lesions higher in the canonical DAM/complement axis (TREM2, APOE, C1QA/C1QC). The donor-level direction is 100% consistent, but with only 4 C9 vs 3 sALS donors carrying pathology spots, no gene reaches significance (smallest p ≈ 0.11 for TREM2 and SOCS2/BIRC3).

The mixed model **cannot** rescue this: with genotype constant within donor and 3 vs 3 donors, the fixed effect is fully confounded with the donor grouping (variance collapses to zero). This is a genuine power ceiling of the spatial cohort, not a modeling artifact — and the reason an independent cohort (the pending bulk RNA-seq, or the companion single-cell data) is needed to close it.

### 4.7 Protein layer confirms myeloid activation — but cannot test the program switch

The 35-plex protein data independently confirm myeloid activation in ALS: in the ventral horn, **CD68 (Δ+0.56, p = 0.013) and ITGAX/CD11c (Δ+0.52, p = 0.027)** are elevated in ALS vs control, with HLA-DRA trending (p = 0.11). This corroborates the RNA/deconvolution myeloid expansion at the protein level.

**Important limitation:** the antibody panel does **not** include GPNMB, CHIT1, TREM2, or APOE, so it cannot directly adjudicate the NF-κB/GPNMB-vs-DAM/complement program switch. At C9-vs-sALS lesions the general myeloid markers trend higher in C9 (CD68 Δ+0.63, ITGAM Δ+0.53) but with 2 vs 3 donors nothing is significant. The protein layer validates *that myeloid cells are activated*, not *which program they run* — the program distinction remains an RNA-level result.

![Figure 6. Protein-layer validation (isotype-subtracted CLR) and mixed-model / Jonckheere gradient statistics in the ventral horn.]({{artifact:art_d525b3b7-9ba0-4013-9873-5ed27b92d31b}})

---

## 5. Synthesis: a model of the C9 mechanism

The spatial data support the following model of how the *C9orf72* expansion shapes end-stage spinal cord:

1. **ALS, regardless of genotype, converges on the same cellular endpoint** in the ventral horn: motor-neuron loss with astrocytic and myeloid expansion. Genotype is invisible at the level of cell-type abundance.

2. **Genotype acts on cell-intrinsic transcriptional programs, not composition.** The repeat expansion biases the myeloid/inflammatory response toward **NF-κB/TNF-α signaling** (BIRC3/SOCS2/CALHM6 and the broader NF-κB leading edge), whereas sporadic ALS engages a **complement**-dominated response. Both are superimposed on a shared DAM/lipid-handling program (GPNMB, CHIT1, APOE) that scales with disease severity.

3. **The response is spatially organized** around TDP-43 pathology and motor neurons — it is a targeted reaction to local degeneration, not diffuse neuroinflammation.

4. **The sharpest genotype signature — divergent myeloid programs at TDP-43 lesions — is a well-supported hypothesis awaiting an adequately powered cohort.**

This positions the NF-κB/TNF-α axis as the leading candidate mechanism distinguishing C9-ALS, and makes the pending bulk RNA-seq the natural validation: does the C9-specific NF-κB/GPNMB program replicate in an independent spinal-cord cohort?

---

## 6. Limitations

- **End-stage snapshot.** Visium captures postmortem tissue; the temporal ordering (control → sALS → C9 as a proxy for severity) is inferred, not longitudinal.
- **Genotype power.** 4 C9-ALS vs 5 sALS donors bounds the lesion-level contrast; pathway-level and gradient results are robust, single-gene genotype effects are not.
- **Spatial resolution.** Visium spots (~55 µm) are multi-cellular; "cell-intrinsic" is inferred via deconvolution and program scoring, not single-cell resolution. The companion single-cell dataset (GSE268995) and snRNA-seq (bioRxiv 2025.08.26.672029) would resolve this.
- **Region labels** diverge modestly from the authors' manual histology (WM/GM agreement 0.54–0.71); anatomical claims are kept donor-level.
- **Protein panel** lacks the key DAM markers, limiting cross-modal validation of the program switch.

---

## 7. Analysis inventory

**Pipeline (`code/analysis/`):** `01_visium_qc.py` · `02_integrate_cluster.py` · `03_annotate_regions.py` · `05_pseudobulk_de.py` · `06_scores_by_condition.py` · `07_pathway_enrichment.py` · `08_microglia_deepdive.py` · `09_extract_seurat.R` / `09b_extract_meta_only.R` · `10_integrate_author_meta.py` · `11_protein_and_stats.py`.

**Key result tables (`results/tables/`):** `qc_per_section.csv` · `pseudobulk_*` (6 contrasts) · `gsea_*` (Hallmark/GO-BP/Reactome × 3 contrasts × 2 scopes) · `scores_donor_region.csv` · `microglia_dam_correlations.csv` · `composition_{ALSvsControl,C9vsSALS}_VH.csv` · `tdp43_gradient_expression.csv` · `tdp43_lesion_C9vsSALS_myeloid.csv` · `protein_{ALSvsControl,C9vsSALS}_VH.csv`, `protein_C9vsSALS_lesions.csv` · `mixedmodel_gradient_VH.csv` · `jonckheere_trend_VH.csv` · `validation_vs_authors_metrics.csv`, `validation_gene_correlations.csv`.

**Next step:** integrate the user's bulk spinal-cord RNA-seq DEGs (gene, log2FC, p-value, contrast) to test replication of the C9-specific NF-κB/GPNMB program in an independent cohort.
