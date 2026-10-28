# Spatial and single-cell dissection of the C9orf72 mechanism in ALS

**Project:** Decode the mechanism by which the *C9orf72* hexanucleotide repeat expansion drives ALS progression.
**Datasets:** GSE288365 — 10x Visium spatial transcriptomics + 35-plex antibody-capture protein, human spinal cord (postmortem). GSE268995 — PBMC scRNA-seq/CITE-seq (unpaired cohort).
**Primary mechanistic contrast:** C9-ALS vs sporadic ALS (sALS), which isolates the repeat-expansion effect from ALS-in-general. Control vs ALS establishes the disease baseline.
**Status:** Spatial arm complete through stage 20; PBMC arm complete. Bulk RNA-seq integration pending user's DEG table.

> ### Revision notice (session 9)
> This report was substantially revised. The earlier version stratified every spatial
> contrast by **our marker-score region labels**, which turned out to be badly wrong:
> only **15%** of spots we called `ventral_horn_MN` fall inside the authors' manually
> annotated **anterior horns** (33% were white matter, 21% meninges — see
> `region_annotation_disagreement.csv`). The failure mode was circular: the marker score
> used myeloid/DAM genes, which are highest in white matter, so it selected inflamed
> white matter and labelled it ventral horn.
> **All region-stratified analyses have been recomputed on the authors' manual histology
> (`auth_manual`)** in stages 18–20. Claims that did not survive that correction are
> listed in §7 and have been removed from the body of this report.

---

## 1. Executive summary

Using end-stage human spinal cord profiled by spatial transcriptomics, with all regions defined by the original authors' manual histology, four findings define the mechanistic picture:

1. **ALS expands the myeloid compartment of the anterior horn; genotype does not change composition.** In the anterior horns, proliferating microglia (+2.25 pp, p = 0.006), microglia (+1.08 pp, p = 0.013) and macrophages (+1.27 pp, p = 0.022) all increase in ALS vs control. **No cell type differs between C9-ALS and sALS** (smallest p = 0.063, astrocytes). The genotype effect is not compositional — it is a difference in what cells express.

2. **The C9-specific signal is NF-κB/TNF-α signalling, and it is localized to grey matter.** *TNF-α signalling via NF-κB* is the top pathway separating C9-ALS from sALS in the anterior horn (NES = 1.86, FDR = 0.014) and in other grey matter (NES = 2.23, FDR < 0.001), accompanied by inflammatory response (FDR = 0.009) and IL-6/JAK/STAT3 (FDR = 0.016). It is **absent in both white-matter compartments** (FDR = 0.53, 0.20). The cleanest form of the result: in other grey matter, C9 engages NF-κB (vs control NES +2.41) where **sALS does not** (vs control NES −1.09).

3. **An ordered contrast (control → sALS → C9) through the myeloid/lipid axis is significant in the true anterior horn** — though see the caveat below, this is driven almost entirely by the control-versus-ALS step: CHIT1 (ρ = 0.76, p < 0.001), DAM score (ρ = 0.72, p = 0.001), GPNMB (ρ = 0.71, p = 0.001), C1QC (ρ = 0.69, p = 0.002), complement score (ρ = 0.66, p = 0.003). The homeostatic-microglia score does not track it (p > 0.05).

4. **The DAM/GPNMB/CHIT1 program is anatomically targeted to motor tracts, not to the anterior horn.** Its peak is **anterolateral white matter** — higher than every other region in 9/9 ALS donors — while the **dorsal columns are spared** (GPNMB and CHIT1 not elevated vs control there, p = 0.38, 0.33). This motor-tract selectivity is consistent with corticospinal degeneration.

**What is *not* supported.** The DAM/GPNMB/CHIT1 program shows **no C9-vs-sALS difference in any region** (p = 0.19–0.91). Neither does myeloid abundance. So the myeloid/lipid axis is a shared ALS marker, not a C9-specific mechanism; the C9-specific signal is confined to NF-κB-family signalling in grey matter.

All results were validated against the original authors' processed Seurat objects (per-spot QC concordance r = 1.00; key-gene Spearman 0.96–1.00).

![Corrected summary. (a) Anterior-horn composition, ALS vs control — myeloid expansion, with neuron and astrocyte changes not significant by deconvolution. (b) NF-κB enrichment by region and contrast — C9-specific in grey matter. (c) Severity gradient in the true anterior horn.]({{artifact:art_9a738732-f96a-4b6e-9222-bc6e2857717f}})

---

## 2. Data and cohort

### 2.1 Dataset structure
GSE288365 comprises **20 Visium sections from 18 donors** (two donors contributed paired top/bottom captures of the same block). Each spot carries two modalities in a single matrix: 18,085 gene-expression features and 35 antibody-capture protein features. Protein counts dominate raw totals (~99.7%), so the modalities were split at load: RNA placed in `.X`, protein in `obsm['protein']`. Tissue imaging is **immunofluorescence** (IBA1, MAP2, TDP-43), not H&E.

### 2.2 Cohort composition
| Condition | Donors | Sections |
|-----------|:------:|:--------:|
| Control   | 9      | 10       |
| C9-ALS    | 4      | 4        |
| sALS      | 5      | 6        |

RNA depth is modest, as expected for postmortem cord on a protein-focused assay: median ~850 counts and ~630 genes per spot (range across sections 377–3,240 counts). Section-level QC is in `qc_per_section.csv` and Figure 1.

![Figure 1. Per-section QC: raw vs QC-passing spot counts, median counts, genes, and mitochondrial fraction across all 20 sections.]({{artifact:art_c4cddaa1-8867-4021-806b-d23f082e0c69}})

### 2.3 Companion PBMC cohort (GSE268995)
40 PBMC samples (18 control, 9 sALS-slow, 9 sALS-fast, 4 C9-ALS) profiled by scRNA-seq with 228 imputed surface proteins. **These are unpaired cohorts** — there is no donor overlap with GSE288365 and no crosswalk between the two deposits, so blood/cord comparisons are population-level, not within-donor.

---

## 3. Methods

### 3.1 Preprocessing and QC (stage 01)
Each Visium section was loaded and the RNA/protein modalities separated. RNA-only QC thresholds: minimum 200 counts, minimum 100 genes, maximum 30% mitochondrial reads per spot. Counts were normalized to median library size and log1p-transformed. QC retained 132,781 spots across the 20 sections.

### 3.2 Region annotation — authors' manual histology (stages 18–20)
**All region-stratified results in this report use `auth_manual`, the original authors' manual histological annotation**, with five classes: Anterior_horns, Grey_matter (other grey), White_matter (anterolateral funiculi), Posterior_white_matter (dorsal columns), Meninges.

Our own marker-score region labels (stage 03) are retained in the object for provenance but are **not used for any claim**. Their concordance with the authors' histology is poor (table below, row-percent):

| Our label | Anterior horns | Grey matter | Anterolat. WM | Posterior WM | Meninges | Border |
|---|:--:|:--:|:--:|:--:|:--:|:--:|
| ventral_horn_MN | **15.0** | 5.2 | 33.4 | 15.4 | 20.6 | 10.4 |
| gray_neuropil | 25.4 | 29.7 | 21.9 | 13.8 | 6.1 | 3.1 |
| white_matter | 7.1 | 3.5 | 52.1 | 29.9 | 2.2 | 5.2 |
| dorsal_horn | 8.0 | 8.2 | 37.4 | 19.5 | 17.9 | 9.1 |

### 3.3 Pseudobulk differential expression (stage 20, supersedes stage 05)
To avoid pseudoreplication, differential expression was computed on **donor × author-region pseudobulk**: raw spot counts summed per donor within each author region (minimum 50 spots per pseudobulk, minimum 3 donors per group), modeled with DESeq2 (pyDESeq2 0.5.4). Three contrasts (C9 vs sALS, C9 vs control, sALS vs control) × four regions = 12 models.

### 3.4 Program scores and pathway enrichment (stage 20, supersedes stages 06–08)
Curated gene-set scores were computed as the mean z-score across member genes for: a disease-associated microglia (DAM) program (*GPNMB, CHIT1, TREM2, TYROBP, APOE, ITGAX, CST7, LPL, CD9, CD63, SPP1, CTSB, CTSD, LGALS3, FABP5, MSR1, CD68*), complement (*C1QA, C1QB, C1QC, C3*), homeostatic microglia (*P2RY12, TMEM119, CX3CR1, CSF1R, SALL1*), and a candidate C9 NF-κB program (*BIRC3, SOCS2, CALHM6*). For the DAM-coupling analyses (stage 18) GPNMB and CHIT1 were **excluded from the DAM score** so that correlating them against it is not circular. Pathway enrichment used GSEA preranked (gseapy 1.3.0) on DESeq2 statistics against MSigDB Hallmark.

### 3.5 Cell-type deconvolution — the authors', not ours
Per-spot cell-type estimates are the **authors' Cell2location (C2L) deconvolution**, imported with their metadata. Their method, per the paper: C2L v0.1.3 trained on a human spinal-cord snRNA-seq reference (Yadav et al.) with the motor-neuron population substituted from Gautier et al.; mitochondrial and haemoglobin genes removed; `N_cells_per_location = 7`, `detection_alpha = 20`, 5,000 epochs with centered CDR and UMI as covariates. The exported 5% posterior quantile gives per-spot cell numbers (`auth_<CellType>`, summing to `auth_Total_Cells`, median ~6.3 cells/spot); `auth_norm_ct_<CellType>` expresses each as a **percent of inferred cells in that spot** (sums to exactly 100).

Note that C2L carries **two** microglial populations — `Microglia` and `ProliferatingMicroglia` — and the latter is more than half the microglial compartment; both are included wherever "microglia" is reported. We did **not** run any deconvolution ourselves; our `ct_*` columns are marker-signature scores with no abundance interpretation.

### 3.6 Validation against authors' processed objects (stages 09–10)
The authors' processed Seurat objects were downloaded and their metadata, reductions, and key-gene matrices extracted with a slot-based reader. Per-spot QC and key-gene expression were barcode-matched for one control and one C9-ALS donor. Full author annotations were merged across all 20 sections (98.7% spot match): C2L deconvolution, TDP-43 pathology grade (`annotation_general`), motor-neuron distance, and IBA1/MAP2/TDP-43 imaging features.

### 3.7 Protein layer (stage 20, supersedes stage 11A)
The 35-plex antibody-capture matrix was CLR-normalized and **isotype-control background subtracted** (mean CLR of mouse IgG2a, IgG1κ, IgG2bκ, rat IgG2a). PTPRC appears twice (two CD45 clones) and the duplicates are averaged. Donor-level Mann–Whitney U contrasts were run within each author region.

### 3.8 PBMC arm (stages 12–16)
40 samples QC'd (296,567 cells; median ~1,250 genes, ~2,000 UMIs, mitochondrial ~0.6%). A systematic **barcode-label shift** in the author RDS deposit was identified and corrected: author file labels are cyclically shifted relative to the GEO matrices (`{A1→C1, B1→D1, C1→E1, D1→F1, E1→A1, F1→B1}`); after remapping on barcode content, all 279,856 author-retained cells match. GEO series-matrix conditions are treated as authoritative; all four C9 samples (C2, D3, F3, H5) are concordant. Azimuth L1/L2 labels were merged, counts aggregated to sample × cell-type pseudobulk, and DESeq2 run within cell type.

### 3.9 Reproducibility and technical notes
- OpenMP thread guards set in every Python entry point before importing numpy/scanpy; scikit-learn avoided (OpenMP abort in this environment).
- Memory ceiling 7 GB: obs-level work reads categoricals directly via h5py; matrix operations stream per-sample from backed h5ad.
- Code in `code/analysis/` (stages 01–20); tables in `results/tables/`; figures in `results/figures/`. Repository: github.com/rikku1983/C9orf72-ALS.

---

## 4. Results

### 4.1 Validation: our pipeline reproduces the authors' object

Barcode-matched comparison showed **perfect QC concordance**: per-spot total counts, gene counts, and mitochondrial fraction vs the authors' metrics gave Pearson r = 1.000 and median ratio = 1.000 for both donors tested. This is the decisive check that the RNA/protein modality split — the main technical hazard in this dataset — was done correctly.

Key-gene per-spot correlations (our log-normalized vs authors' SCT) were 0.96–1.00 across all 16 genes tested (GPNMB 0.98, CHIT1 1.00, C1QA 0.96, CHAT 1.00). Spot-set agreement (Jaccard) was 0.84 (control) and 0.99 (C9-ALS). **Region labels were the one step that diverged, and the divergence was larger than this validation suggested** — see §3.2; this is why all regions are now taken from the authors.

![Figure 2. Validation against authors' Seurat objects: per-spot QC concordance, key-gene correlations, and region-label overlap.]({{artifact:art_92faf353-a7a7-4c3b-a5c3-42ed13e21823}})

### 4.2 ALS expands the anterior-horn myeloid compartment; genotype does not

Using the authors' C2L deconvolution in the authors' anterior horns (donor-level Mann–Whitney):

| Cell type | Control % | ALS % | Δ pp | p |
|-----------|:---------:|:-----:|:----:|:---:|
| Proliferating microglia | 1.23 | 3.48 | +2.25 | 0.006 |
| Microglia | 1.41 | 2.49 | +1.08 | 0.013 |
| Macrophages | 1.00 | 2.27 | +1.27 | 0.022 |
| Lymphocytes | 0.78 | 1.20 | +0.42 | 0.052 |
| Oligodendrocytes | 50.67 | 46.17 | −4.49 | 0.158 |
| Neurons | 9.03 | 6.88 | −2.15 | 0.791 |
| Astrocytes | 26.26 | 27.16 | +0.90 | 0.791 |

Myeloid expansion is the robust compositional signature. **Neuron loss and astrogliosis are not significant by deconvolution in the authors' anterior horns** — a correction from the earlier version of this report, which found both significant in the mislabelled region. Neuronal loss is nonetheless clear at gene level: *NEFL* is the single most significant gene in the C9-vs-control anterior-horn contrast (padj = 1.3 × 10⁻¹²). The deconvolution's neuron estimates are the least reliable column, since C2L's reference had to borrow its motor-neuron population from a second dataset.

**No cell type differs between C9-ALS and sALS** in the anterior horn (smallest p = 0.063, astrocytes −6.16 pp; endothelial p = 0.111). Genotype does not change which cells are present.

### 4.3 The C9-specific signal is NF-κB/TNF-α, localized to grey matter

Hallmark GSEA on donor-level pseudobulk, **C9-ALS vs sALS**, by author region (NES; FDR in parentheses):

| Pathway | Anterior horn | Other grey | Anterolat. WM | Posterior WM |
|---|:--:|:--:|:--:|:--:|
| **TNF-α signalling via NF-κB** | **+1.86 (0.014)** | **+2.23 (<0.001)** | +1.25 (0.53) | +1.41 (0.20) |
| **Inflammatory response** | **+1.82 (0.009)** | +1.51 (0.087) | +1.29 (0.53) | +1.10 (0.73) |
| **IL-6/JAK/STAT3** | **+1.73 (0.016)** | +1.33 (0.18) | +1.23 (0.53) | +1.13 (0.67) |
| Interferon-γ response | +1.58 (0.063) | +1.22 (0.27) | +0.96 (0.75) | −1.01 (0.60) |
| Complement | +1.33 (0.18) | −0.92 (0.63) | +0.70 (0.98) | −1.09 (0.50) |

The C9-specific inflammatory signal sits in **grey matter** and disappears in white matter. Its cleanest form is the three-way comparison of NF-κB enrichment:

| Contrast | Anterior horn | Other grey | Anterolat. WM | Posterior WM |
|---|:--:|:--:|:--:|:--:|
| C9 vs control | +2.55 | +2.41 | +2.81 | +2.67 |
| sALS vs control | +2.00 | **−1.09** | +2.64 | +2.42 |
| C9 vs sALS | +1.86 | **+2.23** | +1.25 | +1.41 |

**In other grey matter, C9 engages NF-κB and sALS does not.** That dissociation — rather than a simple quantitative difference — is the strongest C9-specific finding in the dataset.

**No gene survives FDR correction in the C9-vs-sALS anterior-horn contrast** (0 of 11,072; best is *MAN2A1*, padj = 0.059). The signal is pathway-level only, which is the expected ceiling at 4 vs 5 donors. For reference, the same regions yield 219 (C9 vs control) and 325 (sALS vs control) FDR-significant genes, confirming the pipeline detects real effects when they are large.

### 4.4 An ordered condition contrast in the myeloid/lipid axis

Ordinal trend (control < sALS < C9), donor-level Spearman in the authors' anterior horns (n = 18 donors):

| Target | ρ | p |
|--------|:-:|:-:|
| CHIT1 | 0.758 | <0.001 |
| DAM score | 0.724 | 0.001 |
| GPNMB | 0.713 | 0.001 |
| C1QC | 0.688 | 0.002 |
| Complement score | 0.661 | 0.003 |
| C1QA | 0.661 | 0.003 |
| Proliferating microglia | 0.606 | 0.008 |
| ITGAX protein | 0.559 | 0.016 |
| APOE | 0.556 | 0.017 |
| TREM2 | 0.544 | 0.019 |
| Homeostatic score | 0.398 | 0.102 |

The contrast is robust and **now established in the authors' anterior horns**, which the earlier version of this report could not legitimately claim.

> **What 'ordinal trend' does and does not mean here.** The ordering control < sALS < C9
> is an *assumption* encoded in the coding (0/1/2), not a measured severity variable.
> **This cohort records no survival, disease duration, age at onset, ALSFRS or pathological
> stage** — only age, sex and subdiagnosis. Decomposing the trend shows it is almost entirely
> the control-versus-ALS step: for M5 in anterolateral WM the ordinal coding gives rho = 0.754
> versus rho = 0.739 for a plain binary ALS/control indicator (difference +0.015), and the
> sALS -> C9 step alone is not significant (rho = 0.346, p = 0.36). These results should be
> read as *disease versus control*, not as a dose-response in severity. The homeostatic-microglia score does not track it, consistent with a shift from homeostatic to activated myeloid states rather than a uniform increase in myeloid activity.

### 4.5 The DAM/GPNMB/CHIT1 program is targeted to motor tracts (stage 18)

Mapping DAM core, GPNMB and CHIT1 on the authors' histology shows the program peaks in **anterolateral white matter** — the compartment containing the corticospinal tracts — exceeding every other region in 9/9 ALS donors (p = 0.004 for most comparisons; GPNMB vs posterior WM 8/9, p = 0.020).

The **dorsal columns are spared**: GPNMB and CHIT1 there are not elevated over control (p = 0.38, 0.33), whereas in anterolateral WM they are (p = 0.010, 0.013). Anterior horns are also elevated over control, at lower magnitude (DAM p = 0.003, GPNMB p = 0.005, CHIT1 p = 0.002). This sensory-sparing, motor-tract-selective pattern is what corticospinal degeneration predicts.

GPNMB co-expresses with the complement and myeloid programs in ALS spots, most strongly in anterolateral WM (ρ = +0.254 complement, +0.189 myeloid) and more weakly in the anterior horn (+0.165, +0.078). Within-section GPNMB–DAM coupling rises with disease (median ρ 0.08 control → 0.13 sALS → 0.23 C9; ALS vs control p = 0.045), and across donors in anterolateral WM reaches ρ = 0.86 (0.82 within ALS donors alone, p = 0.007).

**Critically, none of DAM, GPNMB or CHIT1 differs between C9 and sALS in any region** (p = 0.19–0.91). This program is a shared ALS axis, not a genotype discriminator.

![Figure 3. Spatial maps of histological region, DAM core score, GPNMB and CHIT1 in one representative section per condition (selected as the section closest to its group's median DAM core).]({{artifact:art_2f4b5b40-f48a-4119-bf5a-840822f657df}})

![Figure 4. Donor-level quantification by author region: DAM core, GPNMB, CHIT1, macrophages, microglia, total myeloid, GPNMB–DAM coupling per section, and the donor-level GPNMB vs DAM relationship in anterolateral white matter.]({{artifact:art_869b18dc-fafd-4402-9a04-ab199d574877}})

### 4.6 TDP-43 pathology does not couple to the myeloid program (stage 17)

Testing the authors' per-spot TDP-43 pathology annotation (Ring / adjacent / distant / Not_selected) against the myeloid, complement, DAM and NF-κB programs produced an **honest negative**. At spot level a striking genotype divergence appeared — sALS apparently recruiting programs to pathology (+0.12 to +0.19) while C9 did not — with p < 0.001. That was **pseudoreplication**: only 3 donors per genotype carry both pathology-bearing and unaffected spots, and donor-level deltas are sign-inconsistent (C9 myeloid: −0.26, −0.12, +0.03; sALS: −0.20, −0.09, +0.06). No reproducible coupling survives.

Two caveats bound this test: `auth_TDP43_mean` is an immunofluorescence intensity whose disease direction is ambiguous (TDP-43 is normally nuclear and is cleared in pathology), and pathology-bearing spots are genuinely sparse (417 Ring, 472 adjacent across all ALS).

![Figure 5. TDP-43 pathology coupling: spot-level (misleading) versus donor-level (honest) comparison.]({{artifact:art_79cf0856-9a4e-4b82-9583-238d3cc9b8b6}})

### 4.7 Protein layer confirms myeloid activation

In the authors' anterior horns, **ITGAX/CD11c is elevated in ALS vs control (Δ +0.80, p = 0.013)**; CD68 (Δ +0.43, p = 0.133) and HLA-DRA (Δ +0.57, p = 0.185) trend in the same direction but are not significant. ITGAX protein also tracks the ordered condition contrast (ρ = 0.559, p = 0.016). *This is a correction: the earlier version reported CD68 as significant (p = 0.013) in the mislabelled region; on the authors' anterior horns ITGAX is the significant marker and CD68 is not.*

No protein marker separates C9 from sALS (smallest p = 0.111, ITGAM).

**Limitation:** the antibody panel does not include GPNMB, CHIT1, TREM2 or APOE, so it cannot adjudicate DAM-program composition at protein level. It validates *that* myeloid cells are activated, not *which* program they run.

### 4.8 Peripheral blood: a weak myeloid-activation echo, no composition shift

In the unpaired PBMC cohort, **no L1 or L2 cell-type proportion differs by condition after FDR correction** (4 C9 donors; largest raw effect CD8 T +7 pp in C9 vs sALS, p = 0.39). Monocyte gene-level DE yields nothing surviving FDR.

A directionally consistent but non-significant myeloid-activation trend appears in C9 monocytes: **CD68 surface protein up (p = 0.033 vs sALS, FDR 0.43)** with CD9, CD63 and ITGAX mRNA up (CD9 p = 0.033–0.042). Importantly, the **CNS DAM lipid core (GPNMB, APOE, TREM2, SPP1) is not expressed in peripheral monocytes** — it filtered out as too low — so the peripheral correlate is the activation-marker arm only, not the full DAM program. Complement does not go peripheral.

Because the cohorts are unpaired and the compartments differ, this is a population-level convergence, not evidence that the same individuals show both.

![Figure 6. Reconciliation of bulk-style and spatially-resolved myeloid estimates: pseudo-bulk collapse, region-resolved contributions, and the regional composition of the cord.]({{artifact:art_df4208f9-2d86-453c-bd6d-0ab736ccebdf}})

---

## 5. Synthesis: a model of the C9 mechanism

1. **ALS converges on a shared cellular endpoint** in the anterior horn: myeloid expansion (microglia, proliferating microglia, macrophages) with neuronal transcript loss. Genotype is invisible at the level of cell-type abundance.

2. **A shared disease axis separates ALS from control** in the myeloid/lipid/complement program (DAM, GPNMB, CHIT1, C1QA/C1QC, APOE), anatomically targeted to **anterolateral white matter** with sparing of the dorsal columns — the signature of motor-tract degeneration and debris clearance. This axis does **not** distinguish C9 from sALS, and because no severity variable was recorded we cannot say whether it scales with disease severity — only that it separates ALS from control.

3. **Genotype acts on inflammatory signalling in grey matter.** The repeat expansion is associated with NF-κB/TNF-α, inflammatory-response and IL-6/JAK/STAT3 enrichment in the anterior horn and other grey matter — and in other grey matter this is a true dissociation, with C9 engaging NF-κB where sALS does not. This is the leading candidate mechanism distinguishing C9-ALS.

4. **Two compartments, two processes.** The myeloid/DAM abundance signal (white matter) and the C9-specific inflammatory signalling (grey matter) are anatomically distinct. The earlier version of this report blended them, because its region labels mixed the two.

The pending bulk RNA-seq cohort is the natural validation: does NF-κB/inflammatory signalling separate C9 from sALS in grey-matter-containing cord tissue with more donors? Preliminary comparison with the user's bulk deconvolution already agrees on the shared finding — myeloid expansion is an ALS feature that is not demonstrably greater in C9 than sALS.

---

## 6. Limitations

- **No severity or survival data.** The cohort records only age, sex and subdiagnosis — no survival, disease duration, onset age, ALSFRS or stage. The control → sALS → C9 ordering is an assumption, not a measurement, and the apparent trend is driven by the control-versus-ALS step (see §4.4). Visium also captures only an end-stage postmortem snapshot.
- **Confounds.** C9 donors with recorded age are younger (51, 55) than the control mean (~63); two sALS donors carry co-pathology (Alzheimer's, FTD).
- **Genotype power.** 4 C9-ALS vs 5 sALS donors. Pathway-level and gradient results are robust; **no single gene survives FDR in the C9-vs-sALS anterior-horn contrast**.
- **Spatial resolution.** Visium spots (~55 µm) are multi-cellular; cell-intrinsic claims are inferred via deconvolution and program scoring. Macrophages and microglia are transcriptionally close and cannot be cleanly separated at this resolution.
- **Deconvolution is prior-anchored.** C2L was run with `N_cells_per_location = 7`; absolute per-spot cell counts reflect that prior, so only relative comparisons are sound. Motor-neuron estimates are the least reliable column.
- **Unpaired cohorts.** GSE288365 (cord) and GSE268995 (blood) share no donors; cross-compartment statements are population-level.
- **Cord level unverified.** Section identifiers hint at cervical for some donors but spinal level is not a recorded field, so level-matching to external cohorts is unverified.
- **Protein panel** lacks the key DAM markers.
- **Region labels** are the authors' manual histology; our marker-score labels proved unusable (§3.2) and all claims derived from them have been withdrawn.

---

## 7. Claims withdrawn in the session-9 revision

For transparency, the following statements from the earlier version did **not** survive recomputation on the authors' histology:

| Withdrawn claim | Old value | On author anatomy |
|---|---|---|
| Unfolded protein response separates C9 from sALS | NES 1.88, FDR <0.001 | NES 1.23, **FDR 0.21** (anterior horn) |
| mTORC1 separates C9 from sALS | NES 1.88, FDR <0.001 | NES 1.22, **FDR 0.21** (anterior horn) |
| *CANX* up in C9 ventral horn | padj 7 × 10⁻⁴ | **padj 0.16** (anterior horn) |
| *CALHM6 / BIRC3 / SOCS2* are C9-specific markers | p = 1.8 × 10⁻⁴ / 0.023 / 0.030 | CALHM6 **p = 0.47**; BIRC3, SOCS2 below expression filter |
| *ITGAX* down in C9 ventral horn | p = 3.7 × 10⁻⁵ | p = 0.055 |
| Neuron loss in ventral horn (deconvolution) | −2.57 pp, p = 0.004 | −2.15 pp, **p = 0.79** |
| Astrogliosis in ventral horn | +6.21 pp, p = 0.005 | +0.90 pp, **p = 0.79** |
| CD68 protein up in ALS ventral horn | Δ +0.56, p = 0.013 | Δ +0.43, **p = 0.133** (ITGAX is the significant marker instead) |
| "Two routes": C9 → NF-κB, sALS → complement | NES 2.30 vs 1.96 | **Not supported.** Complement is the *top* C9-vs-control pathway in the anterior horn (NES 2.73) and is enriched in both genotypes |
| C9-vs-sALS NF-κB "NES up to 2.30" | 2.30 | Maximum \|NES\| for this contrast is **2.17**; 2.53 belongs to C9 vs control |
| Divergent myeloid programs at TDP-43 lesions | directionally consistent hypothesis | **Negative** — spot-level effect was pseudoreplication (§4.6) |

Claims that **did** survive: myeloid expansion in ALS, absence of any compositional genotype difference, the control → sALS → C9 ordered contrast (now correctly described as a control-vs-ALS effect), the NF-κB/TNF-α C9-vs-sALS enrichment in grey matter, and the validation against the authors' objects.

---

## 8. Analysis inventory

**Pipeline (`code/analysis/`):** `01_visium_qc.py` · `02_integrate_cluster.py` · `03_annotate_regions.py` *(superseded for region use)* · `05_pseudobulk_de.py` *(superseded by 20)* · `06_scores_by_condition.py` *(superseded by 20)* · `07_pathway_enrichment.py` *(superseded by 20)* · `08_microglia_deepdive.py` *(superseded by 18/20)* · `09_extract_seurat.R` / `09b_extract_meta_only.R` · `10_integrate_author_meta.py` · `11_protein_and_stats.py` *(superseded by 20)* · `12_scrna_qc.py` · `13_extract_scrna_author_meta.R` · `14_merge_author_meta.py` · `15_pseudobulk_aggregate.py` · `16_myeloid_markers_ADT.py` · `17_tdp43_program_coupling.py` · `18_spatial_DAM_GPNMB_CHIT1.py` · `19_pseudobulk_authoranat.py` / `19b_gsea_authoranat.py` · `20_recompute_authoranat.py`.

**Current result tables (`results/tables/`):** `region_annotation_disagreement.csv` · `de_authoranat_*` (12 contrasts) + `de_authoranat_summary.csv` · `gsea_authoranat_all.csv` · `composition_authoranat_allcontrasts.csv` · `protein_authoranat_allcontrasts.csv` · `programs_authoranat_allcontrasts.csv` · `trend_authoranat.csv` · `programscores_authoranat_donor_region.csv` · `spatial_DAMcore_GPNMB_CHIT1_*.csv` · `GSE288365_TDP43_coupling_*.csv` · `C9_vs_sALS_evidence_ranked.csv` · `validation_vs_authors_metrics.csv`, `validation_gene_correlations.csv` · PBMC: `GSE268995_*`.

**Superseded tables** (retained for provenance, not cited): `pseudobulk_*_ventral_horn_MN.csv`, `gsea_*_ventral_horn_MN_*.csv`, `scores_donor_region.csv`, `microglia_*.csv`, `composition_*_VH.csv`, `protein_*_VH.csv`, `jonckheere_trend_VH.csv`, `mixedmodel_gradient_VH.csv`.

**Next step:** integrate the user's bulk spinal-cord RNA-seq DEGs to test whether NF-κB/inflammatory signalling separates C9 from sALS in an independent, better-powered cohort.

---

## 9. Phase II: mapping the bulk causal-mediator genes onto cell types and regions

**Input.** 118 genes nominated by an independent bulk-RNA causal-mediation analysis as
potential mediators of the C9orf72 → shorter-survival effect (strong C9 → gene and
gene → survival paths). The bulk data itself was not available to this analysis, so
Phase II **characterises** the gene list in the spatial cohort rather than re-testing
the causal model.

### 9.1 Panel coverage

108 of 118 genes are present in the Visium probe panel; 104 reach ≥1% of spots.
**All nine ribosomal genes (RPL23, RPL24P4, RPL26, RPL27, RPL7A, RPS10, RPS15A, RPS17)
plus GPX1 and AC006011.2 are absent from the panel by design** — the probe set contains
zero ribosomal protein genes. They therefore cannot be evaluated here; whether they
represent an RNA-integrity axis in the bulk data remains an open question for that cohort.
ALOX15B, CHRNA1, DPEP2 and OTOA are in the panel but below detection.

### 9.2 The list is five modules, not 118 independent mediators

Co-expression clustering on donor x region pseudobulk (region-centred, average linkage)
resolves five modules, each attributable to a cell type by correlation with the authors'
Cell2location abundances:

| Module | n | Cell-type correlation | Peak region |
|---|:-:|---|---|
| M1 Oligodendrocyte/myelin | 19 | Oligodendrocytes +0.43 | Anterolateral WM |
| M2 Astrocyte | 25 | Astrocytes +0.29 | Anterolateral WM |
| M3 Mixed, oligodendrocyte-leaning | 8 | no coherent type (6 types across 8 genes) | Anterolateral WM |
| M4 Endothelial/vascular (sparse) | 4 | Endothelial +0.43 (4/4 genes, but all sparsely detected) | Anterior horns |
| **M5 Myeloid/DAM** | **48** | **Prolif. microglia +0.60, macrophages +0.57, microglia +0.56, oligodendrocytes −0.44** | **Anterolateral WM** |

M5 contains 46% of the analysable genes and is by far the most internally coherent
(TREM2, TYROBP, CD68, GPNMB, CHIT1, LGALS3, CD63, APOC1, C1QC, CTSS, CTSD, GRN, LAPTM5,
ITGB2, CD53, SPI1, SYK, HCK, NCF2, CYBA, PYCARD, C3AR1, IFI30 and others). Gene-level
mediation across such a collinear block repeatedly re-tests one latent variable; the
effective number of independent mediators in this list is approximately five.

**Module annotation was subsequently audited** (coherence metrics + independent CellGuide marker
overlap; `mediator_module_annotation_summary.csv`). M1 (high confidence: 17/19 genes oligodendrocyte,
within-module r = 0.42, canonical markers CARNS1/OPALIN p = 0.032) and M5 (high confidence: all 48 genes
immune, r = 0.52, macrophage markers C1QC/CYBA/LAPTM5/TYROBP p = 0.042) are solid. **M3 was renamed from
"Neuronal" — it is not a neuronal module**: only 2/8 genes best-match neurons, its strongest external
evidence is oligodendrocyte (EDIL3, SLC24A2; p = 0.005), and it fails a direct sanity check — CNTN1,
SLC24A2 and MAPK8IP2 are flat in ALS anterior horn (p = 0.72-1.00) while NEFL falls sharply in the same
spots (p = 0.001). M4 is unanimous for endothelium but rests on four sparsely detected, functionally
unrelated genes. M2 (astrocyte) is moderate: half its genes map to the Cell2location 'ependymal' column,
a very small population likely absorbing astrocyte signal, though the set contains textbook astrocyte
genes (GJA1, NDRG2, HEPACAM, GPR37L1, DAAM2, TMEM47).

Note on the suspected muscle-contamination subset: in spinal cord, **CARNS1 and LDB3
track oligodendrocyte abundance** (r = 0.70, 0.73) and PYGM tracks astrocytes (r = 0.30),
so they behave as CNS-expressed genes here rather than as skeletal-muscle carryover.
CHRNA1 is below detection.

![Figure 7. Mediator-gene module structure. (a) Pseudobulk co-expression correlation, hierarchically ordered, with the five modules outlined. (b) Module identity assigned by correlation with Cell2location cell-type abundance.]({{artifact:art_86ca08cf-2624-487c-bf8c-5c43e1b24637}})

### 9.3 Regional localisation

Four of five modules peak in anterolateral white matter. M5 has by far the strongest
regional preference — higher than the anterior horn, other grey matter and posterior
white matter in **9/9 ALS donors** (all p = 0.004) — matching the DAM/GPNMB/CHIT1
localisation established in §4.5. M1, M2 and M3 also favour anterolateral WM but with
much smaller margins, and M4 (vascular) shows no significant regional preference.

![Figure 8. Module scores by donor across the four author histological regions.]({{artifact:art_5f6a240c-9a14-4cc5-a770-b9b0eebdcbe0}})

![Figure 9. Spatial maps of each module score in one representative section per condition.]({{artifact:art_6aae6dac-c194-4249-b9ce-23d6b93de128}})

### 9.4 Only the myeloid module tracks disease — and it is not C9-specific

Donor-level contrasts per module per region:

- **ALS vs control:** only **M5** is significant, in the anterior horn (p = 0.005),
  other grey matter (p = 0.013) and anterolateral WM (p = 0.003). M1–M4 reach
  significance nowhere (smallest p = 0.077, M1 in anterolateral WM).
- **Severity trend (control < sALS < C9):** only **M5**, and in all four regions
  (rho = 0.48–0.75, p = 0.041 to <0.001).
- **C9-ALS vs sALS: no module is significant in any region** (smallest p = 0.063,
  M4 in other grey matter).

![Figure 10. Module-level condition contrasts across regions; asterisks mark p < 0.05.]({{artifact:art_3de3a96c-64d2-4b23-abe2-64e0ffaf5959}})

### 9.5 Peripheral detectability

Of the mediator genes present in the PBMC monocyte pseudobulk, 77 are detectable in blood.
Seven M5 members are **CNS-restricted** — **TREM2, GPNMB, CHIT1, APOC1, PLAU, ADORA3,
SIGLEC8** — consistent with §4.8: the lipid/DAM core of the myeloid program does not
appear in circulating monocytes, so these specific mediators are not accessible as blood
biomarkers.

### 9.6 Interpretation and caveats

The dominant mediator module is, gene for gene, the myeloid/DAM program this project
characterised independently from the spatial data — an encouraging convergence between
two unrelated analyses.

Three caveats bound what Phase II can conclude:

1. **Direction of causation is not testable here — and neither is severity.** If the bulk
   expression is postmortem, it is measured *after* survival is determined, so "gene shortens
   survival" and "gene marks the terminal state of an aggressive course" are observationally
   equivalent. **This spatial cohort cannot adjudicate it: no survival, duration, onset or
   functional-score data were recorded.** What we can say is narrower than first stated — M5 is
   elevated in ALS versus control and shows **no C9-versus-sALS difference in any region**.
   An earlier draft described M5 as scaling with "severity"; that wording was unsupported,
   since the control < sALS < C9 ordering is an assumption rather than a measurement
   (see the note in section 4.4).
2. **Collinearity inflates the apparent mediator count.** Five modules, not 118 mediators.
3. **Composition versus regulation is unresolved.** M5 correlates strongly with
   deconvolved myeloid abundance, so much of its bulk signal plausibly reflects *how much*
   myeloid infiltration a cord contains rather than per-cell regulatory change. Visium
   spots are multicellular and cannot settle this; single-nucleus data would.

A testable prediction follows: if these genes mediate a **C9-specific** effect, the
C9-versus-sALS contrast should be non-null in a better-powered cohort. In this cohort
(4 C9 vs 5 sALS donors) it is null for every module, whereas the shared ALS-versus-control
and severity-trend contrasts are strongly positive. The spatial data therefore fail to support M5 as a **genotype-specific** mediator; whether it is a
severity mediator is untestable here for want of any severity variable — leaving the
NF-kB/grey-matter axis of section 4.3 as the better candidate for C9-specific mediation.


### 9.7 Spatial context: motor-neuron proximity, MAP2 imaging and TDP-43 (stage 22)

Three author annotations were used, with designs chosen to escape the 4-vs-5-donor ceiling.
Motor-neuron distance and MAP2 are **within-donor** gradients: each donor contributes its own
slope or correlation, and those are then tested across all donors.

**Confound control.** 91% of near-motor-neuron spots lie in the anterior horn, so the distance
gradient was computed **within the anterior horn only** (16 of 18 donors have >=15 spots in >=3 of
the four distance bins: 9 control, 4 sALS, 3 C9). Without this restriction the gradient would
simply re-measure the grey/white axis. MAP2 passes its own sanity check (mean 0.89 in anterior horn
versus 0.48 in white matter).

**M1 oligodendrocyte gradient differs between ALS and control** — the strongest well-powered
result in Phase II. In controls the myelin module *rises* with distance from motor neurons
(mean slope +0.017 per bin); in ALS that gradient is flat or reversed (sALS -0.006, C9 +0.006);
ALS versus control p = 0.0079 (BH q = 0.039 across the five module tests). All four sALS donors
have slopes at or below zero, versus 9/9 controls positive. The simplest reading is loss of the
normal myelin enrichment in peri-anterior-horn white matter, consistent with the overall depletion
of M1 in ALS white matter (section 9.4).

**M2 tracks neuronal content; M3 does not.** Per-donor Spearman against MAP2 immunofluorescence
within the anterior horn: M2 rho = +0.127 (14/18 donors positive, p = 0.0028, **BH q = 0.042** —
the only test surviving correction across all 15), consistent with astrocytes occupying neuropil
alongside neurons. **M3 shows no relationship with MAP2 (rho = +0.010, p = 0.50)** — independent
protein-imaging confirmation, orthogonal to the transcriptome, that M3 is not a neuronal module
and that its renaming (section 9.2) is correct.

**The myeloid module shows no motor-neuron proximity gradient** (M5 slope p = 0.21; MAP2
rho = -0.017, p = 1.00). Combined with its anterolateral white-matter peak, this locates the
myeloid program to degenerating tracts rather than to the immediate motor-neuron
microenvironment — an argument against it being a proximate mediator of motor-neuron death.

**TDP-43 pathology: no detectable effect for any module.** The contrast is Ring + TDP43_adjacent
spots (pathology) versus the authors' **TDP43_distant** spots, which serve as a spatially matched
internal control; per-donor deltas were tested by Wilcoxon signed-rank. Mean deltas -0.021 to
+0.036, all p >= 0.156 (M1 1.00, M3 1.00, M4 0.22, M5 0.44, M2 0.16).

*Design note.* An initial version of this test compared all three selected categories (including
TDP43_distant) against Not_selected. That comparison is regionally confounded: Not_selected spots
are 67% white matter whereas Ring/adjacent/distant are ~70-75% anterior horn, so it substantially
measured the grey/white axis rather than pathology. It produced spuriously larger negative deltas
(e.g. M5 -0.094 versus -0.016 in the matched contrast). The matched-control result above supersedes it.

*How strong is this negative?* Weaker than it may appear. Only 6 of 9 ALS donors carry any Ring or
adjacent spots at all (C5_GW24, S2_GW37 and S3_GW18 have none), and two of the six contribute fewer
than 25 pathology spots. With n = 6 the smallest achievable two-sided Wilcoxon p is 0.031, so this
is an underpowered test that detected nothing, not a demonstration that no effect exists. It is
consistent with the stage-17 negative on the original programs, which hit the same ceiling.

![Figure 11. Spatial-context tests. (a) M1 across motor-neuron distance bins, per donor with condition means. (b) Per-donor M1 slopes by condition. (c) Per-donor correlation of each module with MAP2 neuronal immunofluorescence. (d) Motor-neuron distance slopes for all modules.]({{artifact:art_5b828a47-c654-450e-ae4d-2f1b5cfe4a3f}})

**Tables:** `mediator_MNdistance_gradient.csv` - `mediator_MAP2_correlation.csv` -
`mediator_TDP43_donorlevel.csv` (initial, confounded) · `mediator_TDP43_corrected.csv` (matched-control, authoritative). **Code:** `code/analysis/22_mediator_spatial_context.py`.


### 9.8 Composition adjustment and differential co-expression (stage 23)

Two further tests of whether a C9-versus-sALS effect is hidden by the mean-comparison design.

**Composition adjustment.** Each module score was residualised **within donor** against the
Cell2location abundance of its matched cell type (M1 oligodendrocytes, M2 astrocytes + ependymal,
M3 oligodendrocytes + neurons, M4 endothelium, M5 total myeloid), then donor-level contrasts were
repeated on the residuals. **No C9-versus-sALS contrast becomes significant after adjustment**
(all p >= 0.063; smallest is M4 in other grey matter, unchanged at 0.063). Composition was not
masking a genotype effect.

The adjustment is informative about the *disease* effect, however. For M5, ALS-versus-control
survives adjustment in anterolateral white matter (p 0.003 -> 0.004) but is abolished in the
anterior horn (0.005 -> 0.133) and other grey matter (0.013 -> 0.064). Read directly: in grey
matter the myeloid signal is **compositional** -- more myeloid cells, similar per-cell expression --
whereas in anterolateral white matter a **cell-intrinsic** component remains, consistent with a
genuinely activated myeloid state in degenerating tracts rather than simple infiltration.

**Differential co-expression failed as a method.** Per-donor module coherence (mean within-module
gene-gene correlation across that donor's spots) was intended as a better-powered within-donor
statistic. It is not interpretable here: coherence correlates strongly with sequencing depth
(Spearman rho = +0.73 to +0.82 across modules in white matter), so the statistic largely reports
how deeply a donor was sequenced. All condition contrasts were null (C9 vs sALS p >= 0.19;
ALS vs control p >= 0.052) but no weight should be placed on that given the confound.

**Multiple testing.** Across all 60 tests in this stage nothing survives Benjamini-Hochberg
correction (smallest q = 0.107, for the two anterolateral-white-matter composition-adjusted
ALS-versus-control results). The composition/intrinsic distinction above should be read as a
coherent pattern worth reporting, not an established result.

![Figure 12. (a) M5 ALS-versus-control significance before and after adjusting for myeloid abundance, by region. (b) C9-versus-sALS p-values, raw versus composition-adjusted, all modules and regions. (c) Module coherence against per-donor sequencing depth, showing the confound that invalidates the co-expression test.]({{artifact:art_a064754c-ec4d-4d91-b3b4-a402c5ec948e}})

**Tables:** `mediator_composition_adjusted.csv` - `mediator_coherence_bydonor.csv` -
`mediator_coherence_perdonor_raw.csv`.

**Phase II tables:** `mediator_genes_annotated.csv` (per gene: module, cell type, peak
region, module contrast statistics, blood detectability, panel/detection status) ·
`mediator_modules.csv` · `mediator_celltype_correlation.csv` · `mediator_region_profile.csv` ·
`mediator_module_contrasts.csv`. **Code:** `code/analysis/21_mediator_modules.py`.
