# C9orf72-ALS: Single-cell & Spatial Transcriptomics of ALS Spinal Cord

Decoding the mechanism of the *C9orf72* hexanucleotide repeat expansion in ALS
progression, extending bulk RNA-seq findings with single-nucleus and spatial
transcriptomics of **human spinal cord** (C9-ALS, sporadic ALS, and controls).

## Project layout

```
C9orf72-ALS/
├── code/
│   ├── download/     # scripts that fetch public data (GEO, etc.) into data/raw/
│   ├── analysis/     # analysis scripts / notebooks (numbered by stage)
│   └── utils/        # shared helper functions
├── data/
│   ├── raw/          # downloaded data — git-ignored, reproduced via code/download/
│   ├── processed/    # intermediate objects — git-ignored
│   └── metadata/     # small sample sheets / dataset inventory — COMMITTED
├── results/
│   ├── figures/      # committed figures
│   └── tables/       # committed result tables
├── env/              # environment specs (conda / renv) for reproducibility
└── docs/             # notes, methods, data-availability tracking
```

## Reproducibility model

Raw data is **not** stored in git (files are multi-GB; GitHub caps at 100 MB).
Instead, everything needed to regenerate it is committed:

1. Create the environment — see `env/`.
2. Run the download scripts in `code/download/` to populate `data/raw/`.
3. Run the analysis scripts in `code/analysis/` (numbered in execution order).

## Datasets

See `data/metadata/c9als_spinalcord_dataset_inventory.csv` for the full,
verified inventory. Headline sources:

| Accession | Modality | Tissue | Design |
|-----------|----------|--------|--------|
| **GSE288365** | 10x Visium spatial | Human cervical spinal cord | 9 control + 4 sporadic + 5 C9orf72 (paired RNA + protein) |
| GSE272626 | Bulk RNA-seq | Human spinal cord (NYGC ALS Consortium) | 266 samples, ALS / ALS+FTD / CTL — validation cohort |
| GSE271156 | snRNA-seq | Human motor cortex | sporadic ALS / control (cortical comparison) |
| bioRxiv 2025.08.26.672029 | snRNA-seq | Human spinal cord | C9 / SOD1 / sALS — **data embargoed until publication** |

## Status

- [x] Dataset landscape surveyed (spinal cord focus)
- [ ] GSE288365 Visium downloaded & QC'd
- [ ] Bulk DEG list integrated (user-provided)
- [ ] Spatial enrichment of bulk DEGs by anatomical region
