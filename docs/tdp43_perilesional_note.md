# TDP-43 perilesional myeloid halo — parked finding

**Status:** complete, reproducible, **not used in the SfN talk** (decision Nov 2026).
Picked up again? Start here, then run `code/analysis/23_tdp43_perilesional.py`.

## What changed relative to the earlier negative

Stages 17 and 22 reported TDP-43 coupling as an honest negative. Both tested a
**monotone** relationship over the whole anterior horn — a grade ladder in 17, a
proximal-vs-matched-distant contrast in 22. Neither can detect a peak at an
intermediate radius, and both were dominated by far-field spots (median distance
from a lesion 361–1,101 µm depending on section).

The negative stands for what it tested. It was the wrong shape of test.

## The annotation is geometry, not severity

`auth_tdp43_a2_d2` levels are concentric dilation rings around manually marked
TDP-43⁺ spots — identical to the micron in every section, zero variance:

| level | radius | spots (cohort) |
|---|--:|--:|
| `TDP43_positive` | 0 µm (the spot itself) | 90 |
| `TDP43_adjacent` | 112 µm (ring 1) | 382 |
| `Ring` | 200 µm (ring 2) | 417 |
| `TDP43_distant` | 330 µm (ring 3) | 891 |
| `Not_selected` | ~2,600 µm median | 129,310 |

`auth_annotation_general` merges positive + adjacent into a single
`TDP43_adjacent` level (472 spots) = the paper's binary "TDP-43+" class.

**`Ring` is a dilation ring, not a severity tier.** Stage 17 encoded the grade
ladder as `Not_selected < TDP43_distant < TDP43_adjacent < Ring`, treating Ring as
most severe. That ordering is wrong; Ring sits *between* adjacent and distant.
Stage 17's conclusion is unaffected (it rested on donor-level sign inconsistency)
but its docstring and REPORT §4.6 wording should be corrected if revisited.

## The finding

Anterior horn, 4 lesion-richest sections, values centred within section:

| distance | neurons | **myeloid** | oligodendrocytes | astrocytes |
|---|--:|--:|--:|--:|
| 0 (TDP-43⁺ spot) | **+6.37** | −1.24 | −7.65 | +2.76 |
| 2 (200 µm) | −1.23 | +0.60 | −1.77 | +1.18 |
| **3 (300 µm)** | −1.99 | **+2.03** | −0.88 | +0.64 |
| 4 (400 µm) | −0.65 | +1.12 | −0.13 | −0.73 |
| 9–11 | +0.93 | **−1.61** | +1.39 | −0.90 |

The Myeloid/DAM score peaks at ring 3 (+0.030) and GPNMB at ring 4 (+0.143).

**Specificity is the load-bearing evidence.** At ring 3 neurons fall 1.99 pp and
myeloid gains 2.03 pp — essentially the whole displaced fraction — while
oligodendrocytes still *lose* 0.88 pp and astrocytes gain only 0.64 pp. Generic
neuron displacement would lift every non-neuronal type proportionally. It doesn't.

Per section: S1_GW50 +4.31 pp, C4_N280 +2.70, C3_AU87 +1.90, **S4_GW35 −4.33**.
S4_GW35 has only 33% of its lesion spots in the anterior horn (63% posterior
horn), so its anterior-horn distance field is anchored on few unrepresentative
centres. Report as 3 of 4, not 4 of 4.

## Caveats to carry

- Spot-level and pseudoreplicated. The evidence is **section consistency**, not
  the p-values, which are inflated by spatial autocorrelation.
- Peak at 300–400 µm is **further than microglial process reach** (~50–100 µm):
  a local reactive field, not direct perineuronal ensheathment.
- The rightmost bin (17+ rings) rises for **both** myeloid and oligodendrocytes —
  that is the anterior-horn periphery abutting white matter, an anatomical edge
  effect. Not a second peak.
- Ring labels are nearest-lesion, not unique-lesion; shells of neighbouring
  lesions merge in dense clusters (C3_AU87, S1_GW50).
- Only 90 TDP-43⁺ spots cohort-wide. This constrains effect size; it does not
  establish absence of other structure.

## If resumed, next steps worth taking

1. **Donor-level test.** Collapse each section to one halo statistic (e.g.
   rings 2–4 minus rings 7–11) and test across donors. Removes pseudoreplication;
   n = 4 sections currently, so this needs more lesion-annotated tissue.
2. **Lesion-level rather than spot-level.** Treat each connected TDP-43⁺ cluster
   as the unit, build a radial profile per lesion, then average. Avoids the
   merged-shell problem.
3. **Per-genotype halo amplitude.** C9 vs sALS on halo height is the question the
   project actually cares about; currently 2 sections per genotype.
4. **Check against IBA1 immunofluorescence**, which is in the imaging channels
   (`auth_cleaned_IBA1_binary`, `auth_IBA1_positive`) and is an independent
   measurement of myeloid presence at the same coordinates.
5. Re-examine whether the authors' Fig. 6i result (inflammatory genes up in
   TDP-43-proximal neuron-enriched spots) is the rising limb of this halo; their
   contrast pooled rings 0–1 against everything else, mixing the myeloid-poor
   lesion spot with the start of the rise.

## Files

**Script:** `code/analysis/23_tdp43_perilesional.py` (self-contained; regenerates
every table and the figure from the two h5ad files plus the author obs CSV).

**Tables:** `tdp43_perilesional_profile.csv` (bin means/SEM/n, 4 cell types +
2 expression features) · `tdp43_halo_persection.csv` (per-section rings 2–4 vs
7–11) · `tdp43_rings03_gradient.csv` (authors' categorical levels, ordinal ρ) ·
`tdp43_grade_geometry.csv` (the dilation radii) ·
`tdp43_lesion_vs_unaffected_anteriorhorn.csv`.

**Figures:** `tdp43_perilesional_profile.png` (the radial profile) ·
`tdp43_lesion_maps.png` (where the lesions are, rings coloured, 4 sections) ·
`tdp43_ring_gradient.png` (the earlier wide-bin version — superseded, kept for
the record of why the monotone test missed this).
