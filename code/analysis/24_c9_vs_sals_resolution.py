#!/usr/bin/env python
"""
Stage 24 - Resolving the C9-ALS vs sALS direction of the 118 bulk mediator genes.

Session-11 analyses, consolidated and made reproducible. Five questions:

  A. Three-group resolution. Is the apparent conflict with the bulk causal model
     (TREM2/APOE "driven up by C9" there, lower in C9 here) an artefact of the
     comparator?  ->  key_genes_three_way.csv
  B. Can Cell2location abundance be regressed out of expression?  Tests state
     leakage into the deconvolution and the collinearity of the proposed
     covariate, plus an attempted IBA1-immunofluorescence validation.
     ->  deconvolution_state_confound.csv
  C. Why does the per-spot mean-z score disagree with pseudobulk DESeq2 on
     direction in the posterior regions?  ->  score_vs_pseudobulk_decomposition.csv
  D. Split the Myeloid/DAM cluster by direction and test whether the two arms
     are biologically coherent.  ->  myeloid_two_arms_genes.csv, myeloid_two_arms.png
  E. Does whole-section (bulk-like) aggregation reconcile the direction?
     ->  aggregation_reconciliation.csv, myeloid_wholesection_vs_region.csv

Inputs (all produced by earlier stages):
  data/processed/GSE288365_annotated.h5ad          log-normalised X, stage 02
  results/tables/author_meta_merged_obs.csv        merged obs, stage 14
  results/tables/pseudobulk_authoranat_counts.csv  raw counts per donor x region, stage 19
  results/tables/pseudobulk_authoranat_groups.csv
  results/tables/de_authoranat_<contrast>_<region>.csv   pyDESeq2 results, stage 19
  results/tables/mediator_gene_annotation_k3.csv   118-gene k=3 clustering, stage 21

NOTE ON REGION NAMES. The h5ad categories are the authors' labels; the display
names used throughout the report are:
    Anterior_horns         -> anterior horn
    Grey_matter            -> posterior horn
    White_matter           -> anterior/lateral WM
    Posterior_white_matter -> posterior WM

NOTE ON CATEGORICALS. auth_manual stores 3,121 missing values as code -1.  We read
the merged CSV, where they are NaN, rather than decoding codes from the h5ad.

Run from the repository root:
    python code/analysis/24_c9_vs_sals_resolution.py
"""
import os

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
           "NUMBA_NUM_THREADS"):
    os.environ[_v] = "1"
os.environ["KMP_AFFINITY"] = "disabled"
os.environ["OMP_PROC_BIND"] = "FALSE"

import h5py
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy.sparse as ss
from scipy import stats
from statsmodels.stats.multitest import multipletests

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
TAB = os.path.join(ROOT, "results", "tables")
FIG = os.path.join(ROOT, "results", "figures")
H5 = os.path.join(ROOT, "data", "processed", "GSE288365_annotated.h5ad")

REGIONS = ["Anterior_horns", "Grey_matter", "White_matter", "Posterior_white_matter"]
LABEL = dict(zip(REGIONS, ["anterior horn", "posterior horn",
                           "anterior/lateral WM", "posterior WM"]))
MYELOID3 = ["Microglia", "ProliferatingMicroglia", "Macrophages"]

# Genes used in the state-leakage test (B).  Lineage markers are expressed by
# myeloid cells regardless of activation; activation markers are induced.
LINEAGE = ["P2RY12", "TMEM119", "CX3CR1", "CSF1R", "SPI1", "AIF1",
           "C1QA", "C1QB", "FCER1G", "LAPTM5"]
ACTIVATION = ["GPNMB", "CHIT1", "SPP1", "CD68", "APOE", "TREM2", "LGALS3",
              "APOC1", "CTSD", "CTSB"]

# Functional modules for the two-arm coherence test (D).  Curated from canonical
# myeloid biology, assigned before looking at the direction of each gene.
MODULES = {
    "TREM2-DAP12-Syk ITAM module": ["TREM2", "TYROBP", "SYK", "NCKAP1L", "SPI1"],
    "NADPH oxidase / respiratory burst": ["NCF2", "CYBA", "NCF1", "CYBB"],
    "inflammasome / complement receptor": ["PYCARD", "C3AR1", "ADORA3", "TNFRSF14"],
    "lysosomal hydrolase / degradative": ["GRN", "IFI30", "CTSD", "CTSS", "CD63", "LAPTM5"],
    "lipid storage / esterification": ["SOAT1", "GPNMB", "APOC1", "APOBR", "PLA2G7"],
    "osmolyte / amino-acid transport": ["SLC5A3", "SLC38A2", "SLC37A2"],
    "membrane remodelling / exocytosis": ["ANO6", "STX3", "NBEAL2", "ANXA2",
                                          "S100A11", "EMP3", "DSP"],
    "myeloid lineage TF": ["SPI1", "CEBPA"],
    "oxidative-stress defence": ["SOD2"],
}

KEY_GENES = ["TREM2", "APOE", "GPNMB", "CHIT1", "APOC1", "SPI1", "CD68", "TYROBP"]


# --------------------------------------------------------------------------- IO
def load_obs():
    """Merged per-spot metadata.  Missing auth_manual values are already NaN here."""
    obs = pd.read_csv(os.path.join(TAB, "author_meta_merged_obs.csv"), low_memory=False)
    return obs


def load_de(contrast, region):
    """One pyDESeq2 result table from stage 19, indexed by gene symbol."""
    df = pd.read_csv(os.path.join(TAB, f"de_authoranat_{contrast}_{region}.csv"))
    gene_col = [c for c in df.columns
                if c.lower() in ("gene", "symbol", "index", "unnamed: 0")][0]
    return df.rename(columns={gene_col: "gene"}).set_index("gene")


def load_panel():
    """The 118 mediator genes with their k=3 cluster labels."""
    ann = pd.read_csv(os.path.join(TAB, "mediator_gene_annotation_k3.csv"))
    gene_col = [c for c in ann.columns
                if c.lower() in ("gene", "gene_name", "symbol")][0]
    return ann, gene_col, dict(zip(ann[gene_col], ann.cluster_k3))


def pseudobulk_cpm(counts):
    """log2 CPM from a samples x genes raw-count matrix."""
    return np.log2(counts.div(counts.sum(axis=1), axis=0) * 1e6 + 1)


# ------------------------------------------------------ A. three-group ordering
def three_way(cluster_of):
    """
    Pseudobulk log-CPM for KEY_GENES in each region, split control / sALS / C9.

    The bulk causal model compares C9 against sALS.  Our earlier summaries quoted
    only that contrast, which made TREM2 and APOE look like they fall in C9.  Both
    are in fact raised above control; they are simply raised less than in sALS.
    Reporting all three levels is what makes the two statements compatible.
    """
    counts = pd.read_csv(os.path.join(TAB, "pseudobulk_authoranat_counts.csv"), index_col=0)
    groups = pd.read_csv(os.path.join(TAB, "pseudobulk_authoranat_groups.csv"))
    groups["key"] = groups.donor_id + "|" + groups.anat
    groups = groups.set_index("key").reindex(counts.index)

    rows = []
    for region in REGIONS:
        mask = groups.anat.values == region
        expr = pseudobulk_cpm(counts[mask])
        cond = groups.condition[mask].values
        for gene in KEY_GENES:
            if gene not in expr.columns:
                continue
            by = {c: expr.loc[cond == c, gene].values for c in ("control", "sALS", "C9-ALS")}
            rows.append(dict(
                gene=gene, region=LABEL[region],
                control=np.mean(by["control"]), sALS=np.mean(by["sALS"]),
                C9=np.mean(by["C9-ALS"]),
                p_C9_vs_control=stats.mannwhitneyu(by["C9-ALS"], by["control"])[1],
                p_C9_vs_sALS=stats.mannwhitneyu(by["C9-ALS"], by["sALS"])[1],
                p_sALS_vs_control=stats.mannwhitneyu(by["sALS"], by["control"])[1]))
    out = pd.DataFrame(rows)
    out.to_csv(os.path.join(TAB, "key_genes_three_way.csv"), index=False)
    return out


# ------------------------------------- B. can deconvolution be regressed out?
def deconvolution_confound(obs):
    """
    Three questions about using Cell2location abundance as a covariate.

    1.  State leakage.  If the myeloid abundance estimate were largely reading the
        activation programme rather than cell number, it would track ACTIVATION
        markers much better than LINEAGE markers.  Compared at donor x region
        level, region-centred so that between-region differences do not drive it.
    2.  Collinearity.  How much variance in the DAM signal would survive
        regressing out the myeloid fraction.
    3.  External validation.  IBA1 immunofluorescence is image-derived and so
        independent of the RNA.  Compared within section, because staining
        intensity is not comparable across slides.
    """
    counts = pd.read_csv(os.path.join(TAB, "pseudobulk_authoranat_counts.csv"), index_col=0)
    groups = pd.read_csv(os.path.join(TAB, "pseudobulk_authoranat_groups.csv"))
    groups["key"] = groups.donor_id + "|" + groups.anat
    groups = groups.set_index("key").reindex(counts.index)
    expr = pseudobulk_cpm(counts)
    expr["donor"], expr["anat"] = groups.donor_id.values, groups.anat.values

    ct_cols = [c for c in obs.columns if c.startswith("auth_norm_ct_")]
    cell = obs.groupby(["donor_id", "auth_manual"], observed=True)[ct_cols].mean()
    cell.columns = [c.replace("auth_norm_ct_", "") for c in cell.columns]
    cell = cell.reset_index().rename(columns={"donor_id": "donor", "auth_manual": "anat"})
    joined = expr.merge(cell, on=["donor", "anat"])
    joined["MyeloidTotal"] = joined[MYELOID3].sum(axis=1)

    def centred(col):
        return joined[col] - joined.groupby("anat", observed=True)[col].transform("mean")

    rows = []
    for label, genes in (("lineage/homeostatic", LINEAGE), ("activation/DAM", ACTIVATION)):
        for gene in [g for g in genes if g in joined.columns]:
            rows.append(dict(
                test="marker_vs_C2L", marker_class=label, item=gene,
                rho_vs_Microglia=stats.spearmanr(centred(gene), centred("Microglia"))[0],
                rho_vs_totalMyeloid=stats.spearmanr(centred(gene), centred("MyeloidTotal"))[0]))

    present = [g for g in ACTIVATION if g in joined.columns]
    for region in REGIONS:
        sub = joined[joined.anat == region]
        dam = sub[present].mean(axis=1)
        rho = stats.spearmanr(sub.MyeloidTotal, dam)[0]
        rows.append(dict(
            test="collinearity", marker_class="DAM_expr_vs_myeloid_fraction",
            item=LABEL[region], rho_vs_Microglia=rho, rho_vs_totalMyeloid=rho,
            residual_variance=1 - np.corrcoef(sub.MyeloidTotal, dam)[0, 1] ** 2))

    # IBA1 immunofluorescence, within section.  auth_cleaned_IBA1_binary is
    # misnamed: it is a continuous intensity, not a binary call.
    iba1 = pd.to_numeric(obs.get("auth_cleaned_IBA1_binary"), errors="coerce")
    if iba1 is not None and iba1.notna().any():
        img = obs[["section_key", "donor_id", "auth_manual"]].copy()
        img["iba1"] = iba1.values
        img = img.dropna(subset=["auth_manual", "iba1"])
        per = (img.groupby(["section_key", "donor_id", "auth_manual"], observed=True)
                  .iba1.mean().reset_index()
                  .rename(columns={"donor_id": "donor", "auth_manual": "anat"}))
        merged = joined.merge(per, on=["donor", "anat"])

        def sec_centred(col):
            return merged[col] - merged.groupby("section_key", observed=True)[col].transform("mean")

        for col, name in (("Microglia", "C2L_Microglia"),
                          ("MyeloidTotal", "C2L_totalMyeloid"),
                          ("Neurons", "C2L_Neurons")):
            rows.append(dict(
                test="IF_validation", marker_class="IBA1_intensity_within_section",
                item=name,
                rho_vs_Microglia=stats.spearmanr(sec_centred("iba1"), sec_centred(col))[0]))

    out = pd.DataFrame(rows)
    out.to_csv(os.path.join(TAB, "deconvolution_state_confound.csv"), index=False)
    return out


# --------------------------------- C. why the score and pseudobulk disagree
def score_vs_pseudobulk(obs, cluster_of, panel_genes):
    """
    The per-spot mean-z score and the pseudobulk median log2FC can disagree on
    sign in the posterior regions.  They do not disagree about individual genes
    (per-gene correlation +0.48 to +0.71); they disagree because the score is a
    MEAN over genes that a handful of widely-detected or heavy-tailed genes
    dominate, while the median log2FC gives every gene one vote.

    Correction to an earlier claim: z-scoring equalises each gene's contribution
    to the score LEVEL (1/n each), but NOT its contribution to the between-group
    DIFFERENCE, which here ranges from ~36% for a single gene down to zero.
    """
    with h5py.File(H5, "r") as handle:
        gene_names = np.array([g.decode() if isinstance(g, bytes) else g
                               for g in handle["var/_index"][:]])
        n_spots = handle["obs/_index"].shape[0]
        keep = [g for g in panel_genes if g in set(gene_names)]
        idx = np.array([int(np.where(gene_names == g)[0][0]) for g in keep])
        dense = np.zeros((n_spots, len(idx)), dtype=np.float32)
        indptr, indices, data = handle["X/indptr"], handle["X/indices"], handle["X/data"]
        for start in range(0, n_spots, 20000):           # 7 GB RAM ceiling
            stop = min(start + 20000, n_spots)
            ptr = indptr[start:stop + 1]
            block = ss.csr_matrix(
                (data[ptr[0]:ptr[-1]], indices[ptr[0]:ptr[-1]], ptr - ptr[0]),
                shape=(stop - start, len(gene_names)))
            dense[start:stop] = block[:, idx].toarray()

    # Global z per gene across all spots - the same reference used by stage 20.
    zed = (dense - dense.mean(axis=0)) / dense.std(axis=0)
    myeloid_pos = [i for i, g in enumerate(keep) if cluster_of.get(g) == "Myeloid/DAM"]
    myeloid_genes = [keep[i] for i in myeloid_pos]

    frame = pd.DataFrame(zed[:, myeloid_pos], columns=myeloid_genes)
    frame["donor"] = obs.donor_id.values
    frame["cond"] = obs.condition.values
    frame["anat"] = obs.auth_manual.values
    per_donor = (frame.dropna(subset=["anat"])
                      .groupby(["donor", "cond", "anat"], observed=True)
                      .mean(numeric_only=True).reset_index())

    parts = []
    for region in REGIONS:
        sub = per_donor[per_donor.anat == region]
        z_diff = (sub[sub.cond == "C9-ALS"][myeloid_genes].median()
                  - sub[sub.cond == "sALS"][myeloid_genes].median())
        de = load_de("C9ALS_vs_sALS", region)
        detected = pd.Series(
            100 * (dense[obs.auth_manual.values == region][:, myeloid_pos] > 0).mean(axis=0),
            index=myeloid_genes)
        part = pd.DataFrame({
            "z_diff": z_diff,
            "LFC": de.log2FoldChange.reindex(myeloid_genes),
            "baseMean": de.baseMean.reindex(myeloid_genes),
            "detect_pct": detected})
        part["share_pct"] = 100 * part.z_diff / (part.z_diff.mean() * len(part))
        part["region"] = LABEL[region]
        parts.append(part)

    out = pd.concat(parts)
    out.index.name = "gene"
    out.reset_index().to_csv(
        os.path.join(TAB, "score_vs_pseudobulk_decomposition.csv"), index=False)
    return out


# ------------------------------------------- D. the two arms of the myeloid set
def two_arms(cluster_of, myeloid_genes):
    """
    Classify each myeloid gene by the direction of C9-vs-sALS across regions, then
    ask whether the two arms are functionally coherent and whether each is
    disease-induced relative to control.

    A gene is assigned to an arm if it moves the same way in at least 3 of 4
    regions.  Direction consistency, not significance: only TREM2 and NCF2 clear
    genome-wide FDR individually.
    """
    frames = []
    for region in REGIONS:
        vs_sals = load_de("C9ALS_vs_sALS", region)
        vs_ctrl = load_de("C9ALS_vs_control", region)
        sals_ctrl = load_de("sALS_vs_control", region)
        frames.append(pd.DataFrame({
            "gene": myeloid_genes, "region": LABEL[region],
            "LFC_C9vsS": vs_sals.log2FoldChange.reindex(myeloid_genes).values,
            "p_C9vsS": vs_sals.pvalue.reindex(myeloid_genes).values,
            "LFC_C9vsC": vs_ctrl.log2FoldChange.reindex(myeloid_genes).values,
            "LFC_SvsC": sals_ctrl.log2FoldChange.reindex(myeloid_genes).values,
            "baseMean": vs_sals.baseMean.reindex(myeloid_genes).values}))
    long = pd.concat(frames)

    agg = long.groupby("gene").agg(
        n_up=("LFC_C9vsS", lambda s: (s > 0).sum()),
        n_tested=("LFC_C9vsS", "count"),
        mean_LFC=("LFC_C9vsS", "mean"),
        mean_C9vsC=("LFC_C9vsC", "mean"),
        mean_SvsC=("LFC_SvsC", "mean"),
        baseMean=("baseMean", "mean"),
        min_p=("p_C9vsS", "min"))
    agg["arm"] = np.where(agg.n_up >= 3, "UP in C9",
                          np.where(agg.n_up <= 1, "DOWN in C9", "mixed"))
    module_of = {g: m for m, genes in MODULES.items() for g in genes}
    agg["module"] = pd.Series(module_of).reindex(agg.index)

    wide = long.pivot_table(index="gene", columns="region", values="LFC_C9vsS")
    wide = wide[[LABEL[r] for r in REGIONS]]

    # Crossover: each arm is induced in both genotypes, but each genotype favours
    # one.  Genes are co-regulated, so this p is descriptive rather than calibrated.
    up, down = agg[agg.arm == "UP in C9"], agg[agg.arm == "DOWN in C9"]
    crossover_p = stats.mannwhitneyu(up.mean_C9vsC - up.mean_SvsC,
                                     down.mean_C9vsC - down.mean_SvsC)[1]
    # Guard against a detection-level explanation for the split.
    basemean_p = stats.mannwhitneyu(up.baseMean.dropna(), down.baseMean.dropna())[1]

    order = pd.concat([
        up.sort_values("mean_LFC", ascending=False),
        agg[agg.arm == "mixed"].sort_values("mean_LFC", ascending=False),
        down.sort_values("mean_LFC")])
    order.reset_index().to_csv(os.path.join(TAB, "myeloid_two_arms_genes.csv"), index=False)
    _plot_two_arms(order, wide.reindex(order.index), agg, crossover_p)
    return agg, crossover_p, basemean_p


def _plot_two_arms(order, wide, agg, crossover_p):
    plt.rcParams.update({"font.size": 8, "axes.spines.top": False,
                         "axes.spines.right": False, "figure.dpi": 150})
    n_up = int((order.arm == "UP in C9").sum())
    n_mixed = int((order.arm == "mixed").sum())
    n_down = len(order) - n_up - n_mixed
    vmax = np.nanmax(np.abs(wide.values))

    fig = plt.figure(figsize=(7.4, 7.2))
    grid = fig.add_gridspec(1, 2, width_ratios=[1, 1], wspace=0.70,
                            left=0.135, right=0.93, top=0.885, bottom=0.095)

    ax = fig.add_subplot(grid[0, 0])
    img = ax.imshow(wide.values, cmap="RdBu_r", vmin=-vmax, vmax=vmax, aspect="auto")
    ax.set_xticks(range(wide.shape[1]))
    ax.set_xticklabels(wide.columns, rotation=38, ha="right")
    ax.set_yticks(range(len(wide)))
    ax.set_yticklabels(wide.index, fontsize=5.5)
    plt.setp(ax.get_yticklabels(), style="italic")
    for y in (n_up - 0.5, n_up + n_mixed - 0.5):
        ax.axhline(y, color="k", lw=1.2)
    bands = [((0, n_up), f"higher in C9-ALS\n(n={n_up})", "#b2182b"),
             ((n_up, n_up + n_mixed), f"mixed\n(n={n_mixed})", "#777777"),
             ((n_up + n_mixed, len(wide)), f"higher in sALS\n(n={n_down})", "#2166ac")]
    for (lo, hi), text, colour in bands:
        ax.add_patch(mpl.patches.Rectangle((wide.shape[1] - 0.4, lo - 0.5), 0.22, hi - lo,
                                           color=colour, clip_on=False))
        ax.text(wide.shape[1] + 0.05, (lo + hi) / 2 - 0.5, text, rotation=270,
                va="center", ha="left", fontsize=6.6, color=colour, fontweight="bold")
    ax.set_title("Mediator myeloid genes split into two\nopposing arms", loc="left", pad=7)
    ax.text(-0.33, 1.045, "a", transform=ax.transAxes, fontsize=11, fontweight="bold")
    cax = fig.add_axes([0.135, -0.045, 0.26, 0.016])
    bar = fig.colorbar(img, cax=cax, orientation="horizontal")
    bar.set_label("log$_2$FC, C9-ALS vs sALS", fontsize=6)
    bar.ax.tick_params(labelsize=5.5)

    ax2 = fig.add_subplot(grid[0, 1])
    for arm, colour in (("UP in C9", "#b2182b"), ("DOWN in C9", "#2166ac")):
        sub = agg[agg.arm == arm]
        for _, row in sub.iterrows():
            ax2.plot([0, 1], [row.mean_SvsC, row.mean_C9vsC],
                     color=colour, lw=0.55, alpha=0.28, zorder=1)
        ax2.plot([0, 1], [sub.mean_SvsC.mean(), sub.mean_C9vsC.mean()],
                 color=colour, lw=2.8, marker="o", ms=7.5, zorder=3)
        label = "higher in\nC9-ALS" if arm == "UP in C9" else "higher in\nsALS"
        ax2.text(1.06, sub.mean_C9vsC.mean(), f"{label}\n(n={len(sub)})",
                 color=colour, va="center", fontsize=6.8, fontweight="bold")
    ax2.axhline(0, color="#999999", lw=0.8, ls="--")
    ax2.text(-0.10, 0.02, "control level", fontsize=6, color="#999999", va="bottom")
    ax2.set_xticks([0, 1])
    ax2.set_xticklabels(["sALS", "C9-ALS"])
    ax2.set_xlim(-0.15, 1.52)
    ax2.set_ylabel("log$_2$FC vs control (mean of 4 regions)")
    ax2.set_title("Both arms are induced in disease;\neach genotype favours one",
                  loc="left", pad=7)
    ax2.text(0.02, 0.03, f"crossover $p$ = {crossover_p:.1e}",
             transform=ax2.transAxes, fontsize=6.5, color="#333333")
    ax2.text(-0.26, 1.045, "b", transform=ax2.transAxes, fontsize=11, fontweight="bold")

    fig.savefig(os.path.join(FIG, "myeloid_two_arms.png"), dpi=300, bbox_inches="tight")
    plt.close(fig)


# ------------------------------------------- E. does bulk-like pooling reconcile?
def aggregation_reconciliation(cluster_of, panel_genes):
    """
    The bulk cohort is whole cervical cord; this analysis is region-stratified.
    Pooling all four regions into one pseudobulk per donor mimics bulk RNA-seq.

    A shift-share decomposition separates the part of the whole-section difference
    caused by a different regional mix (composition) from the part caused by
    different expression within regions.  C9 sections carry ~8 pp more posterior
    WM and ~6 pp less anterior/lateral WM than sALS, but composition accounts for
    only 10-14% of the signal; most of the shift comes from the fact that log-CPM
    on pooled RNA weights by transcript abundance.
    """
    counts = pd.read_csv(os.path.join(TAB, "pseudobulk_authoranat_counts.csv"), index_col=0)
    groups = pd.read_csv(os.path.join(TAB, "pseudobulk_authoranat_groups.csv"))
    groups["key"] = groups.donor_id + "|" + groups.anat
    groups = groups.set_index("key").reindex(counts.index)

    whole = counts.groupby(groups.donor_id.values).sum()
    cond = groups.groupby("donor_id").condition.first().reindex(whole.index)
    whole_expr = pseudobulk_cpm(whole)
    als = whole.index[cond.isin(["C9-ALS", "sALS"])]

    by_region, weights = {}, {}
    for region in REGIONS:
        sub = counts[groups.anat.values == region].copy()
        sub.index = groups.loc[sub.index, "donor_id"].values
        by_region[region] = pseudobulk_cpm(sub)
    for donor in whole.index:
        total = groups[groups.donor_id == donor].n_spots.sum()
        weights[donor] = {r: groups[(groups.donor_id == donor) & (groups.anat == r)]
                             .n_spots.sum() / total for r in REGIONS}
    weights = pd.DataFrame(weights).T

    c9 = [d for d in als if cond[d] == "C9-ALS"]
    sals = [d for d in als if cond[d] == "sALS"]

    def contrast(expr, genes):
        return expr.reindex(c9)[genes].mean() - expr.reindex(sals)[genes].mean()

    rows, per_gene = [], {}
    for cluster in ("Myeloid/DAM", "Astrocyte", "Oligodendrocyte/myelin"):
        genes = [g for g in panel_genes if cluster_of.get(g) == cluster
                 and g in counts.columns]
        pooled = contrast(whole_expr, genes)
        stratified = pd.concat([contrast(by_region[r], genes) for r in REGIONS],
                               axis=1).mean(axis=1)
        w9, ws = weights.loc[c9].mean(), weights.loc[sals].mean()
        e9 = pd.DataFrame({r: by_region[r].reindex(c9)[genes].mean() for r in REGIONS})
        es = pd.DataFrame({r: by_region[r].reindex(sals)[genes].mean() for r in REGIONS})
        composition = ((w9 - ws) * ((e9 + es) / 2)).sum(axis=1)
        expression = (((w9 + ws) / 2) * (e9 - es)).sum(axis=1)
        rows.append(dict(cluster=cluster, n=len(genes),
                         region_avg=stratified.median(), whole_section=pooled.median(),
                         pct_up_regionavg=100 * (stratified > 0).mean(),
                         pct_up_wholesection=100 * (pooled > 0).mean(),
                         composition_term=composition.median(),
                         expression_term=expression.median()))
        if cluster == "Myeloid/DAM":
            per_gene = pd.DataFrame({"whole_section": pooled,
                                     "region_stratified": stratified})

    out = pd.DataFrame(rows)
    out.to_csv(os.path.join(TAB, "aggregation_reconciliation.csv"), index=False)
    per_gene.to_csv(os.path.join(TAB, "myeloid_wholesection_vs_region.csv"))
    shift = 100 * (weights.loc[c9].mean() - weights.loc[sals].mean())
    return out, shift


# -------------------------------------------- panel-restricted FDR, for reference
def panel_restricted_fdr(cluster_of, panel_genes):
    """
    BH within the panel rather than genome-wide.  The more generous test: 97-106
    genes instead of ~13,000.  Ten gene x region results survive, nine of them
    lower in C9.
    """
    cluster_rows, gene_rows = [], []
    for region in REGIONS:
        de = load_de("C9ALS_vs_sALS", region).reset_index()
        de["cluster"] = de.gene.map(cluster_of)
        sub = de[de.gene.isin(set(panel_genes))].dropna(subset=["pvalue"]).copy()
        sub["padj_panel"] = multipletests(sub.pvalue, method="fdr_bh")[1]
        genome_median = de.log2FoldChange.dropna().median()
        for _, row in sub[sub.padj_panel < 0.05].sort_values("pvalue").iterrows():
            gene_rows.append(dict(region=LABEL[region], gene=row.gene,
                                  cluster=row.cluster, log2FC=row.log2FoldChange,
                                  p=row.pvalue, padj_panel=row.padj_panel,
                                  padj_genomewide=row.get("padj", np.nan)))
        for cluster in ("Myeloid/DAM", "Astrocyte", "Oligodendrocyte/myelin",
                        "ALL 118 (detected)"):
            part = sub if cluster.startswith("ALL") else sub[sub.cluster == cluster]
            values = part.log2FoldChange.dropna()
            if len(values) <= 3:
                continue
            cluster_rows.append(dict(
                region=LABEL[region], cluster=cluster, n=len(values),
                median_LFC=values.median(),
                vs_genomewide=values.median() - genome_median,
                p_signrank=stats.wilcoxon(values - genome_median)[1],
                n_nominal_p05=int((part.pvalue < 0.05).sum())))
    clusters = pd.DataFrame(cluster_rows)
    genes = pd.DataFrame(gene_rows)
    clusters.to_csv(os.path.join(TAB, "panel118_C9vsSALS_cluster_summary.csv"), index=False)
    genes.to_csv(os.path.join(TAB, "panel118_C9vsSALS_significant_genes.csv"), index=False)
    return clusters, genes


def main():
    os.makedirs(TAB, exist_ok=True)
    os.makedirs(FIG, exist_ok=True)
    obs = load_obs()
    ann, gene_col, cluster_of = load_panel()
    panel_genes = list(ann[gene_col])

    three = three_way(cluster_of)
    print(f"[A] three-group ordering written for {three.gene.nunique()} genes "
          f"x {three.region.nunique()} regions")

    confound = deconvolution_confound(obs)
    marker = confound[confound.test == "marker_vs_C2L"]
    print("[B] C2L vs markers (mean rho, total myeloid): " + ", ".join(
        f"{k} {v:+.2f}" for k, v in
        marker.groupby("marker_class").rho_vs_totalMyeloid.mean().items()))
    coll = confound[confound.test == "collinearity"]
    print(f"    collinearity rho range {coll.rho_vs_totalMyeloid.min():+.2f} to "
          f"{coll.rho_vs_totalMyeloid.max():+.2f}")

    decomposition = score_vs_pseudobulk(obs, cluster_of, panel_genes)
    for region, part in decomposition.groupby("region", sort=False):
        top5 = part.z_diff.sort_values(ascending=False).head(5).sum() / part.z_diff.sum()
        print(f"[C] {region:22s} score diff {part.z_diff.mean():+.3f}  "
              f"median LFC {part.LFC.median():+.3f}  top-5 share {100 * top5:3.0f}%  "
              f"corr {part[['z_diff', 'LFC']].corr().iloc[0, 1]:+.2f}")

    myeloid_genes = [g for g in panel_genes if cluster_of.get(g) == "Myeloid/DAM"]
    agg, crossover_p, basemean_p = two_arms(cluster_of, myeloid_genes)
    print(f"[D] arms: {(agg.arm == 'UP in C9').sum()} up in C9, "
          f"{(agg.arm == 'DOWN in C9').sum()} down, {(agg.arm == 'mixed').sum()} mixed; "
          f"crossover p={crossover_p:.1e}; baseMean difference p={basemean_p:.3f}")

    recon, shift = aggregation_reconciliation(cluster_of, panel_genes)
    for _, row in recon.iterrows():
        print(f"[E] {row.cluster:24s} region-avg {row.region_avg:+.3f} -> "
              f"whole-section {row.whole_section:+.3f}")
    print("    spot-share shift C9 minus sALS (pp): " + ", ".join(
        f"{LABEL[r]} {shift[r]:+.1f}" for r in REGIONS))

    clusters, genes = panel_restricted_fdr(cluster_of, panel_genes)
    print(f"[F] panel-restricted FDR: {len(genes)} gene x region results survive "
          f"({(genes.log2FC < 0).sum()} lower in C9)")


if __name__ == "__main__":
    main()
