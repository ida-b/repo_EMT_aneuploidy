"""
Milan Multiome - EMT / aneuploidy analysis
==========================================
Single-cell multiome (snRNA + snATAC), 6 samples:
    DMSO and REV (reversine washout) at day 0, 3, 7.

Section 1  Load, QC, integrate RNA; score EMT programs; load epiAneufinder+
           CNV calls; compute per-cell aneuploidy burden vs the cell-line
           baseline karyotype (DMSO_D0 consensus).

Section 2  [next] copyKat (RNA) cross-check of the high-aneuploidy population.
"""

import os
import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
import muon as mu
import matplotlib.pyplot as plt
import matplotlib.patches as patches

# ----------------------------------------------------------------------
# Config
# ----------------------------------------------------------------------
DATA_ROOT       = "/store24/project24/ladcol_023/Milan_Multiome/data"
RESULTS_ROOT    = "/store24/project24/ladcol_023/Milan_Multiome/results"
ANNOT_DIR       = os.path.join(DATA_ROOT, "annotations")
EPIANEU_TABLE   = "result_table_holmes"      # epiAneufinder+ per-cell somy table
BASELINE_SAMPLE = "DMSO_D0"                   # consensus karyotype of the cell line

SAMPLES = ["DMSO_D0", "REV_D0", "DMSO_D3", "REV_D3", "DMSO_D7", "REV_D7"]

# RNA QC thresholds (from visual inspection of the original notebooks)
QC = dict(min_genes=1000, max_genes=10000, max_counts=75000, max_pct_mt=40)

# Aneuploidy population thresholds (fraction of genome deviating from baseline)
HIGH_BURDEN = 0.30
LOW_BURDEN  = 0.05

# copyKat (RNA) cross-check: per-sample runs with an external diploid reference
COPYKAT_REF_TAG  = "_TSP4_Mammary_NA_10X_1_"   # tags the reference (diploid) cells
COPYKAT_BASELINE = "DMSO_D0"                    # cell-line consensus reference
COPYKAT_PATHS = {
    "DMSO_D0": os.path.join(RESULTS_ROOT, "S60989_DMSO_D0",
                            "copykat_results", "S60989_DMSO_D0_copykat_CNA_results.txt"),
    "REV_D0":  os.path.join(RESULTS_ROOT, "S60990_REV_D0",
                            "copykat_results", "S60990_REV_D0_copykat_CNA_results.txt"),
    "DMSO_D3": os.path.join(RESULTS_ROOT, "S60991_DMSO_D3",
                            "copykat_results", "S60991_DMSO_D3_copykat_CNA_results.txt"),
    "REV_D3":  os.path.join(RESULTS_ROOT, "S60992_REV_D3",
                            "copykat_results", "S60992_REV_D3_copykat_CNA_results.txt"),
    "DMSO_D7": os.path.join(RESULTS_ROOT, "S60993_DMSO_D7",
                            "copykat_results", "S60993_DMSO_D7_copykat_CNA_results.txt"),
    "REV_D7":  os.path.join(RESULTS_ROOT, "S60994_REV_D7",
                            "copykat_results", "S60994_REV_D7_copykat_CNA_results.txt"),
}
REV_TIMEPOINTS = ["REV_D0", "REV_D3", "REV_D7"]


# ----------------------------------------------------------------------
# Helpers
# ----------------------------------------------------------------------
def _chrom_key(seq):
    """Sort key for chromosomes: 1..22, then X, Y, M."""
    c = str(seq).replace("chr", "")
    return {"X": 23, "Y": 24, "M": 25, "MT": 25}.get(c, int(c) if c.isdigit() else 99)


def tag_barcodes(names, sample):
    """Replace the cellranger '-1' suffix with '-<sample>' (anchored, so an
    internal '-1' is never touched); otherwise append '-<sample>'."""
    if names.str.endswith("-1").all():
        return names.str.replace(r"-1$", f"-{sample}", regex=True)
    return names + f"-{sample}"


# ----------------------------------------------------------------------
# 1. Load + QC each RNA modality
# ----------------------------------------------------------------------
def load_rna(sample):
    path = os.path.join(DATA_ROOT, sample, "outs", "filtered_feature_bc_matrix")
    mdata = mu.read_10x_mtx(path)
    mdata.var_names_make_unique()
    rna = mdata.mod["rna"]
    sc.pp.filter_genes(rna, min_cells=1)
    rna.var["mt"] = rna.var_names.str.startswith("MT-")
    sc.pp.calculate_qc_metrics(rna, qc_vars=["mt"], inplace=True, log1p=False)
    return rna


def qc_filter(rna):
    sc.pp.filter_cells(rna, min_genes=QC["min_genes"])
    sc.pp.filter_cells(rna, max_genes=QC["max_genes"])
    sc.pp.filter_cells(rna, max_counts=QC["max_counts"])
    return rna[rna.obs["pct_counts_mt"] < QC["max_pct_mt"]].copy()


def plot_qc(rnas):
    metrics = ["n_genes_by_counts", "total_counts", "pct_counts_mt"]
    fig, axes = plt.subplots(len(metrics), len(rnas), figsize=(3 * len(rnas), 9),
                             sharey="row", sharex="col")
    for j, (label, rna) in enumerate(rnas.items()):
        for i, m in enumerate(metrics):
            ax = axes[i, j]
            ax.violinplot(rna.obs[m].dropna())
            ax.set_title(label if i == 0 else "")
            ax.set_ylabel(m if j == 0 else "")
            ax.set_xticks([]); ax.grid(True)
    fig.tight_layout()
    return fig


# ----------------------------------------------------------------------
# 2. Integrate RNA
# ----------------------------------------------------------------------
def build_rna(samples):
    rnas = {}
    for s in samples:
        r = qc_filter(load_rna(s))
        r.obs_names = tag_barcodes(r.obs_names, s)
        rnas[s] = r

    rna = ad.concat(rnas, join="outer", label="sample", index_unique=None)
    rna.obs["sample"] = rna.obs["sample"].astype("category")
    rna.obs["condition"] = rna.obs["sample"].astype(str).str.split("_").str[0]
    rna.obs["day"] = rna.obs["sample"].astype(str).str.split("_").str[1]
    rna.layers["counts"] = rna.X.copy()

    sc.pp.normalize_total(rna, target_sum=1e4)
    sc.pp.log1p(rna)
    sc.pp.highly_variable_genes(rna, n_top_genes=2000, batch_key="sample")
    sc.pp.pca(rna)
    sc.pp.neighbors(rna)
    sc.tl.umap(rna)
    return rna


# ----------------------------------------------------------------------
# 3. EMT program scoring
# ----------------------------------------------------------------------
def score_emt(rna):
    emt = pd.read_csv(os.path.join(
        ANNOT_DIR, "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION.v2025.1.Hs.tsv"), sep="\t")
    hallmark = emt.loc[emt["STANDARD_NAME"] == "GENE_SYMBOLS",
                       "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION"].values[0].split(",")

    gs = pd.read_csv(os.path.join(ANNOT_DIR, "1-s2.0-S2589004224013415-mmc2.csv"), sep=";")
    col = lambda name: [x for x in gs[name].values if isinstance(x, str)]

    signatures = {
        "epi_score":         ["CDH1", "EPCAM"],
        "mesen_score":       ["VIM", "ZEB1", "CDH2", "ITGA5", "FOXC2", "TWIST1"],
        "EMT_score":         hallmark,
        "epithelial_score":  col("Epithelial genes"),
        "mesenchymal_score": col("Mesenchymal genes"),
        "luminal_score":     col("Luminal genes"),
        "basal_score":       col("Basal genes"),
        "pEMT_score":        col("pEMT genes"),
        "EMT_down_score":    col("EMT_down"),
        "EMT_partial_score": col("EMT_partial"),
        "EMT_up_score":      col("EMT_up"),
    }
    for name, genes in signatures.items():
        sc.tl.score_genes(rna, gene_list=genes, score_name=name)
    return rna


# ----------------------------------------------------------------------
# 4. epiAneufinder+ CNV calls -> aligned bins x cells matrix
# ----------------------------------------------------------------------
def load_epianeufinder(sample):
    """DataFrame indexed by (seq, start, end), columns = cells, values = somy
    (0 loss / 1 base / 2 gain). Barcodes tagged with sample."""
    path = os.path.join(RESULTS_ROOT, sample, "pyepiAneufinder", "outs",
                        f"{EPIANEU_TABLE}.tsv.gz")
    df = pd.read_csv(path, sep="\t", index_col=0).set_index(["seq", "start", "end"])
    df.columns = tag_barcodes(df.columns, sample)
    return df


def build_cnv_matrix(samples, keep_cells):
    somies = {s: load_epianeufinder(s) for s in samples}

    common = None
    for df in somies.values():
        common = df.index if common is None else common.intersection(df.index)
    common = sorted(common, key=lambda b: (_chrom_key(b[0]), b[1]))
    common = pd.MultiIndex.from_tuples(common, names=["seq", "start", "end"])

    parts = []
    for df in somies.values():
        df = df.reindex(common)
        parts.append(df[[c for c in df.columns if c in keep_cells]])
    return pd.concat(parts, axis=1)          # bins x cells, shared bin index


def add_cnv_burden(rna, cnv):
    """cnv: bins x cells. Baseline = rounded mean somy over BASELINE_SAMPLE.
    Burden = number of bins differing from baseline."""
    shared = [c for c in rna.obs_names if c in cnv.columns]
    rna = rna[shared].copy()
    cnv = cnv[shared].T                       # cells x bins
    rna.obsm["cnv"] = cnv.values

    base_cells = [c for c in cnv.index if c.endswith(f"-{BASELINE_SAMPLE}")]
    baseline = np.round(cnv.loc[base_cells].mean(axis=0)).values

    rna.obs["cnv_burden"] = (cnv.values != baseline).sum(axis=1)
    rna.obs["cnv_burden_proportion"] = rna.obs["cnv_burden"] / cnv.shape[1]
    rna.uns["cnv_baseline"] = baseline
    rna.uns["cnv_bins"] = cnv.columns.to_frame(index=False)
    return rna


def plot_karyotype(baseline, bins):
    """Scatter the baseline somy call along the genome (loss/base/gain)."""
    counts = bins["seq"].value_counts()
    counts = counts.reindex(sorted(counts.index, key=_chrom_key))

    fig = plt.figure(figsize=(12, 1.2))
    gs = fig.add_gridspec(2, 1, height_ratios=[0.05, 0.1], hspace=0.05)
    ax0, ax1 = fig.add_subplot(gs[0]), None

    start = 0
    for i, (chrom, n) in enumerate(counts.items()):
        ax0.add_patch(patches.Rectangle((start, 0), n, 1,
                      color="lightgrey" if i % 2 == 0 else "grey"))
        ax0.text(start + n / 2, 1.5, str(chrom).replace("chr", ""),
                 ha="center", va="center", fontsize=8)
        start += n
    ax0.set_xlim(0, start); ax0.set_ylim(0, 2); ax0.axis("off")

    ax1 = fig.add_subplot(gs[1], sharex=ax0)
    x = np.arange(len(baseline))
    for val, color in [(1, "#00EE76"), (0, "#9A32CD"), (2, "#CD0000")]:
        idx = np.where(baseline == val)[0]
        ax1.scatter(x[idx], np.repeat(1, len(idx)), color=color, s=2)
    start = 0
    for _, n in counts.items():
        start += n
        ax1.axvline(start, color="black", ls="--", lw=0.5)
    ax1.set_ylim(0.8, 1.2); ax1.set_yticklabels([])
    return fig


# ----------------------------------------------------------------------
# 5. copyKat (RNA) cross-check of the high-aneuploidy population
# ----------------------------------------------------------------------
# copyKat was run per sample with an EXTERNAL diploid reference (Tabula
# Sapiens mammary cells, tagged COPYKAT_REF_TAG), so its CNA values measure
# deviation from TRUE diploid and every cancer-line cell looks aneuploid.
# To make it comparable to the epiAneufinder burden (deviation from the
# DMSO_D0 cell-line consensus), we re-reference copyKat the same way:
# baseline = median copyKat profile over DMSO_D0 cells; score = mean
# |deviation| from it. Join RNA<->ATAC on the bare 16bp barcode.

BARCODE_RE = r"([ACGT]{16})"

def _extract_bc(names):
    return pd.Index(names).str.extract(BARCODE_RE, expand=False)


def load_copykat(path, ref_tag=COPYKAT_REF_TAG):
    """copyKat CNA matrix -> (windows x sample-cells); reference cells dropped.
    Window index = (chrom, chrompos)."""
    ck = pd.read_csv(path, sep="\t")
    coord = ck.columns[:3]                         # chrom, chrompos, abspos
    cells = ck.columns[3:]
    sample_cells = cells[~cells.str.contains(ref_tag)]
    mat = ck[sample_cells].copy()
    mat.index = pd.MultiIndex.from_frame(ck[list(coord[:2])])
    return mat


def copykat_deviation(baseline_path, target_path):
    """Per-cell mean |deviation| of target cells from the copyKat cell-line
    consensus (median over baseline-sample cells), on shared windows.
    Returns a Series indexed by the 16bp cell barcode."""
    base = load_copykat(baseline_path)
    tgt  = load_copykat(target_path)
    common = base.index.intersection(tgt.index)
    base, tgt = base.loc[common], tgt.loc[common]
    consensus = base.median(axis=1)
    dev = tgt.sub(consensus, axis=0).abs().mean(axis=0)
    dev.index = _extract_bc(dev.index)
    return dev[~dev.index.duplicated()]


def add_copykat_crosscheck(rna, sample, baseline=COPYKAT_BASELINE):
    """Attach copyKat RNA-side deviation to rna.obs['copykat_dev'] for `sample`."""
    dev = copykat_deviation(COPYKAT_PATHS[baseline], COPYKAT_PATHS[sample])
    mask = rna.obs["sample"] == sample
    rna_bc = _extract_bc(rna.obs_names[mask])
    mapped = pd.Series(rna_bc, index=rna.obs_names[mask]).map(dev)
    if "copykat_dev" not in rna.obs:
        rna.obs["copykat_dev"] = np.nan
    rna.obs.loc[mapped.index, "copykat_dev"] = mapped.values
    matched = mapped.notna().sum()
    print(f"[{sample}] copyKat joined to {matched}/{mask.sum()} cells")
    return rna


def validate_high_aneuploid(rna, sample, high=HIGH_BURDEN, low=LOW_BURDEN,
                            dev_thr=None, normal_quantile=0.95):
    """Split epiAneufinder high-burden cells by whether copyKat agrees.
    Discordant (high ATAC burden, low RNA deviation) = suspected inference
    failures. Default 'normal' ceiling = `normal_quantile` of copyKat_dev
    over the sample's low-burden (DMSO-like) cells -- i.e. the upper edge of
    the normal range, a stringent test for genuine aneuploidy."""
    sub = rna.obs[(rna.obs["sample"] == sample) & rna.obs["copykat_dev"].notna()]
    if dev_thr is None:
        normal = sub.loc[sub["cnv_burden_proportion"] < low, "copykat_dev"]
        dev_thr = normal.quantile(normal_quantile)
    hi = sub[sub["cnv_burden_proportion"] > high]
    concordant = hi.index[hi["copykat_dev"] >= dev_thr]
    discordant = hi.index[hi["copykat_dev"] <  dev_thr]
    print(f"[{sample}] high-burden: {len(hi)} | copyKat-concordant: {len(concordant)}"
          f" | suspected inference failures: {len(discordant)} (dev_thr={dev_thr:.3f})")
    return concordant, discordant, dev_thr


def get_aneuploid_populations(rna, sample, high=HIGH_BURDEN, low=LOW_BURDEN):
    """Return (validated_high, low) barcode sets for a REV timepoint:
    validated_high = epiAneufinder high-burden cells confirmed by copyKat;
    low = DMSO-like low-burden cells. Feeds analyses 2 and 3."""
    validated_high, _, _ = validate_high_aneuploid(rna, sample, high, low)
    low_cells = rna.obs_names[(rna.obs["sample"] == sample) &
                              (rna.obs["cnv_burden_proportion"] < low)]
    return list(validated_high), list(low_cells)


def plot_crosscheck(rna, sample, high=HIGH_BURDEN, low=LOW_BURDEN, dev_thr=None):
    from scipy.stats import spearmanr
    sub = rna.obs[(rna.obs["sample"] == sample) & rna.obs["copykat_dev"].notna()]
    fig, ax = plt.subplots(figsize=(5, 5))
    ax.scatter(sub["cnv_burden_proportion"], sub["copykat_dev"], s=6, alpha=0.5)
    ax.axvline(high, color="red", ls="--", lw=0.8)
    ax.axvline(low, color="grey", ls="--", lw=0.8)
    if dev_thr is not None:
        ax.axhline(dev_thr, color="black", ls=":", lw=0.8)
    rho = spearmanr(sub["cnv_burden_proportion"], sub["copykat_dev"]).correlation
    ax.text(0.05, 0.95, f"Spearman rho = {rho:.2f}", transform=ax.transAxes, va="top")
    ax.set_xlabel("epiAneufinder+ burden (ATAC)")
    ax.set_ylabel("copyKat deviation from baseline (RNA)")
    ax.set_title(f"{sample}: ATAC vs RNA aneuploidy")
    return fig


# ----------------------------------------------------------------------
# Run
# ----------------------------------------------------------------------
if __name__ == "__main__":
    rna = build_rna(SAMPLES)
    rna = score_emt(rna)
    cnv = build_cnv_matrix(SAMPLES, set(rna.obs_names))
    rna = add_cnv_burden(rna, cnv)

    # copyKat cross-check over every REV timepoint (baseline = DMSO_D0)
    populations = {}
    for tp in REV_TIMEPOINTS:
        rna = add_copykat_crosscheck(rna, tp)
        populations[tp] = get_aneuploid_populations(rna, tp)   # (high_validated, low)
    print(rna)

    # Save
    rna.write(DATA_ROOT + "/rna_processed.h5ad", compression='gzip')