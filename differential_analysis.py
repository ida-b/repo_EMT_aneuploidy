"""
Milan Multiome - differential expression (DEG) and accessibility (DAP) exports
==============================================================================
Per timepoint, CSV per comparison. Only CNV-safe comparisons are run:

DEG (RNA)
  - REV: high vs low aneuploidy
  - low-aneuploidy REV vs DMSO
  (REV high-vs-low EMT dropped: circular, groups defined on expression)

DAP (ATAC, common peak set)
  - REV: high vs low EMT, WITHIN the low-aneuploidy subset (CNV-controlled)
  - low-aneuploidy REV vs DMSO (CNV-matched; residual <5% burden acknowledged,
    not corrected, per design decision)

Deliberately NOT run (confounded by copy number):
  - DAP REV high vs low aneuploidy
  - DAP DMSO vs high-aneuploidy REV

Cell groupings come from the RNA object (cnv_burden_proportion, copykat_dev,
EMT_up_score) and are transferred to ATAC by the shared multiome barcode.
"""

import os
import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad

# ----------------------------------------------------------------------
# Config
# ----------------------------------------------------------------------
PROJECT    = "/store24/project24/ladcol_023/Milan_Multiome"
RNA_H5AD   = os.path.join(PROJECT, "data", "rna_processed.h5ad")   # must carry
                                                                  # copykat_dev + EMT_up_score
ATAC_OUT   = os.path.join(PROJECT, "ATAC_processed")
EXPORT_DIR = os.path.join(PROJECT, "scripts", "repo_EMT_aneuploidy", "results", "differential")

TIMEPOINTS  = ["D0", "D3", "D7"]
HIGH_BURDEN = 0.30
LOW_BURDEN  = 0.05
NORMAL_Q    = 0.95
EMT_SCORE   = "EMT_up_score"
MIN_CELLS   = 10            # warn below this per group
METHOD      = "wilcoxon"
DA_METHODS  = ["lognorm", "tfidf"]   # ATAC normalizations to run per DAP comparison
ATAC_MIN_PCT = 0.05   # test a peak only if accessible in >=5% of cells in a group
                      # (drops near-empty peaks -> power + lighter multiple-testing)

# significance thresholds for the summary printout only (CSVs are unaffected)
SIG_PADJ     = 0.05
SIG_LFC      = 1.0    # DEG (RNA) effect-size cutoff
SIG_LFC_ATAC = 0.25   # DAP (ATAC) effect-size cutoff; per-peak logFCs are
                      # compressed, so the RNA cutoff of 1.0 is too strict here


# ----------------------------------------------------------------------
# Population definitions (from RNA obs)
# ----------------------------------------------------------------------
def aneuploidy_groups(rna, sample, high=HIGH_BURDEN, low=LOW_BURDEN, normal_q=NORMAL_Q):
    """(validated_high, low) barcodes for a REV sample. If copyKat deviation is
    present, the high set is restricted to copyKat-concordant cells."""
    s = rna.obs["sample"] == sample
    burden = rna.obs.loc[s, "cnv_burden_proportion"]
    names = pd.Index(rna.obs_names[s])
    low_cells = names[(burden < low).values]
    if "copykat_dev" in rna.obs:
        dev = rna.obs.loc[s, "copykat_dev"]
        thr = dev[burden < low].quantile(normal_q)
        high_cells = names[((burden > high) & (dev >= thr)).values]
    else:
        print(f"  ! copykat_dev missing: high-aneuploidy set NOT copyKat-validated")
        high_cells = names[(burden > high).values]
    return list(high_cells), list(low_cells)


def emt_quartile_groups(rna, cells, score=EMT_SCORE):
    """(high_EMT, low_EMT) = top vs bottom quartile of `score` within `cells`."""
    v = rna.obs.loc[cells, score]
    q1, q3 = v.quantile(0.25), v.quantile(0.75)
    return list(v.index[v >= q3]), list(v.index[v <= q1])


# ----------------------------------------------------------------------
# Core differential test -> tidy frame
# ----------------------------------------------------------------------
def tfidf(a, scale_factor=1e4):
    """Signac-style log-TF-IDF (the transform underlying LSI). TF = per-cell
    depth normalization; IDF = per-peak specificity weight (ncells / peak_total);
    then log1p(TF * IDF * scale_factor)."""
    import scipy.sparse as sp
    X = sp.csr_matrix(a.X).astype(float)
    cell = np.asarray(X.sum(1)).ravel(); cell[cell == 0] = 1
    peak = np.asarray(X.sum(0)).ravel()
    tf = X.multiply(1.0 / cell[:, None])
    idf = X.shape[0] / (peak + 1)
    M = tf.multiply(idf[None, :]).tocsr()
    M = M.multiply(scale_factor).tocsr()
    M.data = np.log1p(M.data)
    a.X = M
    return a


def _normalize(a, norm):
    """norm: None (already normalized, e.g. RNA) | 'lognorm' (CP10k+log1p) | 'tfidf'."""
    if norm == "lognorm":
        sc.pp.normalize_total(a, target_sum=1e4)
        sc.pp.log1p(a)
    elif norm == "tfidf":
        a = tfidf(a)
    elif norm is not None:
        raise ValueError(f"unknown norm: {norm}")
    return a


def differential(adata, group_cells, ref_cells, g_label, r_label,
                 norm=None, annotation=None, min_pct=0.0):
    for lbl, cs in [(g_label, group_cells), (r_label, ref_cells)]:
        if len(cs) < MIN_CELLS:
            print(f"  ! {lbl}: only {len(cs)} cells (<{MIN_CELLS}) - result is weak")

    a = adata[list(group_cells) + list(ref_cells)].copy()
    a.obs["_grp"] = [g_label] * len(group_cells) + [r_label] * len(ref_cells)

    # detection-rate pre-filter (ATAC): keep a feature only if accessible in
    # >= min_pct of cells in at least one group. Done on raw counts, before
    # normalization, so it both raises power and shrinks the FDR correction.
    if min_pct > 0:
        gX = a[a.obs["_grp"] == g_label].X
        rX = a[a.obs["_grp"] == r_label].X
        pg = np.asarray((gX > 0).mean(0)).ravel()
        pr = np.asarray((rX > 0).mean(0)).ravel()
        keep = np.maximum(pg, pr) >= min_pct
        a = a[:, keep].copy()
        print(f"     features tested: {int(keep.sum())}/{len(keep)} "
              f"(detected in >= {min_pct:.0%} of a group)")

    a = _normalize(a, norm)

    sc.tl.rank_genes_groups(a, "_grp", groups=[g_label], reference=r_label, method=METHOD)
    df = sc.get.rank_genes_groups_df(a, group=g_label)

    gX = a[a.obs["_grp"] == g_label].X
    rX = a[a.obs["_grp"] == r_label].X
    stats = pd.DataFrame({
        "names": a.var_names,
        f"mean_{g_label}": np.asarray(gX.mean(0)).ravel(),
        f"mean_{r_label}": np.asarray(rX.mean(0)).ravel(),
        f"pct_{g_label}":  np.asarray((gX > 0).mean(0)).ravel(),
        f"pct_{r_label}":  np.asarray((rX > 0).mean(0)).ravel(),
    })
    df = df.merge(stats, on="names", how="left")
    df[f"n_{g_label}"], df[f"n_{r_label}"] = len(group_cells), len(ref_cells)

    if annotation is not None:
        df = df.merge(annotation, left_on="names", right_on="peak", how="left") \
               .drop(columns="peak")

    lfc_cut = SIG_LFC_ATAC if norm is not None else SIG_LFC   # ATAC when norm is set
    n_sig = ((df["pvals_adj"] < SIG_PADJ) & (df["logfoldchanges"].abs() > lfc_cut)).sum()
    print(f"  -> {g_label} vs {r_label}: {n_sig} sig "
          f"(padj<{SIG_PADJ}, |lfc|>{lfc_cut}) of {len(df)}")
    return df


def collapse_annotation(ann):
    """One row per peak: nearest gene + distance + peak_type + all genes."""
    a = ann.assign(_abs=ann["distance"].abs()).sort_values("_abs")
    nearest = a.groupby("peak").agg(
        gene=("gene", "first"), distance=("distance", "first"),
        peak_type=("peak_type", "first"),
        all_genes=("gene", lambda x: ";".join(map(str, x.unique()))),
    ).reset_index()
    return nearest


def _numpy2_pickle_shim():
    """Let numpy<2 unpickle arrays written by numpy>=2 (numpy._core -> numpy.core)."""
    import sys, importlib
    for name in ["numpy._core", "numpy._core.numeric", "numpy._core.multiarray",
                 "numpy._core.umath", "numpy._core._multiarray_umath",
                 "numpy._core.overrides"]:
        if name in sys.modules:
            continue
        try:
            sys.modules[name] = importlib.import_module(name.replace("_core", "core"))
        except Exception:
            pass


def load_peak_annotation(tp_dir):
    """Read unified peak annotation as CSV; if only the (numpy-2) pickle exists,
    read it via the shim and cache a CSV for next time."""
    csv = os.path.join(tp_dir, "unified_peak_annotation.csv")
    pkl = os.path.join(tp_dir, "unified_peak_annotation.pkl")
    if os.path.exists(csv):
        return pd.read_csv(csv)
    _numpy2_pickle_shim()
    ann = pd.read_pickle(pkl)
    ann.to_csv(csv, index=False)
    return ann


def attach(atac, rna, cols):
    ref = rna.obs.reindex(atac.obs_names)
    for c in cols:
        atac.obs[c] = ref[c].values
    return atac


# ----------------------------------------------------------------------
# Run
# ----------------------------------------------------------------------
if __name__ == "__main__":
    os.makedirs(EXPORT_DIR, exist_ok=True)
    rna = ad.read_h5ad(RNA_H5AD)
    for col in ["cnv_burden_proportion", EMT_SCORE]:
        assert col in rna.obs, f"RNA obs missing required column: {col}"

    for tp in TIMEPOINTS:
        print(f"\n==================== {tp} ====================")
        dmso_s, rev_s = f"DMSO_{tp}", f"REV_{tp}"

        high_aneu, low_aneu = aneuploidy_groups(rna, rev_s)
        dmso_cells = list(rna.obs_names[rna.obs["sample"] == dmso_s])

        # ---------- DEG (RNA) ----------
        differential(rna, high_aneu, low_aneu,
                     "REVhighAneu", "REVlowAneu").to_csv(
            os.path.join(EXPORT_DIR, f"DEG_REV_aneuploidy_highVSlow_{tp}.csv"), index=False)

        differential(rna, low_aneu, dmso_cells,
                     "REVlowAneu", "DMSO").to_csv(
            os.path.join(EXPORT_DIR, f"DEG_REVlowAneu_vs_DMSO_{tp}.csv"), index=False)

        # ---------- DAP (ATAC) ----------
        atac_dmso = ad.read_h5ad(os.path.join(ATAC_OUT, tp, "DMSO.h5ad"))
        atac_rev  = ad.read_h5ad(os.path.join(ATAC_OUT, tp, "REV.h5ad"))
        ann = collapse_annotation(load_peak_annotation(os.path.join(ATAC_OUT, tp)))

        atac_rev = attach(atac_rev, rna, ["cnv_burden_proportion", EMT_SCORE])

        # EMT high vs low WITHIN the low-aneuploidy REV subset (quartiles in-subset)
        low_aneu_rev = list(atac_rev.obs_names[
            atac_rev.obs["cnv_burden_proportion"] < LOW_BURDEN])
        emt_hi, emt_lo = emt_quartile_groups(rna, low_aneu_rev)
        combined = ad.concat([atac_dmso, atac_rev], join="inner")

        # run each comparison under both ATAC normalizations
        for norm in DA_METHODS:
            print(f"  [{norm}]")
            differential(atac_rev, emt_hi, emt_lo, "EMThigh", "EMTlow",
                         norm=norm, annotation=ann, min_pct=ATAC_MIN_PCT).to_csv(
                os.path.join(EXPORT_DIR,
                             f"DAP_REV_EMT_highVSlow_lowAneu_{tp}_{norm}.csv"), index=False)
            differential(combined, low_aneu_rev, list(atac_dmso.obs_names),
                         "REVlowAneu", "DMSO", norm=norm, annotation=ann,
                         min_pct=ATAC_MIN_PCT).to_csv(
                os.path.join(EXPORT_DIR,
                             f"DAP_REVlowAneu_vs_DMSO_{tp}_{norm}.csv"), index=False)

    print(f"\nCSVs written to {EXPORT_DIR}")