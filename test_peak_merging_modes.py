"""
Milan Multiome - common peak set + fragment requantification (ATAC)
===================================================================
For analysis 3 (DMSO vs REV differential accessibility) peaks must live in a
shared feature space, but CellRanger ARC called peaks per sample. For each
timepoint we build a common peak set from DMSO u REV, requantify fragments
per cell against it, attach a unified peak annotation, and write paired
AnnData objects.

Cell set is inherited from the RNA object (same multiome barcodes), so ATAC
and RNA stay aligned cell-for-cell.
"""

import os
import numpy as np
import pandas as pd
import anndata as ad
import muon as mu
import pybedtools

try:
    from sortedcontainers import SortedList
    _HAVE_SORTED = True
except ImportError:                       # pure-python fallback is fine for our sizes
    _HAVE_SORTED = False

# ----------------------------------------------------------------------
# Config
# ----------------------------------------------------------------------
DATA_ROOT  = "/store24/project24/ladcol_023/Milan_Multiome/data"
ATAC_OUT   = "/store24/project24/ladcol_023/Milan_Multiome/ATAC_requantified"
RNA_H5AD   = os.path.join(DATA_ROOT, "rna_processed.h5ad")

TIMEPOINTS = ["D0", "D3", "D7"]
PEAK_MERGE_METHOD = "iterative_overlap"    # "iterative_overlap" (recommended) | "union"

STANDARD_CHROMS = ["chr%d" % i for i in range(1, 23)] + ["chrX", "chrY"]
BLACKLIST_BED = "/store24/project24/ladcol_023/Milan_Multiome/hg38-blacklist.v2.bed"

BARCODE_RE = r"([ACGT]{16})"


# ----------------------------------------------------------------------
# Helpers
# ----------------------------------------------------------------------
def _extract_bc(names):
    return pd.Index(names).str.extract(BARCODE_RE, expand=False)


def _raw_barcodes(rna, sample):
    """cellranger '-1' barcodes of the RNA cells belonging to `sample`."""
    names = rna.obs_names[rna.obs["sample"] == sample]
    return set(_extract_bc(names) + "-1")


# ----------------------------------------------------------------------
# 1. Load ATAC per sample, inherit the RNA cell set, register fragments
# ----------------------------------------------------------------------
def load_sample_atac(sample, rna):
    """Return the ATAC AnnData for `sample`, subset to the RNA cells, with
    obs_names kept as cellranger '-1' and the fragments file registered."""
    outs = os.path.join(DATA_ROOT, sample, "outs")
    mdata = mu.read_10x_mtx(os.path.join(outs, "filtered_feature_bc_matrix"))
    mdata.var_names_make_unique()

    # register fragments explicitly (lives in outs/, not the matrix dir)
    mu.atac.tl.locate_fragments(
        mdata, fragments=os.path.join(outs, "atac_fragments.tsv.gz"))

    atac = mdata.mod["atac"]
    keep = atac.obs_names.isin(_raw_barcodes(rna, sample))
    atac = atac[keep].copy()
    atac.uns["files"] = mdata.mod["atac"].uns.get("files", {})
    return atac


# ----------------------------------------------------------------------
# 2. Peaks -> scored BED frame; build common peak set
# ----------------------------------------------------------------------
def extract_peaks(atac):
    """Peaks from var_names as a DataFrame (chrom,start,end,score), score =
    total accessibility (used to rank peaks for iterative overlap)."""
    coords = atac.var_names.str.extract(r"(chr[\w]+)[:\-_](\d+)[:\-_](\d+)")
    coords.columns = ["chrom", "start", "end"]
    coords = coords.dropna()
    coords["start"] = coords["start"].astype(int)
    coords["end"] = coords["end"].astype(int)
    score = np.asarray(atac[:, coords.index].X.sum(axis=0)).ravel()
    coords["score"] = score
    return coords.reset_index(drop=True)


def _iterative_overlap(peaks):
    """ArchR-style: keep highest-scoring peak, drop overlappers, repeat.
    Returns a non-overlapping subset of `peaks`."""
    keep = []
    for chrom, grp in peaks.groupby("chrom", sort=False):
        grp = grp.sort_values("score", ascending=False)
        kept = SortedList() if _HAVE_SORTED else []
        for idx, s, e in zip(grp.index, grp["start"], grp["end"]):
            if _HAVE_SORTED:
                pos = kept.bisect_left((s, s))
                nbrs = []
                if pos > 0:
                    nbrs.append(kept[pos - 1])
                if pos < len(kept):
                    nbrs.append(kept[pos])
            else:
                nbrs = kept
            if not any(ps < e and s < pe for ps, pe in nbrs):
                kept.add((s, e)) if _HAVE_SORTED else kept.append((s, e))
                keep.append(idx)
    return peaks.loc[keep].sort_values(["chrom", "start"])


def build_common_peaks(atac_dmso, atac_rev, method=PEAK_MERGE_METHOD):
    peaks = pd.concat([extract_peaks(atac_dmso), extract_peaks(atac_rev)],
                      ignore_index=True)
    peaks = peaks[peaks["chrom"].isin(STANDARD_CHROMS)]

    if method == "union":
        bed = pybedtools.BedTool.from_dataframe(peaks[["chrom", "start", "end"]])
        merged = bed.sort().merge().to_dataframe(names=["chrom", "start", "end"])
    elif method == "iterative_overlap":
        merged = _iterative_overlap(peaks)[["chrom", "start", "end"]]
    else:
        raise ValueError(f"unknown method: {method}")

    if BLACKLIST_BED:
        bl = pybedtools.BedTool(BLACKLIST_BED)
        merged = (pybedtools.BedTool.from_dataframe(merged)
                  .intersect(bl, v=True)
                  .to_dataframe(names=["chrom", "start", "end"]))

    merged = merged.sort_values(["chrom", "start"]).reset_index(drop=True)
    merged["peak_name"] = (merged["chrom"] + ":" +
                           merged["start"].astype(str) + "-" +
                           merged["end"].astype(str))
    return merged


# ----------------------------------------------------------------------
# Run
# ----------------------------------------------------------------------
if __name__ == "__main__":
    # Load preprocessed RNA data
    rna = ad.read_h5ad('/store24/project24/ladcol_023/Milan_Multiome/data/rna_processed.h5ad')
    
    for tp in TIMEPOINTS:
        print(f"\n==================== {tp} ====================")
        dmso_key, rev_key = f"DMSO_{tp}", f"REV_{tp}"
    
        atac_dmso = load_sample_atac(dmso_key, rna)
        atac_rev  = load_sample_atac(rev_key, rna)
        
        pk_d = extract_peaks(atac_dmso)
        pk_r = extract_peaks(atac_rev)
        union = build_common_peaks(atac_dmso, atac_rev, method="union")
        io    = build_common_peaks(atac_dmso, atac_rev, method="iterative_overlap")
        print(f"DMSO input: {len(pk_d)} | REV input: {len(pk_r)} | sum: {len(pk_d)+len(pk_r)}")
        print(f"union: {len(union)} | iterative: {len(io)}")
        print(f"iterative median width: {int(io.eval('end - start').median())}")