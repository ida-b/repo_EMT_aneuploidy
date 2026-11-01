"""
Export per-timepoint ATAC data for chromVAR (run in the `scanalysis` env).
========================================================================
chromVAR runs in R, so we hand it plain files instead of .h5ad:
  chromvar_input/<tp>/counts.mtx   peaks x cells (integer fragment counts)
  chromvar_input/<tp>/peaks.bed    chrom, start, end, peak_name  (BED, 0-based)
  chromvar_input/<tp>/cells.csv    barcode + condition/sample + grouping columns

The R script (chromvar_analysis.R) reads these back.
"""

import os
import numpy as np
import pandas as pd
import anndata as ad
from scipy.io import mmwrite
from scipy.sparse import csr_matrix

from differential_analysis import (
    PROJECT, RNA_H5AD, ATAC_OUT, TIMEPOINTS, EMT_SCORE, attach,
)

OUT = os.path.join(PROJECT, "chromvar_input")


if __name__ == "__main__":
    rna = ad.read_h5ad(RNA_H5AD)

    for tp in TIMEPOINTS:
        d = os.path.join(OUT, tp)
        os.makedirs(d, exist_ok=True)

        atac_dmso = ad.read_h5ad(os.path.join(ATAC_OUT, tp, "DMSO.h5ad"))
        atac_rev  = ad.read_h5ad(os.path.join(ATAC_OUT, tp, "REV.h5ad"))
        combined = ad.concat([atac_dmso, atac_rev], join="inner")
        combined = attach(combined, rna, ["cnv_burden_proportion", EMT_SCORE])

        # peaks.bed  (parse chr:start-end from the common-peak names)
        pk = combined.var_names.str.extract(r"(chr[\w]+)[:\-_](\d+)[:\-_](\d+)")
        pk.columns = ["chrom", "start", "end"]
        pk["name"] = combined.var_names
        pk[["chrom", "start", "end", "name"]].to_csv(
            os.path.join(d, "peaks.bed"), sep="\t", header=False, index=False)

        # counts.mtx  (peaks x cells, integer)
        mmwrite(os.path.join(d, "counts.mtx"),
                csr_matrix(combined.X).T.astype("int32"))

        # cells.csv  (grouping metadata, same order as counts columns)
        cells = combined.obs[["condition", "sample",
                              "cnv_burden_proportion", EMT_SCORE]].copy()
        cells.insert(0, "barcode", combined.obs_names)
        cells.to_csv(os.path.join(d, "cells.csv"), index=False)

        print(f"{tp}: {combined.n_vars} peaks x {combined.n_obs} cells -> {d}")

    print(f"\nchromVAR inputs written under {OUT}")
