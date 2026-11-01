"""
Milan Multiome - gene-activity differential analysis (ATAC, aggregated)
=======================================================================
Per-peak differential accessibility is badly underpowered on sparse multiome
ATAC (~240k peaks, genome-wide correction). This script aggregates peak
accessibility to the GENE level: for each gene, sum the accessibility of the
common peaks annotated to it -> a cells x genes matrix that is tested exactly
like RNA. ~20k features instead of 240k means far more power, and the output is
directly comparable to the DEG tables gene-by-gene.

Only the CNV-safe comparisons are run (gene activity is accessibility-derived,
so copy-number still confounds aneuploidy contrasts):
  - low-aneuploidy REV vs DMSO
  - REV high vs low EMT, within the low-aneuploidy subset

Outputs: results/differential/GA_*.csv
"""

import os
import numpy as np
import pandas as pd
import scipy.sparse as sp
import anndata as ad

# reuse the tested machinery + population definitions
from differential_analysis import (
    differential, emt_quartile_groups, attach,
    PROJECT, RNA_H5AD, ATAC_OUT, EXPORT_DIR,
    TIMEPOINTS, LOW_BURDEN, EMT_SCORE, ATAC_MIN_PCT,
)

# include peaks within this distance (bp) of a gene; None = every annotated peak
GA_MAX_DIST = None


def build_gene_activity(adata, ann, max_dist=GA_MAX_DIST):
    """cells x genes accessibility = peak_counts (cells x peaks) @ P2G (peaks x
    genes), where P2G[p, g] = 1 if peak p is annotated to gene g. A peak near
    several genes contributes to each (standard one-to-many gene activity)."""
    ann = ann[ann["gene"].notna()].copy()
    if max_dist is not None:
        ann = ann[ann["distance"].abs() <= max_dist]

    peaks = pd.Index(adata.var_names)
    genes = pd.Index(sorted(ann["gene"].astype(str).unique()))
    pi = peaks.get_indexer(pd.Index(ann["peak"]))
    gi = genes.get_indexer(pd.Index(ann["gene"].astype(str)))
    m = (pi >= 0) & (gi >= 0)

    p2g = sp.csr_matrix((np.ones(m.sum()), (pi[m], gi[m])),
                        shape=(len(peaks), len(genes)))
    X = sp.csr_matrix(adata.X) @ p2g
    ga = ad.AnnData(X=X, obs=adata.obs.copy(), var=pd.DataFrame(index=genes))
    ga.var_names.name = "gene"
    return ga


if __name__ == "__main__":
    os.makedirs(EXPORT_DIR, exist_ok=True)
    rna = ad.read_h5ad(RNA_H5AD)

    for tp in TIMEPOINTS:
        print(f"\n==================== {tp} ====================")
        atac_dmso = ad.read_h5ad(os.path.join(ATAC_OUT, tp, "DMSO.h5ad"))
        atac_rev  = ad.read_h5ad(os.path.join(ATAC_OUT, tp, "REV.h5ad"))
        ann = pd.read_csv(os.path.join(ATAC_OUT, tp, "unified_peak_annotation.csv"))

        combined = ad.concat([atac_dmso, atac_rev], join="inner")
        ga = build_gene_activity(combined, ann)
        ga = attach(ga, rna, ["cnv_burden_proportion", EMT_SCORE])
        print(f"  gene-activity features: {ga.n_vars} genes")

        low_aneu_rev = list(ga.obs_names[(ga.obs["condition"] == "REV") &
                                         (ga.obs["cnv_burden_proportion"] < LOW_BURDEN)])
        dmso_cells = list(ga.obs_names[ga.obs["condition"] == "DMSO"])
        emt_hi, emt_lo = emt_quartile_groups(rna, low_aneu_rev)

        # low-aneuploidy REV vs DMSO  (compare against DEG_REVlowAneu_vs_DMSO_*)
        differential(ga, low_aneu_rev, dmso_cells, "REVlowAneu", "DMSO",
                     norm="lognorm", min_pct=ATAC_MIN_PCT).to_csv(
            os.path.join(EXPORT_DIR, f"GA_REVlowAneu_vs_DMSO_{tp}.csv"), index=False)

        # EMT high vs low within the low-aneuploidy REV subset
        differential(ga, emt_hi, emt_lo, "EMThigh", "EMTlow",
                     norm="lognorm", min_pct=ATAC_MIN_PCT).to_csv(
            os.path.join(EXPORT_DIR, f"GA_REV_EMT_highVSlow_lowAneu_{tp}.csv"), index=False)

    print(f"\nGene-activity CSVs written to {EXPORT_DIR}")
