"""
Characterize gene-activity hits that are NOT DEGs (REVlow vs DMSO).
==================================================================
Most GA-significant genes have no matched expression change. This script, per
timepoint:
  1. splits GA-significant genes into OPENED (lfc>0) and CLOSED (lfc<0);
  2. buffered-vs-primed diagnostic: for the OPENED genes, is their expression
     (DEG) logFC shifted in the same direction (buffered) or centered at 0
     (decoupled / primed)?
  3. enrichment of the OPENED set against the EMT / epithelial / mesenchymal
     gene sets already used for scoring (hypergeometric, local gene sets);
  4. writes ranked opened / closed gene lists.

Outputs under results/differential/: GAopen_* CSVs, and a printed summary.
"""

import os
import numpy as np
import pandas as pd
from scipy.stats import binomtest, hypergeom

from differential_analysis import EXPORT_DIR, PROJECT, TIMEPOINTS

ANNOT_DIR = os.path.join(PROJECT, "data", "annotations")
PADJ, LFC = 0.05, 0.25


def load_gene_sets():
    emt = pd.read_csv(os.path.join(
        ANNOT_DIR, "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION.v2025.1.Hs.tsv"), sep="\t")
    hallmark = emt.loc[emt["STANDARD_NAME"] == "GENE_SYMBOLS",
                       "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION"].values[0].split(",")
    gs = pd.read_csv(os.path.join(ANNOT_DIR, "1-s2.0-S2589004224013415-mmc2.csv"), sep=";")
    col = lambda n: set(x for x in gs[n].values if isinstance(x, str))
    return {
        "HALLMARK_EMT": set(hallmark),
        "Epithelial":   col("Epithelial genes"),
        "Mesenchymal":  col("Mesenchymal genes"),
        "pEMT":         col("pEMT genes"),
        "EMT_up":       col("EMT_up"),
        "EMT_down":     col("EMT_down"),
    }


if __name__ == "__main__":
    gene_sets = load_gene_sets()
    summ, enr_rows = [], []

    for tp in TIMEPOINTS:
        ga  = pd.read_csv(os.path.join(EXPORT_DIR, f"GA_REVlowAneu_vs_DMSO_{tp}.csv"))
        deg = pd.read_csv(os.path.join(EXPORT_DIR, f"DEG_REVlowAneu_vs_DMSO_{tp}.csv"))
        universe = set(ga["names"])

        ga_sig = ga[(ga["pvals_adj"] < PADJ) & (ga["logfoldchanges"].abs() > LFC)]
        opened = ga_sig[ga_sig["logfoldchanges"] > 0].copy()
        closed = ga_sig[ga_sig["logfoldchanges"] < 0].copy()
        opened.sort_values("logfoldchanges", ascending=False).to_csv(
            os.path.join(EXPORT_DIR, f"GAopen_opened_{tp}.csv"), index=False)
        closed.sort_values("logfoldchanges").to_csv(
            os.path.join(EXPORT_DIR, f"GAopen_closed_{tp}.csv"), index=False)

        # buffered-vs-primed: expression logFC of the OPENED genes
        deg_lfc = deg.set_index("names")["logfoldchanges"]
        opened_deg = deg_lfc.reindex(opened["names"]).dropna()
        n = len(opened_deg)
        same = int((opened_deg > 0).sum())
        med = float(opened_deg.median()) if n else np.nan
        bp = binomtest(same, n, 0.5).pvalue if n else np.nan
        summ.append({"timepoint": tp, "opened": len(opened), "closed": len(closed),
                     "opened_w_DEGlfc": n, "median_DEGlfc": round(med, 3),
                     "pct_DEGlfc_up": round(100 * same / n, 1) if n else np.nan,
                     "sign_p": bp})

        # enrichment of OPENED genes against EMT gene sets (hypergeometric)
        opened_genes = set(opened["names"]) & universe
        N, nn = len(universe), len(opened_genes)
        for name, gset in gene_sets.items():
            K = len(gset & universe)
            k = len(opened_genes & gset)
            p = hypergeom.sf(k - 1, N, K, nn) if (K and nn) else np.nan
            enr_rows.append({"timepoint": tp, "gene_set": name, "set_in_universe": K,
                             "opened_hits": k, "expected": round(nn * K / N, 2) if N else np.nan,
                             "p": p})

    summary = pd.DataFrame(summ)
    enrich = pd.DataFrame(enr_rows)
    summary.to_csv(os.path.join(EXPORT_DIR, "GAopen_summary.csv"), index=False)
    enrich.to_csv(os.path.join(EXPORT_DIR, "GAopen_enrichment.csv"), index=False)
    print("Opened/closed + buffered-vs-primed:\n", summary.to_string(index=False))
    print("\nEMT enrichment of OPENED genes:\n", enrich.to_string(index=False))
