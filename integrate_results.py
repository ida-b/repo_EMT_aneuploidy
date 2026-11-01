"""
Integrate expression (DEG) and gene-level accessibility (GA) results.
=====================================================================
Only REVlowAneu-vs-DMSO has both an RNA table (DEG_) and an accessibility table
(GA_), so this asks, per timepoint: which genes are significant in BOTH, and do
expression and accessibility move in the SAME direction (opened chromatin <->
higher expression)? Direction concordance is what separates a coordinated
regulatory change from accessibility noise.

Writes:
  results/differential/INTEGRATED_REVlowAneu_vs_DMSO_<tp>.csv   shared genes
  results/differential/INTEGRATED_summary.csv                   counts per tp
  results/figures/integrated_<tp>.png                           DEG vs GA logFC
"""

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import spearmanr, binomtest, hypergeom

from differential_analysis import EXPORT_DIR, PROJECT, TIMEPOINTS

FIG_DIR = os.path.join(PROJECT, "scripts", "repo_EMT_aneuploidy", "results", "figures")
os.makedirs(FIG_DIR, exist_ok=True)

PADJ       = 0.05
LFC        = 0.25   # effect-size cut used for the strict double-overlap
STRICT_DEG = 1.0    # RNA cut defining the "focused" DEG set for the sign test


def sig(df, lfc):
    return (df["pvals_adj"] < PADJ) & (df["logfoldchanges"].abs() > lfc)


if __name__ == "__main__":
    rows = []
    for tp in TIMEPOINTS:
        deg = pd.read_csv(os.path.join(EXPORT_DIR, f"DEG_REVlowAneu_vs_DMSO_{tp}.csv"))
        ga  = pd.read_csv(os.path.join(EXPORT_DIR, f"GA_REVlowAneu_vs_DMSO_{tp}.csv"))
        m = deg.merge(ga, on="names", suffixes=("_deg", "_ga"))   # shared gene universe
        N = len(m)

        # ---- (1) strict double-overlap + hypergeometric enrichment ----
        deg_sig = set(deg["names"][sig(deg, LFC)])
        ga_sig  = set(ga["names"][sig(ga, LFC)])
        both = deg_sig & ga_sig
        k = len(both)
        Kd = m["names"].isin(deg_sig).sum()   # DEG-sig within shared universe
        Kg = m["names"].isin(ga_sig).sum()    # GA-sig within shared universe
        hyper_p = hypergeom.sf(k - 1, N, Kg, Kd) if (Kd and Kg) else np.nan
        ov = m[m["names"].isin(both)].copy()
        conc = np.sign(ov["logfoldchanges_deg"]) == np.sign(ov["logfoldchanges_ga"])
        ov["dir_concordant"] = conc.values
        ov.sort_values("pvals_adj_deg").to_csv(
            os.path.join(EXPORT_DIR, f"INTEGRATED_REVlowAneu_vs_DMSO_{tp}.csv"), index=False)

        # ---- (2) genome-wide effect-size coupling (the powered read) ----
        rho, rho_p = spearmanr(m["logfoldchanges_deg"], m["logfoldchanges_ga"])

        # ---- (3) sign test: do strict DEGs shift accessibility the same way? ----
        strict = m[sig(deg.set_index("names").reindex(m["names"]).reset_index(), STRICT_DEG).values]
        same = int((np.sign(strict["logfoldchanges_deg"]) ==
                    np.sign(strict["logfoldchanges_ga"])).sum())
        n_strict = len(strict)
        binom_p = binomtest(same, n_strict, 0.5).pvalue if n_strict else np.nan

        rows.append({
            "timepoint": tp, "DEG_sig": len(deg_sig), "GA_sig": len(ga_sig),
            "overlap": k, "overlap_concordant": int(conc.sum()),
            "overlap_hyper_p": hyper_p,
            "spearman_lfc": round(rho, 3), "spearman_p": rho_p,
            "nDEGstrict": n_strict,
            "pct_GA_same_sign": round(100 * same / n_strict, 1) if n_strict else np.nan,
            "sign_binom_p": binom_p,
        })

        # scatter of ALL shared genes (global coupling), overlap highlighted
        fig, ax = plt.subplots(figsize=(5, 5))
        ax.scatter(m["logfoldchanges_deg"], m["logfoldchanges_ga"], s=3, alpha=0.2,
                   c="#999999")
        if k:
            ax.scatter(ov["logfoldchanges_deg"], ov["logfoldchanges_ga"], s=20,
                       c=np.where(conc, "#4C72B0", "#C44E52"))
        ax.axhline(0, lw=0.5, c="grey"); ax.axvline(0, lw=0.5, c="grey")
        ax.set_xlabel("expression log2FC (DEG)")
        ax.set_ylabel("gene-activity log2FC (GA)")
        ax.set_title(f"{tp} REVlow vs DMSO  (Spearman rho={rho:.2f}, n={N})")
        fig.tight_layout()
        fig.savefig(os.path.join(FIG_DIR, f"integrated_{tp}.png"), dpi=150)
        plt.close(fig)

    summary = pd.DataFrame(rows)
    summary.to_csv(os.path.join(EXPORT_DIR, "INTEGRATED_summary.csv"), index=False)
    print(summary.to_string(index=False))