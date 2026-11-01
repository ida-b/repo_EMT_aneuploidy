"""
Rank chromVAR motifs by effect size (|deltaZ|), not p-value.
============================================================
differentialDeviations tests across cells, so with thousands of cells per group
its p-values are inflated and padj<0.05 flags a large fraction of motifs. deltaZ
(difference in mean motif z-score between groups) does not scale with cell count,
so rank by it. EMT-relevant TFs are flagged.

Writes results/differential/chromVARranked_*.csv and prints the top motifs.
"""

import os
import pandas as pd

from differential_analysis import EXPORT_DIR, TIMEPOINTS

COMPARISONS = ["REVlowAneu_vs_DMSO", "EMT_highVSlow_lowAneu"]
TOP = 20

# EMT-associated TFs (mesenchymal drivers + epithelial factors)
EMT_TFS = {
    "SNAI1", "SNAI2", "ZEB1", "ZEB2", "TWIST1", "TWIST2", "TCF3", "TCF4", "TCF12",
    "PRRX1", "FOXC2", "TEAD1", "TEAD2", "TEAD3", "TEAD4", "RUNX1", "RUNX2",
    "ETS1", "FOS", "FOSL1", "FOSL2", "JUN", "JUNB", "JUND", "SMAD3", "SMAD4",
    "LEF1", "NFKB1", "RELA",
    "GRHL2", "GRHL3", "OVOL1", "OVOL2", "ELF3", "KLF4", "KLF5",
    "TFAP2A", "TFAP2C", "GATA3", "FOXA1",
}


def is_emt(tf):
    parts = {p.strip().upper() for p in str(tf).replace("::", ";").split(";")}
    return bool(parts & EMT_TFS)


if __name__ == "__main__":
    for comp in COMPARISONS:
        for tp in TIMEPOINTS:
            f = os.path.join(EXPORT_DIR, f"chromVAR_{comp}_{tp}.csv")
            if not os.path.exists(f):
                continue
            df = pd.read_csv(f)
            df["abs_deltaZ"] = df["deltaZ"].abs()
            df["EMT_TF"] = df["tf"].map(is_emt)
            df = df.sort_values("abs_deltaZ", ascending=False)
            df.to_csv(os.path.join(EXPORT_DIR,
                      f"chromVARranked_{comp}_{tp}.csv"), index=False)

            top = df.head(TOP)[["tf", "deltaZ", "p_adj", "EMT_TF"]]
            n_emt_top = int(top["EMT_TF"].sum())
            print(f"\n=== {comp}  {tp}  (top {TOP} by |deltaZ|, "
                  f"{n_emt_top} EMT TFs) ===")
            print(top.to_string(index=False))
