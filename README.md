# Multiome analysis: EMT and aneuploidy in a breast epithelial cancer line

Analysis code and results for the single-cell multiome experiment on EMT
(epithelial–mesenchymal transition) and aneuploidy. This file describes the
experiment, the scripts, the comparisons, and how to read the output tables. A
glossary is at the end.

## Experimental design

Six samples: DMSO (control) and reversine (REV) at day 0, 3, and 7. Reversine
inhibits MPS1, which allows chromosome mis-segregation and induces aneuploidy; it
was washed out before the day-0 sample. Each sample is single-cell multiome —
paired snRNA (gene expression) and snATAC (open chromatin) measured in the same
cells.

## Aneuploidy

The REV samples are less aneuploid than expected. Small subpopulations carry
copy-number changes absent from DMSO, and aneuploidy decreases from day 0 to day
7, but the effect is small. The likely explanation is that the most aneuploid
cells were too unstable to pass QC and were not sequenced.

Aneuploidy was called from the ATAC data with epiAneufinder+ and cross-checked
against copyKat run on the RNA. The two callers agree, so the small effect is not
an artifact of the caller. Cells called highly aneuploid are additionally
required to be concordant between the two callers; this removes a small number of
ATAC inference failures, and the number removed increases at day 7 where genuine
aneuploidy is rarest.

## Repository layout

```
scripts/
  processing.py ....................... load and QC the 6 samples, integrate the
                                        RNA, score EMT, call aneuploidy burden,
                                        run the copyKat cross-check. Writes the
                                        main RNA object used downstream.
  fragment_requantification_common_peaks.py
                                        build one ATAC peak set per timepoint and
                                        requantify each cell against it.
  differential_analysis.py ............ expression (DEG) and per-peak accessibility
                                        (DAP) tests. Writes DEG_ and DAP_ CSVs.
  gene_activity_analysis.py ........... accessibility summed per gene, tested like
                                        the RNA. Writes GA_ CSVs.
  integrate_results.py ................ compares expression and gene-activity
                                        results gene by gene. Writes INTEGRATED_.
  gene_activity_opened_analysis.py .... splits the significant gene-activity genes
                                        into opened/closed and tests EMT
                                        enrichment. Writes GAopen_.
  export_for_chromvar.py
  chromvar_analysis.R ................. chromVAR TF-motif analysis. Writes chromVAR_.
  chromvar_rank.py .................... re-ranks the chromVAR motifs by |deltaZ|
                                        and flags EMT TFs. Writes chromVARranked_.
  residual_cnv_check.py ............... checks whether accessibility hits sit in
                                        residual copy-number regions.
  test_peak_merging_modes.py .......... reports peak counts per merge strategy.
  *.sbatch ............................ SLURM wrappers.
  logs/

notebooks/
  overview_plots.ipynb ................ summary figures: aneuploidy per sample,
                                        EMT per sample, UMAPs, aneuploidy vs EMT,
                                        the ATAC-vs-RNA cross-check.

results/
  differential/ ....................... output CSV tables.
  figures/ ............................ figures from the notebook.

environment.yml
```

## Common peak set

CellRanger calls ATAC peaks separately per sample, so the DMSO and REV peak sets
differ and cannot be compared cell-for-cell. `fragment_requantification_common_peaks.py`
builds one peak set per timepoint from the DMSO and REV peaks (iterative-overlap
merge: rank peaks by accessibility, keep the strongest, drop weaker overlapping
peaks) and requantifies every cell's fragments against it. Cells are taken from
the RNA object, so RNA and ATAC stay matched per cell.

## Comparisons

Run per timepoint (day 0, 3, 7).

Expression (RNA):
- REV high- vs low-aneuploidy
- low-aneuploidy REV vs DMSO

Accessibility (ATAC), at three levels — per peak, gene activity, and TF motif
(chromVAR):
- low-aneuploidy REV vs DMSO
- high- vs low-EMT, within the low-aneuploidy REV cells

Accessibility comparisons are restricted to homogeneous groups in terms of aneuploidy.
Copy number changes ATAC signal directly, so comparing groups that differ in
aneuploidy would confound copy number with regulation. High- vs low-aneuploidy is
therefore tested in expression only, and the EMT comparison is run inside the
low-aneuploidy REV cells so that EMT, not aneuploidy, is what differs.

## Reading the output tables

Files are in `results/differential/`. One row per gene (expression, gene
activity) or per peak (per-peak accessibility). Columns:

| column | meaning |
|---|---|
| `names` | gene or peak |
| `logfoldchanges` | log2 fold change; positive = higher in the first-named group |
| `pvals` | raw p-value |
| `pvals_adj` | p-value corrected for multiple testing; use this one |
| `mean_*`, `pct_*`, `n_*` | per-group mean, fraction of cells with signal, cell count |
| `gene`, `distance`, `peak_type` | per-peak tables: nearest gene and distance to it |

Suggested significance filter: `pvals_adj` < 0.05 with `|logfoldchanges|` > 1 for
expression. Per-peak fold changes are smaller, so a lower cutoff (`|logfoldchanges|` > 0.25) is used
for accessibility.

File prefixes:

- `DEG_` — expression. `DAP_` — per-peak accessibility, in two normalizations,
  `_lognorm` and `_tfidf` (TF-IDF is the ATAC standard; lognorm is included as a
  check). Peaks significant in both are the most robust.
- `GA_` — gene activity: open signal summed per gene. More sensitive than
  per-peak testing, and comparable to the `DEG_` tables gene for gene.
- `chromVAR_` — TF motifs. One row per motif. `deltaZ` is the direction and size
  of the difference between groups; `meanZ_*` are the per-group mean motif
  z-scores. **Rank by `|deltaZ|`, not by `p_adj`.** `differentialDeviations` tests
  across cells, so with thousands of cells per group the p-values are inflated by
  cell count — a large fraction of the 692 motifs pass `p_adj` < 0.05 regardless
  of effect size, so that threshold is not informative here. `deltaZ` does not
  scale with cell number. `chromvar_rank.py` writes `chromVARranked_` files sorted
  by `|deltaZ|` with EMT-associated TFs flagged; use those, not the raw p-value
  order.
  Motifs from the same TF family recognise near-identical sequences and move
  together (e.g. AP-1: FOS/JUN/FOSL/BATF; p53: TP53/TP63/TP73; E-box:
  ZEB/SNAI/TCF/bHLH), so interpret by family, not by single TF.
  A `tf` name with a double colon (e.g. `MAFG::NFE2L1`) is a heterodimer motif:
  the site bound by the two factors together, not either alone. A factor may also
  appear on its own in a separate row; the solo and paired entries are different
  motifs.
- `INTEGRATED_` — expression vs gene activity, gene by gene, for low-aneuploidy
  REV vs DMSO. `INTEGRATED_summary.csv` gives, per timepoint, the number of genes
  significant in both, their direction agreement, and the genome-wide correlation
  between the two fold changes. The per-timepoint files list the shared genes.
- `GAopen_` — the gene-activity genes that are not also DEGs. `GAopen_opened_` and
  `GAopen_closed_` list the genes that gain or lose accessibility. `GAopen_summary.csv`
  reports whether the opened genes show any matched expression shift.
  `GAopen_enrichment.csv` tests the opened genes against the EMT gene sets.

To test gene lists for pathway enrichment, paste the significant genes from any
table into e.g. Enrichr (https://maayanlab.cloud/Enrichr/).

## Caveats

- No biological replicates (one culture per condition/timepoint). Tests use
  individual cells as the unit, so p-values are optimistic; rely on effect sizes
  and consistency across timepoints.
- Per-peak accessibility testing is underpowered at this peak count and sparsity.
  Gene activity and chromVAR aggregate the signal and are more sensitive.
- The low-aneuploidy REV cells (burden < 0.05) are close to, not identical to,
  DMSO. The 0.05 cutoff is below the mean DMSO burden against the consensus
  profile (~0.075).

## Reproducing

```bash
conda env create -f environment.yml
conda activate scanalysis
cd scripts
sbatch processing.sbatch                               # data/rna_processed.h5ad
sbatch fragment_requantification_common_peaks.sbatch   # ATAC_processed/<tp>/*.h5ad
sbatch differential_analysis.sbatch                    # DEG_, DAP_
sbatch gene_activity_analysis.sbatch                   # GA_
python integrate_results.py                            # INTEGRATED_
python gene_activity_opened_analysis.py                # GAopen_
sbatch chromvar_analysis.sbatch                        # chromVAR_ (see below)
python chromvar_rank.py                                # chromVARranked_
# figures: notebooks/overview_plots.ipynb
```

Paths are set at the top of each script (`*_ROOT` / `PROJECT`) and point at the
cluster location of the data. epiAneufinder+ and copyKat run upstream and are not
included here.

chromVAR runs in a separate R environment. The Python export step uses
`scanalysis`; the R step needs these packages:

```bash
conda activate caf_atac_r
conda install -c conda-forge -c bioconda \
  bioconductor-chromvar bioconductor-motifmatchr bioconductor-tfbstools \
  bioconductor-jaspar2022 bioconductor-bsgenome.hsapiens.ucsc.hg38
```

`chromvar_analysis.sbatch` runs the export (`export_for_chromvar.py`) and the R
analysis (`chromvar_analysis.R`) in sequence.

## Glossary

- EMT — epithelial-to-mesenchymal transition.
- Aneuploidy / CNV — abnormal chromosome or segment copy number (copy-number
  variation).
- MPS1 / reversine — the kinase inhibited and the drug used to induce aneuploidy.
  REV = reversine; DMSO = control.
- Multiome — paired RNA (expression) and ATAC (open chromatin) in the same cells.
- Peak — an accessible region in ATAC.
- DEG / DAP — differentially expressed gene / differentially accessible peak.
- Gene activity — ATAC signal summed per gene.
- epiAneufinder+ / copyKat — aneuploidy callers, from ATAC and RNA respectively.
- UMAP — 2-D embedding placing similar cells together.
