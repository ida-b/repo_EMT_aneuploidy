#!/usr/bin/env Rscript
# chromVAR TF-motif differential accessibility (run in the `caf_atac_r` env)
# =========================================================================
# Reads the per-timepoint exports from export_for_chromvar.py, computes per-cell
# TF-motif accessibility deviations (JASPAR2022 CORE, human), and tests the two
# CNV-safe comparisons per timepoint:
#   - low-aneuploidy REV vs DMSO
#   - EMT high vs low, within the low-aneuploidy REV subset
# Writes results/differential/chromVAR_*.csv (one row per TF motif).

suppressMessages({
  library(chromVAR); library(motifmatchr)
  library(SummarizedExperiment); library(GenomicRanges)
  library(Matrix); library(BiocParallel)
  library(TFBSTools); library(JASPAR2022)
  library(BSgenome.Hsapiens.UCSC.hg38)
})
ncores <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "4"))
register(MulticoreParam(workers = ncores))
set.seed(2024)

PROJECT      <- "/store24/project24/ladcol_023/Milan_Multiome"
IN_DIR       <- file.path(PROJECT, "chromvar_input")
OUT_DIR      <- file.path(PROJECT, "scripts", "repo_EMT_aneuploidy", "results", "differential")
TIMEPOINTS   <- c("D0", "D3", "D7")
LOW_BURDEN   <- 0.05
MIN_CELLS    <- 10
MIN_PCT_PEAK <- 0.02   # drop peaks accessible in <2% of cells (speed)
N_BG         <- 25     # background iterations for deviations (default 50)
genome       <- BSgenome.Hsapiens.UCSC.hg38
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
cat(sprintf("chromVAR with %d cores\n", ncores))

# JASPAR2022 CORE vertebrate motifs; lookup from motif ID -> TF name
motifs   <- getMatrixSet(JASPAR2022, list(species = 9606, collection = "CORE"))
motif_id <- vapply(motifs, ID,   character(1))
names(motifs) <- motif_id
tf_lookup <- data.frame(motif = motif_id,
                        tf    = vapply(motifs, name, character(1)),
                        stringsAsFactors = FALSE)

run_cmp <- function(dev, cells_a, cells_b, label_a, label_b, fname) {
  if (length(cells_a) < MIN_CELLS || length(cells_b) < MIN_CELLS) {
    cat(sprintf("    skip %s: groups too small (%d vs %d)\n",
                fname, length(cells_a), length(cells_b)))
    return(invisible())
  }
  sub <- dev[, c(cells_a, cells_b)]
  grp <- factor(c(rep(label_a, length(cells_a)), rep(label_b, length(cells_b))),
                levels = c(label_a, label_b))
  res <- differentialDeviations(sub, grp)
  z   <- deviationScores(sub)
  ma  <- rowMeans(z[, grp == label_a, drop = FALSE], na.rm = TRUE)
  mb  <- rowMeans(z[, grp == label_b, drop = FALSE], na.rm = TRUE)

  out <- data.frame(motif = rownames(res), stringsAsFactors = FALSE)
  out$tf <- tf_lookup$tf[match(out$motif, tf_lookup$motif)]
  out[[paste0("meanZ_", label_a)]] <- ma[out$motif]
  out[[paste0("meanZ_", label_b)]] <- mb[out$motif]
  out$deltaZ  <- ma[out$motif] - mb[out$motif]
  out$p_value <- res$p_value
  out$p_adj   <- res$p_value_adjusted
  out <- out[order(out$p_adj), ]
  write.csv(out, file.path(OUT_DIR, fname), row.names = FALSE)
  cat(sprintf("    %s vs %s: %d motifs, %d sig (padj<0.05)\n",
              label_a, label_b, nrow(out), sum(out$p_adj < 0.05, na.rm = TRUE)))
}

for (tp in TIMEPOINTS) {
  cat(sprintf("==================== %s ====================\n", tp))
  d <- file.path(IN_DIR, tp)

  peaks <- read.table(file.path(d, "peaks.bed"), sep = "\t",
                      col.names = c("chrom", "start", "end", "name"))
  gr <- GRanges(peaks$chrom, IRanges(peaks$start + 1L, peaks$end))  # BED -> 1-based
  names(gr) <- peaks$name

  counts <- as(readMM(file.path(d, "counts.mtx")), "CsparseMatrix")
  cells  <- read.csv(file.path(d, "cells.csv"), stringsAsFactors = FALSE)
  rownames(cells) <- cells$barcode
  rownames(counts) <- peaks$name
  colnames(counts) <- cells$barcode

  SE <- SummarizedExperiment(assays = list(counts = counts),
                             rowRanges = gr, colData = DataFrame(cells))
  SE <- addGCBias(SE, genome = genome)
  SE <- filterPeaks(SE, non_overlapping = TRUE)

  # detection pre-filter: fewer peaks -> much faster matchMotifs / deviations
  det <- Matrix::rowMeans(assay(SE, "counts") > 0)
  SE <- SE[det >= MIN_PCT_PEAK, ]
  cat(sprintf("    peaks after filter: %d\n", nrow(SE)))

  motif_ix <- matchMotifs(motifs, SE, genome = genome)
  bg  <- getBackgroundPeaks(SE, niterations = N_BG)
  dev <- computeDeviations(object = SE, annotations = motif_ix,
                           background_peaks = bg)
  cd  <- as.data.frame(colData(dev))

  lowrev <- rownames(cd)[cd$condition == "REV" &
                         cd$cnv_burden_proportion < LOW_BURDEN]
  dmso   <- rownames(cd)[cd$condition == "DMSO"]
  run_cmp(dev, lowrev, dmso, "REVlowAneu", "DMSO",
          sprintf("chromVAR_REVlowAneu_vs_DMSO_%s.csv", tp))

  emt <- cd[lowrev, "EMT_up_score"]
  qs  <- quantile(emt, c(0.25, 0.75), na.rm = TRUE)
  hi  <- lowrev[!is.na(emt) & emt >= qs[2]]
  lo  <- lowrev[!is.na(emt) & emt <= qs[1]]
  run_cmp(dev, hi, lo, "EMThigh", "EMTlow",
          sprintf("chromVAR_EMT_highVSlow_lowAneu_%s.csv", tp))
}

cat(sprintf("\nchromVAR CSVs written to %s\n", OUT_DIR))
