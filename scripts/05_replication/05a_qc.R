#!/usr/bin/env Rscript
# ==============================================================================
# 05a_qc.R
# Replication quality control (median +/- 3 MAD) and doublet scoring.
# Produces: replication cell set
# Inputs  : peaks/replication_combined.rds
# Outputs : qc/04_QCmetric_distributions_with_cutoffs.png
#           qc/01_DensityScatter_ATAC.png
#           qc/02_histogram_nucleosome_signal.png
#           qc/03_TSSplot.png
#           qc/05_VlnPlot_before_filter.png
#           qc/06_VlnPlot_before_after_filter.png
#           qc/07_VlnPlot_after_filter_byCelltype.png
#           ...
# Run     : Rscript scripts/05_replication/05a_qc.R   (from the repository root)
# ==============================================================================
# 05a_qc.R  — scATAC/scRNA quality control, European cohort (GSE244585 pools 9/11/12)
#
# QC STRATEGY: follows the ORIGINAL PAPER's method (Bhat-Nakshatri et al., Nat Med 2024,
# 30:3482-3494 -- the source of GSE244585). Their Methods define QC only qualitatively:
#   snRNA : exclude low-quality cells and cells with "extremely high or low" numbers of
#           detected genes and UMIs  (detected-genes filter NOT applied here -- the
#           CELLxGENE RNA is already QC-filtered upstream by the authors)
#   snATAC: discard cells with low TSS signal enrichment, and "extremely high or low"
#           numbers of reads falling in detected peaks
# They publish NO numeric cutoffs. Two readings of their criteria are COMPUTED and REPORTED
# SIDE BY SIDE, but only ONE IS EVER APPLIED:
#   APPLIED  : median +/- NMADS*MAD (log10 for count metrics, scater isOutlier-style), a
#              literal reading of the paper's "extremely high or low". Selected by
#              QC_SCHEME below, chosen after comparing retention.
#   REPORTED : the fixed thresholds of the 4-donor discovery analysis
#              (03_qc_doublets.R: ATAC 1e3-1e5, RNA <25e3, NS<1.2, TSS>1, blacklist<0.02),
#              with the RNA lower bound excluded (APPLY_RNA_MIN). Computed every run and
#              reported side by side, but not applied.
# Both schemes' pass/fail counts, the fixed-vs-MAD concordance table, and which scheme ran
# are PRINTED to logs/05_QC.log. Per-cell pass flags for both are kept in
# qc/QC_metrics_per_cell.csv (qc_pass, qc_pass_fixed, qc_pass_mad).
#
# WHY THE MAD SCHEME IS APPLIED:
#   (1) RETENTION. Fixed retained 9,713/12,469 (77.9%) even after dropping the RNA floor;
#       MAD retains 10,700 (85.8%).
#   (2) THE EXTRA CELLS ARE CLEAN, NOT NOISY. The 2,025 cells MAD keeps and fixed drops have
#       median TSS enrichment 3.83 -- ABOVE the 3.64 of cells both schemes keep. They are
#       low-DEPTH, not low-signal-to-noise, and TSS enrichment is the metric that separates
#       real nuclei from ambient DNA.
#   (3) The paper's criteria are stated RELATIVELY ("extremely high or low"), which the MAD
#       rule implements literally; the fixed numbers come from a different, deeper cohort.
# TWO CONSEQUENCES THAT DIFFER FROM DISCOVERY, both reported every run:
#   * MAD applies a percent.mt upper bound (~1,338 cells); the discovery thresholds have no
#     MT criterion. This is why immune/adipocyte retention is LOWER under MAD (T-cells 64.3%,
#     Adi-2 68.7%) than under fixed, while BM/LASP improve to 89.1%.
#   * MAD's nCount_ATAC lower bound is ~280 vs the discovery floor of 1,000.
# The fixed scheme is still computed and reported every run for comparison.
#
# Cohort-specific notes:
#   * hg38 (not hg19): blacklist_hg38; EnsDb.Hsapiens.v86 annotation already in the object
#   * autosomes chr1-22 only (chrX/Y dropped upstream in 03c; GWAS SNPs are autosomal)
#   * DOUBLET SCORING (section 1b): scDblFinder on RNA, per-pool
#     simulation (samples=pool -- the 3 pools are separate 10x lanes, so simulating a
#     doublet across lanes would model a droplet that cannot exist).
#     CRITERION: the discovery rule scDblFinder_score <= 0.95 (DOUBLET_SCORE_MAX), matching
#     03_qc_doublets.R and the ../results/singlecell/ outputs. scDblFinder's
#     own binary class call is reported alongside for contrast but is NOT the applied rule;
#     it is considerably stricter.
#     AMULET (ATAC-fragment doublet detection) is not used: it reads the whole 5.0 GB fragment
#     file and did not complete in reasonable time, and the discovery cohort had no ATAC
#     doublet detection, so omitting it keeps the two cohorts comparable.
#     REMOVAL IS OFF BY DEFAULT (DOUBLET_FILTER=FALSE) because the CELLxGENE cells are
#     already the authors' QC'd set, so a second doublet pass filters an already-filtered
#     population. Inspect qc/doublet_scores.csv and decide.
#     Cross-donor doublets cannot be assumed to have been removed by genotype demultiplexing.
#     Inspection of the h5ad obs (03b_extract_cellxgene_embeddings.py) shows
#     donor_id is a POOL-LEVEL string -- "pooled [K107181,K109128,...]", identical for every
#     cell in a pool, 3 distinct values across all 12,469 cells. The deposited object was
#     never genetically demultiplexed, so cross-donor doublets may well remain, and the
#     observed 7.3% class=="doublet" rate (evenly spread across pools) is consistent with
#     that. Treat the doublet decision on its own evidence, not on an assumed demux step.
#   * nCount_ATAC here IS "reads in peaks" (the ATAC assay is the peak x cell matrix), so it
#     maps directly onto the paper's ATAC criterion.
#   * nucleosome_signal, blacklist_fraction and percent.mt are additionally filtered (upper
#     tail only), using the same data-driven MAD rule. These go beyond the paper's stated
#     scheme but match the 4-donor original, which filtered NS and blacklist fraction.
#   * blacklist_fraction is DIAGNOSTIC ONLY here, not a filter. 04 excludes every peak
#     overlapping the ENCODE UNIFIED hg38 blacklist, so no counts remain in blacklist
#     regions and this metric is identically 0 (verified: 0 of the peaks overlap it;
#     both peak and blacklist GRanges are UCSC-style, so the 0 is real, not a join
#     failure). Artifact regions are therefore removed at the PEAK level rather than the
#     CELL level. The 4-donor original did the opposite: its blacklist subsetByOverlaps
#     was commented out, so it kept those peaks and filtered cells on blacklist_fraction
#     < 0.02. Only one of the two can do work; peak-level removal is preferred because it
#     keeps artifact regions out of chromVAR and the TF-GWAS entirely.
#
# ============================= DATA SOURCES (imports) =============================
#  peaks/replication_combined.rds  <- 04_build_combined_object.R
#      (12,469 cells; assays RNA + ATAC; 230,894 autosomal hg38 peaks; NO reductions;
#  combined/combined_fragments.tsv.gz  <- 03c_combine_fragments.sbatch
#      (needed by NucleosomeSignal / TSSEnrichment / FragmentHistogram)
# Outputs -> qc/ : figures 01-07 (PNG only, one file each),
#                      QC_metrics_per_cell.csv, doublet_scores.{csv,rds},
#                      replication_postQC.rds
# QC summaries (derived cutoffs, cells failing each criterion, retention by cell type and pool)
# are printed to the log; the per-cell values are in QC_metrics_per_cell.csv (qc_pass,
# qc_pass_fixed, qc_pass_mad).
# ==================================================================================
## ---- paths (set the environment variables; see README) ----
TNBC_ROOT <- Sys.getenv("TNBC_ROOT", ".")
suppressPackageStartupMessages({
  library(Signac); library(Seurat); library(GenomicRanges); library(GenomeInfoDb)
  library(ggplot2); library(patchwork); library(data.table)
})
setwd(file.path(TNBC_ROOT, "replication"))
dir.create("qc", showWarnings=FALSE)
set.seed(1234)
ts <- function(m) message("[", format(Sys.time(),"%H:%M:%S"), "] ", m)

## ---- QC configuration ----
## BOTH schemes are computed every run and reported side by side; QC_SCHEME selects the one
## APPLIED, and the log records which ran.
##
## WHY "mad":
## The fixed thresholds retained 9,713 / 12,469 (77.9%) even after dropping the RNA floor,
## versus 10,700 (85.8%) under MAD. Inspection of the 2,025 cells MAD keeps and fixed drops
## showed they are SHALLOW BUT CLEAN -- median TSS enrichment 3.83, ABOVE the 3.64 of cells
## both schemes keep -- so they are not low-signal-to-noise cells, just low-depth ones.
## NOTE two consequences of this choice, both of which differ from the discovery analysis:
##   (1) MAD applies a percent.mt upper bound (removes ~1,338 cells); discovery had none.
##   (2) MAD's nCount_ATAC lower bound is ~280, vs the discovery floor of 1,000.
## Both are reported in QC_derived_cutoffs.csv / QC_cutoff_comparison.csv.
QC_SCHEME         <- "mad"  # "mad"   = median +/- NMADS*MAD  (APPLIED)
                            # "fixed" = 4-donor discovery thresholds
NMADS             <- 3      # 3 MADs from the median
FILTER_NUCLEOSOME <- TRUE   # nucleosome signal filtered (as the 4-donor original did)

## fixed thresholds of the 4-donor discovery analysis (03_qc_doublets.R). Kept as named
## constants so the applied numbers appear in one place and can be cited in Methods.
FIX_ATAC_MIN <- 1e3;   FIX_ATAC_MAX <- 1e5
FIX_RNA_MIN  <- 1e3;   FIX_RNA_MAX  <- 25e3

## nCount_RNA LOWER BOUND NOT APPLIED. The discovery analysis used
## nCount_RNA > 1000, but its snRNA is deeper than this cohort's: here median nCount_RNA is
## 2,232 with q10 = 748, so a 1,000-UMI floor removed 2,208 cells -- the single largest
## contributor to the 63.9% retention seen in run 221339 -- and did so SELECTIVELY by cell
## type (T-cells 29.8% retained, Adi-2 53.0%, vs LASP 76.1%, LHS 75.3%), because small
## RNA-poor cell types cannot clear a fixed UMI floor.
## It is dropped because the CELLxGENE RNA is ALREADY QC-filtered by the authors: these
## 12,469 cells are their post-QC set, so a UMI floor here is a second pass over an
## already-filtered population. This is the same reasoning that already excludes the
## detected-genes filter (see the snRNA note at the top of this header).
## The UPPER bound is retained -- it still removes implausibly high-UMI cells.
APPLY_RNA_MIN <- FALSE
FIX_NS_MAX   <- 1.2;   FIX_TSS_MIN  <- 1;      FIX_BL_MAX <- 0.02

## FINAL POST-QC: mitochondrial fraction. The 4-donor discovery analysis applied
## percent.mt < 20 as its LAST filter (05_atac_umap_figures.R), after annotation,
## and that filter is what defines its final_cell_barcodes.csv (16,451 of 16,677 cells).
## Applied here at the equivalent point so this cohort's final cell set is defined the
## same way. Applied under BOTH schemes -- it is a post-QC step, not part of either
## threshold scheme -- so switching QC_SCHEME does not silently drop it.
FIX_MT_MAX   <- 20

## ---- doublet configuration (section 1b) ----
DOUBLET_SCORE  <- TRUE    # compute scDblFinder scores (RNA)
DOUBLET_FILTER <- FALSE   # TRUE only after inspecting qc/doublet_scores.csv (see header)
DBL_CACHE      <- "qc/doublet_scores.rds"

## DISCOVERY-MATCHED CRITERION. The 4-donor discovery analysis filtered on the scDblFinder
## SCORE, not the binary class call:
##     atac_subset <- subset(atac_subset, subset = scDblFinder_score <= 0.95)
## (hence the ../results/singlecell/ output directory). Replication uses the
## identical rule so the two cohorts' doublet handling is comparable -- same argument as
## the fixed QC thresholds. NOTE this is a PERMISSIVE cutoff: it removes only cells scDblFinder is
## ~certain about, and is far weaker than class=="doublet". Both are reported below so the
## difference is visible rather than buried.
DOUBLET_SCORE_MAX <- 0.95

## ---- pool11 structural-zero peaks (diagnostic; see section 0b) ----
## pool11's GEO fragments are truncated beyond chr9:133.35Mb, so peaks there are called
## from pool9+pool12 signal only and ALL pool11 cells read exactly 0 for TECHNICAL, not
## biological, reasons. Such peaks are perfectly pool-discriminating: TF-IDF up-weights
## them precisely BECAUSE they are sparse and pool-specific, so they load onto early LSI
## components and separate pool11 as a batch artifact -- a plausible cause of a muddy
## scATAC UMAP that the scRNA UMAP does not show.
DROP_STRUCTURAL_ZERO <- TRUE   # FALSE = report only, keep the peaks
CHR9_TRUNC_POS       <- 133350000

ts("load combined object (04)")
obj <- readRDS("peaks/replication_combined.rds")
DefaultAssay(obj) <- "ATAC"
ts(paste("loaded:", ncol(obj), "cells;", nrow(obj), "peaks"))
# loaded: 12469 cells; 230894 peaks

## peak-width filter (moved here from 04; same bounds as 03_qc_doublets.R)
PEAK_WMIN <- 20; PEAK_WMAX <- 10000
pw   <- width(granges(obj))
keep_pk <- pw > PEAK_WMIN & pw < PEAK_WMAX
ts(sprintf("peak width filter (%d < w < %d): %d of %d peaks removed",
           PEAK_WMIN, PEAK_WMAX, sum(!keep_pk), length(keep_pk)))
if (any(!keep_pk)) obj[["ATAC"]] <- subset(obj[["ATAC"]], features = rownames(obj[["ATAC"]])[keep_pk])
ts(paste("pre-QC:", ncol(obj), "cells;", nrow(obj), "peaks"))
# peak width filter (20 < w < 10000): 0 of 230894 peaks removed
# pre-QC: 12469 cells; 230894 peaks

######################################################
###  (0b) pool11 structural-zero peaks (batch)     ###
######################################################
## Identify peaks where pool11 reads EXACTLY zero while pool9+pool12 carry signal, and
## which sit in pool11's truncated region (chr9 > CHR9_TRUNC_POS). These are technical
## zeros, not biology. See the config block for why they matter to the ATAC embedding.
ts("structural-zero peak scan (pool11 truncation)")
## NOTE: GetAssayData(<Seurat>, slot=) is DEFUNCT in SeuratObject 5.4.0 (not merely
## deprecated -- it errors). The ChromatinAssay method still accepts slot=, so go through
## the assay object rather than the Seurat object.
atac_ct  <- GetAssayData(obj[["ATAC"]], slot="counts")
is_p11   <- obj$pool == "pool11"
n_p11    <- Matrix::rowSums(atac_ct[, is_p11,  drop=FALSE])
n_other  <- Matrix::rowSums(atac_ct[, !is_p11, drop=FALSE])
gr_pk    <- granges(obj[["ATAC"]])
in_trunc <- as.character(seqnames(gr_pk)) == "chr9" & start(gr_pk) > CHR9_TRUNC_POS

struct_zero <- (n_p11 == 0) & (n_other > 0)
ts(sprintf("  peaks with pool11 == 0 and pool9/12 > 0 : %d", sum(struct_zero)))
ts(sprintf("  of those, inside chr9 > %.2f Mb        : %d",
           CHR9_TRUNC_POS/1e6, sum(struct_zero & in_trunc)))
ts(sprintf("  of those, OUTSIDE the truncated region : %d  (biological or low-depth)",
           sum(struct_zero & !in_trunc)))

## Only peaks that are BOTH pool11-zero AND inside the truncated region are dropped.
## A pool11-zero peak elsewhere may be genuine biology or ordinary sparsity, so removing
## it would discard real signal.
drop_pk <- struct_zero & in_trunc
if (DROP_STRUCTURAL_ZERO && any(drop_pk)) {
  obj[["ATAC"]] <- subset(obj[["ATAC"]], features = rownames(obj[["ATAC"]])[!drop_pk])
  ts(paste("  DROPPED", sum(drop_pk), "structural-zero peaks ->", nrow(obj), "peaks remain"))
} else {
  ts(paste("  DROP_STRUCTURAL_ZERO=FALSE - reporting only;", sum(drop_pk), "peaks retained"))
}
rm(atac_ct, n_p11, n_other); invisible(gc())
# [16:08:43] structural-zero peak scan (pool11 truncation)
# [16:08:44]   peaks with pool11 == 0 and pool9/12 > 0 : 699
# [16:08:44]   of those, inside chr9 > 133.35 Mb        : 699
# [16:08:44]   of those, OUTSIDE the truncated region : 0  (biological or low-depth)
# [16:08:47]   DROPPED 699 structural-zero peaks -> 230195 peaks remain

######################################################
###            (1) compute QC metrics              ###
######################################################
ts("TSSEnrichment (fast=FALSE -> enables TSSPlot)")
obj <- TSSEnrichment(obj, fast=FALSE)
ts("NucleosomeSignal")
obj <- NucleosomeSignal(obj)
ts("blacklist_fraction (FractionCountsInRegion, hg38)")
obj$blacklist_fraction <- FractionCountsInRegion(object=obj, assay="ATAC", regions=blacklist_hg38)
mt <- grep("^MT-", rownames(obj[["RNA"]]), value=TRUE)
obj$percent.mt <- if (length(mt)) PercentageFeatureSet(obj, features=mt, assay="RNA") else 0
ts(paste("MT genes found in RNA assay:", length(mt)))
# [16:14:26] MT genes found in RNA assay: 13

######################################################
###        (1b) DOUBLET SCORING                    ###
######################################################
## scDblFinder on RNA only (see header). Scored here so the values are available to the
## distribution table, the per-cell CSV and the figures; removal is gated on DOUBLET_FILTER.
obj$scDblFinder_score <- NA_real_; obj$scDblFinder_class <- NA_character_

if (DOUBLET_SCORE) {
  if (file.exists(DBL_CACHE)) {
    dbl <- readRDS(DBL_CACHE)
    ## validate the cache before trusting it (same discipline as 04's FeatureMatrix cache,
    ## which silently produced a near-empty matrix once already)
    if (setequal(dbl$cell, colnames(obj))) {
      i <- match(colnames(obj), dbl$cell)
      obj$scDblFinder_score <- dbl$scDblFinder_score[i]
      obj$scDblFinder_class <- dbl$scDblFinder_class[i]
      ts("doublet: cache VALID")
    } else { ts("doublet: cache STALE (cells differ) - recomputing"); file.remove(DBL_CACHE) }
  }

  if (!file.exists(DBL_CACHE)) {
    ## RNA -- scDblFinder. samples=pool is REQUIRED: the 3 pools are separate 10x lanes,
    ## so simulating doublets across lanes models a droplet that cannot exist.
    ts("doublet: RNA (scDblFinder), per-pool simulation")
    if (requireNamespace("scDblFinder", quietly=TRUE) &&
        requireNamespace("SingleCellExperiment", quietly=TRUE)) {
      set.seed(1234)
      sce <- SingleCellExperiment::SingleCellExperiment(
               ## RNA is an Assay5; LayerData is the v5 accessor (slot= is defunct)
               list(counts = LayerData(obj, assay="RNA", layer="counts")))
      sce <- scDblFinder::scDblFinder(sce, samples = obj$pool)
      i <- match(colnames(obj), colnames(sce))
      obj$scDblFinder_score <- sce$scDblFinder.score[i]
      obj$scDblFinder_class <- as.character(sce$scDblFinder.class)[i]
      rm(sce); invisible(gc())
      print(table(class=obj$scDblFinder_class, pool=obj$pool, useNA="ifany"))
    } else ts("  !! scDblFinder/SingleCellExperiment not installed - RNA scoring SKIPPED")

    saveRDS(data.frame(cell=colnames(obj),
                       scDblFinder_score=obj$scDblFinder_score,
                       scDblFinder_class=obj$scDblFinder_class,
                       stringsAsFactors=FALSE), DBL_CACHE)
    ts(paste("doublet: cached ->", DBL_CACHE))
  }
}

## Two criteria, both recorded. NA-safe: a missing score must not silently become a
## singlet, so NA is treated as FAILING the test it cannot be judged by.
##  - doublet_score : DISCOVERY-MATCHED and APPLIED. score > 0.95
##                    (i.e. fails the discovery rule scDblFinder_score <= 0.95)
##  - doublet_class : scDblFinder's own binary call. Reported for contrast only; STRICTER
##                    than the discovery rule, so the two counts differ substantially.
obj$doublet_score <- is.na(obj$scDblFinder_score) | obj$scDblFinder_score > DOUBLET_SCORE_MAX
obj$doublet_class <- !is.na(obj$scDblFinder_class) & obj$scDblFinder_class == "doublet"
obj$doublet_rna   <- obj$doublet_score   # back-compat alias; the APPLIED criterion
obj$doublet_call  <- ifelse(obj$doublet_score, "doublet", "singlet")

ts("doublet: summary")
cat("  DISCOVERY RULE  scDblFinder_score >", DOUBLET_SCORE_MAX, "(APPLIED):",
      sum(obj$doublet_score), sprintf("(%.2f%%)", 100*mean(obj$doublet_score)), "\n",
    " scDblFinder class=='doublet' (reported only) :",
      sum(obj$doublet_class), sprintf("(%.2f%%)", 100*mean(obj$doublet_class)), "\n")
## how much stricter is the class call than the discovery score rule?
print(table(discovery_score_rule = obj$doublet_score, scDblFinder_class = obj$doublet_class))
if (any(!is.na(obj$scDblFinder_score))) {
  ts("scDblFinder_score distribution:")
  print(round(quantile(obj$scDblFinder_score,
                       c(0.5,0.75,0.9,0.95,0.99,0.995,0.999,1), na.rm=TRUE), 4))
}
print(table(call=obj$doublet_call, pool=obj$pool))
print(table(call=obj$doublet_call, celltype=obj$celltype))

# [16:14:26] doublet: summary
# DISCOVERY RULE  scDblFinder_score > 0.95 (APPLIED): 676 (5.42%)
# scDblFinder class=='doublet' (reported only) : 915 (7.34%)
# scDblFinder_class
# discovery_score_rule FALSE  TRUE
# FALSE 11554   239
# TRUE      0   676
# [16:14:26] scDblFinder_score distribution:
#   50%    75%    90%    95%    99%  99.5%  99.9%   100%
# pool
# call      pool11 pool12 pool9
# doublet    178    246   252
# singlet   3230   4398  4165
# celltype
# call      Adi-1 Adi-2   BM Endo-1 Endo-2 Fibroblasts LASP  LHS Macrophages
# doublet    14    13  170     48     46          70  163  116          22
# singlet   212    70 4250    439    479        1015 3165 1455         235
# celltype
# call      T-cells
# doublet      14
# singlet     473

fwrite(data.table(cell=colnames(obj), pool=obj$pool,
                  celltype=as.character(obj$celltype),
                  scDblFinder_score=obj$scDblFinder_score,
                  scDblFinder_class=obj$scDblFinder_class,
                  doublet_call=obj$doublet_call),
       "qc/doublet_scores.csv")

######################################################
###   (2) DISTRIBUTION of QC metrics (inspection)  ###
######################################################
metrics <- c("nCount_RNA","nFeature_RNA","nCount_ATAC","nFeature_ATAC",
             "TSS.enrichment","nucleosome_signal","blacklist_fraction","percent.mt",
             "scDblFinder_score")
metrics <- metrics[metrics %in% colnames(obj@meta.data)]
qs <- c(0.001,0.005,0.01,0.05,0.10,0.25,0.50,0.75,0.90,0.95,0.99,0.995,0.999)

dist_tab <- rbindlist(lapply(metrics, function(m) {
  x <- obj@meta.data[[m]]
  c(list(metric=m, n=sum(!is.na(x)), mean=mean(x,na.rm=TRUE), sd=sd(x,na.rm=TRUE),
         MAD=mad(x,na.rm=TRUE), min=min(x,na.rm=TRUE), max=max(x,na.rm=TRUE)),
    setNames(as.list(round(quantile(x, qs, na.rm=TRUE),4)), paste0("q",qs*100)))
}))
ts("QC metric distribution:"); print(as.data.frame(dist_tab), digits=4)

######################################################
###   (3) DERIVE cutoffs from the distributions    ###
######################################################
# median +/- NMADS*MAD, on log10 scale for count metrics (their "extremely high or low").
# direction: "both" (high and low are both outliers), "lower" (only a low value is bad),
#            "higher" (only a high value is bad).
mad_bounds <- function(x, nmads=NMADS, use_log=TRUE, direction="both") {
  v   <- if (use_log) log10(x + 1) else x
  med <- median(v, na.rm=TRUE); m <- mad(v, na.rm=TRUE)
  back <- function(z) if (use_log) 10^z - 1 else z
  lo <- if (direction == "higher") -Inf else back(med - nmads*m)
  hi <- if (direction == "lower")   Inf else back(med + nmads*m)
  c(lower=lo, upper=hi)
}

spec <- list(
  # paper's snRNA criterion: UMIs, extremely high or low.
  # nFeature_RNA is NOT filtered: the CELLxGENE RNA matrix is already QC-filtered by the
  # authors, so a detected-gene filter here would be a redundant second pass. Reported only.
  nCount_RNA        = list(log=TRUE,  dir="both",   apply=TRUE),
  nFeature_RNA      = list(log=TRUE,  dir="both",   apply=FALSE),
  # paper's snATAC criterion: reads in peaks extremely high or low; TSS enrichment low
  nCount_ATAC       = list(log=TRUE,  dir="both",   apply=TRUE),
  TSS.enrichment    = list(log=FALSE, dir="lower",  apply=TRUE),
  # additional ATAC/RNA quality metrics (beyond the paper's stated scheme; the 4-donor
  # original 03_qc_doublets.R also filtered on nucleosome signal and blacklist fraction).
  # Same data-driven rule: only "extremely high" values are removed.
  # NOTE: if a metric has MAD == 0 (e.g. blacklist_fraction is identically 0 because 04
  # already removed blacklist-overlapping peaks) the bounds collapse to +/-Inf and the
  # filter cannot remove any cell. That is reported explicitly in QC_derived_cutoffs.csv.
  nucleosome_signal  = list(log=FALSE, dir="higher", apply=FILTER_NUCLEOSOME),
  blacklist_fraction = list(log=FALSE, dir="higher", apply=FALSE),
  percent.mt         = list(log=FALSE, dir="higher", apply=TRUE)
)
spec <- spec[names(spec) %in% colnames(obj@meta.data)]

cut_tab <- rbindlist(lapply(names(spec), function(m) {
  s <- spec[[m]]; x <- obj@meta.data[[m]]
  b <- if (mad(x, na.rm=TRUE) == 0) c(lower=-Inf, upper=Inf) else
       mad_bounds(x, use_log=s$log, direction=s$dir)
  n_lo <- sum(x < b["lower"], na.rm=TRUE); n_hi <- sum(x > b["upper"], na.rm=TRUE)
  data.table(metric=m, scale=ifelse(s$log,"log10","linear"), direction=s$dir,
             nmads=NMADS, lower=round(b["lower"],4), upper=round(b["upper"],4),
             n_below=n_lo, n_above=n_hi, applied=s$apply)
}))
ts(paste0("derived cutoffs (median +/- ", NMADS, " MAD):")); print(as.data.frame(cut_tab))

## cutoff-label formatting, shared by figures 04 and 06
CUTLAB_SIZE <- 2.6
fmt_cut <- function(v) {
  if (!is.finite(v)) return("")
  if (abs(v) >= 1e4 || (abs(v) > 0 && abs(v) < 1e-2)) format(v, digits = 3, scientific = TRUE)
  else format(round(v, 3), big.mark = ",", trim = TRUE)
}

# distribution figure with the derived cutoffs drawn on
md <- obj@meta.data
plots <- lapply(names(spec), function(m) {
  b  <- cut_tab[metric==m]
  lg <- spec[[m]]$log
  p  <- ggplot(md, aes(x=.data[[m]])) +
        geom_histogram(bins=100, fill="grey35") +
        labs(title=paste0(m, ifelse(spec[[m]]$apply, "  [FILTER]", "  [diagnostic]")),
             x=NULL, y=NULL) + theme_bw(base_size=9)
  if (lg) p <- p + scale_x_log10()
  # cutoff lines PLUS the exact numeric value, so the figure is readable without
  # cross-referencing the log
  if (is.finite(b$lower)) p <- p +
    geom_vline(xintercept=b$lower, color="red4", linetype=2) +
    annotate("text", x=b$lower, y=Inf, label=fmt_cut(b$lower), colour="red4",
             size=CUTLAB_SIZE, hjust=-0.08, vjust=1.6)
  if (is.finite(b$upper)) p <- p +
    geom_vline(xintercept=b$upper, color="red4", linetype=2) +
    annotate("text", x=b$upper, y=Inf, label=fmt_cut(b$upper), colour="red4",
             size=CUTLAB_SIZE, hjust=1.08, vjust=1.6)
  p
})
png("qc/04_QCmetric_distributions_with_cutoffs.png", width=14, height=8, units="in", res=200)
print(wrap_plots(plots, ncol=3) +
      plot_annotation(title=paste0("QC metric distributions with data-derived cutoffs (median +/- ",
                                   NMADS, " MAD; red dashed)")))
dev.off()

######################################################
###          (4) other QC figures (pre-filter)     ###
######################################################
ts("QC figures -> qc/")
png("qc/01_DensityScatter_ATAC.png", width=9, height=7, units="in", res=200)
print(DensityScatter(obj, x="nCount_ATAC", y="TSS.enrichment", log_x=TRUE, quantiles=TRUE))
dev.off()
# The canonical Signac ATAC QC plot: sequencing depth (log₁₀ fragments in peaks)
# against TSS enrichment. Real nuclei form a dense cloud at high depth and high TSS. The bottom-left corner is ambient/debris — few fragments, no promoter enrichment. Quantile lines mark where the cutoffs bite.

# What to say: cells cluster in the high-depth/high-TSS region (median 10,600 fragments,
# TSS 3.67), confirming the library is dominated by real nuclei rather than ambient DNA.

## cutoff lines on figures 02/03 follow the applied QC scheme, so the plots agree with the
## filter that was used.
ns_hi <- if (QC_SCHEME == "fixed") FIX_NS_MAX else cut_tab[metric=="nucleosome_signal", upper]
obj$nucleosome_group <- if (is.finite(ns_hi))
    ifelse(obj$nucleosome_signal > ns_hi, "above cutoff", "below cutoff") else "all cells"
png("qc/02_histogram_nucleosome_signal.png", width=12, height=7, units="in", res=200)
print(FragmentHistogram(object=obj, group.by="nucleosome_group"))
dev.off()

tss_lo <- if (QC_SCHEME == "fixed") FIX_TSS_MIN else cut_tab[metric=="TSS.enrichment", lower]
png("qc/03_TSSplot.png", width=9, height=6, units="in", res=200)
obj$high.tss <- ifelse(obj$TSS.enrichment > tss_lo, "above cutoff", "below cutoff")
print(TSSPlot(obj, group.by="high.tss") + NoLegend())
dev.off()

## Violin features. blacklist_fraction is DELIBERATELY EXCLUDED: 04 removes blacklist peaks
## at the peak level so it is identically 0 here -- plotting it yields a flat line and four
## "All cells have the same value" warnings. percent.mt takes its slot instead: under the MAD
## scheme it is the largest single filter (1,338 cells).
qc_feats <- c("nCount_RNA","nCount_ATAC","TSS.enrichment","nucleosome_signal","percent.mt")
png("qc/05_VlnPlot_before_filter.png", width=14, height=6, units="in", res=200)
print(VlnPlot(obj, features=qc_feats, pt.size=0, ncol=5) +
      plot_annotation(title=sprintf("Pre-QC QC metrics by cell type (%d cells)", ncol(obj))))
dev.off()

######################################################
###      (5) apply the derived cutoffs             ###
######################################################
## BOTH schemes are evaluated on every run; QC_SCHEME decides which one is APPLIED.
applied <- cut_tab[applied == TRUE]
ts(paste("evaluating", nrow(applied), "data-derived (MAD) filters:",
         paste(applied$metric, collapse=", ")))
keep_mad <- rep(TRUE, ncol(obj))
for (i in seq_len(nrow(applied))) {
  m <- applied$metric[i]; x <- obj@meta.data[[m]]
  keep_mad <- keep_mad & !is.na(x) & x >= applied$lower[i] & x <= applied$upper[i]
}

## fixed thresholds of the 4-donor discovery analysis (03_qc_doublets.R).
## NOTE blacklist_fraction is identically 0 here (04 removes blacklist-overlapping peaks
## at the PEAK level), so that clause is always TRUE and does no work -- kept only so the
## expression matches the discovery script literally.
ts("evaluating fixed 4-donor discovery thresholds")
keep_fixed <- with(obj@meta.data,
  nCount_ATAC > FIX_ATAC_MIN & nCount_ATAC < FIX_ATAC_MAX &
  (if (APPLY_RNA_MIN) nCount_RNA > FIX_RNA_MIN else TRUE) & nCount_RNA < FIX_RNA_MAX &
  nucleosome_signal < FIX_NS_MAX & TSS.enrichment > FIX_TSS_MIN &
  blacklist_fraction < FIX_BL_MAX)
keep_fixed[is.na(keep_fixed)] <- FALSE

## APPLIED FILTER: selected by QC_SCHEME. The other scheme is still computed and reported.
keep <- switch(QC_SCHEME,
               mad   = keep_mad,
               fixed = keep_fixed,
               stop("QC_SCHEME must be 'mad' or 'fixed', got: ", QC_SCHEME))

## FINAL POST-QC (applied to whichever scheme was selected): mitochondrial fraction.
## See FIX_MT_MAX above for why this sits outside the scheme switch.
mt_ok  <- !is.na(obj$percent.mt) & obj$percent.mt < FIX_MT_MAX
n_mt   <- sum(keep & !mt_ok)
keep   <- keep & mt_ok
ts(paste0("FINAL POST-QC percent.mt < ", FIX_MT_MAX,
          " -> additionally removing ", n_mt, " cells"))

## doublets join the filter only if explicitly enabled (see header for why not by default)
if (DOUBLET_FILTER) {
  n_dbl <- sum(keep & obj$doublet_call == "doublet")
  keep  <- keep & obj$doublet_call == "singlet"
  ts(paste("DOUBLET_FILTER=TRUE - additionally removing", n_dbl, "called doublets"))
} else ts("DOUBLET_FILTER=FALSE - doublet scores recorded, no cells removed on that basis")

obj$qc_pass       <- keep
obj$qc_pass_mad   <- keep_mad
obj$qc_pass_fixed <- keep_fixed
ts(paste0("APPLIED SCHEME: '", QC_SCHEME, "' -> ", sum(keep), " of ", ncol(obj), " cells pass"))

# what each filter removes on its own
fail <- rbindlist(lapply(seq_len(nrow(applied)), function(i) {
  m <- applied$metric[i]; x <- obj@meta.data[[m]]
  data.table(criterion=paste0(m, " outside [", signif(applied$lower[i],4), ", ",
                              signif(applied$upper[i],4), "]"),
             n_failing=sum(is.na(x) | x < applied$lower[i] | x > applied$upper[i]))
}))
print(as.data.frame(fail))

# per-criterion failures under the FIXED thresholds (the applied scheme by default)
fixed_spec <- list(
  "nCount_ATAC outside (1e3, 1e5)"   = with(obj@meta.data, nCount_ATAC > FIX_ATAC_MIN & nCount_ATAC < FIX_ATAC_MAX),
  "nCount_RNA >= 25e3 (lower bound NOT applied)" = with(obj@meta.data, nCount_RNA < FIX_RNA_MAX),
  "[not applied] nCount_RNA <= 1e3"  = with(obj@meta.data, nCount_RNA  > FIX_RNA_MIN),
  "nucleosome_signal >= 1.2"         = with(obj@meta.data, nucleosome_signal < FIX_NS_MAX),
  "TSS.enrichment <= 1"              = with(obj@meta.data, TSS.enrichment > FIX_TSS_MIN),
  "blacklist_fraction >= 0.02"       = with(obj@meta.data, blacklist_fraction < FIX_BL_MAX),
  "percent.mt >= 20 (FINAL post-QC)" = with(obj@meta.data, percent.mt < FIX_MT_MAX))
fail_fixed <- rbindlist(lapply(names(fixed_spec), function(k)
  data.table(criterion=k, n_failing=sum(!fixed_spec[[k]] | is.na(fixed_spec[[k]])))))
ts("failing per FIXED criterion:"); print(as.data.frame(fail_fixed))

# side-by-side comparison of BOTH schemes; 'applied' marks the one actually used
cmp <- data.table(
  scheme   = c("fixed (4-donor discovery)", "data-derived MAD"),
  applied  = c(QC_SCHEME == "fixed", QC_SCHEME == "mad"),
  n_pass   = c(sum(keep_fixed), sum(keep_mad)),
  n_fail   = c(sum(!keep_fixed), sum(!keep_mad)),
  pct_kept = round(100*c(sum(keep_fixed), sum(keep_mad))/ncol(obj), 1))
ts("cutoff scheme comparison:"); print(as.data.frame(cmp))

# concordance between the two schemes -- large off-diagonals mean the choice matters
ts("scheme concordance (rows=fixed, cols=MAD):")
print(table(fixed=keep_fixed, mad=keep_mad))

obj_qc <- subset(obj, subset = qc_pass)
ts(paste("post-QC:", ncol(obj_qc), "cells retained of", ncol(obj),
         sprintf("(%.1f%%)", 100*ncol(obj_qc)/ncol(obj))))
print(table(obj_qc$celltype, useNA="ifany")); print(table(obj_qc$pool))
# [16:14:36] post-QC: 10700 cells retained of 12469 (85.8%)
#
# Adi-1       Adi-2          BM      Endo-1      Endo-2 Fibroblasts
# LASP         LHS Macrophages     T-cells
#
# pool11 pool12  pool9
# celltype pre_QC post_QC pct_retained
# 1        Adi-1    226     190         84.1
# 2        Adi-2     83      57         68.7
# 3           BM   4420    3937         89.1
# 4       Endo-1    487     392         80.5
# 5       Endo-2    525     392         74.7
# 6  Fibroblasts   1085     906         83.5
# 7         LASP   3328    2964         89.1
# 8          LHS   1571    1350         85.9
# 9  Macrophages    257     199         77.4
# 10     T-cells    487     313         64.3
# pool pre_QC post_QC pct_retained
# 1 pool11   3408    2880         84.5
# 2 pool12   4644    3914         84.3
# 3  pool9   4417    3906         88.4

ret <- merge(as.data.frame(table(celltype=obj$celltype),    responseName="pre_QC"),
             as.data.frame(table(celltype=obj_qc$celltype), responseName="post_QC"), all=TRUE)
ret$pct_retained <- round(100*ret$post_QC/ret$pre_QC, 1)
print(ret)
retp <- merge(as.data.frame(table(pool=obj$pool),    responseName="pre_QC"),
              as.data.frame(table(pool=obj_qc$pool), responseName="post_QC"), all=TRUE)
retp$pct_retained <- round(100*retp$post_QC/retp$pre_QC, 1)
print(retp)

fwrite(data.table(cell=colnames(obj), pool=obj$pool, celltype=as.character(obj$celltype),
                  nCount_RNA=obj$nCount_RNA, nFeature_RNA=obj$nFeature_RNA,
                  nCount_ATAC=obj$nCount_ATAC,
                  TSS.enrichment=obj$TSS.enrichment, nucleosome_signal=obj$nucleosome_signal,
                  blacklist_fraction=obj$blacklist_fraction, percent.mt=obj$percent.mt,
                  scDblFinder_score=obj$scDblFinder_score,
                  scDblFinder_class=obj$scDblFinder_class,
                  doublet_call=obj$doublet_call,
                  qc_pass=obj$qc_pass,
                  qc_pass_fixed=obj$qc_pass_fixed, qc_pass_mad=obj$qc_pass_mad),
       "qc/QC_metrics_per_cell.csv")

######################################################
###        (6) before/after filter figures         ###
######################################################
ts("before/after violin panel")
## Cutoff lines follow the APPLIED scheme, both bounds where finite. Previously only
## TSS.enrichment got a line, so panels for metrics that were actively filtering
## (nucleosome_signal, percent.mt under MAD) showed no threshold at all.
cut_lo <- function(m) {
  if (QC_SCHEME == "fixed")
    return(switch(m, TSS.enrichment = FIX_TSS_MIN, nCount_ATAC = FIX_ATAC_MIN,
                  nCount_RNA = if (APPLY_RNA_MIN) FIX_RNA_MIN else NA_real_, NA_real_))
  b <- cut_tab[metric == m]
  if (nrow(b) && is.finite(b$lower) && b$applied) b$lower else NA_real_
}
cut_hi <- function(m) {
  if (QC_SCHEME == "fixed")
    return(switch(m, nucleosome_signal = FIX_NS_MAX, nCount_ATAC = FIX_ATAC_MAX,
                  nCount_RNA = FIX_RNA_MAX, NA_real_))
  b <- cut_tab[metric == m]
  if (nrow(b) && is.finite(b$upper) && b$applied) b$upper else NA_real_
}
add_cuts <- function(p, m) {
  lo <- cut_lo(m); hi <- cut_hi(m)
  if (!is.na(lo)) p <- p +
    geom_hline(yintercept = lo, color = "red4", linetype = 2) +
    annotate("text", x = Inf, y = lo, label = fmt_cut(lo), colour = "red4",
             size = CUTLAB_SIZE, hjust = 1.05, vjust = -0.5)
  if (!is.na(hi)) p <- p +
    geom_hline(yintercept = hi, color = "red4", linetype = 2) +
    annotate("text", x = Inf, y = hi, label = fmt_cut(hi), colour = "red4",
             size = CUTLAB_SIZE, hjust = 1.05, vjust = 1.4)
  p
}

## SHARED Y-LIMITS for the before/after comparison.
## Limits are taken from the PRE-filter object, which is a superset of the post-filter
## cells, so nothing is clipped -- and the same metric is drawn on the same scale in both
## rows. Previously the bottom row was a single VlnPlot() with its own per-feature scaling
## and no log transform, so before/after panels were not comparable.
ylim_for <- function(m, logy) {
  v <- obj@meta.data[[m]]
  v <- v[is.finite(v)]
  if (logy) v <- v[v > 0]
  if (!length(v)) return(NULL)
  rng <- range(v)
  if (logy) c(rng[1] * 0.9, rng[2] * 1.1) else {
    pad <- diff(rng) * 0.04
    c(max(0, rng[1] - pad), rng[2] + pad)
  }
}

## one violin panel; `o` selects pre- or post-filter object, scale args are identical
vln_panel <- function(o, m, logy = FALSE, ylim = NULL, tag = "") {
  p <- VlnPlot(o, features = m, ncol = 1, pt.size = 0) +
       theme(legend.position = "none") +
       labs(title = paste0(m, tag))
  p <- if (logy) p + scale_y_log10(limits = ylim) else p + scale_y_continuous(limits = ylim)
  add_cuts(p, m)
}

png("qc/06_VlnPlot_before_after_filter.png", width=16, height=10, units="in", res=200)
## log-scale the two count metrics; identical settings used for both rows
vln_logy <- c(nCount_RNA = TRUE, nCount_ATAC = TRUE, TSS.enrichment = FALSE,
              nucleosome_signal = FALSE, percent.mt = FALSE)
vln_ylim <- lapply(qc_feats, function(m) ylim_for(m, vln_logy[[m]]))
names(vln_ylim) <- qc_feats

pre_row  <- lapply(qc_feats, function(m)
  vln_panel(obj,    m, vln_logy[[m]], vln_ylim[[m]], "  (pre-QC)"))
post_row <- lapply(qc_feats, function(m)
  vln_panel(obj_qc, m, vln_logy[[m]], vln_ylim[[m]], "  (post-QC)"))
print(wrap_plots(c(pre_row, post_row), ncol = 5, nrow = 2) +
      plot_annotation(title=sprintf(
        "Pre-QC (top) vs post-QC (bottom), IDENTICAL y-scales per metric; red dashed = applied %s cutoffs: %d -> %d cells",
        ifelse(QC_SCHEME=="mad", paste0(NMADS," MAD"), "fixed"), ncol(obj), ncol(obj_qc))))
dev.off()

png("qc/07_VlnPlot_after_filter_byCelltype.png", width=18, height=6, units="in", res=200)
print(VlnPlot(obj_qc, features=qc_feats, group.by="celltype", pt.size=0, ncol=5) +
      plot_annotation(title=sprintf("Post-QC QC metrics by cell type (%d cells, %s scheme)",
                                    ncol(obj_qc), QC_SCHEME)))
dev.off()

######################################################
###   (7) save post-QC object                      ###
######################################################
## No embeddings are computed here. RNA PCA/UMAP, ATAC LSI/UMAP and clustering are done in
## 05b_annotate_rna.R, which
## follows the discovery cohort's 04_annotate_cell_types.R (SCTransform -> PCA -> elbow
## -> clustering -> marker-based annotation), so both cohorts use one embedding recipe.
## qc/replication_postQC.rds therefore carries no reductions.
saveRDS(obj_qc, "qc/replication_postQC.rds")
ts(paste("DONE -> qc/replication_postQC.rds  (", ncol(obj_qc), "cells; no reductions )"))
ts("Next: 05b_annotate_rna.R (embeddings + clustering + RNA cell-type annotation)")
