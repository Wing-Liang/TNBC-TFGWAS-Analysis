#!/usr/bin/env Rscript
# ==============================================================================
# 04_build_combined_object.R
# Build the replication Seurat object (RNA and ATAC peaks).
# Produces: replication object
# Inputs  : data/cellxgene/rna_replication/
#           celltypes/cellxgene_celltypes_replication.csv
#           peaks/replication_combined.rds
# Outputs : peaks/peaks_combined.rds
#           peaks/peaks_combined.bed
#           peaks/replication_combined.rds
# Run     : Rscript scripts/05_replication/04_build_combined_object.R   (from the repository root)
# ==============================================================================
# 04_build_combined_object.R  — combined European cohort (GSE244585 pools 9/11/12)
# Mirrors 02_call_peaks_macs2.R: build ONE combined Seurat object (RNA + ATAC), then call
# peaks on the combined object's cells (matches manuscript "peak calling on all cells").
# hg38. RNA + cell types from CELLxGENE (12,469 QC-d cells). Combined fragments from 03c.
# Outputs: peaks/replication_combined.rds (RNA+ATAC, no reductions),
#          peaks_combined.{rds,bed}.  UMAPs are generated in 05a_qc.R (pre- and post-QC).
#
# ============================= DATA SOURCES (imports) =============================
#  data/cellxgene/rna_replication/                  <- 03a_extract_cellxgene_data.py
#      (CELLxGENE h5ad /raw/X; dataset 61af564d; data/cellxgene/breast_atlas.h5ad)
#  celltypes/cellxgene_celltypes_replication.csv <- 03a_extract_cellxgene_data.py
#      (CELLxGENE h5ad obs$author_cell_type)
#  combined/combined_fragments.tsv.gz       <- 03c_combine_fragments.sbatch
#      (GEO fragments GSM7821199 + data/fragments_hg38/pool11[from 02] + GSM7821201,
#       pool-prefixed, filtered to data/cellxgene/pool*_cells.txt [12,469 cells])
#  data/cellxgene/pool*_cells.txt                <- 03a_extract_cellxgene_data.py (used by 03c)
#  BSgenome.Hsapiens.UCSC.hg38, EnsDb.Hsapiens.v86 <- Bioconductor (conda: signac_env)
# ==================================================================================
## ---- paths (set the environment variables; see README) ----
TNBC_ROOT <- Sys.getenv("TNBC_ROOT", ".")
suppressPackageStartupMessages({
  library(Signac); library(Seurat); library(GenomicRanges); library(GenomeInfoDb)
  library(BSgenome.Hsapiens.UCSC.hg38); library(EnsDb.Hsapiens.v86)
  library(Matrix); library(data.table)
})
setwd(file.path(TNBC_ROOT, "replication"))
macs2 <- Sys.which("macs2")
dir.create("peaks", showWarnings=FALSE); dir.create("peaks/macs2", showWarnings=FALSE)
ts <- function(m) message("[", format(Sys.time(),"%H:%M:%S"), "] ", m)

## (1) Combined RNA (CELLxGENE; cells pool9_/pool11_/pool12_ ; 12,469)
ts("load combined RNA (CELLxGENE)")
rna_mat <- Read10X("data/cellxgene/rna_replication/")            # SOURCE: 03a_extract_cellxgene_data.py (CELLxGENE /raw/X)
obj <- CreateSeuratObject(counts=rna_mat, assay="RNA", project="European")
ct <- fread("celltypes/cellxgene_celltypes_replication.csv") # SOURCE: 03a_extract_cellxgene_data.py (obs$author_cell_type)
obj$celltype <- setNames(ct$author_cell_type, ct$cell)[colnames(obj)]
obj$pool     <- sub("_.*$","",colnames(obj))

ts(paste("RNA object:", nrow(obj), "genes x", ncol(obj), "cells"))
print(table(obj$celltype, useNA="ifany")); print(table(obj$pool))
# Adi-1       Adi-2          BM      Endo-1      Endo-2 Fibroblasts        LASP         LHS
# Macrophages     T-cells

# pool11 pool12  pool9

Idents(obj) <- "celltype"
obj_ct <- subset(obj, cells = colnames(obj)[!is.na(obj$celltype)])
Idents(obj)
# pool9_ATGGCCGGTAGCCTCC-1 pool9_ATGGCTAGTAAATTGC-1 pool9_ATGGCTAGTACAAAGA-1 pool9_ATGGCTAGTCTAGCTT-1
# BM                       BM                      LHS              Fibroblasts
# pool9_ATGGCTAGTGCGCAAT-1 pool9_ATGGCTAGTGGTTATG-1 pool9_ATGGCTAGTTTAGCTC-1 pool9_ATGGCTGTCAGCCAAT-1
# BM                   Endo-1                  T-cells              Fibroblasts
# [ reached 'max' / getOption("max.print") -- omitted 11469 entries ]

## (2) Combined ATAC fragments (03c: one file, pool-prefixed, valid CELLxGENE cells)
ts("attach combined fragment object")
frag <- CreateFragmentObject("combined/combined_fragments.tsv.gz",   # SOURCE: 03c_combine_fragments.sbatch
                             cells=colnames(obj),
                             validate.fragments=FALSE)
head(frag)
# chrom start   end                   barcode readCount
# 1  chr1 10000 10097  pool9_TCCATGCTCTCACATT-1         1
# 2  chr1 10001 10204 pool11_GCGTGCTAGCTAAGTC-1         1
# 3  chr1 10005 10013 pool11_CCGTGCTGTCACCAAA-1         1
# 4  chr1 10006 10030 pool11_CGTTTGGAGCTTGCTC-1         1
# 5  chr1 10006 10036  pool9_ATGACGAAGCTTATGA-1         1
# 6  chr1 10006 10049 pool11_TTGACATCAGCCAGAA-1         1

## (3) Peak calling on the combined cells (MACS2, Signac params)
## The MACS2 call is EXPENSIVE and is cached as raw, unfiltered peaks. The downstream
## filters are CHEAP and are re-applied on every run, so changing a filter takes effect
## without re-calling peaks (caching the *filtered* set would silently ignore such changes).
RAW_PEAKS <- "peaks/peaks_macs2_raw.rds"
if (file.exists(RAW_PEAKS)) {
  ts("load cached RAW MACS2 peaks")
  peaks_raw <- readRDS(RAW_PEAKS)
} else if (file.exists("peaks/macs2/combined_peaks.narrowPeak")) {
  ts("import RAW MACS2 peaks from peaks/macs2/combined_peaks.narrowPeak")
  xc <- c(signalValue="numeric", pValue="numeric", qValue="numeric", peak="integer")
  peaks_raw <- rtracklayer::import("peaks/macs2/combined_peaks.narrowPeak", format="BED", extraCols=xc)
  saveRDS(peaks_raw, RAW_PEAKS)
} else {
  ts("CallPeaks (MACS2) on combined cells")
  peaks_raw <- CallPeaks(frag, macs2.path=macs2, outdir="peaks/macs2", name="combined", cleanup=FALSE)
  saveRDS(peaks_raw, RAW_PEAKS)
}
ts(paste("raw MACS2 peaks:", length(peaks_raw)))

## filters (re-applied every run)
peaks <- keepSeqlevels(peaks_raw, paste0("chr",1:22), pruning.mode="coarse")
ts(paste("  after autosome chr1-22        :", length(peaks)))
# ENCODE UNIFIED hg38 blacklist (Signac data). The unified list is the ENCODE-recommended
# set for hg38 and is broader than blacklist_hg38 (910 vs 636 regions).
peaks <- subsetByOverlaps(peaks, blacklist_hg38_unified, invert=TRUE)
ts(paste("  after ENCODE unified blacklist:", length(peaks)))
# NOTE: the peak-width filter (>20 & <10000 bp) lives in 05a_qc.R
saveRDS(peaks, "peaks/peaks_combined.rds")
rtracklayer::export.bed(peaks, "peaks/peaks_combined.bed")

## (4) Quantify peaks x cells -> ChromatinAssay
# FeatureMatrix is the expensive step (~13 min). Cache it, and validate the cache against
# the current peaks and cells before using it (a stale cache would give a wrong matrix).
FM_CACHE <- "peaks/featurematrix_counts.rds"
counts <- NULL
if (file.exists(FM_CACHE)) {
  ts("found cached FeatureMatrix counts - validating against current peaks/cells")
  cc <- readRDS(FM_CACHE)
  ok <- identical(colnames(cc), colnames(obj)) &&
        nrow(cc) == length(peaks) &&
        identical(rownames(cc), GRangesToString(peaks))
  if (ok) { counts <- cc; ts(paste("cache VALID:", nrow(cc), "peaks x", ncol(cc), "cells")) }
  else    { ts("cache STALE (peaks/cells do not match) - recomputing"); rm(cc) }
}
if (is.null(counts)) {
  ts("FeatureMatrix (peaks x cells)")
  counts <- FeatureMatrix(fragments=frag,
                          features=peaks,
                          cells=colnames(obj),
                          process_n=2000)#Signac's chunk size for the fragment scan — how many peaks FeatureMatrix pulls into memory at once before counting.
  saveRDS(counts, FM_CACHE)
  ts(paste("cached ->", FM_CACHE))
}
annot <- GetGRangesFromEnsDb(EnsDb.Hsapiens.v86)   # EnsDb.Hsapiens.v86: Bioconductor
seqlevelsStyle(annot) <- "UCSC"                     # Ensembl "1" -> UCSC "chr1"
genome(annot) <- "hg38"                             # must match ChromatinAssay genome="hg38"
chrom <- CreateChromatinAssay(counts=counts,
                              fragments=frag,
                              genome="hg38",
                              annotation=annot)

## (5) Combined object: RNA + ATAC (matched cells)
obj[["ATAC"]] <- chrom
ts(paste("combined object:", ncol(obj), "cells; assays:", paste(Assays(obj), collapse=",")))
# combined object: 12469 cells; assays: RNA,ATAC
DefaultAssay(obj) <- "ATAC"

## (6) Save combined object (RNA + ATAC). Dimensionality reduction / UMAP is NOT done
## here: 04 only builds the object. UMAPs are produced in 05a_qc.R, before and after QC.
DefaultAssay(obj) <- "RNA"
saveRDS(obj, "peaks/replication_combined.rds")
ts("DONE -> peaks/replication_combined.rds  (+ peaks_combined.{rds,bed})")

obj<-readRDS("peaks/replication_combined.rds")

identical(colnames(obj[["RNA"]]), colnames(obj[["ATAC"]]))

library(Seurat); library(Signac); library(Matrix); library(data.table)
rna_counts  <- LayerData(obj, assay = "RNA", layer = "counts")
atac_counts <- obj[["ATAC"]]@counts        # ChromatinAssay: direct slot is version-safe

cells       <- colnames(obj)
rna_counts  <- LayerData(obj, assay = "RNA", layer = "counts")[, cells, drop = FALSE]
atac_counts <- LayerData(obj, assay = "ATAC", layer = "counts")[, cells, drop = FALSE]
identical(colnames(rna_counts), colnames(atac_counts))

qc <- data.table(
  cell      = colnames(obj),
  pool      = obj$pool,
  celltype  = as.character(obj$celltype),
  rna_umi   = Matrix::colSums(rna_counts),
  rna_genes = Matrix::colSums(rna_counts > 0),
  atac_frag = Matrix::colSums(atac_counts),
  atac_peak = Matrix::colSums(atac_counts > 0)
)

cat("cells total      :", nrow(qc), "\n")
# cells total      : 12469
cat("ATAC zero signal :", sum(qc$atac_frag == 0), "\n")
# ATAC zero signal : 0
cat("RNA  zero signal :", sum(qc$rna_umi   == 0), "\n")
# ATAC zero signal : 0
cat("both modalities  :", sum(qc$atac_frag > 0 & qc$rna_umi > 0), "\n")
# both modalities  : 12469

qc[, .(n = .N, atac_zero = sum(atac_frag == 0), rna_zero = sum(rna_umi == 0)), by = celltype][order(-n)]
# celltype     n atac_zero rna_zero
# <char> <int>     <int>    <int>
#   1:          BM  4420         0        0
# 2:        LASP  3328         0        0
# 3:         LHS  1571         0        0
# 4: Fibroblasts  1085         0        0
# 5:      Endo-2   525         0        0
# 6:     T-cells   487         0        0
# 7:      Endo-1   487         0        0
# 8: Macrophages   257         0        0
# 9:       Adi-1   226         0        0
# 10:       Adi-2    83         0        0

qc[, .(n = .N, atac_zero = sum(atac_frag == 0)), by = pool]
# pool     n atac_zero
# <char> <int>     <int>
#   1:  pool9  4417         0
# 2: pool11  3408         0
# 3: pool12  4644         0
