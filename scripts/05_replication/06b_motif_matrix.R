# ==============================================================================
# 06b_motif_matrix.R
# JASPAR 2022 motif annotation and the replication variant x peak x motif matrix.
# Produces: input to 07b_tfgwas_model.R
# Inputs  : qc/replication_postQC.rds
#           celltypes/cellxgene_celltypes_replication.csv
#           tfgwas/replication_variant_peak_motif_GWAS_ALL.csv.gz
# Outputs : tfgwas/replication_jaspar2022_motif_list.csv
#           tfgwas/replication_motifs_missing.csv
#           tfgwas/replication_chromvar.rds
#           tfgwas/replication_motif_list.csv
#           tfgwas/replication_matrix_QC.csv
#           tfgwas/replication_LDblock_summary.csv
#           tfgwas/replication_motif_summary_combined.csv
# Run     : Rscript scripts/05_replication/06b_motif_matrix.R   (from the repository root)
# ==============================================================================
# 06b_motif_matrix.R
# European replication (GSE244585): add motifs, run chromVAR, build GWAS×peak×motif matrix
#
# ============================= DATA SOURCES =====================================
#  qc/replication_postQC.rds                <- 05a_qc.R  (10,700 cells; 230,894 peaks)
#  ../results/gwas/bcac_gwas_ldblock_hg38filled.csv <- BCAC GWAS + LD blocks,
#      with pos38 gaps recovered by 06a_liftover_gwas_hg38.R (falls back to the original)
# Outputs -> tfgwas/
#  replication_chromvar.rds              Seurat object with motifs + chromvar assay
#  replication_motif_list.csv            motif id → JASPAR2022 TF name
#  replication_variant_peak_motif_GWAS_ALL.csv.gz all-variant × all-motif matrix (no filter)
#  replication_matrix_QC.csv             per-motif hit frequency (chr1)
#  replication_LDblock_summary.csv       LD block extents + SNP counts (lmer random-effect QC)
#  replication_motif_summary_combined.csv  per-motif GWAS variant counts + curated list + family
# ===============================================================================
## ---- paths (set the environment variables; see README) ----
TNBC_ROOT <- Sys.getenv("TNBC_ROOT", ".")
if (!require("BiocManager", quietly = TRUE))

suppressPackageStartupMessages({
  library(Seurat)
  library(Signac)
  library(GenomicRanges)
  library(BSgenome.Hsapiens.UCSC.hg38)
  library(JASPAR2022)
  library(TFBSTools)
  library(motifmatchr)
  library(chromVAR)
  library(BiocParallel)
  library(data.table)
  library(Matrix)
})

ts <- function(msg) message("[", format(Sys.time(), "%H:%M:%S"), "] ", msg)

## Print head() of each dataset as it is created, so the actual values can be eyeballed in
## the log rather than inferred from summary statistics. Cheap, and it catches build /
## chromosome-naming / column-order problems at the point they enter the pipeline.
peek <- function(x, label, n = 6) {
  cat(sprintf("
---- head: %s ----
", label))
  print(utils::head(x, n))
  if (is.data.frame(x)) cat(sprintf("   dim: %d x %d
", nrow(x), ncol(x)))
  cat("
")
}
setwd(file.path(TNBC_ROOT, "replication"))
set.seed(42)
dir.create("tfgwas", showWarnings = FALSE)

# ─────────────────────────────────────────────
#  1. Load post-QC object + cell types
# ─────────────────────────────────────────────
ts("Loading replication_postQC.rds")
obj <- readRDS("qc/replication_postQC.rds")
ts(paste("Object:", ncol(obj), "cells;", nrow(obj[["ATAC"]]), "peaks"))

if (!"celltype" %in% colnames(obj@meta.data) || all(is.na(obj$celltype))) {
  library(data.table)
  ct <- fread("celltypes/cellxgene_celltypes_replication.csv", header = TRUE)
  obj$celltype <- setNames(ct$author_cell_type, ct$cell)[colnames(obj)]
}
obj$celltype <- as.character(obj$celltype)
ts(paste("Cell types:", sum(!is.na(obj$celltype)), "of", ncol(obj), "labelled"))
print(table(obj$celltype, useNA = "ifany"))

# [16:08:29] Object: 10700 cells; 230195 peaks
# [16:08:29] Cell types: 10700 of 10700 labelled
#
# Adi-1       Adi-2          BM      Endo-1      Endo-2 Fibroblasts
# LASP         LHS Macrophages     T-cells

# ─────────────────────────────────────────────
#  2. JASPAR2022 motifs + chromVAR
# ─────────────────────────────────────────────
## MOTIF SET — MUST MATCH THE DISCOVERY COHORT
## The discovery cohort's 57 TFs were derived from this motif universe, so the replication
## must use the same one to be comparable.
##
## Discovery used:
##     getMatrixSet(JASPAR2022, opts=list(all_versions=TRUE, species='Homo sapiens'))
##   i.e. NO collection="CORE" restriction, ALL versions, species (not tax_id) filter,
##   then subset to a curated ID list that resolves duplicate base motifs to ONE version
##   and retains heterodimers.
##
## Note: collection="CORE", tax_id=9606, all_versions=FALSE gives a different, smaller
## universe and must not be used here.
MOTIF_KEEP <- "../results/singlecell/jaspar2022_motifs_727.csv"

ts("Fetching JASPAR2022 human motifs (discovery-matched: all_versions=TRUE, species filter)")
pfm_all <- getMatrixSet(JASPAR2022,
                        opts = list(all_versions = TRUE, species = "Homo sapiens"))
ts(paste("JASPAR2022 human motifs (unfiltered):", length(pfm_all)))

## Capture id / name / FAMILY for every motif before filtering. The TF-family tag is used to
## group motifs (Table S8) and is only available from the PFM tags.
motif_data_all <- data.table(
  motif      = vapply(pfm_all@listData, function(x) x@ID,   character(1)),
  motif.name = vapply(pfm_all@listData, function(x) x@name, character(1)),
  family     = vapply(pfm_all@listData, function(x)
                        if (is.null(x@tags$family)) NA_character_
                        else paste(x@tags$family, collapse = ";"), character(1)))
fwrite(motif_data_all, "tfgwas/replication_jaspar2022_motif_list.csv")
ts(paste("Full motif annotation written:", nrow(motif_data_all), "motifs"))

## Subset to the discovery cohort's curated ID list (one version per base motif, heterodimers kept)
if (!file.exists(MOTIF_KEEP))
  stop("Discovery motif list not found: ", MOTIF_KEEP,
       "\n  Required so the replication uses the SAME motif universe as the discovery run.")
motif_keep <- fread(MOTIF_KEEP, header = TRUE)
pfm <- pfm_all[names(pfm_all) %in% motif_keep$motif]
ts(sprintf("Motif set after discovery filter: %d of %d requested (%d in JASPAR2022)",
           length(pfm), nrow(motif_keep), length(pfm_all)))
missing_motifs <- setdiff(motif_keep$motif, names(pfm_all))
if (length(missing_motifs)) {
  ts(paste("  WARNING:", length(missing_motifs),
           "discovery motif IDs absent from this JASPAR2022 build -- listed in tfgwas/"))
  fwrite(data.table(motif = missing_motifs), "tfgwas/replication_motifs_missing.csv")
}
rm(pfm_all); gc()

DefaultAssay(obj) <- "ATAC"
ts("AddMotifs (hg38)")
obj <- AddMotifs(object = obj,
                 genome = BSgenome.Hsapiens.UCSC.hg38,
                 pfm    = pfm,
                 assay  = "ATAC")

ts("RunChromVAR (4 cores)")
register(MulticoreParam(4))
obj <- RunChromVAR(object = obj,
                   genome = BSgenome.Hsapiens.UCSC.hg38,
                   assay  = "ATAC")

ts("Saving replication_chromvar.rds")
saveRDS(obj, "tfgwas/replication_chromvar.rds")
ts("Saved tfgwas/replication_chromvar.rds")

# ─────────────────────────────────────────────
#  3. Extract motif info, then free the large obj from memory
# ─────────────────────────────────────────────
ts("Extracting motif info from object")
motif_mat <- GetMotifData(object = obj[["ATAC"]], slot = "data")   # sparse peaks × motifs
motif_ids <- colnames(motif_mat)
peaks_gr  <- granges(obj)
ts(paste("Motifs:", length(motif_ids), "| Peaks:", length(peaks_gr)))
# [16:50:46] Motifs: 727 | Peaks: 230195
#
# ---- head: peaks_gr (granges from the ATAC assay) ----
#   GRanges object with 6 ranges and 0 metadata columns:
#   seqnames        ranges strand
# <Rle>     <IRanges>  <Rle>
#   -------
#   seqinfo: 22 sequences from an unspecified genome; no seqlengths
#
#
# ---- head: motif_mat (peaks x motifs, first 6x6) ----
#                     MA0002.1 MA0030.1 MA0031.1 MA0051.1 MA0059.1 MA0065.1
# chr1-817072-817380    FALSE    FALSE    FALSE    FALSE    FALSE    FALSE
# chr1-817920-818121    FALSE    FALSE    FALSE    FALSE    FALSE    FALSE
# chr1-819907-820631    FALSE    FALSE    FALSE    FALSE    FALSE    FALSE
# chr1-820813-821104    FALSE    FALSE    FALSE    FALSE    FALSE    FALSE
# chr1-821674-822053    FALSE    FALSE    FALSE    FALSE    FALSE    FALSE
# chr1-826665-827758    FALSE    FALSE    FALSE    FALSE    FALSE    FALSE

peek(peaks_gr, "peaks_gr (granges from the ATAC assay)")
cat("
---- head: motif_mat (peaks x motifs, first 6x6) ----
")
print(as.matrix(motif_mat[1:6, 1:6])); cat("
")

motif_tf_names <- sapply(motif_ids, function(id) {
  m <- pfm[[id]]; if (is.null(m)) NA_character_ else m@name
})
## carry the family tag through, so downstream PWM-redundancy clustering can group
## motifs without re-querying JASPAR
motif_list <- data.table(motif = motif_ids, TF_name = motif_tf_names)
motif_list <- merge(motif_list, motif_data_all[, .(motif, family)],
                    by = "motif", all.x = TRUE, sort = FALSE)
setcolorder(motif_list, c("motif", "TF_name", "family"))
fwrite(motif_list, "tfgwas/replication_motif_list.csv")
ts(paste("Motif list saved:", nrow(motif_list), "entries;",
         sum(!is.na(motif_list$family)), "with a family tag"))
# [16:50:46] Motif list saved: 727 entries; 727 with a family tag
#
# ---- head: motif_list (motif -> TF_name -> family) ----
#   motif     TF_name
# <char>      <char>
#   1: MA0002.1       RUNX1
# 2: MA0030.1       FOXF2
# 3: MA0031.1       FOXD1
# 4: MA0051.1        IRF2
# 5: MA0059.1    MAX::MYC
# 6: MA0065.1 PPARG::RXRA
# family
# <char>
#   1:                                                       Runt-related factors
# 2:                                                                        FOX
# 3:                                                                        FOX
# 4:                                              Interferon-regulatory factors
# 5:                                                                   bHLH-ZIP
# 6: Thyroid hormone receptor-related factors (NR1);RXR-related receptors (NR2)
# dim: 727 x 3

peek(motif_list, "motif_list (motif -> TF_name -> family)")

# Free the large Seurat object before building the matrix
rm(obj, pfm); gc()
ts("Seurat object freed from memory")

# ─────────────────────────────────────────────
#  4. Load BCAC GWAS data
# ─────────────────────────────────────────────
## hg38-FILLED GWAS TABLE (06a_liftover_gwas_hg38.R)
## The original bcac_gwas_ldblock.csv has pos38 = NA for 173,263 of 9,576,981 variants
## (1.81%), which this script silently dropped. Those NAs came from a left_join against
## ld_blocks_hg38_ALL.csv, not from the variants being intrinsically unmappable, so
## 06a_liftover_gwas_hg38.R re-lifted them directly: 164,596 recovered (95.0%), 0 ambiguous,
## 8,667 genuinely absent from hg38.
##
## That script also VERIFIED the pre-existing pos38 values -- 200,000 sampled variants
## re-lifted independently, 100.0000% agreement, 0 disagreements -- so the recovered
## coordinates sit on the same footing as the originals.
##
## pos38_final = original where present, recovered otherwise; pos38_source records which.
BCAC_FILE <- "../results/gwas/bcac_gwas_ldblock_hg38filled.csv"
if (!file.exists(BCAC_FILE)) {
  BCAC_FILE <- "../results/gwas/bcac_gwas_ldblock.csv"
  ts(paste("hg38-filled table absent - falling back to", BCAC_FILE))
}
ts(paste("Loading BCAC GWAS (", basename(BCAC_FILE), ")"))
# [16:50:48] Loading BCAC GWAS ( bcac_gwas_ldblock_hg38filled.csv )

BCAC <- fread(BCAC_FILE, header = TRUE)
pos_col <- if ("pos38_final" %in% names(BCAC)) "pos38_final" else "pos38"
ts(paste("  hg38 coordinate column:", pos_col))
# [16:51:18]   hg38 coordinate column: pos38_final
# liftOver_recovered           original           unmapped

if ("pos38_source" %in% names(BCAC)) print(table(BCAC$pos38_source))
BCAC <- BCAC[!is.na(get(pos_col)) & !is.na(block)]
BCAC[, Pos38 := as.integer(get(pos_col))]
ts(paste("BCAC variants:", nrow(BCAC), "| chromosomes:", length(unique(BCAC$Chr))))
peek(BCAC, "BCAC (GWAS input)")
cat("   Chr values:", paste(head(sort(unique(BCAC$Chr)), 4), collapse=", "),
    "...  build check: Pos (hg37) vs Pos38 (hg38) both present

")

# ─────────────────────────────────────────────
#  5. Build GWAS × peak × motif matrix (chr by chr)
#     Write each chromosome directly to the gz file to avoid accumulating in RAM.
#     All motifs included — no pre-filter; as in the discovery matrix.
# ─────────────────────────────────────────────
ts("Building variant × motif matrix chromosome by chromosome (streaming to gz)")
out_gz   <- "tfgwas/replication_variant_peak_motif_GWAS_ALL.csv.gz"
chroms   <- paste0("chr", 1:22)
qc_rows  <- list()   # lightweight per-chr summary for QC
first_ch <- TRUE

for (ch in chroms) {
  bcac_ch <- BCAC[Chr == ch]
  if (nrow(bcac_ch) == 0) next

  vars_gr     <- GRanges(seqnames = ch,
                         ranges   = IRanges(start = bcac_ch$Pos38, end = bcac_ch$Pos38))
  peak_idx_ch <- which(as.character(seqnames(peaks_gr)) == ch)
  if (length(peak_idx_ch) == 0) { ts(paste(ch, ": no peaks — skip")); next }

  peaks_ch <- peaks_gr[peak_idx_ch]
  mm_ch    <- motif_mat[peak_idx_ch, , drop = FALSE]

  ov       <- findOverlaps(vars_gr, peaks_ch)
  var_peak <- integer(nrow(bcac_ch))
  var_peak[queryHits(ov)] <- subjectHits(ov)
  in_peak  <- var_peak > 0

  ## Column set matches 01_build_variant_motif_matrix.R: both the p-value and the
  ## chi-square are carried for every subtype (the motif-summary QC in section 8 reports
  ## per-motif p-value distributions).
  ## NB LuminalB's p-value column is named LuminalB_P_meta upstream, not LuminalB_pval.
  base <- data.table(
    Chr                        = bcac_ch$Chr,
    Pos                        = bcac_ch$Pos,
    Pos38                      = bcac_ch$Pos38,
    chrpos37                   = bcac_ch$chrpos37,
    block                      = bcac_ch$block,
    TNBC_pval                  = bcac_ch$TNBC_pval,
    TNBC_chisquare             = bcac_ch$TNBC_chisquare,
    LuminalA_pval              = bcac_ch$LuminalA_pval,
    LuminalA_chisquare         = bcac_ch$LuminalA_chisquare,
    LuminalB_P_meta            = bcac_ch$LuminalB_P_meta,
    LuminalB_chisquare         = bcac_ch$LuminalB_chisquare,
    LuminalB_HER2Neg_pval      = bcac_ch$LuminalB_HER2Neg_pval,
    LuminalB_HER2Neg_chisquare = bcac_ch$LuminalB_HER2Neg_chisquare,
    HER2_pval                  = bcac_ch$HER2_pval,
    HER2_chisquare             = bcac_ch$HER2_chisquare,
    Chromatin                  = as.integer(in_peak)
  )
  N_ANNOT_COLS <- ncol(base)   # used by the header check in section 7

  mm_dense  <- as.matrix(mm_ch)
  motif_bin <- matrix(0L, nrow = nrow(bcac_ch), ncol = length(motif_ids))
  colnames(motif_bin) <- motif_ids
  if (any(in_peak))
    motif_bin[in_peak, ] <- mm_dense[var_peak[in_peak], , drop = FALSE]
  rm(mm_dense); gc()

  chr_dt <- cbind(base, as.data.table(motif_bin))
  if (first_ch) {
    ## first chromosome only -- show what is actually being written, annotation columns
    ## plus a few motifs, and a row that landed inside a peak
    peek(chr_dt[, c(seq_len(N_ANNOT_COLS), N_ANNOT_COLS + 1:4), with = FALSE],
         paste0("chr_dt (", ch, ", annotation + first 4 motifs)"))
    inpk <- which(chr_dt$Chromatin == 1)
    if (length(inpk))
      peek(chr_dt[head(inpk, 3), c("Chr","Pos","Pos38","block","Chromatin",
                                   head(motif_ids, 5)), with = FALSE],
           "chr_dt rows with Chromatin == 1")
  }
  fwrite(chr_dt, out_gz, compress = "gzip", append = !first_ch, col.names = first_ch)
  first_ch <- FALSE

  ## cache this chromosome's per-motif variant counts while motif_bin is still in memory.
  ## Section 8b sums these; without it the only way to get genome-wide counts is to read the
  ## whole gz back, which is the single most expensive thing this script could do.
  saveRDS(setNames(colSums(motif_bin), motif_ids),
          sprintf("tfgwas/.motifcount_%s.rds", ch))

  qc_rows[[ch]] <- data.table(Chr = ch, n_variants = nrow(bcac_ch),
                               n_in_peak = sum(in_peak),
                               pct_in_peak = round(100 * sum(in_peak) / nrow(bcac_ch), 1))
  ts(sprintf("%s: %d variants | %d in peak (%.1f%%)",
             ch, nrow(bcac_ch), sum(in_peak), 100 * sum(in_peak) / nrow(bcac_ch)))
  rm(chr_dt, motif_bin, base); gc()
}
ts(paste("Matrix written ->", out_gz))
# [16:51:50] chr1: 739826 variants | 38987 in peak (5.3%)
# [16:52:20] chr2: 798543 variants | 38620 in peak (4.8%)
# [16:52:40] chr3: 676766 variants | 30474 in peak (4.5%)
# [16:53:04] chr4: 702519 variants | 30113 in peak (4.3%)
# [16:53:27] chr5: 605831 variants | 27709 in peak (4.6%)
# [16:53:54] chr6: 653004 variants | 36784 in peak (5.6%)
# [16:54:12] chr7: 563793 variants | 28008 in peak (5.0%)
# [16:54:33] chr8: 523019 variants | 24406 in peak (4.7%)
# [16:54:48] chr9: 406193 variants | 18364 in peak (4.5%)
# [16:55:08] chr10: 487880 variants | 25075 in peak (5.1%)
# [16:55:27] chr11: 473503 variants | 25033 in peak (5.3%)
# [16:55:46] chr12: 464019 variants | 22882 in peak (4.9%)
# [16:56:00] chr13: 353699 variants | 15903 in peak (4.5%)
# [16:56:13] chr14: 315114 variants | 15473 in peak (4.9%)
# [16:56:25] chr15: 271191 variants | 14259 in peak (5.3%)
# [16:56:37] chr16: 295935 variants | 14696 in peak (5.0%)
# [16:56:49] chr17: 259451 variants | 16723 in peak (6.4%)
# [16:57:01] chr18: 272032 variants | 13825 in peak (5.1%)
# [16:57:12] chr19: 227304 variants | 16729 in peak (7.4%)
# [16:57:22] chr20: 210171 variants | 12484 in peak (5.9%)
# [16:57:31] chr21: 132453 variants | 7329 in peak (5.5%)
# [16:57:40] chr22: 136068 variants | 9421 in peak (6.9%)

# [16:57:43] Matrix written -> tfgwas/replication_variant_peak_motif_GWAS_ALL.csv.gz

# ─────────────────────────────────────────────
#  6. QC (summary only — spot-check needs object reload, skipped here)
# ─────────────────────────────────────────────
ts("=== MATRIX QC ===")
qc_summary <- rbindlist(qc_rows)
n_total  <- sum(qc_summary$n_variants)
n_inpeak <- sum(qc_summary$n_in_peak)
ts(sprintf("Variants: %d total | %d in peak (%.1f%%)", n_total, n_inpeak, 100*n_inpeak/n_total))
ts("Per-chromosome summary:")
print(qc_summary)

# Per-motif frequency (read back one chr to estimate, then scale)
ts("Per-motif QC: reading back chr1 to compute hit frequency")
## BUGFIX: this awk had no -F, so it split on WHITESPACE, not commas. In a CSV
## that makes $1 the entire line, which never equals "chr1", so only the header survived and
## every motif count came back 0 -- the log read "Motifs with >=1 hit in chr1: 0 / 727" while
## section 8 correctly reported 727/727 from the in-memory counts. -F, fixes the split.
chr1_dt    <- fread(cmd = paste("zcat", out_gz, "| awk -F, 'NR==1 || $1==\"chr1\"'"),
                    header = TRUE, showProgress = FALSE)
motif_hits <- colSums(chr1_dt[, ..motif_ids])
motif_qc   <- data.table(motif    = motif_ids,
                          TF_name  = motif_tf_names,
                          n_hits_chr1 = as.integer(motif_hits))[order(-n_hits_chr1)]
ts(sprintf("Motifs with ≥1 hit in chr1: %d / %d",
           sum(motif_qc$n_hits_chr1 > 0), nrow(motif_qc)))
fwrite(motif_qc, "tfgwas/replication_matrix_QC.csv")
ts("QC -> tfgwas/replication_matrix_QC.csv")

# ─────────────────────────────────────────────
#  7. Finalise
# ─────────────────────────────────────────────
## The matrix was streamed to out_gz chromosome by chromosome in section 5 (to avoid holding
## it in memory); verify what was written.
if (!file.exists(out_gz)) stop("Expected streamed matrix not found: ", out_gz)
ts(sprintf("Matrix on disk: %s (%.1f MB)", out_gz, file.size(out_gz) / 1024^2))
hdr        <- fread(cmd = paste("zcat", out_gz, "| head -1"), header = TRUE, nrows = 0)
motif_cols <- grep("^MA[0-9]", names(hdr), value = TRUE)   # discovery's motif-column pattern
ts(sprintf("Columns: %d = %d annotation + %d motif (^MA[0-9])",
           ncol(hdr), ncol(hdr) - length(motif_cols), length(motif_cols)))
if (length(motif_cols) != length(motif_ids))
  ts(sprintf("  !! motif column count %d != motif set size %d - investigate",
             length(motif_cols), length(motif_ids)))
ts(paste("Rows written (variants):", n_total))

full_dataset<-fread("tfgwas/replication_variant_peak_motif_GWAS_ALL.csv.gz")

# ─────────────────────────────────────────────
#  8. LD-block summary + per-motif summary
# ─────────────────────────────────────────────
## Ports the summary sections of 01_build_variant_motif_matrix.R.
##
## NOT PORTED: that script's DA-motif integration, i.e.
## FindAllMarkers on the chromvar assay to get differential motif ACTIVITY per cell type
## (discovery's Motif_DiffAct_general_cluster.csv, 01_singlecell_annotation/12_motif_annotation_chromvar.R) and the
## intersection of those DA motifs with GWAS-variant motif counts.
## If it is wanted later: run FindAllMarkers(assay="chromvar", only.pos=TRUE,
## mean.fxn=rowMeans, test.use="LR", fc.name="avg_diff", latent.vars="nCount_ATAC") on
## replication_chromvar.rds, filter p_val_adj < 0.05 & avg_diff > 0, and merge on `motif` with
## replication_motif_summary_combined.csv. It needs no rerun of this script.
## (a) LD block extents. Discovery writes this as results/tfgwas/LDblock_summary.csv.
##     The lmer model in 07b_tfgwas_model.R uses `block` as the random effect, so block size/extent is a
##     direct diagnostic on that model: a block with one SNP contributes no within-block
##     variance, and a pathologically wide block signals a liftover problem.
ts("LD block summary")
ld_summary <- BCAC[, .(min_pos37 = min(Pos,   na.rm = TRUE),
                       max_pos37 = max(Pos,   na.rm = TRUE),
                       min_pos38 = min(Pos38, na.rm = TRUE),
                       max_pos38 = max(Pos38, na.rm = TRUE),
                       snp_count = .N), by = .(Chr, block)]
peek(ld_summary, "ld_summary (LD block extents)")
fwrite(ld_summary, "tfgwas/replication_LDblock_summary.csv")
ts(sprintf("  blocks: %d | SNPs/block median %.0f, min %d, max %d | singleton blocks: %d",
           nrow(ld_summary), median(ld_summary$snp_count),
           min(ld_summary$snp_count), max(ld_summary$snp_count),
           sum(ld_summary$snp_count == 1)))

## (b) Per-motif summary joined to the curated motif list, as discovery's
##     motif_summary_combined.csv. Counts how many GWAS variants carry each
##     motif, genome-wide rather than chr1-only.
ts("Per-motif variant counts (all chromosomes)")
motif_counts <- integer(length(motif_ids)); names(motif_counts) <- motif_ids
for (ch in chroms) {
  f <- sprintf("tfgwas/.motifcount_%s.rds", ch)
  if (file.exists(f)) motif_counts <- motif_counts + readRDS(f)
}
if (all(motif_counts == 0)) {
  ts("  no per-chromosome counts cached - falling back to a single streamed pass")
  con <- gzfile(out_gz, "r"); hdr_line <- strsplit(readLines(con, 1), ",")[[1]]
  idx <- match(motif_ids, hdr_line); close(con)
  full <- fread(out_gz, select = idx, showProgress = FALSE)
  motif_counts <- colSums(full); names(motif_counts) <- motif_ids
  rm(full); gc()
}
motif_summary <- data.table(motif = motif_ids, Count = as.integer(motif_counts))
motif_summary <- merge(motif_keep, motif_summary, by = "motif", all.x = TRUE, sort = FALSE)
motif_summary <- merge(motif_summary, motif_data_all[, .(motif, family)],
                       by = "motif", all.x = TRUE, sort = FALSE)
setorder(motif_summary, -Count, na.last = TRUE)
peek(motif_summary, "motif_summary (per-motif GWAS variant counts)")
fwrite(motif_summary, "tfgwas/replication_motif_summary_combined.csv")
ts(sprintf("  motifs with >=1 GWAS variant: %d / %d",
           sum(motif_summary$Count > 0, na.rm = TRUE), nrow(motif_summary)))

# [16:57:43] Per-motif QC: reading back chr1 to compute hit frequency
# [17:00:15] Motifs with ≥1 hit in chr1: 727 / 727
# [17:00:15] QC -> tfgwas/replication_matrix_QC.csv
# [17:00:15] Matrix on disk: tfgwas/replication_variant_peak_motif_GWAS_ALL.csv.gz (1155.5 MB)
# [17:00:15] Columns: 743 = 16 annotation + 727 motif (^MA[0-9])
# [17:00:15] Rows written (variants): 9568314
# [17:00:15] LD block summary
# ---- head: ld_summary (LD block extents) ----
#   Chr    block min_pos37 max_pos37 min_pos38 max_pos38 snp_count
# <char>   <char>     <int>     <int>     <int>     <int>     <int>
#   1:  chr10  chr10_7   8775169  10544178   8733206  10502215      6721
# 2:  chr10 chr10_50  97066807 100082385  95307050  98322628      9508
# 3:  chr10 chr10_51 100082534 102090924  98322777 100331167      6427
# 4:  chr10  chr10_1     60684   1747614     14744   1705420      6915
# 5:  chr10 chr10_52 102091470 106139294 100331713 104379536     10398
# 6:  chr10  chr10_8  10544501  12859007  10502538  12817008      9476
# dim: 1336 x 7
#
# [17:00:15]   blocks: 1336 | SNPs/block median 6914, min 612, max 21528 | singleton blocks: 0
# [17:00:15] Per-motif variant counts (all chromosomes)
#
# ---- head: motif_summary (per-motif GWAS variant counts) ----
#   motif motif.name base_motif_ids motif_version  Count
# <char>     <char>         <char>         <int>  <int>
#   1: MA0130.1    ZNF354C         MA0130             1 273060
# 2: MA1723.1      PRDM9         MA1723             1 211184
# 3: MA0149.1 EWSR1-FLI1         MA0149             1 176724
# 4: MA1596.1     ZNF460         MA1596             1 176236
# 5: MA1630.2     ZNF281         MA1630             2 162709
# 6: MA0080.5       SPI1         MA0080             5 162204
# family
# <char>
#   1:            More than 3 adjacent zinc fingers
# 2: Factors with multiple dispersed zinc fingers
# 3:                                  Ets-related
# 4:            More than 3 adjacent zinc fingers
# 5:            More than 3 adjacent zinc fingers
# 6:                                  Ets-related
# dim: 727 x 6
#
# [17:00:15]   motifs with >=1 GWAS variant: 727 / 727

## clean up the per-chromosome count caches
unlink(sprintf("tfgwas/.motifcount_%s.rds", chroms))

ts("DONE -> tfgwas/replication_variant_peak_motif_GWAS_ALL.csv.gz")
ts("Outputs: replication_chromvar.rds, replication_motif_list.csv,")
ts("         replication_jaspar2022_motif_list.csv, replication_variant_peak_motif_GWAS_ALL.csv.gz,")
ts("         replication_matrix_QC.csv, replication_LDblock_summary.csv,")
ts("         replication_motif_summary_combined.csv")
ts("=== 06b_motif_matrix.R complete ===")

#
# [17:00:15] DONE -> tfgwas/replication_variant_peak_motif_GWAS_ALL.csv.gz
# [17:00:15] Outputs: replication_chromvar.rds, replication_motif_list.csv,
# [17:00:15]          replication_jaspar2022_motif_list.csv, replication_variant_peak_motif_GWAS_ALL.csv.gz,
# [17:00:15]          replication_matrix_QC.csv, replication_LDblock_summary.csv,
# [17:00:15]          replication_motif_summary_combined.csv
