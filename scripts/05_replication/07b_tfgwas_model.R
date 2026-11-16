#!/usr/bin/env Rscript
# ==============================================================================
# 07b_tfgwas_model.R
# Replication TF-GWAS model, one batch of motifs per array task.
# Produces: Table S7
# Run     : Rscript scripts/05_replication/07b_tfgwas_model.R   (from the repository root)
# ==============================================================================
# 07b_tfgwas_model.R — TF-GWAS linear mixed models, batched BY MOTIF
#                     European replication cohort (GSE244585)
#
# Model (per motif, per subtype):   {subtype}_chisquare ~ TF + (1 | block)
# SLURM array: SLURM_ARRAY_TASK_ID = motif batch index (1..n_batch)
#
# ============ WHY BY MOTIF AND NOT BY CHROMOSOME ============
# Each motif is fitted once on all variants genome-wide, rather than per chromosome with an
# inverse-variance-weighted combination. Batching by motif has three advantages:
#
#  1. CORRECTNESS. Fitting per chromosome and IVW-combining is an APPROXIMATION to the
#     pooled fit. It assumes the 22 chromosome-specific effects estimate one common
#     parameter with independent errors -- fine asymptotically, but it is not the model
#     the manuscript describes. Fitting each motif once on ALL variants genome-wide is the
#     model as written, and removes the meta-analysis step entirely. LD blocks are nested
#     within chromosomes, so (1 | block) is unaffected by pooling.
#  2. LOAD BALANCE. Chromosomes differ ~7-fold in variant count, so the by-chr array was
#     bottlenecked by chr1/chr2 (observed: 1-00:48 to 1-09:07 across tasks). Motifs are
#     near-interchangeable in cost, so batches finish together.
#  3. IT MATCHES DISCOVERY. The discovery model (03_tf_gwas/02_tfgwas_model.R) fits batches of
#     20 motifs against the combined genome-wide matrix (variant_peak_motif_subtypes_ALL.csv),
#     not chromosomes. The matrix itself is built per chromosome and then combined
#     (06b_motif_matrix.R), as in the discovery cohort.
#
# CONSEQUENCE FOR THE COMBINE STEP: 08_tfgwas_combine.R does NOT meta-analyse. Estimates are
# already genome-wide, so it only rbinds the batches and applies BH-FDR across motifs.
#
# ============================= DATA SOURCES (imports) =============================
#  tfgwas/lmer_batches/{annot,motifs_batch{N}}.rds  <- 07a_split_motif_batches.R  (preferred)
#  tfgwas/replication_variant_peak_motif_GWAS_ALL.csv.gz        <- 06b_motif_matrix.R (fallback if not split)
#  tfgwas/replication_motif_list.csv                   <- 06b_motif_matrix.R (motif -> TF name)
# Outputs -> tfgwas/lmer_bymotif/tfgwas_batch{N}.csv
# ==================================================================================

## ---- paths (set the environment variables; see README) ----
TNBC_ROOT <- Sys.getenv("TNBC_ROOT", ".")
suppressPackageStartupMessages({
  library(data.table); library(lme4); library(lmerTest)
})
setwd(file.path(TNBC_ROOT, "replication"))
ts <- function(m) message("[", format(Sys.time(), "%H:%M:%S"), "] ", m)
dir.create("tfgwas/lmer_bymotif", showWarnings = FALSE, recursive = TRUE)

BATCH_SIZE <- 20
BATCH_DIR  <- "tfgwas/lmer_batches"
IN_GZ      <- "tfgwas/replication_variant_peak_motif_GWAS_ALL.csv.gz"
ANNOT_COLS <- c("block",
                "TNBC_chisquare", "LuminalA_chisquare", "LuminalB_chisquare",
                "LuminalB_HER2Neg_chisquare", "HER2_chisquare")
subtypes <- c(
  TNBC             = "TNBC_chisquare",
  LuminalA         = "LuminalA_chisquare",
  LuminalB         = "LuminalB_chisquare",
  LuminalB_HER2Neg = "LuminalB_HER2Neg_chisquare",
  HER2             = "HER2_chisquare"
)

# ─────────────────────────────────────────────
#  Resolve batch index
# ─────────────────────────────────────────────
batch <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID", unset = "0"))
args  <- commandArgs(trailingOnly = TRUE)
if (length(args) >= 1) batch <- as.integer(args[1])
if (is.na(batch) || batch < 1)
  stop("No batch index. Pass SLURM_ARRAY_TASK_ID or an integer argument (1..n_batch).")

out_file <- sprintf("tfgwas/lmer_bymotif/tfgwas_batch%d.csv", batch)
ts(sprintf("=== 07b_tfgwas_model.R  BATCH=%d ===", batch))

# ─────────────────────────────────────────────
#  Load this batch's data
# ─────────────────────────────────────────────
split_ok <- file.exists(file.path(BATCH_DIR, "annot.rds")) &&
            file.exists(file.path(BATCH_DIR, sprintf("motifs_batch%d.rds", batch)))

if (split_ok) {
  ts("loading pre-split batch (07a_split_motif_batches.R output)")
  annot     <- readRDS(file.path(BATCH_DIR, "annot.rds"))
  motif_mat <- readRDS(file.path(BATCH_DIR, sprintf("motifs_batch%d.rds", batch)))
  batch_motifs <- colnames(motif_mat)
  dat <- cbind(annot, as.data.table(motif_mat))
  rm(annot, motif_mat)
} else {
  ## Fallback: read straight from the gz. Correct but slow -- fread must decompress and
  ## parse the entire file to reach this batch's ~20 columns, and every array task repeats
  ## that. Run 07a_split_motif_batches.R once to avoid it.
  ts("pre-split files absent - reading from gz (SLOW; consider running 07a_split_motif_batches.R)")
  if (!file.exists(IN_GZ)) stop("Not found: ", IN_GZ, "\n  Run 06b_motif_matrix.R first.")
  hdr        <- names(fread(cmd = paste("zcat", IN_GZ, "| head -1"), header = TRUE, nrows = 0))
  motif_cols <- grep("^MA[0-9]", hdr, value = TRUE)
  lo <- (batch - 1) * BATCH_SIZE + 1
  if (lo > length(motif_cols))
    stop("Batch ", batch, " is beyond the motif count (", length(motif_cols), ")")
  batch_motifs <- motif_cols[lo:min(batch * BATCH_SIZE, length(motif_cols))]
  dat <- fread(IN_GZ, select = c(ANNOT_COLS, batch_motifs), showProgress = FALSE)
  dat <- dat[!is.na(block)]
}
ts(sprintf("  %d variants | %d motifs in this batch", nrow(dat), length(batch_motifs)))
ts(sprintf("  LD blocks: %d", uniqueN(dat$block)))

## motif -> TF name from the discovery curated motif list, joined by motif as in
## 03_tf_gwas/02_tfgwas_model.R. 06b_motif_matrix.R filters the PFM set with the same file, so
## names and motif universe cannot drift apart.
MOTIF_ID_NAME <- "../results/singlecell/jaspar2022_motifs_727.csv"
motif_id_name <- fread(MOTIF_ID_NAME, header = TRUE)
tf_name_map   <- setNames(motif_id_name$motif.name, motif_id_name$motif)

# ─────────────────────────────────────────────
#  Fit
# ─────────────────────────────────────────────
## Failures are recorded as NA rather than aborting the batch: a motif carried by too few
## variants can produce a singular or non-converging fit, and losing the other 19 motifs in
## the batch to one bad model would be worse than an NA row that the combine step reports.
run_lmer <- function(outcome_col, tf_col, data) {
  tryCatch({
    m  <- lmer(as.formula(paste0("`", outcome_col, "` ~ `", tf_col, "` + (1 | block)")),
               data = data, control = lmerControl(optimizer = "bobyqa"))
    sm <- summary(m)$coefficients
    list(B = sm[2, 1], SE = sm[2, 2], p = sm[2, 5])
  }, error = function(e) list(B = NA_real_, SE = NA_real_, p = NA_real_))
}

## Each motif's row is appended to out_file as soon as it completes (~7 min per motif), so a
## job killed part-way keeps its finished fits and progress can be inspected.
## Start from a clean file so a re-run cannot append to a previous run's rows.
if (file.exists(out_file)) {
  ts(paste("removing existing", out_file, "(fresh run)"))
  file.remove(out_file)
}

ts(sprintf("fitting %d motifs x %d subtypes (genome-wide)",
           length(batch_motifs), length(subtypes)))
results <- lapply(seq_along(batch_motifs), function(i) {
  tf <- batch_motifs[i]
  n1 <- sum(dat[[tf]], na.rm = TRUE)
  ts(sprintf("  [%d/%d] %s  (n_motif1=%d)", i, length(batch_motifs), tf, n1))
  row <- list(batch = batch, motif = tf, TF_name = unname(tf_name_map[tf]),
              n_variants = nrow(dat), n_motif1 = n1)
  ## a motif present on no variant (or all) has no contrast -- skip the fit, record NA
  if (n1 == 0 || n1 == nrow(dat)) {
    for (sub_name in names(subtypes)) {
      row[[paste0(sub_name, "_B")]]  <- NA_real_
      row[[paste0(sub_name, "_SE")]] <- NA_real_
      row[[paste0(sub_name, "_p")]]  <- NA_real_
    }
    fwrite(as.data.table(row), out_file, append = (i > 1))
    return(row)
  }
  for (sub_name in names(subtypes)) {
    fit <- run_lmer(subtypes[[sub_name]], tf, dat)
    row[[paste0(sub_name, "_B")]]  <- fit$B
    row[[paste0(sub_name, "_SE")]] <- fit$SE
    row[[paste0(sub_name, "_p")]]  <- fit$p
  }
  ## append = (i > 1): FALSE on the first row (writes the header), TRUE after.
  ## Deliberately NOT file.exists(out_file). Motifs are ~15 min apart, which is longer
  ## than the NFS attribute-cache lifetime, so file.exists() came back FALSE on a loaded
  ## mount even though the file was there -- and fwrite then TRUNCATED instead of
  ## appending, leaving only the most recent row. The file is unlinked once above, so the
  ## loop index is an exact and filesystem-independent test for "first write".
  fwrite(as.data.table(row), out_file, append = (i > 1))
  row
})

## Rows are already on disk. Read back and verify rather than rewriting -- a rewrite here
## would mask a mismatch between what was streamed and what the loop returned.
res <- rbindlist(results, fill = TRUE)
if (!file.exists(out_file)) stop("No rows were written to ", out_file)
chk <- fread(out_file)
if (nrow(chk) != length(batch_motifs))
  stop(sprintf("wrote %d rows but batch has %d motifs -- incomplete output",
               nrow(chk), length(batch_motifs)))
if (!setequal(chk$motif, batch_motifs))
  stop("motifs on disk do not match this batch's motif set")
ts(sprintf("wrote %s (%d rows, verified)", out_file, nrow(chk)))
ts(sprintf("  fits returning NA (TNBC): %d / %d", sum(is.na(res$TNBC_p)), nrow(res)))
ts("=== 07b_tfgwas_model.R complete ===")
