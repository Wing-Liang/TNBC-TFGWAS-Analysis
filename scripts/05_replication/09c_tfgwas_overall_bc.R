#!/usr/bin/env Rscript
# ==============================================================================
# 09c_tfgwas_overall_bc.R
# TF-GWAS model for overall breast cancer (replication; discovery with OVERALLBC_COHORT=discovery).
# Produces: Tables S4 and S7 (overall BC)
# Run     : Rscript scripts/05_replication/09c_tfgwas_overall_bc.R   (from the repository root)
# ==============================================================================
# 09c_tfgwas_overall_bc.R — overall breast cancer (overall BC) TF-GWAS, batched BY MOTIF.
#                    European replication cohort (GSE244585).
#
# Model (per motif):  OverallBC_chisquare ~ TF + (1 | block)
# SLURM array: SLURM_ARRAY_TASK_ID = motif batch index (1..37)
#
# BY MOTIF, NOT BY CHROMOSOME — same reasoning as 07b_tfgwas_model.R: fitting per chromosome
# and inverse-variance-combining is only an approximation to the pooled fit, is badly load
# balanced, and does not match the discovery design (02_tfgwas_model.R fits batches of 20 motifs
# against the combined genome-wide matrix). Each motif is fitted ONCE on all variants.
#
# REUSES the existing motifs_batch{N}.rds unchanged, so overall BC costs one outcome rather
# than re-fitting the five subtypes. Only the annotation differs (annot_OverallBC.rds), whose
# row order is asserted identical to annot.rds by 09a_overall_bc_annotation.R.
#
# COVERAGE: overall BC is available for 84.4% of variants (the BCAC overall-BC file is hg38 and
# multi-allelic positions are excluded as unresolvable). Rows without it are dropped for
# this outcome only; all 1,336 LD blocks still carry overall BC variants, so the random effect
# is unaffected. Subtype estimates are untouched and still use the full variant set, so
# overall BC standard errors are mildly larger than the subtypes' -- note when comparing.
#
# ============================= DATA SOURCES (imports) =============================
#  tfgwas/lmer_batches/annot_OverallBC.rds       <- 09a_overall_bc_annotation.R (block + overall BC chisq)
#  tfgwas/lmer_batches/motifs_batch{N}.rds   <- 07a_split_motif_batches.R (unchanged)
#  ../results/singlecell/jaspar2022_motifs_727.csv
# Outputs -> tfgwas/lmer_OverallBC/tfgwas_OverallBC_batch{N}.csv
# ==================================================================================
## ---- paths (set the environment variables; see README) ----
TNBC_ROOT <- Sys.getenv("TNBC_ROOT", ".")
suppressPackageStartupMessages({ library(data.table); library(lme4); library(lmerTest) })
setwd(file.path(TNBC_ROOT, "replication"))
ts <- function(m) message("[", format(Sys.time(), "%H:%M:%S"), "] ", m)
dir.create("tfgwas/lmer_OverallBC", showWarnings = FALSE, recursive = TRUE)

## COHORT-parameterised so the SAME model/code runs on both datasets (env OVERALLBC_COHORT).
##   european  (default) -> batches from 07a_split_motif_batches.R + annot from 09a_overall_bc_annotation.R
##   discovery           -> batches + annot from 09b_split_discovery_overall_bc.R
COHORT    <- Sys.getenv("OVERALLBC_COHORT", unset = "european")
BATCH_DIR <- if (COHORT == "discovery") "tfgwas/discovery_OverallBC_batches" else "tfgwas/lmer_batches"
OUT_DIR   <- if (COHORT == "discovery") "tfgwas/lmer_OverallBC_discovery"    else "tfgwas/lmer_OverallBC"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
OUTCOME   <- "OverallBC_chisquare"

batch <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID", unset = "0"))
args  <- commandArgs(trailingOnly = TRUE)
if (length(args) >= 1) batch <- as.integer(args[1])
if (is.na(batch) || batch < 1) stop("No batch index (SLURM_ARRAY_TASK_ID or argument).")
out_file <- sprintf("%s/tfgwas_OverallBC_batch%d.csv", OUT_DIR, batch)
ts(sprintf("=== 09c_tfgwas_overall_bc.R  COHORT=%s  BATCH=%d ===", COHORT, batch))

annot <- readRDS(file.path(BATCH_DIR, "annot_OverallBC.rds"))
mm    <- readRDS(file.path(BATCH_DIR, sprintf("motifs_batch%d.rds", batch)))
if (nrow(annot) != nrow(mm))
  stop("Row mismatch: annot_OverallBC.rds ", nrow(annot), " vs motif batch ", nrow(mm))
batch_motifs <- colnames(mm)
dat <- cbind(annot[, .(block, outcome = get(OUTCOME))], as.data.table(mm))

n_all <- nrow(dat)
dat   <- dat[!is.na(outcome)]                    # overall BC missing -> unusable for this model
ts(sprintf("  variants: %s of %s with overall BC (%.1f%%) | LD blocks: %d | motifs: %d",
           format(nrow(dat), big.mark = ","), format(n_all, big.mark = ","),
           100 * nrow(dat) / n_all, uniqueN(dat$block), length(batch_motifs)))

MOTIF_ID_NAME <- "../results/singlecell/jaspar2022_motifs_727.csv"
tf_name_map <- { x <- fread(MOTIF_ID_NAME, header = TRUE); setNames(x$motif.name, x$motif) }

## failures recorded as NA rather than aborting: one singular fit must not cost the
## other 19 motifs in the batch (same policy as 07b_tfgwas_model.R)
run_lmer <- function(tf_col, data) {
  tryCatch({
    m  <- lmer(as.formula(paste0("outcome ~ `", tf_col, "` + (1 | block)")),
               data = data, control = lmerControl(optimizer = "bobyqa"))
    sm <- summary(m)$coefficients
    list(B = sm[2, 1], SE = sm[2, 2], p = sm[2, 5])
  }, error = function(e) list(B = NA_real_, SE = NA_real_, p = NA_real_))
}

## fresh file, then append per completed motif, so a job killed part-way keeps its finished
## fits (as in 07b_tfgwas_model.R)
if (file.exists(out_file)) { ts(paste("removing existing", out_file)); file.remove(out_file) }

for (i in seq_along(batch_motifs)) {
  tf <- batch_motifs[i]
  n1 <- sum(dat[[tf]], na.rm = TRUE)
  ts(sprintf("  [%d/%d] %s  (n_motif1=%d)", i, length(batch_motifs), tf, n1))
  row <- list(batch = batch, motif = tf, TF_name = unname(tf_name_map[tf]),
              n_variants = nrow(dat), n_motif1 = n1,
              OverallBC_B = NA_real_, OverallBC_SE = NA_real_, OverallBC_p = NA_real_)
  if (n1 > 0 && n1 < nrow(dat)) {          # no contrast -> leave NA
    fit <- run_lmer(tf, dat)
    row$OverallBC_B <- fit$B; row$OverallBC_SE <- fit$SE; row$OverallBC_p <- fit$p
  }
  fwrite(as.data.table(row), out_file, append = (i > 1))
}

chk <- fread(out_file)
if (nrow(chk) != length(batch_motifs))
  stop("Wrote ", nrow(chk), " rows but expected ", length(batch_motifs))
ts(sprintf("DONE -> %s (%d motifs)", out_file, nrow(chk)))
