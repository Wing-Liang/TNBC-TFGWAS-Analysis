## 02_tfgwas_model.R
## TF-GWAS linear mixed model, discovery cohort (Table S4, Fig 3A).
##
## For each of the 727 JASPAR 2022 motifs, and one breast-cancer subtype at a time:
##
##     <subtype>_chisquare ~ motif_in_peak + (1 | LD block)          (lme4 / lmerTest)
##
## where motif_in_peak = 1 if the variant lies in an accessible peak carrying the motif.
## Motifs can be fitted in batches of 20 (e.g. as a SLURM array) and merged afterwards.
##
## Usage (from the repository root):
##   Rscript scripts/03_tf_gwas/02_tfgwas_model.R TNBC                 all motifs, one run
##   Rscript scripts/03_tf_gwas/02_tfgwas_model.R TNBC 3               batch 3 (motifs 41-60)
##   Rscript scripts/03_tf_gwas/02_tfgwas_model.R TNBC combine         merge the batch files
## Subtypes: TNBC, LumA, LumB, LumB_HER2Neg, HER2
##
## Inputs   results/tfgwas/variant_peak_motif_subtypes_ALL.csv                  (01_build_variant_motif_matrix.R)
##          results/singlecell/jaspar2022_motifs_727.csv
## Outputs  results/tfgwas/tfgwas_results_<subtype>_<batch>.csv                 (batch mode)
##          results/tfgwas/tfgwas_results_<subtype>_ALL.csv                     (all motifs / combine)
##
## Note: this script consolidates the per-batch scripts used for the published analysis; the
## model, the extracted coefficients and the output columns are unchanged.

TNBC_ROOT <- Sys.getenv("TNBC_ROOT", ".")

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(lme4)
  library(lmerTest)
})

PHENOTYPE <- c(TNBC = "TNBC_chisquare", LumA = "LuminalA_chisquare", LumB = "LuminalB_chisquare",
               LumB_HER2Neg = "LuminalB_HER2Neg_chisquare", HER2 = "HER2_chisquare")
BATCH_SIZE <- 20

args    <- commandArgs(trailingOnly = TRUE)
subtype <- if (length(args) >= 1) args[1] else stop("give a subtype: ", paste(names(PHENOTYPE), collapse = ", "))
if (!subtype %in% names(PHENOTYPE)) stop("unknown subtype '", subtype, "'")
mode    <- if (length(args) >= 2) args[2] else "all"

setwd(TNBC_ROOT)
OUT <- "results/tfgwas"
out_file <- function(tag) file.path(OUT, sprintf("tfgwas_results_%s_%s.csv", subtype, tag))

## ---- merge batch files ----
if (mode == "combine") {
  files <- list.files(OUT, pattern = sprintf("^tfgwas_results_%s_[0-9]+\\.csv$", subtype), full.names = TRUE)
  if (!length(files)) stop("no batch files for ", subtype)
  res <- rbindlist(lapply(files, fread), use.names = TRUE)
  res <- unique(res, by = "TF")
  fwrite(res, out_file("ALL"))
  cat("combined", length(files), "batches ->", nrow(res), "motifs ->", out_file("ALL"), "\n")
  quit(save = "no")
}

## ---- data ----
motif_id_name <- fread("results/singlecell/jaspar2022_motifs_727.csv", header = TRUE)
dataset <- fread("results/tfgwas/variant_peak_motif_subtypes_ALL.csv", header = TRUE)

tf_cols <- grep("MA", names(dataset), value = TRUE)
tf_cols <- tf_cols[!grepl("Chromatin", tf_cols)]          # the accessibility indicator is not a motif

if (mode == "all") {
  todo <- tf_cols; tag <- "ALL"
} else {
  b <- as.integer(mode)
  idx <- ((b - 1) * BATCH_SIZE + 1):min(b * BATCH_SIZE, length(tf_cols))
  if (idx[1] > length(tf_cols)) stop("batch ", b, " is beyond the ", length(tf_cols), " motifs")
  todo <- tf_cols[idx]; tag <- as.character(b)
}
response <- PHENOTYPE[[subtype]]
cat(sprintf("%s: fitting %d motif(s), response %s\n", subtype, length(todo), response))

## ---- model per motif ----
model1_results <- rbindlist(lapply(todo, function(tf) {
  fit <- lmer(as.formula(sprintf("%s ~ get(tf) + (1 | block)", response)), data = dataset)
  coef_table <- summary(fit)$coefficients
  data.table(
    TF                = tf,
    motif             = sub("_.*", "", tf),
    Intercept_B       = coef_table[1, 1],
    Intercept_SE      = coef_table[1, 2],
    Intercept_p_value = coef_table[1, 5],
    TF_B              = coef_table[2, 1],
    TF_SE             = coef_table[2, 2],
    TF_p_value        = coef_table[2, 5]          # Pr(>|t|), Satterthwaite (lmerTest)
  )
}))
model1_results <- left_join(model1_results, motif_id_name, by = "motif")

dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
fwrite(model1_results, out_file(tag))
cat("wrote", out_file(tag), "\n")
