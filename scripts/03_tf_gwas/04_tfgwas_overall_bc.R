# ==============================================================================
# 04_tfgwas_overall_bc.R
# TF-GWAS model for overall breast cancer, discovery cohort.
# Produces: Table S4 (overall BC)
# Run     : Rscript scripts/03_tf_gwas/04_tfgwas_overall_bc.R   (from the repository root)
# ==============================================================================
# 04_tfgwas_overall_bc.R
# Overall breast cancer (overall BC) TF-GWAS model, discovery cohort, using the BCAC overall
# breast-cancer association statistics (Zhang et al. 2020, Nat Genet). Table S4, overall BC columns.
#
# All coordinates: hg19 (original ATAC peaks + GWAS are both hg19 — no liftover).
#
# Inputs:
#   results/tfgwas/variant_peak_motif_subtypes_ALL.csv
#   $BCAC_DIR/Overall_BC_Zhang2020_clean.txt
#   results/singlecell/jaspar2022_motifs_727.csv
#   results/tfgwas/tfgwas_significant_TNBC_unique.csv
#
# Outputs (results/tfgwas/):
#   OverallBC_variant_peak_motif_dataset.csv.gz
#   OverallBC_tfgwas_allMotifs_results.csv
#   OverallBC_tfgwas_significant.csv
#   OverallBC_vs_TNBC_concordance.csv

## ---- paths (set the environment variables; see README) ----
BCAC_DIR <- Sys.getenv("BCAC_DIR", "data/BCAC")
TNBC_ROOT <- Sys.getenv("TNBC_ROOT", ".")
suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(lme4)
  library(lmerTest)
})

ts <- function(msg) message("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] ", msg)
setwd(TNBC_ROOT)
outdir <- "results/singlecell/GWAS"

ts("Loading peak-motif-GWAS dataset")
dataset <- fread(file.path(outdir, "variant_peak_motif_subtypes_ALL.csv"),
                 header = TRUE)
ts(paste("Dataset rows:", nrow(dataset), "cols:", ncol(dataset)))

ts("Loading overall BC GWAS (Zhang 2020 overall, hg19)")
overall_bc_path <- path.expand(
  file.path(BCAC_DIR, "Overall_BC_Zhang2020_clean.txt"))
overall_bc <- fread(overall_bc_path, header = TRUE, nThread = 4)
ts(paste("Overall BC rows:", nrow(overall_bc)))

overall_bc[, chrpos37 := paste0("chr", chr.iCOGs, ":", Position.iCOGs)]
overall_bc_select <- unique(overall_bc[!is.na(chi.meta),
                              .(chrpos37, OverallBC_chisquare = chi.meta,
                                OverallBC_pval = p.meta)])
ts(paste("Overall BC unique positions:", nrow(overall_bc_select)))

ts("Merging with peak-motif table")
merged <- left_join(as.data.frame(dataset),
                    as.data.frame(overall_bc_select), by = "chrpos37")
ts(paste("Variants matched to overall BC:", sum(!is.na(merged$OverallBC_chisquare))))
fwrite(merged, file.path(outdir, "OverallBC_variant_peak_motif_dataset.csv.gz"),
       compress = "gzip")

motif_cols <- grep("^MA", names(merged), value = TRUE)
ts(paste("Motif columns:", length(motif_cols)))

motif_id_map <- fread(
  "results/singlecell/jaspar2022_motifs_727.csv",
  header = TRUE)
clean_name <- function(x) toupper(sub("::.*", "", sub("\\.[0-9]+$", "", x)))
motif_id_map$TF_clean <- clean_name(motif_id_map$motif.name)

ts("Running overall BC TF-GWAS lmer model")
dat_model <- merged[!is.na(merged$block) & !is.na(merged$OverallBC_chisquare), ]
ts(paste("Model dataset rows:", nrow(dat_model)))

overall_bc_results <- lapply(motif_cols, function(tf) {
  tryCatch({
    m  <- lmer(OverallBC_chisquare ~ get(tf) + (1 | block),
               data = dat_model,
               control = lmerControl(optimizer = "bobyqa"))
    sm <- summary(m)$coefficients
    data.frame(motif = tf, TF_B = sm[2,1], TF_SE = sm[2,2],
               TF_t = sm[2,4], TF_p = sm[2,5], stringsAsFactors = FALSE)
  }, error = function(e) {
    data.frame(motif = tf, TF_B = NA, TF_SE = NA, TF_t = NA, TF_p = NA,
               stringsAsFactors = FALSE)
  })
})
overall_bc_results <- rbindlist(overall_bc_results)

strip_ver <- function(x) sub("\\.[0-9]+$", "", x)
overall_bc_results$motif_base <- strip_ver(overall_bc_results$motif)
motif_meta <- unique(motif_id_map[, .(motif_base = strip_ver(motif),
                                      motif.name, TF_clean)])
overall_bc_results <- left_join(as.data.frame(overall_bc_results),
                           as.data.frame(motif_meta), by = "motif_base")
overall_bc_results$FDR <- p.adjust(overall_bc_results$TF_p, method = "BH")
overall_bc_results <- overall_bc_results[order(overall_bc_results$TF_p, na.last = TRUE), ]
fwrite(overall_bc_results, file.path(outdir, "OverallBC_tfgwas_allMotifs_results.csv"))

sig_overall_bc <- overall_bc_results[!is.na(overall_bc_results$FDR) & overall_bc_results$FDR < 0.05, ]
ts(paste("Overall BC significant TFs (FDR<0.05):", nrow(sig_overall_bc)))
fwrite(sig_overall_bc, file.path(outdir, "OverallBC_tfgwas_significant.csv"))

tnbc57 <- fread(file.path(outdir, "tfgwas_significant_TNBC_unique.csv"), header = TRUE)
tnbc57$TF_clean <- clean_name(tnbc57$motif.name)
concordance <- left_join(
  tnbc57[, .(TF_clean, motif.name, TNBC_p = TF_p_value)],
  overall_bc_results[, .(TF_clean, motif, OverallBC_B = TF_B, OverallBC_p = TF_p, OverallBC_FDR = FDR)],
  by = "TF_clean")
concordance$in_OverallBC_nominal <- !is.na(concordance$OverallBC_p) & concordance$OverallBC_p < 0.05
concordance$in_OverallBC_FDR     <- !is.na(concordance$OverallBC_FDR) & concordance$OverallBC_FDR < 0.05
fwrite(concordance, file.path(outdir, "OverallBC_vs_TNBC_concordance.csv"))

ts("=== OVERALL BC SUMMARY ===")
ts(paste("Overall BC significant (FDR<0.05):", nrow(sig_overall_bc)))
ts(paste("TNBC 57 TFs in overall BC (nominal p<0.05):",
         sum(concordance$in_OverallBC_nominal, na.rm = TRUE), "of 57"))
ts(paste("TNBC 57 TFs in overall BC (FDR<0.05):",
         sum(concordance$in_OverallBC_FDR, na.rm = TRUE), "of 57"))
ts(paste("G2 gate: if overall BC signal recovers TNBC TFs -> restructure Results;",
         "if weaker -> keep TNBC-first, retract restructure claim"))
ts("=== DONE ===")
