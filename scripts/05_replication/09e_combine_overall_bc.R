#!/usr/bin/env Rscript
# 09e_combine_overall_bc.R
# Combine the overall breast-cancer TF-GWAS batch results (09c_tfgwas_overall_bc.R) into one
# table per cohort: one row per motif with the motif coefficient (OverallBC_B, OverallBC_SE, OverallBC_p)
# and the Benjamini-Hochberg FDR across motifs.
#
# Inputs : tfgwas/lmer_OverallBC_discovery/tfgwas_OverallBC_batch<N>.csv   (09d_tfgwas_overall_bc_discovery.sbatch)
#          tfgwas/lmer_OverallBC/tfgwas_OverallBC_batch<N>.csv             (09c_tfgwas_overall_bc.sbatch)
# Outputs: tfgwas/Discovery_OverallBC_ALL.csv
#          tfgwas/replication_OverallBC_ALL.csv
# Run    : Rscript scripts/05_replication/09e_combine_overall_bc.R   (from the repository root)

TNBC_ROOT <- Sys.getenv("TNBC_ROOT", ".")
suppressPackageStartupMessages(library(data.table))
setwd(file.path(TNBC_ROOT, "replication"))

COHORTS <- c(Discovery = "tfgwas/lmer_OverallBC_discovery", replication = "tfgwas/lmer_OverallBC")
for (co in names(COHORTS)) {
  files <- list.files(COHORTS[[co]], pattern = "^tfgwas_OverallBC_batch[0-9]+\\.csv$", full.names = TRUE)
  if (!length(files)) { message(co, ": no batch files in ", COHORTS[[co]], " -- skipped"); next }
  res <- rbindlist(lapply(files, fread))[order(motif)]
  res[, OverallBC_FDR := p.adjust(OverallBC_p, method = "BH")]
  out <- sprintf("tfgwas/%s_OverallBC_ALL.csv", co)
  fwrite(res, out)
  message(sprintf("%s: %d motifs from %d batches -> %s (Bonferroni-significant: %d)",
                  co, nrow(res), length(files), out, sum(res$OverallBC_p < 0.05 / 727, na.rm = TRUE)))
}
