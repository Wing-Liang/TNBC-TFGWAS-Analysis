#!/usr/bin/env Rscript
# ==============================================================================
# 15_cohort_matrix_check.R
# Structural check of the replication matrix against discovery.
# Produces: data check
# Run     : Rscript scripts/05_replication/15_cohort_matrix_check.R   (from the repository root)
# ==============================================================================
# Structural check of the replication matrix against discovery.
# Key question: are the motif columns aligned with their headers, and do the
# same motifs show comparable prevalence in both cohorts?
suppressPackageStartupMessages(library(data.table))
setDTthreads(4)
ts <- function(...) message("[", format(Sys.time(), "%H:%M:%S"), "] ", ...)

DISC <- "results/tfgwas/variant_peak_motif_subtypes_ALL.csv"
REPL <- "replication/tfgwas/replication_variant_peak_motif_GWAS_ALL.csv.gz"

grab <- function(f, lab) {
  hdr <- names(fread(f, nrows = 0))
  mc  <- setdiff(grep("^MA", hdr, value = TRUE), "Chromatin")
  ts(lab, ": ", length(hdr), " cols, ", length(mc), " motif cols")
  d <- fread(f, select = c("Chr", "block", "Chromatin", mc), showProgress = FALSE)
  s <- vapply(mc, function(cc) sum(d[[cc]]), numeric(1))
  list(n = nrow(d), chrom = sum(d$Chromatin == 1), sums = s, mc = mc,
       # how many motifs is a typical accessible variant assigned?
       per_var = sum(s) / max(sum(d$Chromatin == 1), 1))
}

D <- grab(DISC, "discovery")
R <- grab(REPL, "replication")

cat("\n================ STRUCTURE ================\n")
cat(sprintf("%-28s %14s %14s\n", "", "discovery", "replication"))
cat(sprintf("%-28s %14s %14s\n", "rows", format(D$n, big.mark=","), format(R$n, big.mark=",")))
cat(sprintf("%-28s %14s %14s\n", "Chromatin==1",
            format(D$chrom, big.mark=","), format(R$chrom, big.mark=",")))
cat(sprintf("%-28s %13.2f%% %13.2f%%\n", "Chromatin rate",
            100*D$chrom/D$n, 100*R$chrom/R$n))
cat(sprintf("%-28s %14s %14s\n", "total motif hits",
            format(sum(D$sums), big.mark=","), format(sum(R$sums), big.mark=",")))
cat(sprintf("%-28s %14.1f %14.1f\n", "motifs per accessible variant", D$per_var, R$per_var))

stopifnot(identical(D$mc, R$mc))
cat("\nmotif column names identical and in the same order: TRUE\n")

cat("\n========== PER-MOTIF PREVALENCE ==========\n")
cs <- cor(D$sums, R$sums)
cr <- cor(rank(D$sums), rank(R$sums))
cat(sprintf("Pearson r (per-motif hit counts, discovery vs replication): %.4f\n", cs))
cat(sprintf("Spearman r                                                : %.4f\n", cr))
cat(sprintf("ratio replication/discovery: median %.3f, IQR %.3f-%.3f\n",
            median(R$sums/D$sums), quantile(R$sums/D$sums, .25), quantile(R$sums/D$sums, .75)))
cat(sprintf("motifs with ZERO hits: discovery %d, replication %d\n",
            sum(D$sums == 0), sum(R$sums == 0)))

o <- order(-D$sums)[1:10]
cat("\ntop 10 motifs by discovery prevalence:\n")
cat(sprintf("  %-11s %12s %12s %8s\n", "motif", "discovery", "replication", "ratio"))
for (i in o) cat(sprintf("  %-11s %12s %12s %8.2f\n", D$mc[i],
                         format(D$sums[i], big.mark=","), format(R$sums[i], big.mark=","),
                         R$sums[i]/D$sums[i]))
o2 <- order(R$sums/D$sums)[1:6]
cat("\n6 motifs most depleted in replication:\n")
for (i in o2) cat(sprintf("  %-11s disc %10s  repl %10s  ratio %.3f\n", D$mc[i],
                          format(D$sums[i], big.mark=","), format(R$sums[i], big.mark=","),
                          R$sums[i]/D$sums[i]))
