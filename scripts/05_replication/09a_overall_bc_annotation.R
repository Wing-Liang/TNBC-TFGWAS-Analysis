#!/usr/bin/env Rscript
# 09a_overall_bc_annotation.R
# Attach the overall breast-cancer chi-square to the replication motif batches.
#
# Writes tfgwas/lmer_batches/annot_OverallBC.rds: the rows of annot.rds, in the same order, plus
# OverallBC_chisquare; the motif batches are reused unchanged. Because a misaligned row would
# silently pair each variant with another variant's motif indicators, the script rebuilds the
# block vector with the same filter as 07a_split_motif_batches.R and stops unless it is
# identical to annot.rds.
#
# Inputs : tfgwas/replication_variant_peak_motif_GWAS_ALL.csv.gz   (06b_motif_matrix.R)
#          tfgwas/lmer_batches/annot.rds                            (07a_split_motif_batches.R)
#          tfgwas/OverallBC_chisq_hg38.rds                              (03_tf_gwas/03_build_overall_bc_gwas.R)
# Outputs: tfgwas/lmer_batches/annot_OverallBC.rds
# Run    : Rscript scripts/05_replication/09a_overall_bc_annotation.R   (from the repository root)
## ---- paths (set the environment variables; see README) ----
TNBC_ROOT <- Sys.getenv("TNBC_ROOT", ".")
suppressPackageStartupMessages({ library(data.table) })
setwd(file.path(TNBC_ROOT, "replication"))
ts <- function(m) message("[", format(Sys.time(), "%H:%M:%S"), "] ", m)

IN_GZ <- "tfgwas/replication_variant_peak_motif_GWAS_ALL.csv.gz"
ts("reading key columns from the combined matrix")
hdr    <- names(fread(cmd = paste("zcat", IN_GZ, "| head -1"), header = TRUE, nrows = 0))
pos_col <- intersect(c("Pos38", "pos38"), hdr)[1]
if (is.na(pos_col)) stop("No Pos38/pos38 column in ", IN_GZ)
d <- fread(IN_GZ, select = c("Chr", pos_col, "block"), showProgress = FALSE)
setnames(d, pos_col, "pos38")
d[, pos38 := as.integer(pos38)]
ts(sprintf("  %s rows", format(nrow(d), big.mark = ",")))

d <- d[!is.na(block)]                       # identical filter to 07a_split_motif_batches.R
ts(sprintf("  after dropping NA block: %s rows", format(nrow(d), big.mark = ",")))

ts("verifying row alignment against annot.rds")
annot <- readRDS("tfgwas/lmer_batches/annot.rds")
if (nrow(annot) != nrow(d))
  stop("Row count mismatch: annot.rds ", nrow(annot), " vs reconstructed ", nrow(d))
if (!identical(as.character(annot$block), as.character(d$block)))
  stop("block vectors differ - row order does NOT match the motif batches. Aborting.")
ts("  PASS: row order identical to annot.rds")

ts("attaching overall BC chi-square")
overall <- readRDS("tfgwas/OverallBC_chisq_hg38.rds")
d[, n38 := .N, by = .(Chr, pos38)]          # matrix-side duplicated positions are unusable
## explicit match() on a composite key: returns exactly one index per row of d, so the
## assignment cannot silently mis-length (a join with which=TRUE did, and data.table
## caught it only because the lengths happened to differ).
i <- match(paste(d$Chr, d$pos38), paste(overall$Chr, overall$pos38))
stopifnot(length(i) == nrow(d))
d[, OverallBC_chisquare := overall$OverallBC_chisquare[i]]
d[n38 > 1, OverallBC_chisquare := NA_real_]   # ambiguous position -> unusable
n_ok <- sum(!is.na(d$OverallBC_chisquare))
ts(sprintf("  overall BC available for %s of %s rows (%.1f%%)",
           format(n_ok, big.mark = ","), format(nrow(d), big.mark = ","), 100 * n_ok / nrow(d)))
ts(sprintf("  LD blocks with >=1 overall BC variant: %s of %s",
           format(uniqueN(d[!is.na(OverallBC_chisquare)]$block), big.mark = ","),
           format(uniqueN(d$block), big.mark = ",")))

out <- cbind(annot, OverallBC_chisquare = d$OverallBC_chisquare)
saveRDS(out, "tfgwas/lmer_batches/annot_OverallBC.rds")
ts("saved -> tfgwas/lmer_batches/annot_OverallBC.rds")
ts("=== 09a_overall_bc_annotation.R complete ===")
