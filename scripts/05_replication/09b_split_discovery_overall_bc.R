#!/usr/bin/env Rscript
# 09b_split_discovery_overall_bc.R
# Split the discovery variant x motif matrix into batches of 20 motifs and attach the overall
# breast-cancer chi-square, so that the overall-BC model is fitted the same way in both cohorts.
#
# The discovery matrix is keyed on GRCh37 positions, while the overall-BC statistics are GRCh38;
# they are joined on the matrix's own GRCh38 position column. This join reproduces the
# GRCh37-derived TNBC p-values for 99.9997% of variants once multi-allelic positions are
# excluded (03_tf_gwas/03_build_overall_bc_gwas.R); those positions are dropped here on both
# sides, since neither table stores alleles. Needs ~40 GB of memory.
#
# Inputs : results/tfgwas/variant_peak_motif_subtypes_ALL.csv      (03_tf_gwas/01_build_variant_motif_matrix.R)
#          tfgwas/OverallBC_chisq_hg38.rds                              (03_tf_gwas/03_build_overall_bc_gwas.R)
# Outputs: tfgwas/discovery_OverallBC_batches/annot_OverallBC.rds, motifs_batch<N>.rds, batch_index.csv, n_batch.txt
# Run    : Rscript scripts/05_replication/09b_split_discovery_overall_bc.R   (from the repository root)
## ---- paths (set the environment variables; see README) ----
TNBC_ROOT <- Sys.getenv("TNBC_ROOT", ".")
suppressPackageStartupMessages({ library(data.table) })
setwd(file.path(TNBC_ROOT, "replication"))
ts <- function(m) message("[", format(Sys.time(), "%H:%M:%S"), "] ", m)

BATCH_SIZE <- 20                       # same batch size as 03_tf_gwas/02_tfgwas_model.R
IN_CSV  <- "../results/tfgwas/variant_peak_motif_subtypes_ALL.csv"
OUT_DIR <- "tfgwas/discovery_OverallBC_batches"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

ts("reading header")
hdr        <- names(fread(cmd = paste("head -1", shQuote(IN_CSV)), header = TRUE, nrows = 0))
motif_cols <- grep("^MA[0-9]", hdr, value = TRUE)
pos_col    <- intersect(c("pos38", "Pos38"), hdr)[1]
need       <- c("Chr", pos_col, "block")
if (anyNA(c(pos_col))) stop("No pos38 column")
if (length(setdiff(need, hdr))) stop("Missing: ", paste(setdiff(need, hdr), collapse = ", "))
ts(sprintf("  %d columns: %d motifs", length(hdr), length(motif_cols)))
n_batch <- ceiling(length(motif_cols) / BATCH_SIZE)
ts(sprintf("  -> %d batches of %d", n_batch, BATCH_SIZE))

ts("reading matrix (single pass; expensive)")
d <- fread(IN_CSV, select = c(need, motif_cols), showProgress = FALSE)
setnames(d, pos_col, "pos38")
d[, pos38 := as.integer(pos38)]
ts(sprintf("  loaded: %s variants x %d cols", format(nrow(d), big.mark = ","), ncol(d)))

n0 <- nrow(d); d <- d[!is.na(block)]
ts(sprintf("  dropped %s with NA block -> %s remain",
           format(n0 - nrow(d), big.mark = ","), format(nrow(d), big.mark = ",")))

ts("attaching overall BC chi-square (multi-allelic positions excluded on both sides)")
overall <- readRDS("tfgwas/OverallBC_chisq_hg38.rds")
d[, n38 := .N, by = .(Chr, pos38)]
i <- match(paste(d$Chr, d$pos38), paste(overall$Chr, overall$pos38))
stopifnot(length(i) == nrow(d))
d[, OverallBC_chisquare := overall$OverallBC_chisquare[i]]
d[n38 > 1, OverallBC_chisquare := NA_real_]
n_ok <- sum(!is.na(d$OverallBC_chisquare))
ts(sprintf("  overall BC on %s of %s rows (%.1f%%) | LD blocks with overall BC: %s of %s",
           format(n_ok, big.mark = ","), format(nrow(d), big.mark = ","), 100 * n_ok / nrow(d),
           format(uniqueN(d[!is.na(OverallBC_chisquare)]$block), big.mark = ","),
           format(uniqueN(d$block), big.mark = ",")))

ts("writing shared annotation")
saveRDS(d[, .(block, OverallBC_chisquare)], file.path(OUT_DIR, "annot_OverallBC.rds"))

idx <- data.table(batch = integer(), motif = character())
for (b in seq_len(n_batch)) {
  cols <- motif_cols[((b - 1) * BATCH_SIZE + 1):min(b * BATCH_SIZE, length(motif_cols))]
  m <- as.matrix(d[, ..cols]); storage.mode(m) <- "integer"
  saveRDS(m, file.path(OUT_DIR, sprintf("motifs_batch%d.rds", b)))
  idx <- rbind(idx, data.table(batch = b, motif = cols))
  if (b %% 5 == 0 || b == n_batch) ts(sprintf("  batch %d / %d", b, n_batch))
}
fwrite(idx, file.path(OUT_DIR, "batch_index.csv"))
writeLines(as.character(n_batch), file.path(OUT_DIR, "n_batch.txt"))
ts(sprintf("DONE -> %s  (submit array 1-%d)", OUT_DIR, n_batch))
