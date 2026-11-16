#!/usr/bin/env Rscript
# 07a_split_motif_batches.R
# Split the replication variant x motif matrix into batches of 20 motifs (one RDS per batch,
# plus one RDS with the LD block and the five subtype chi-squares), so that each task of the
# model array (07b_tfgwas_model.sbatch) reads only its own motifs instead of re-parsing the full
# compressed matrix. Rows with no LD block are dropped. Optional: 07b_tfgwas_model.R reads the
# full matrix when the batch files are absent. Needs ~35 GB of memory.
#
# Inputs : tfgwas/replication_variant_peak_motif_GWAS_ALL.csv.gz   (06b_motif_matrix.R)
# Outputs: tfgwas/lmer_batches/annot.rds, motifs_batch<N>.rds, batch_index.csv, n_batch.txt
# Run    : Rscript scripts/05_replication/07a_split_motif_batches.R   (from the repository root)

## ---- paths (set the environment variables; see README) ----
TNBC_ROOT <- Sys.getenv("TNBC_ROOT", ".")
suppressPackageStartupMessages({ library(data.table) })
setwd(file.path(TNBC_ROOT, "replication"))
ts <- function(m) message("[", format(Sys.time(), "%H:%M:%S"), "] ", m)

BATCH_SIZE <- 20   # same batch size as 03_tf_gwas/02_tfgwas_model.R
IN_GZ      <- "tfgwas/replication_variant_peak_motif_GWAS_ALL.csv.gz"
OUT_DIR    <- "tfgwas/lmer_batches"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

## the model needs `block` (random effect) and the 5 subtype outcomes; nothing else
ANNOT_COLS <- c("block",
                "TNBC_chisquare", "LuminalA_chisquare", "LuminalB_chisquare",
                "LuminalB_HER2Neg_chisquare", "HER2_chisquare")

if (!file.exists(IN_GZ)) stop("Not found: ", IN_GZ, "\n  Run 06b_motif_matrix.R first.")

ts("reading header")
hdr        <- names(fread(cmd = paste("zcat", IN_GZ, "| head -1"), header = TRUE, nrows = 0))
motif_cols <- grep("^MA[0-9]", hdr, value = TRUE)   # discovery's motif-column pattern
missing    <- setdiff(ANNOT_COLS, hdr)
if (length(missing)) stop("Matrix is missing required column(s): ", paste(missing, collapse = ", "))
ts(sprintf("  %d columns: %d annotation + %d motifs", length(hdr),
           length(hdr) - length(motif_cols), length(motif_cols)))

n_batch <- ceiling(length(motif_cols) / BATCH_SIZE)
ts(sprintf("  -> %d batches of %d", n_batch, BATCH_SIZE))

ts("reading full matrix (single pass; this is the expensive step)")
d <- fread(IN_GZ, select = c(ANNOT_COLS, motif_cols), showProgress = FALSE)
ts(sprintf("  loaded: %d variants x %d columns", nrow(d), ncol(d)))

## drop variants the model cannot use: no LD block => no random effect
n0 <- nrow(d)
d  <- d[!is.na(block)]
ts(sprintf("  dropped %d variants with NA block -> %d remain", n0 - nrow(d), nrow(d)))

ts("writing shared annotation")
saveRDS(d[, ..ANNOT_COLS], file.path(OUT_DIR, "annot.rds"))

idx <- data.table(batch = integer(), motif = character())
for (b in seq_len(n_batch)) {
  cols <- motif_cols[((b - 1) * BATCH_SIZE + 1):min(b * BATCH_SIZE, length(motif_cols))]
  m    <- as.matrix(d[, ..cols])
  storage.mode(m) <- "integer"          # motif indicators are 0/1
  saveRDS(m, file.path(OUT_DIR, sprintf("motifs_batch%d.rds", b)))
  idx  <- rbind(idx, data.table(batch = b, motif = cols))
  if (b %% 5 == 0 || b == n_batch) ts(sprintf("  batch %d / %d", b, n_batch))
}
fwrite(idx, file.path(OUT_DIR, "batch_index.csv"))
writeLines(as.character(n_batch), file.path(OUT_DIR, "n_batch.txt"))

ts(sprintf("DONE -> %s  (%d batches; submit array 1-%d)", OUT_DIR, n_batch, n_batch))
ts("=== 07a_split_motif_batches.R complete ===")
