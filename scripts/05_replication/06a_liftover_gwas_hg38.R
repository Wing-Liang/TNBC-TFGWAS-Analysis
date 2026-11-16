#!/usr/bin/env Rscript
# ==============================================================================
# 06a_liftover_gwas_hg38.R
# Fill missing GRCh38 coordinates in the BCAC variant table.
# Produces: input to 06b_motif_matrix.R
# Run     : Rscript scripts/05_replication/06a_liftover_gwas_hg38.R   (from the repository root)
# ==============================================================================
# 06a_liftover_gwas_hg38.R -- fill in the missing hg38 coordinates in the BCAC GWAS table.
#
# WHY
# bcac_gwas_ldblock.csv carries pos38, but 173,263 of 9,576,981 variants (1.81%) have it
# as NA and are therefore dropped by 06b_motif_matrix.R. `block` is 100% complete, so the
# only thing missing is the hg38 coordinate. pos38 originally came from a left_join against
# ld_blocks_hg38_ALL.csv (02_gwas_ld_blocks/02_assign_ld_blocks.R), so the NAs are variants that failed
# liftOver when THAT file was built -- not variants that are intrinsically unmappable.
# Re-lifting them directly may recover most.
#
# THIS SCRIPT DOES TWO THINGS
#  1. VERIFIES the existing pos38 values by independently re-lifting a random sample and
#     comparing. If the existing column disagrees with a fresh liftOver, everything 06b_motif_matrix.R has
#     produced is suspect and that must surface before anything is added.
#  2. RECOVERS the missing ones and writes an augmented table. Variants that map to more
#     than one hg38 locus, or to none, stay NA -- an ambiguous coordinate is worse than a
#     missing one.
#
# The source file is NOT overwritten: it is shared with the discovery cohort's analyses.
#
# INPUT   ../results/gwas/bcac_gwas_ldblock.csv
#         $LIFTOVER_DIR/hg19ToHg38.over.chain
# OUTPUT  ../results/gwas/bcac_gwas_ldblock_hg38filled.csv
#         qc/liftover_hg38_recovery_report.csv
## ---- paths (set the environment variables; see README) ----
TNBC_ROOT <- Sys.getenv("TNBC_ROOT", ".")
LIFTOVER_DIR <- Sys.getenv("LIFTOVER_DIR", "data/liftover")
suppressPackageStartupMessages({
  library(data.table); library(GenomicRanges); library(GenomeInfoDb); library(rtracklayer)
})
setwd(file.path(TNBC_ROOT, "replication"))
ts <- function(m) message("[", format(Sys.time(), "%H:%M:%S"), "] ", m)
set.seed(1234)

SRC   <- "../results/gwas/bcac_gwas_ldblock.csv"
CHAIN <- file.path(LIFTOVER_DIR, "hg19ToHg38.over.chain")
OUT   <- "../results/gwas/bcac_gwas_ldblock_hg38filled.csv"

ts("reading BCAC")
b <- fread(SRC, header = TRUE, showProgress = FALSE)
ts(sprintf("  %d variants | pos38 NA: %d (%.2f%%)",
           nrow(b), sum(is.na(b$pos38)), 100 * mean(is.na(b$pos38))))

ts("importing chain")
ch <- import.chain(CHAIN)

lift_one <- function(chr, pos) {
  gr <- GRanges(chr, IRanges(pos, width = 1))
  lf <- liftOver(gr, ch)
  n  <- elementNROWS(lf)
  out <- rep(NA_integer_, length(gr))
  ok  <- n == 1                       # exactly one target: unambiguous
  if (any(ok)) out[ok] <- start(unlist(lf[ok]))
  list(pos = out, n_multi = sum(n > 1), n_none = sum(n == 0))
}

# ---------------------------------------------
#  1. VERIFY the existing pos38 against a fresh liftOver
# ---------------------------------------------
ts("verifying existing pos38 on a 200k random sample")
have <- b[!is.na(pos38)]
s    <- have[sample(.N, min(200000L, .N))]
r    <- lift_one(s$Chr, s$Pos)
cmp  <- data.table(existing = s$pos38, relifted = r$pos)[!is.na(relifted)]
agree <- sum(cmp$existing == cmp$relifted)
cat(sprintf("\n=== VERIFY existing pos38 ===\n"))
cat(sprintf("  sampled                    : %d\n", nrow(s)))
cat(sprintf("  re-lifted unambiguously    : %d\n", nrow(cmp)))
cat(sprintf("  AGREE with existing pos38  : %d (%.4f%%)\n", agree, 100*agree/nrow(cmp)))
cat(sprintf("  DISAGREE                   : %d\n", nrow(cmp) - agree))
if (nrow(cmp) - agree > 0) {
  cat("  -- examples of disagreement --\n")
  print(head(cmp[existing != relifted], 5))
  cat("  !! existing pos38 is NOT reproducible by liftOver; investigate before trusting 06\n")
}

# ---------------------------------------------
#  2. RECOVER the missing ones
# ---------------------------------------------
miss <- which(is.na(b$pos38))
ts(sprintf("re-lifting the %d variants with NA pos38", length(miss)))
r2 <- lift_one(b$Chr[miss], b$Pos[miss])
b[miss, pos38_lifted := r2$pos]
n_rec <- sum(!is.na(b$pos38_lifted[miss]))
cat(sprintf("\n=== RECOVER missing pos38 ===\n"))
cat(sprintf("  missing before             : %d\n", length(miss)))
cat(sprintf("  recovered (1:1 mapping)    : %d (%.1f%%)\n", n_rec, 100*n_rec/length(miss)))
cat(sprintf("  mapped to >1 locus (kept NA): %d\n", r2$n_multi))
cat(sprintf("  failed to map              : %d\n", r2$n_none))

## merge: keep the original value wherever it exists, fill only the gaps
b[, pos38_final := fifelse(is.na(pos38), pos38_lifted, pos38)]
b[, pos38_source := fifelse(is.na(pos38) & !is.na(pos38_lifted), "liftOver_recovered",
                     fifelse(is.na(pos38), "unmapped", "original"))]
cat(sprintf("\n  usable variants BEFORE     : %d\n", sum(!is.na(b$pos38) & !is.na(b$block))))
cat(sprintf("  usable variants AFTER      : %d  (+%d, +%.2f%%)\n",
            sum(!is.na(b$pos38_final) & !is.na(b$block)),
            sum(!is.na(b$pos38_final) & !is.na(b$block)) - sum(!is.na(b$pos38) & !is.na(b$block)),
            100*(sum(!is.na(b$pos38_final) & !is.na(b$block)) /
                 sum(!is.na(b$pos38) & !is.na(b$block)) - 1)))
print(table(b$pos38_source))

fwrite(b, OUT)
ts(paste("wrote", OUT))
fwrite(data.table(metric = c("total","na_pos38_before","recovered","multi_mapped",
                             "unmapped","usable_before","usable_after"),
                  value  = c(nrow(b), length(miss), n_rec, r2$n_multi, r2$n_none,
                             sum(!is.na(b$pos38) & !is.na(b$block)),
                             sum(!is.na(b$pos38_final) & !is.na(b$block)))),
       "qc/liftover_hg38_recovery_report.csv")
ts("=== 06a_liftover_gwas_hg38.R complete ===")
