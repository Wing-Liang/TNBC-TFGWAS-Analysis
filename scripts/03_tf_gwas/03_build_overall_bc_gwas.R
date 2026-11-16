#!/usr/bin/env Rscript
# ==============================================================================
# 03_build_overall_bc_gwas.R
# Join the BCAC overall breast-cancer GWAS to the variant table (GRCh38).
# Produces: Table S4 (overall BC)
# Inputs  : ../results/gwas/bcac_gwas_ldblock.csv
# Outputs : tfgwas/OverallBC_chisq_hg38.rds
# Run     : Rscript scripts/03_tf_gwas/03_build_overall_bc_gwas.R   (from the repository root)
# ==============================================================================
# 03_build_overall_bc_gwas.R — build the overall breast cancer (overall BC) chi-square annotation
#                     shared by the DISCOVERY (hg19) and EUROPEAN (hg38) TF-GWAS analyses,
#                     and PROVE the join is variant-correct before writing it.
#
# GENOME BUILD — the reason this script is careful
# The subtype chi-squares already in both matrices were joined on hg19: the raw BCAC files
# key on chr.iCOGs / Position.iCOGs (verified: Position.iCOGs 10001753 -> matrix Pos
# 10001753 with identical p). The only Overall-BC summary stats available are hg38 (_b38).
# The DISCOVERY single-cell data is hg19 and the EUROPEAN is hg38, but that does not matter
# for attaching a GWAS statistic: the statistic belongs to a VARIANT, and every matrix row
# carries both Pos (hg19) and pos38. What matters is that the join lands on the RIGHT
# variant, which is what this script tests rather than assumes.
#
# MULTI-ALLELIC POSITIONS — the correctness trap
# A first attempt joined on (Chr, pos38) after unique(by = position), and 59,762 of 8.06M
# variants (0.74%) came back with the WRONG statistic (e.g. p = 0.997 vs 0.0003). Cause:
# 50,043 positions carry >1 variant in the hg38 file (60,386 in the hg19 file), and
# unique() silently kept an arbitrary allele. Neither the matrix nor the BCAC key stores
# alleles, so such positions CANNOT be disambiguated -- they are DROPPED, not guessed.
#
# ROBUSTNESS: p is recomputed from beta/se rather than read by column name. The _spredixcan
# files do not share a column order -- in Triple_Neg the header says col8=pvalue, col9=freq
# but the values are reversed (verified against beta/se). Deriving p is immune to that.
#
# CHI-SQUARE DEFINITION: (beta/se)^2, i.e. z^2 -- the 1-df chi-square.
# This is the SAME quantity the subtypes hold (verified: TNBC_chisquare ==
# computed stably. qchisq(1 - p, df = 1) OVERFLOWS TO Inf for genome-wide-significant
# variants: when p is tiny, 1 - p rounds to exactly 1 and qchisq(1, 1) = Inf. That produced
# 4,832 Inf values, and lmer fails with "NA/NaN/Inf in y", silently NA-ing all 727 fits.
# Overall BC is far better powered than any subtype, so it reaches those p-values where the
# subtype columns (max chi-square 218-1470, zero Inf) never did.
#
# ============================= DATA SOURCES (imports) =============================
#  $BCAC_DIR/GWAS_QC_Overall_BC_spredixcan.txt  <- BCAC overall BC (hg38)
#  $BCAC_DIR/GWAS_QC_Triple_Neg_spredixcan.txt  <- ground truth for the
#      join test: its hg38 p must reproduce the hg19-derived TNBC_pval already in the matrix
#  ../results/gwas/bcac_gwas_ldblock.csv                     <- the variant/LD-block table
# Outputs -> tfgwas/OverallBC_chisq_hg38.rds  (Chr, pos38, OverallBC_chisquare)
#            written ONLY if the join test passes
# ==================================================================================
## ---- paths (set the environment variables; see README) ----
BCAC_DIR <- Sys.getenv("BCAC_DIR", "data/BCAC")
TNBC_ROOT <- Sys.getenv("TNBC_ROOT", ".")
suppressPackageStartupMessages({ library(data.table) })
setwd(file.path(TNBC_ROOT, "replication"))
ts <- function(m) message("[", format(Sys.time(), "%H:%M:%S"), "] ", m)
dir.create("tfgwas", showWarnings = FALSE)
GW <- BCAC_DIR
p_from_beta_se <- function(beta, se) 2 * pnorm(-abs(beta / se))

## unambiguous = exactly one variant at that (Chr, pos38) in the source file
read_unambiguous <- function(f, value_name) {
  d <- fread(f, select = c("chr", "pos", "beta", "se"), showProgress = FALSE)
  d[, `:=`(Chr = paste0("chr", sub("^chr", "", chr)), pos38 = as.integer(pos))]
  d[, n_at_pos := .N, by = .(Chr, pos38)]
  n_multi <- sum(d$n_at_pos > 1)
  d <- d[n_at_pos == 1 & is.finite(beta) & is.finite(se) & se > 0]
  d[, (value_name) := p_from_beta_se(beta, se)]
  d[, chisq_z2 := (beta / se)^2]        # 1-df chi-square, computed stably (see header)
  ts(sprintf("  %s: kept %s unambiguous; dropped %s at multi-allelic positions",
             basename(f), format(nrow(d), big.mark = ","), format(n_multi, big.mark = ",")))
  d[, c("Chr", "pos38", value_name, "chisq_z2"), with = FALSE]
}

ts("reading BCAC LD-block table")
m <- fread("../results/gwas/bcac_gwas_ldblock.csv",
           select = c("Chr", "Pos", "pos38", "TNBC_pval"), showProgress = FALSE)
m[, pos38 := as.integer(pos38)]
m[, n38 := .N, by = .(Chr, pos38)]                 # positions duplicated on the MATRIX side
ts(sprintf("  %s rows; %s at duplicated pos38 (also unusable)",
           format(nrow(m), big.mark = ","), format(sum(m$n38 > 1), big.mark = ",")))

## ---- 1. JOIN TEST on Triple_Neg (ground truth) ----
ts("join test: Triple_Neg hg38 must reproduce the hg19-derived TNBC_pval")
tn <- read_unambiguous(file.path(GW, "GWAS_QC_Triple_Neg_spredixcan.txt"), "p_hg38")
j  <- merge(m[n38 == 1], tn, by = c("Chr", "pos38"), all.x = TRUE, sort = FALSE)
ok <- j[!is.na(p_hg38)]
agree_rate <- mean(abs(ok$p_hg38 - ok$TNBC_pval) < 1e-6)
ts(sprintf("  matched %s | AGREE %.4f%% | DISAGREE %s",
           format(nrow(ok), big.mark = ","), 100 * agree_rate,
           format(sum(abs(ok$p_hg38 - ok$TNBC_pval) >= 1e-6), big.mark = ",")))
if (agree_rate < 0.999)
  stop("Join test FAILED (", round(100*agree_rate, 3), "% agreement). overall BC NOT written.")
ts("  PASS -> the (Chr,pos38) join on unambiguous positions is variant-correct")

## ---- 2. build overall BC using the identical rule ----
ts("building overall BC chi-square")
overall <- read_unambiguous(file.path(GW, "GWAS_QC_Overall_BC_spredixcan.txt"), "OverallBC_p")
overall[, OverallBC_chisquare := chisq_z2]
n_inf <- sum(!is.finite(overall$OverallBC_chisquare))
if (n_inf) ts(sprintf("  dropping %s non-finite chi-squares", format(n_inf, big.mark = ",")))
overall <- overall[is.finite(OverallBC_chisquare), .(Chr, pos38, OverallBC_chisquare)]
ts(sprintf("  chi-square range: %.3g to %.1f (subtypes max ~1470; Inf would break lmer)",
           min(overall$OverallBC_chisquare), max(overall$OverallBC_chisquare)))
setkey(overall, Chr, pos38)
saveRDS(overall, "tfgwas/OverallBC_chisq_hg38.rds")

cov <- m[n38 == 1][overall, on = .(Chr, pos38), nomatch = NULL]
ts(sprintf("  overall BC variants: %s | usable matrix rows covered: %s of %s (%.1f%%)",
           format(nrow(overall), big.mark = ","), format(nrow(cov), big.mark = ","),
           format(nrow(m), big.mark = ","), 100 * nrow(cov) / nrow(m)))
ts("saved -> tfgwas/OverallBC_chisq_hg38.rds")
ts("=== 03_build_overall_bc_gwas.R complete ===")
