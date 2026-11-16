#!/usr/bin/env Rscript
# ==============================================================================
# 14_tf_family_models.R
# TF-GWAS at motif, TF-family and JASPAR-cluster level.
# Produces: Table S8
# Run     : Rscript scripts/05_replication/14_tf_family_models.R   (from the repository root)
# ==============================================================================
# =============================================================================
# 14_tf_family_models.R -- TF-GWAS at motif, TF-family and JASPAR-cluster level (Table S8)
#
# One pass over a cohort's matrix; for every response (5 subtypes + overall BC) and
# every level (motif 727 / JASPAR122 / TF family 94) computes the observed
# t statistic, p value and Bonferroni significance.
#
# Model, identical to the permutation pipeline:
#     <RESP> ~ group_indicator + (1 | LD block)
#   group indicator = 1 if the variant lies in an accessible peak containing ANY
#   member motif of that group (union / logical OR). At motif level each motif is
#   its own group, reproducing the published scan.
#
# Bonferroni denominator is the number of tests AT THAT LEVEL (727 / 122 / 94).
#
# Build: rotation is irrelevant here (no permutation), but rows are still ordered
# on the cohort's own overlap build for reproducibility -- discovery hg19 `Pos`,
# replication hg38 `Pos38`.
# =============================================================================

suppressPackageStartupMessages({ library(data.table); library(Matrix); library(lme4) })

CSV     <- Sys.getenv("OBS_CSV",    "results/tfgwas/variant_peak_motif_subtypes_ALL.csv")
POSCOL  <- Sys.getenv("OBS_POSCOL", "Pos")
COHORT  <- Sys.getenv("OBS_COHORT", "discovery")
OUTDIR  <- Sys.getenv("OBS_OUT",    "replication/tfgwas_grouped")
OVERALL_FILE <- "replication/tfgwas/OverallBC_chisq_hg38.rds"
MAPS    <- c(motif = "reference/group_motif.csv",
             JASPAR122 = "reference/group_JASPAR122.csv",
             TFfamily  = "reference/group_TFfamily.csv")
ALPHA   <- 0.05
setDTthreads(as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "4")))
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)
ts <- function(...) message("[", format(Sys.time(), "%H:%M:%S"), "] ", ...)
ts("cohort=", COHORT, " poscol=", POSCOL)

maps <- lapply(MAPS, fread)
motif_order <- maps$motif$motif_id

hdr <- names(fread(CSV, nrows = 0))
motif_cols <- setdiff(grep("^MA", hdr, value = TRUE), "Chromatin")
stopifnot(setequal(motif_cols, motif_order))
pcol <- intersect(c("pos38", "Pos38"), hdr)[1]
resp_cols <- c("TNBC_chisquare","LuminalA_chisquare","LuminalB_chisquare",
               "LuminalB_HER2Neg_chisquare","HER2_chisquare")
stopifnot(all(resp_cols %in% hdr))

keep <- unique(c("Chr", POSCOL, pcol, "block", "Chromatin", resp_cols, motif_order))
keep <- keep[!is.na(keep)]
ts("reading ", CSV, " (", length(keep), " cols) ...")
dat <- fread(CSV, select = keep,
             colClasses = setNames(rep("integer", length(motif_order)+1L),
                                   c(motif_order, "Chromatin")), showProgress = FALSE)
ts("rows: ", format(nrow(dat), big.mark = ","))

# overall BC response, joined on (Chr, hg38 pos); multi-allelic hg38 positions excluded
if (!is.na(pcol)) {
  setnames(dat, pcol, "hg38pos")
  # In replication the ordering column IS the hg38 column, so it has just been
  # renamed out from under POSCOL. Repoint it or the later setorderv fails.
  if (identical(POSCOL, pcol)) POSCOL <- "hg38pos"
  dat[, hg38pos := as.integer(hg38pos)]
  overall <- readRDS(OVERALL_FILE)
  dat[, n38 := .N, by = .(Chr, hg38pos)]
  i <- match(paste(dat$Chr, dat$hg38pos), paste(overall$Chr, overall$pos38))
  dat[, OverallBC_chisquare := overall$OverallBC_chisquare[i]]
  dat[n38 > 1, OverallBC_chisquare := NA_real_]
  ts(sprintf("overall BC joined on %s of %s rows (%.1f%%)",
             format(sum(!is.na(dat$OverallBC_chisquare)), big.mark=","),
             format(nrow(dat), big.mark=","), 100*mean(!is.na(dat$OverallBC_chisquare))))
  rm(overall); gc()
}
dat <- dat[!is.na(block) & nzchar(as.character(block))]
setorderv(dat, c("Chr", POSCOL))

# ---- motif design, then the collapsed designs ------------------------------
ts("building designs ...")
cols <- lapply(motif_order, function(cc) { v <- dat[[cc]]
  if (!all(v %in% c(0L,1L))) stop("non-binary: ", cc); which(v == 1L) })
n <- nrow(dat)
Xm <- sparseMatrix(i = unlist(cols, use.names=FALSE),
                   j = rep.int(seq_along(cols), lengths(cols)),
                   x = 1, dims = c(n, length(motif_order)))
for (cc in motif_order) set(dat, j = cc, value = NULL); rm(cols); gc()

designs <- list()
for (lv in names(MAPS)) {
  g <- factor(maps[[lv]]$group[match(motif_order, maps[[lv]]$motif_id)])
  M <- sparseMatrix(i = seq_along(g), j = as.integer(g), x = 1,
                    dims = c(length(g), nlevels(g)))
  X <- Xm %*% M; X@x[] <- 1; X <- drop0(X); colnames(X) <- levels(g)
  designs[[lv]] <- X
  ts("  ", lv, ": ", ncol(X), " tests, nnz = ", format(length(X@x), big.mark=","))
}
rm(Xm); gc()

blockchr <- as.character(dat$block)

RESPS <- c(TNBC="TNBC_chisquare", LuminalA="LuminalA_chisquare",
           LuminalB="LuminalB_chisquare", LuminalB_HER2Neg="LuminalB_HER2Neg_chisquare",
           HER2="HER2_chisquare", OverallBC="OverallBC_chisquare")
summ <- list()
for (rn in names(RESPS)) {
  rc <- RESPS[[rn]]
  if (!rc %in% names(dat)) { ts("skip ", rn); next }
  ok <- !is.na(dat[[rc]])
  ts(rn, ": subsetting to ", format(sum(ok), big.mark=","), " complete rows and fitting null")
  # SUBSET, do not zero-weight: block sizes and block sums must be computed on
  # the analysed rows only. Zero-weighting inflates n_b and the block sums for
  # any response with missing values (overall BC loses ~17% to the hg38 join).
  y  <- as.numeric(dat[[rc]])[ok]
  gf <- factor(blockchr[ok]); gidx <- as.integer(gf); nb <- nlevels(gf)
  nb_vec <- as.numeric(tabulate(gidx, nbins = nb))
  nsub <- sum(ok)
  B <- sparseMatrix(i = seq_len(nsub), j = gidx, x = 1, dims = c(nsub, nb))
  fit0 <- lmer(y ~ 1 + (1|gf), REML=TRUE, control=lmerControl(calc.derivs=FALSE))
  vc <- as.data.frame(VarCorr(fit0))
  s2b <- vc$vcov[vc$grp=="gf"]; s2e <- vc$vcov[vc$grp=="Residual"]; s_e <- sqrt(s2e)
  lam <- 1 - s_e/sqrt(s2e + nb_vec*s2b); k2 <- 2*lam - lam^2
  yt <- (y - lam[gidx]*as.numeric(tapply(y, gidx, mean))[gidx])/s_e
  one <- (1-lam[gidx])/s_e
  p_one <- sum(one^2); one_y <- sum(one*yt); yt2 <- sum(yt^2)
  Ytb <- as.numeric(tapply(yt, gidx, sum))
  nobs <- nsub
  for (lv in names(designs)) {
    X <- designs[[lv]][ok, , drop = FALSE]; m <- ncol(X)
    S <- as.matrix(Matrix::crossprod(B, X))
    XtY <- as.numeric(Matrix::crossprod(X, yt))
    nX <- (colSums(S) - colSums(k2*S^2/nb_vec))/s2e
    v1 <- colSums((1-lam)^2*S)/s2e
    xy <- (XtY - colSums(lam*S*Ytb/nb_vec))/s_e
    Sxx <- nX - v1^2/p_one
    Sxy <- xy - v1*(one_y/p_one)
    ssey <- yt2 - one_y^2/p_one
    df <- nobs - 2
    beta <- Sxy/Sxx; sse <- ssey - beta*Sxy
    tt <- ifelse(Sxx > .Machine$double.eps,
                 beta/sqrt(pmax(sse,.Machine$double.eps)/df/Sxx), NA_real_)
    pp <- 2*pt(abs(tt), df, lower.tail=FALSE)
    bonf <- ALPHA/m
    per <- data.table(cohort=COHORT, response=rn, level=lv, group=colnames(X),
                      t_observed=tt, p_observed=pp,
                      bonferroni_significant = !is.na(pp) & pp < bonf)
    fwrite(per, file.path(OUTDIR, sprintf("obs_%s_%s_%s.csv", COHORT, lv, rn)))
    summ[[length(summ)+1]] <- data.table(cohort=COHORT, response=rn, level=lv,
      n_variants=nobs, n_tests=m, bonferroni=bonf,
      n_significant=sum(per$bonferroni_significant))
    ts(sprintf("   %-10s %3d/%3d significant (p < %.2e)", lv,
               sum(per$bonferroni_significant), m, bonf))
  }
}
out <- rbindlist(summ)
fwrite(out, file.path(OUTDIR, sprintf("obs_summary_%s.csv", COHORT)))
print(as.data.frame(out))
ts("done -> ", normalizePath(OUTDIR))
