#!/usr/bin/env Rscript
# ==============================================================================
# 08_tfgwas_combine.R
# Combine replication batches and compare with the discovery cohort.
# Produces: Table S7
# Inputs  : tfgwas/replication_motif_list.csv
#           tfgwas/replication_chromvar.rds
# Outputs : tfgwas/tfgwas_results_ALL.csv
#           tfgwas/concordance
#           tfgwas/replication_replication_enrichment.csv
#           tfgwas/replication_replication_summary.csv
# Run     : Rscript scripts/05_replication/08_tfgwas_combine.R   (from the repository root)
# ==============================================================================
# 08_tfgwas_combine.R — combine motif batches, apply FDR, score concordance vs the discovery
#                     cohort, and compute per-cell-type chromVAR enrichment.
#
# 07b_tfgwas_model.R fits each motif once on all variants genome-wide, as in the discovery
# cohort, so there is nothing to combine statistically: this script binds the batch files
# and applies BH-FDR across motifs, per subtype.
#
# ============================= DATA SOURCES (imports) =============================
#  tfgwas/lmer_bymotif/tfgwas_batch{N}.csv   <- 07b_tfgwas_model.R
#  tfgwas/replication_chromvar.rds              <- 06b_motif_matrix.R
#  tfgwas/replication_motif_list.csv            <- 06b_motif_matrix.R (motif, TF_name, family)
#  ../results/tfgwas/tfgwas_significant_*_unique.csv  <- discovery significant TFs
# Outputs -> tfgwas/
#  tfgwas_results_ALL.csv                                one row per motif, all 5 subtypes + FDR
#  concordance/replication_{subtype}_concordance.csv
#  chromvar_enrichment/replication_{subtype}_celltype_chromVAR.csv
# ==================================================================================

## ---- paths (set the environment variables; see README) ----
TNBC_ROOT <- Sys.getenv("TNBC_ROOT", ".")
suppressPackageStartupMessages({
  library(Seurat); library(Signac); library(data.table); library(dplyr)
})
setwd(file.path(TNBC_ROOT, "replication"))
ts <- function(m) message("[", format(Sys.time(), "%H:%M:%S"), "] ", m)
clean_name <- function(x) toupper(gsub("[^A-Za-z0-9]", "", x))

subtypes <- c("TNBC", "LuminalA", "LuminalB", "LuminalB_HER2Neg", "HER2")

# ─────────────────────────────────────────────
#  1. Collect batches
# ─────────────────────────────────────────────
batch_files <- list.files("tfgwas/lmer_bymotif", pattern = "^tfgwas_batch[0-9]+\\.csv$",
                          full.names = TRUE)
if (!length(batch_files)) stop("No batch files in tfgwas/lmer_bymotif/ - run 07b_tfgwas_model.R")

## verify the array ran to completion rather than silently combining a partial set
nb_file <- "tfgwas/lmer_batches/n_batch.txt"
if (file.exists(nb_file)) {
  n_expected <- as.integer(readLines(nb_file)[1])
  got <- sort(as.integer(gsub(".*batch([0-9]+)\\.csv$", "\\1", batch_files)))
  miss <- setdiff(seq_len(n_expected), got)
  if (length(miss))
    stop("Missing batch file(s): ", paste(miss, collapse = ", "),
         "\n  Re-run those array tasks before combining; FDR over a partial set is wrong.")
  ts(sprintf("all %d batches present", n_expected))
} else ts("n_batch.txt absent - cannot verify completeness (check batch count yourself)")

res <- rbindlist(lapply(batch_files, fread), use.names = TRUE, fill = TRUE)
ts(sprintf("loaded %d rows | %d unique motifs", nrow(res), uniqueN(res$motif)))
if (anyDuplicated(res$motif))
  stop("Duplicate motifs across batches - batches overlap; check 07a_split_motif_batches.R")

# ─────────────────────────────────────────────
#  2. FDR across motifs, per subtype
# ─────────────────────────────────────────────
for (sub in subtypes)
  res[, paste0(sub, "_FDR") := p.adjust(get(paste0(sub, "_p")), method = "BH")]

## motif -> TF name from the discovery curated motif list, joined by motif as in 03_tf_gwas/02_tfgwas_model.R
MOTIF_ID_NAME <- "../results/singlecell/jaspar2022_motifs_727.csv"
motif_id_name <- fread(MOTIF_ID_NAME, header = TRUE)
res <- left_join(res, motif_id_name, by = "motif")
setDT(res)
## family tag (for the PWM-redundancy clustering) is added by 06b_motif_matrix.R, not in the discovery list
ml <- fread("tfgwas/replication_motif_list.csv", header = TRUE)
if ("family" %in% names(ml)) res <- merge(res, ml[, .(motif, family)], by = "motif",
                                          all.x = TRUE, sort = FALSE)
res[, TF_clean := clean_name(motif.name)]   # discovery cleans motif.name likewise
setorder(res, TNBC_p, na.last = TRUE)
fwrite(res, "tfgwas/tfgwas_results_ALL.csv")
ts(paste("tfgwas_results_ALL.csv:", nrow(res), "motifs"))
for (sub in subtypes)
  ts(sprintf("  %-16s nominal p<0.05: %4d | FDR<0.05: %4d | NA fits: %4d", sub,
             sum(res[[paste0(sub, "_p")]]   < 0.05, na.rm = TRUE),
             sum(res[[paste0(sub, "_FDR")]] < 0.05, na.rm = TRUE),
             sum(is.na(res[[paste0(sub, "_p")]]))))

# ─────────────────────────────────────────────
#  3. Concordance vs discovery + cell-type chromVAR enrichment
# ─────────────────────────────────────────────
dir.create("tfgwas/concordance",         showWarnings = FALSE)
dir.create("tfgwas/chromvar_enrichment", showWarnings = FALSE)

subtype_cfg <- list(
  TNBC             = list(sig_file = "tfgwas_significant_TNBC_unique.csv",         pre = "TNBC"),
  LuminalA         = list(sig_file = "tfgwas_significant_LumA_unique.csv",         pre = "LuminalA"),
  LuminalB         = list(sig_file = "tfgwas_significant_LumB_unique.csv",         pre = "LuminalB"),
  LuminalB_HER2Neg = list(sig_file = "tfgwas_significant_LumB_HER2Neg_unique.csv", pre = "LuminalB_HER2Neg"),
  HER2             = list(sig_file = "tfgwas_significant_HER2_unique.csv",         pre = "HER2")
)

ts("loading replication_chromvar.rds for cell-type enrichment")
obj <- readRDS("tfgwas/replication_chromvar.rds")
## SeuratObject 5 made the `slot` argument defunct (it errors, it does not warn):
##   GetAssayData(..., slot="data") -> deprecate_stop -> "Execution halted"
## This killed the combine AFTER tfgwas_results_ALL.csv had been written, so the main result
## survived but the cell-type enrichment block never ran. Use `layer` instead.
cv_scores <- GetAssayData(obj[["chromvar"]], layer = "data")
cv_ids    <- rownames(cv_scores)

## ─────────────────────────────────────────────
##  MATCHING IS BY MOTIF ID, NOT TF NAME
## ─────────────────────────────────────────────
## The previous join used TF_clean = toupper(gsub("[^A-Za-z0-9]","",motif.name)). That is
## lossy: one TF name maps to several JASPAR versions -- CEBPG, NFIC and FOSL2::JUND each
## appear TWICE in the discovery TNBC list -- so a 57-motif discovery set expanded to 64
## joined rows carrying different coefficients, and the TNBC enrichment read 1.54x
## (p=6.3e-4) instead of the correct 1.68x (p=6.5e-5).
## Every discovery motif ID is present in the European 727 (verified: 0 missing across all
## five subtypes), so the ID join is exact and lossless. The chromvar assay is itself keyed
## by motif ID, so cell-type enrichment now indexes it directly instead of via cleaned names.
if (!all(res$motif %in% cv_ids))
  ts(sprintf("  note: %d result motifs absent from the chromvar assay",
             sum(!res$motif %in% cv_ids)))

meta_ct    <- data.frame(cell = colnames(obj), celltype = obj$celltype, stringsAsFactors = FALSE)
meta_ct    <- meta_ct[!is.na(meta_ct$celltype), ]
cell_types <- sort(unique(meta_ct$celltype))
ts(paste("cell types:", paste(cell_types, collapse = ", ")))

N_MOTIF     <- nrow(res)      # 727 -- the universe for every enrichment test below
## Bonferroni is the significance rule THROUGHOUT, in both cohorts.
## Defined here because section 3 needs it; section 4 re-derives the same value.
BONF        <- 0.05 / N_MOTIF
enrich_rows <- list()

for (sub_name in names(subtype_cfg)) {
  cfg      <- subtype_cfg[[sub_name]]
  sig_path <- file.path("../results/singlecell/GWAS", cfg$sig_file)
  if (!file.exists(sig_path)) { ts(paste("WARNING: not found:", sig_path)); next }

  b_col <- paste0(cfg$pre, "_B");  se_col  <- paste0(cfg$pre, "_SE")
  p_col <- paste0(cfg$pre, "_p");  fdr_col <- paste0(cfg$pre, "_FDR")

  orig <- unique(fread(sig_path, header = TRUE)[, .(motif, motif.name, TF_p_value)])
  euro <- res[, c("motif", b_col, se_col, p_col, fdr_col), with = FALSE]
  conc <- merge(orig, euro, by = "motif", all.x = TRUE, sort = FALSE)   # <- BY MOTIF ID
  n_absent <- sum(!orig$motif %in% res$motif)
  if (n_absent)
    ts(sprintf("  !! %d discovery motif IDs absent from the European 727 - check 06 motif set",
               n_absent))

  ## Bonferroni in BOTH cohorts is the primary criterion; the discovery set was defined the
  ## same way, so the comparison is symmetric. Nominal and FDR are retained as extra columns
  ## but must not be quoted without stating the threshold.
  conc$replicated_bonf    <- !is.na(conc[[p_col]])   & conc[[p_col]]   < BONF
  conc$replicated_nominal <- !is.na(conc[[p_col]])   & conc[[p_col]]   < 0.05
  conc$replicated_FDR     <- !is.na(conc[[fdr_col]]) & conc[[fdr_col]] < 0.05
  conc$same_direction     <- !is.na(conc[[b_col]])   & conc[[b_col]]   > 0
  fwrite(conc, file.path("tfgwas/concordance",
                         paste0("replication_", sub_name, "_concordance.csv")))

  ## ENRICHMENT OVER CHANCE.  "x of y discovery TFs replicated" is uninterpretable alone:
  ## 27%-51% of ALL 727 motifs reach nominal p<0.05 in this cohort, so a high replication
  ## rate can be pure background. The hypergeometric test asks whether the discovery set
  ## replicates MORE often than a random set of the same size drawn from the same 727.
  K <- sum(res[[p_col]] < BONF, na.rm = TRUE)          # Bonferroni-significant, full European set
  n <- sum(!is.na(conc[[p_col]]))                      # discovery motifs testable here
  k <- sum(conc$replicated_bonf, na.rm = TRUE)         # of those, replicated at Bonferroni
  expect  <- n * K / N_MOTIF
  hyper_p <- if (n > 0) phyper(k - 1, K, N_MOTIF - K, n, lower.tail = FALSE) else NA_real_

  enrich_rows[[sub_name]] <- data.table(
    subtype = sub_name, analysis = "discovery_set_replication",
    discovery_n = nrow(orig), testable = n, background_sig = K,
    replicated = k, expected = round(expect, 2),
    ratio = round(k / expect, 2), same_direction = sum(conc$same_direction, na.rm = TRUE),
    replicated_nominal = sum(conc$replicated_nominal, na.rm = TRUE),
    replicated_FDR = sum(conc$replicated_FDR, na.rm = TRUE), hyper_p = signif(hyper_p, 3))

  ts(sprintf("%-16s %3d discovery | Bonf %3d repl (exp %.1f, %.2fx, p=%.3g) | %3d same-dir  [ref: nominal %3d, FDR %3d]",
             sub_name, nrow(orig), k, expect, k / expect, hyper_p,
             sum(conc$same_direction, na.rm = TRUE),
             sum(conc$replicated_nominal, na.rm = TRUE),
             sum(conc$replicated_FDR, na.rm = TRUE)))

  ## cell-type chromVAR enrichment, indexed BY MOTIF ID
  conc_valid <- as.data.table(conc)[!is.na(get(p_col))]
  if (!nrow(conc_valid)) { ts("  no testable TFs - chromVAR enrichment skipped"); next }
  conc_valid[, cv_row := match(motif, cv_ids)]
  enrichment <- sapply(cell_types, function(ctype) {
    cells <- intersect(meta_ct$cell[meta_ct$celltype == ctype], colnames(cv_scores))
    rows  <- conc_valid$cv_row; valid <- !is.na(rows)
    out   <- rep(NA_real_, nrow(conc_valid))
    if (length(cells) > 0)
      out[valid] <- rowMeans(cv_scores[rows[valid], cells, drop = FALSE])
    out
  })
  enrichment <- as.data.frame(enrichment); colnames(enrichment) <- cell_types
  fwrite(data.frame(motif = conc_valid$motif, TF = conc_valid$motif.name, enrichment),
         file.path("tfgwas/chromvar_enrichment",
                   paste0("replication_", sub_name, "_celltype_chromVAR.csv")))
}

# ─────────────────────────────────────────────
#  4. Subtype-UNIQUE motifs, by discovery's own rule
# ─────────────────────────────────────────────
## 04_celltype_enrichment/03_subtype_specific_motifs.R defines "unique" as:
##   Bonferroni-significant in THIS subtype (raw p < 0.05/727), AND not Bonferroni-significant
##   in ANY of the other four. Note it is a raw-p Bonferroni cut, NOT FDR.
## Applying the identical rule to the European results makes the two cohorts' unique sets
## directly comparable. It is far more specific than nominal p<0.05 because it removes the
## ubiquitous open-chromatin motifs that dominate the nominal background -- which is why the
## TNBC enrichment RISES from 1.68x to 3.92x under this definition rather than falling.
ts("=== SUBTYPE-UNIQUE MOTIFS (discovery rule: p < 0.05/727, not sig in any other subtype) ===")
ts(sprintf("Bonferroni threshold: %.3e", BONF))   # defined in section 3

sig_sets <- lapply(subtypes, function(s) {
  p <- res[[paste0(s, "_p")]]; res$motif[!is.na(p) & p < BONF] })
names(sig_sets) <- subtypes
uniq_sets <- lapply(subtypes, function(s)
  setdiff(sig_sets[[s]], unlist(sig_sets[setdiff(subtypes, s)])))
names(uniq_sets) <- subtypes

for (s in subtypes)
  res[, paste0(s, "_unique") := motif %in% uniq_sets[[s]]]
fwrite(res, "tfgwas/tfgwas_results_ALL.csv")   # rewrite with the *_unique flags

for (s in subtypes)
  ts(sprintf("  %-16s Bonferroni-sig: %3d | UNIQUE: %3d",
             s, length(sig_sets[[s]]), length(uniq_sets[[s]])))
ts(sprintf("  motifs significant in >=1 subtype: %d / %d",
           length(unique(unlist(sig_sets))), N_MOTIF))

## cross-cohort overlap of the two unique sets
for (sub_name in names(subtype_cfg)) {
  sig_path <- file.path("../results/singlecell/GWAS", subtype_cfg[[sub_name]]$sig_file)
  if (!file.exists(sig_path)) next
  d      <- unique(fread(sig_path, header = TRUE)$motif)
  e      <- uniq_sets[[sub_name]]
  shared <- intersect(d, e)
  expect <- length(e) * length(d) / N_MOTIF
  hyper_p <- phyper(length(shared) - 1, length(d), N_MOTIF - length(d),
                    length(e), lower.tail = FALSE)
  enrich_rows[[paste0(sub_name, "_unique")]] <- data.table(
    subtype = sub_name, analysis = "unique_vs_unique",
    discovery_n = length(d), testable = length(e), background_sig = NA_integer_,
    replicated = length(shared), expected = round(expect, 2),
    ratio = round(length(shared) / expect, 2), same_direction = NA_integer_,
    replicated_FDR = NA_integer_, hyper_p = signif(hyper_p, 3))
  ts(sprintf("  %-16s disc %3d | euro %3d | SHARED %3d (exp %.2f, %.2fx, p=%.3g)",
             sub_name, length(d), length(e), length(shared), expect,
             length(shared) / max(expect, 1e-12), hyper_p))
  if (length(shared))
    fwrite(res[motif %in% shared,
               c("motif", "motif.name", paste0(sub_name, c("_B", "_SE", "_p", "_FDR"))),
               with = FALSE][order(get(paste0(sub_name, "_p")))],
           file.path("tfgwas/concordance",
                     paste0("replication_", sub_name, "_unique_shared.csv")))
}

enrich <- rbindlist(enrich_rows, fill = TRUE)
fwrite(enrich, "tfgwas/replication_replication_enrichment.csv")
ts("enrichment summary -> tfgwas/replication_replication_enrichment.csv")
print(as.data.frame(enrich), row.names = FALSE)

# ─────────────────────────────────────────────
#  5. REPLICATION SUMMARY
# ─────────────────────────────────────────────
## For each subtype, how many of the DISCOVERY motifs replicate in this cohort, at three
## increasingly strict criteria:
##   p05   nominal p < 0.05                       <- MAIN REPLICATION RESULT
##   bonf  Bonferroni  p < 0.05/727
##   uniq  uniquely significant (Bonferroni here, AND not Bonferroni-significant in any
##         other subtype) -- discovery's own definition (03_subtype_specific_motifs.R)
ts("=== REPLICATION SUMMARY (main result: Bonferroni 0.05/727 in BOTH cohorts, and uniquely significant) ===")

summary_rows <- list()
for (sub_name in names(subtype_cfg)) {
  sig_path <- file.path("../results/singlecell/GWAS", subtype_cfg[[sub_name]]$sig_file)
  if (!file.exists(sig_path)) next
  d  <- unique(fread(sig_path, header = TRUE)$motif)
  p  <- res[[paste0(sub_name, "_p")]]; names(p) <- res$motif
  pd <- p[d]

  ## THRESHOLDS. BOTH cohorts use the SAME correction:
  ## Bonferroni p < 0.05/727 = 6.878e-05, exactly as the discovery set is defined and as the
  ## manuscript reports. The reported tiers are therefore:
  ##     bonf_727  Bonferroni 0.05/727 in both cohorts   <- PRIMARY
  ##     uniq      also uniquely significant in both      <- PRIMARY
  ## Two looser tiers are retained as columns only, for reference, and should NOT be quoted
  ## alongside the primary ones without stating their threshold:
  ##     p05       nominal 0.05 in the replication (asymmetric with discovery)
  ##     bonf_set  Bonferroni within the tested set, 0.05/n (corrects for tests performed
  ##               rather than the 727-motif universe; TNBC 22/57, enrichment P = 4.5e-07)
  BONF_SET <- 0.05 / length(d)
  k1 <- sum(pd < 0.05, na.rm = TRUE)                              # nominal
  k2 <- sum(pd < BONF_SET, na.rm = TRUE)                          # Bonferroni within set (reference)
  k4 <- sum(pd < BONF, na.rm = TRUE)                              # Bonferroni 0.05/727    <- PRIMARY
  k3 <- length(intersect(d, uniq_sets[[sub_name]]))               # uniquely significant

  summary_rows[[sub_name]] <- data.table(
    subtype       = sub_name,
    discovery_n   = length(d),
    thr_bonf_set  = signif(BONF_SET, 3),
    bonf_727      = k4, bonf_727_pct = round(100 * k4 / length(d)),
    uniq          = k3, uniq_pct     = round(100 * k3 / length(d)),
    p05           = k1, p05_pct      = round(100 * k1 / length(d)),
    bonf_set      = k2, bonf_set_pct = round(100 * k2 / length(d)))
}
repl <- rbindlist(summary_rows, fill = TRUE)
## Fixed manuscript order, NOT sorted by any metric: sorting by a percentage put whichever
## subtype happened to lead at the top, which reads as a ranking it is not (LuminalA leads on
## bonf_727_pct at 28% vs TNBC 26%, yet its enrichment is P = 0.44 against TNBC's 6.6e-05).
.ord <- c("TNBC", "LuminalA", "LuminalB", "LuminalB_HER2Neg", "HER2")
repl <- repl[order(match(subtype, .ord))]
fwrite(repl, "tfgwas/replication_replication_summary.csv")
ts("replication summary -> tfgwas/replication_replication_summary.csv")
print(as.data.frame(repl), row.names = FALSE)

for (i in seq_len(nrow(repl)))
  ts(sprintf("%-16s n=%3d | Bonf 0.05/727: %3d (%2d%%) | unique: %3d (%2d%%)   [ref: p<0.05 %3d (%2d%%), Bonf 0.05/n %3d (%2d%%)]",
             repl$subtype[i], repl$discovery_n[i],
             repl$bonf_727[i], repl$bonf_727_pct[i],
             repl$uniq[i],     repl$uniq_pct[i],
             repl$p05[i],      repl$p05_pct[i],
             repl$bonf_set[i], repl$bonf_set_pct[i]))

ts("=== 08_tfgwas_combine.R complete ===")
