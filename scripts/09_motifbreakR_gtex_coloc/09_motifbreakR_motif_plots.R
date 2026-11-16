# ==============================================================================
# 09_motifbreakR_motif_plots.R
# motifbreakR sequence-logo plots for promoter variants.
# Produces: Fig 5D-H
# Inputs  : results/motifbreakR/motifbreakR_TNBC_all.csv
#           results/tfgwas/tfgwas_significant_TNBC_unique.csv
# Run     : Rscript scripts/09_motifbreakR_gtex_coloc/09_motifbreakR_motif_plots.R   (from the repository root)
# ==============================================================================
## 09_motifbreakR_motif_plots.R — Fig 5D-H: motifbreakR sequence-logo plots (plotMB) for the
## promoter SNPs that alter TNBC-associated motifs.
##
## WHY A RE-RUN. plotMB() needs the full motifbreakR GRanges with the PWMs attached; the original
## run (03_motifbreakR_tnbc.R) saved only flattened CSVs. This script re-runs
## motifbreakR on the promoter SNPs with the ORIGINAL settings, verifies that the re-run reproduces
## motifbreakR_TNBC_all.csv (alleleDiff per SNP x motif x allele), then plots.
##   original settings: SNPlocs.Hsapiens.dbSNP155.GRCh37, BSgenome.Hsapiens.UCSC.hg19,
##   MotifDb dataSource "jaspar2022" (Hsapiens), threshold = 1e-4, method = "ic",
##   bkg = uniform, filterp = TRUE, then calculatePvalue(granularity = 1e-6).
## PWMs are scored independently, so restricting pwmList to the motifs of interest gives the same
## hits for those motifs as the full-JASPAR run.
##
## SNP SET (env var MOTIF_SET):
##   "specific"    (default) TNBC-specific motifs: 18 SNPs / 20 motifs / 10 genes
##   "significant" any TNBC Bonferroni-significant motif: 57 SNPs / 106 motifs / 20 genes
## Each plot shows only the SNP x motif pairs in the chosen set.
##
## Packages: motifbreakR, MotifDb and SNPlocs.Hsapiens.dbSNP155.GRCh37 (Bioconductor).
## ---- paths (set the environment variables; see README) ----
TNBC_ROOT <- Sys.getenv("TNBC_ROOT", ".")
suppressMessages({library(data.table); library(motifbreakR); library(MotifDb); library(BiocParallel)
                  library(SNPlocs.Hsapiens.dbSNP155.GRCh37); library(BSgenome.Hsapiens.UCSC.hg19)})
setwd(TNBC_ROOT)
OUT <- "results/target_genes/"; FIG <- "results/target_genes/figures/plotMB"
dir.create(OUT, showWarnings = FALSE, recursive = TRUE); dir.create(FIG, showWarnings = FALSE, recursive = TRUE)
SET <- Sys.getenv("MOTIF_SET", "specific"); stopifnot(SET %in% c("specific", "significant"))

## ---- SNP x motif pairs to plot ----
mb  <- fread("results/motifbreakR/motifbreakR_TNBC_all.csv")
u57 <- unique(fread("results/tfgwas/tfgwas_significant_TNBC_unique.csv")$motif)
p   <- mb[gsub('"', '', TNBC_sig) == "Y" & grepl("Promoter", annotation_category)]
if (SET == "specific") p <- p[providerName %in% u57]
## Intronic snoRNAs excluded: SNORD53 / SNORD92 lie inside WDR43, so their
## "promoter" SNPs are promoter only relative to the snoRNA start (intronic to WDR43), and snoRNAs
## are not reliably quantified by poly-A RNA-seq. The cached full motifbreakR object may still
## contain them; they are dropped here, before verification and plotting.
SNORNA <- c("SNORD53", "SNORD92")
p <- p[!SYMBOL %in% SNORNA]
## 1q32 (MDM4 locus) restricted to the selected SNPs: the 7 MDM4 promoter SNPs;
## the PIK3C2B promoter SNPs at the same locus are not plotted.
KEEP_1Q32 <- c("rs4245736", "rs4245737", "rs4252685", "rs4252687", "rs4951078", "rs4951391", "rs61817958")
p <- p[!SYMBOL %in% c("MDM4", "PIK3C2B") | SNP_id %in% KEEP_1Q32]
pairs <- unique(p[, .(SNP_id, gene = SYMBOL, REF, ALT, providerName, geneSymbol, effect,
                      alleleDiff_csv = alleleDiff, distanceToTSS)])
pairs[, key := paste(SNP_id, providerName, ALT)]
cat(sprintf("MOTIF_SET = %s | %d SNPs, %d motifs, %d genes, %d SNP x motif x allele pairs\n",
            SET, uniqueN(pairs$SNP_id), uniqueN(pairs$providerName), uniqueN(pairs$gene), nrow(pairs)))

## ---- re-run motifbreakR (cached) ----
rds_full <- file.path(OUT, sprintf("motifbreakR_promoterSNPs_%s_full.rds", SET))
if (file.exists(rds_full)) {
  results <- readRDS(rds_full)
} else {
  pwm <- subset(MotifDb, dataSource %in% c("jaspar2022") & mcols(MotifDb)$organism == "Hsapiens")
  pwm <- pwm[mcols(pwm)$providerName %in% unique(pairs$providerName)]
  cat("PWMs used:", length(pwm), "\n")
  snps.mb <- snps.from.rsid(rsid = unique(pairs$SNP_id), dbSNP = SNPlocs.Hsapiens.dbSNP155.GRCh37,
                            search.genome = BSgenome.Hsapiens.UCSC.hg19)
  results <- motifbreakR(snpList = snps.mb, filterp = TRUE, pwmList = pwm, threshold = 1e-4,
                         method = "ic", bkg = c(A = 0.25, C = 0.25, G = 0.25, T = 0.25),
                         BPPARAM = SerialParam())
  results <- calculatePvalue(results, granularity = 1e-6)
  saveRDS(results, rds_full)
}

## ---- verify the re-run reproduces the published motifbreakR table ----
R <- data.table(SNP_id = results$SNP_id, providerName = results$providerName,
                ALT = as.character(results$ALT), alleleDiff_rerun = results$alleleDiff,
                effect_rerun = results$effect)
R[, key := paste(SNP_id, providerName, ALT)]
V <- merge(pairs, R[, .(key, alleleDiff_rerun, effect_rerun)], by = "key", all.x = TRUE)
V[, match := !is.na(alleleDiff_rerun) & abs(alleleDiff_rerun - alleleDiff_csv) < 1e-3]
fwrite(V[, !"key"], file.path(OUT, sprintf("plotMB_pairs_verification_%s.csv", SET)))
cat(sprintf("verification: %d/%d pairs reproduced (alleleDiff within 1e-3); %d missing from re-run; %d differ\n",
            sum(V$match), nrow(V), sum(is.na(V$alleleDiff_rerun)),
            sum(!is.na(V$alleleDiff_rerun) & !V$match)))
if (any(!V$match)) print(as.data.frame(V[match == FALSE, .(SNP_id, gene, providerName, geneSymbol, ALT,
                                                          alleleDiff_csv, alleleDiff_rerun)]), row.names = FALSE)

## ---- keep only the pairs in the set; carry the motifbreakR attributes plotMB needs ----
sub <- results[paste(results$SNP_id, results$providerName, as.character(results$ALT)) %in% pairs$key]
for (a in setdiff(names(attributes(results)), names(attributes(sub)))) attr(sub, a) <- attr(results, a)
saveRDS(sub, file.path(OUT, sprintf("motifbreakR_promoterSNPs_%s_plotted.rds", SET)))

## ---- one plot per SNP x alternative allele ----
todo <- unique(pairs[key %in% paste(sub$SNP_id, sub$providerName, as.character(sub$ALT)),
                     .(SNP_id, gene, REF, ALT)])
setorder(todo, gene, SNP_id)
## per-SNP figure width (inches); default 7 (four MDM4 plots made narrower)
WIDTH_SNP <- setNames(rep(5.5, 7), c("rs4245736", "rs4245737", "rs4252685", "rs4252687",
                                     "rs4951078", "rs4951391", "rs61817958"))   # all 7 MDM4 SNPs, uniform
## rs4951391's ZNF140 motif is 21 bp, so its sequence window is longer: at 5.5 in Gviz drops the
## sequence letters for coloured boxes. It needs more width to keep the letters.
WIDTH_SNP[["rs4951391"]] <- 6.5
## optional: redraw only these rsIDs, e.g. PLOTMB_ONLY=rs4951391,rs61817958
only <- strsplit(Sys.getenv("PLOTMB_ONLY", ""), ",")[[1]]
if (length(only)) todo <- todo[SNP_id %in% only]
for (i in seq_len(nrow(todo))) {
  r <- todo$SNP_id[i]; alt <- todo$ALT[i]
  s <- sub[sub$SNP_id == r & as.character(sub$ALT) == alt]
  n <- length(s); eff <- unique(s$effect)
  f <- file.path(FIG, sprintf("plotMB_%s_%s_%s-%s.png", todo$gene[i], r, todo$REF[i], alt))
  w <- if (r %in% names(WIDTH_SNP)) WIDTH_SNP[[r]] else 7
  png(f, width = w, height = 2.2 + 1.0 * n, units = "in", res = 300)
  ok <- tryCatch({ plotMB(results = sub, rsid = r, effect = eff, altAllele = alt); TRUE },
                 error = function(e) { message(r, ": ", conditionMessage(e)); FALSE })
  dev.off()
  if (!ok) unlink(f)
  cat(sprintf("%-8s %-12s %s>%s  motifs = %2d  %s\n", todo$gene[i], r, todo$REF[i], alt, n,
              if (ok) basename(f) else "FAILED"))

  ## optional (PLOTMB_PER_MOTIF=1): one figure per SNP x motif, in figures/plotMB/per_motif/
  if (Sys.getenv("PLOTMB_PER_MOTIF", "") == "1") {
    FIG_M <- file.path(FIG, "per_motif"); dir.create(FIG_M, showWarnings = FALSE)
    for (pn in unique(s$providerName)) {
      s1 <- s[s$providerName == pn]
      for (a in setdiff(names(attributes(sub)), names(attributes(s1)))) attr(s1, a) <- attr(sub, a)
      tf <- unique(s1$geneSymbol)[1]
      f1 <- file.path(FIG_M, sprintf("plotMB_%s_%s_%s-%s_%s_%s.png", todo$gene[i], r, todo$REF[i], alt,
                                     gsub("[^A-Za-z0-9]+", "-", tf), pn))
      png(f1, width = w, height = 3.2, units = "in", res = 300)
      ok1 <- tryCatch({ plotMB(results = s1, rsid = r, effect = unique(s1$effect), altAllele = alt); TRUE },
                      error = function(e) { message(r, " ", pn, ": ", conditionMessage(e)); FALSE })
      dev.off()
      if (!ok1) unlink(f1)
      cat(sprintf("    per motif: %-8s %s\n", tf, if (ok1) basename(f1) else "FAILED"))
    }
  }
}
cat("\nwrote plots to", FIG, "and motifbreakR objects to", OUT, "\n")
