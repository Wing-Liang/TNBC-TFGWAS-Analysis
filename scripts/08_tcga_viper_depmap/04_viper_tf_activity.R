# ==============================================================================
# 04_viper_tf_activity.R
# msVIPER activity of the 61 TFs (ARACNe-BRCA regulons), TCGA-BRCA.
# Produces: Fig 4A-B, Table S10
# Inputs  : results/tfgwas/tfgwas_significant_TNBC_unique.csv
#           results/tcga/TCGA_BRCA_data_all.rds
#           results/tcga/TCGA_BRCA_PAM50_metadata.rds
# Run     : Rscript scripts/08_tcga_viper_depmap/04_viper_tf_activity.R   (from the repository root)
# ==============================================================================
## 04_viper_tf_activity.R — Fig 4A-B: msVIPER regulon activity (ARACNe-BRCA) of the 61 TF genes that
## encode the 57 TNBC-specific motifs, TCGA-BRCA, on DESeq2 VARIANCE-STABILISED counts.
##
## WHY VST (rather than log2(TPM+1)):
##   - VIPER itself only needs a continuous expression matrix; nothing in the package or paper
##     prescribes TPM (viper vignette: "a gene expression matrix or ExpressionSet").
##   - The ARACNe-BRCA network was reverse-engineered from TCGA RNA-seq whose raw counts were
##     library-size normalised and variance-stabilised with DESeq (Alvarez et al. 2016, Nat Genet
##     48:838, Online Methods, doi:10.1038/ng.3593). Scoring it on DESeq2-VST data matches the network.
##   - The same raw counts, samples and pre-filter as the DESeq2 analysis (02_tcga_deseq2_tf_genes.R), so the TCGA
##     expression panel and the VIPER panel share one normalisation.
##
## PIPELINE (each step is checkpointed in results/tf_activity/viper_vst/steps/):
##   1 samples     the same 836: Primary Tumor + Solid Tissue Normal (NAT), PAM50 Basal/Her2/
##                 LumA/LumB, female, white (TCGA_BRCA_PAM50_metadata.rds; NAT already recoded)
##   2 VST         DESeqDataSetFromMatrix(STAR "unstranded" counts, design ~ group) -> DESeq2 pre-filter
##                 (>= 1 count in >= 5 samples) -> estimateSizeFactors -> vst(blind = FALSE)
##                 blind = FALSE: dispersion trend fitted with the group design, as DESeq2's vignette
##                 advises when the transformed data are used for downstream group comparisons
##                 (blind = TRUE would treat real subtype differences as noise and over-shrink).
##   3 genes       Ensembl -> gene_name (GENCODE v36 rowData); _PAR_Y rows dropped; symbols upper-case;
##                 symbol with > 1 Ensembl row -> the row with the highest mean VST kept
##   4 filter DESeq2-vignette pre-filter applied to the genes entering
##                 VIPER: raw count >= COUNT_MIN (10) in >= N_MIN samples, N_MIN = smallest group (Her2, 39)
##                 [DESeq2 vignette: rowSums(counts(dds) >= 10) >= smallestGroupSize]. 56,633 -> 28,821 genes;
##                 16,992 of the 19,307 ARACNe-BRCA interactome genes kept. It removes every gene that is
##                 constant in both groups of a contrast or of a permutation (undefined t-test = NaN): a
##                 no-filter run gave NaN z for 188 genes (Her2 signature) and 4,271 genes (LumB null).
##                 The VST itself is computed on the DESeq2 gene set (step 2) and is not affected.
##                 msviper's default ges.filter = TRUE then restricts the signature to interactome genes.
##   5 signature   viper::rowTtest(Basal, comparator) on VST -> z = qnorm(p/2) * sign(t)
##   6 null        viper::ttestNull(per = 1000, repos = TRUE), seed 42 per contrast
##   7 regulon     ARACNe-BRCA (Entrez -> symbol), built by 03_build_aracne_regulon.R;
##                 normalisation-independent, reused from results/tf_activity/viper/
##   8 msviper     minsize = 5; every TF also flagged `reliable_size` = tested_size >= 25, the size below
##                 which Alvarez et al. 2016 report accuracy drops sharply
##   9 outputs     tables (NES, p, FDR_61, FDR_network, tested_size, reliable_size), leading edges,
##                 viper's plot(msviper) per contrast (61 TFs, top N_PLOT by p), summary dot plot, and a
##                 VST-vs-TPM comparison (NES correlation, agreement of p < 0.05 calls)
##  10 extras      leading-edge (vignette 5.5), shadow (6.1) and synergy (6.2) analyses -- see STEP 10
## Outputs go to results/tf_activity/{data,figures}/viper_vst/min<MINSIZE>/ (shared caches in data/viper_vst/steps/).
## Contrasts: Basal vs LumA / LumB / Her2 / NAT.
## Caveat for Methods: ARACNe-BRCA was inferred from TCGA-BRCA, the cohort scored here.
## Run with Rscript (DESeq2 1.48.2, viper 1.44.0). Expect ~45 min (4 null models).
## ---- paths (set the environment variables; see README) ----
TNBC_ROOT <- Sys.getenv("TNBC_ROOT", ".")
suppressMessages({
  library(SummarizedExperiment); library(DESeq2); library(viper)
  library(data.table); library(ggplot2)
})
setwd(TNBC_ROOT)
MS_TAG <- sprintf("f10n39_min%s", Sys.getenv("VIPER_MINSIZE", "5"))   # filter + minsize tag (caches, outputs)
OUT  <- file.path("results/tf_activity/viper_vst", MS_TAG); STEP <- "results/tf_activity/viper_vst/steps"
FIG  <- file.path("results/tf_activity/figures/viper_vst", MS_TAG); FIGP <- file.path(FIG, "msviper_plot")
for (d in c(OUT, STEP, FIG, FIGP)) dir.create(d, showWarnings = FALSE, recursive = TRUE)
F_REG   <- "results/tf_activity/viper/aracne_brca_regulon_V2_symbol.rds"          # step 7 input
F_TPM   <- "results/tf_activity/viper/VIPER_msviper_61TF_4contrasts_ARACNe_BRCA.csv" # optional log2(TPM+1) sensitivity run; comparison skipped if absent
## minsize: TFs with fewer regulon targets in the signature are not tested. Set at run time:
##   VIPER_MINSIZE=25 (viper default; Alvarez et al. 2016) or 5 (max coverage, flagged). Caches are per value.
MINSIZE <- as.integer(Sys.getenv("VIPER_MINSIZE", "5")); MIN_RELIABLE <- 25
SH_TOP  <- 25     # shadow: genome-wide top regulators added to the 61-set TFs with p < 0.05 (vignette uses 25)
COUNT_MIN <- 10   # step-4 filter: raw count >= COUNT_MIN ...
F_TAG   <- "f10n39"  # ... in >= N_MIN samples (N_MIN = smallest group, checked = 39 below); tags sig/null caches
N_PLOT  <- as.integer(Sys.getenv("VIPER_NPLOT", "15"))
LABEL <- c(LumA = "Basal vs LumA", LumB = "Basal vs LumB", Her2 = "Basal vs Her2", NAT = "Basal vs NAT")
t0 <- Sys.time()
ts <- function(...) cat(sprintf("[+%.1fm] ", as.numeric(difftime(Sys.time(), t0, units = "mins"))), ..., "\n", sep = "")
ck <- function(f, fn) { if (file.exists(f)) { ts("  cached: ", basename(f)); return(readRDS(f)) }; v <- fn(); saveRDS(v, f); v }

## ---- the 61 TF genes ----
motifs <- unique(fread("results/tfgwas/tfgwas_significant_TNBC_unique.csv")$motif.name)
TF61   <- toupper(unique(unlist(strsplit(motifs, "::"))))
stopifnot(length(motifs) == 57, length(TF61) == 61)

## ---- STEP 1: samples (identical to 02_tcga_deseq2_tf_genes.R) ----
data <- readRDS("results/tcga/TCGA_BRCA_data_all.rds")
metadata <- readRDS("results/tcga/TCGA_BRCA_PAM50_metadata.rds")
colData(data) <- DataFrame(metadata, row.names = rownames(metadata))
keep <- metadata$sample_type %in% c("Primary Tumor", "Solid Tissue Normal") &
  metadata$paper_BRCA_Subtype_PAM50 %in% c("Basal", "Her2", "LumA", "LumB", "NAT") &
  metadata$race %in% "white" & metadata$gender %in% "female"
se  <- data[, keep]; rm(data)
grp <- factor(se$paper_BRCA_Subtype_PAM50, levels = c("Basal", "LumA", "LumB", "Her2", "NAT"))
ts("STEP 1 samples: ", paste(names(table(grp)), table(grp), sep = " = ", collapse = ", "))
stopifnot(sum(keep) == 836, as.vector(table(grp)) == c(116, 446, 131, 39, 104))

## ---- STEP 2: DESeq2 VST ----
V <- ck(file.path(STEP, "vst_836_blindFALSE.rds"), function() {
  dds <- DESeqDataSetFromMatrix(countData = assay(se, "unstranded"),
                                colData = data.frame(group = grp, row.names = colnames(se)), design = ~ group)
  dds <- dds[rowSums(counts(dds) >= 1) >= 5, ]          # DESeq2 pre-filter (the only gene filter)
  dds <- estimateSizeFactors(dds)
  vs  <- vst(dds, blind = FALSE)
  list(vst = assay(vs), size_factors = sizeFactors(dds))
})
ts("STEP 2 VST: ", nrow(V$vst), " Ensembl genes x ", ncol(V$vst), " samples | size factors ",
   sprintf("%.2f-%.2f", min(V$size_factors), max(V$size_factors)))

## ---- STEP 4: count filter (DESeq2 vignette), applied to the genes entering VIPER ----
N_MIN <- min(table(grp)); stopifnot(N_MIN == 39)
pass  <- rowSums(assay(se, "unstranded")[rownames(V$vst), ] >= COUNT_MIN) >= N_MIN
ts("STEP 4 filter (>= ", COUNT_MIN, " counts in >= ", N_MIN, " samples): ", sum(pass), " of ", length(pass), " genes kept")

## ---- STEP 3: Ensembl -> symbol, one row per symbol (on the filtered genes) ----
rd  <- rowData(se)[rownames(V$vst)[pass], c("gene_id", "gene_name", "gene_type")]
ann <- data.table(ens = rownames(V$vst)[pass], sym = toupper(rd$gene_name), mean_vst = rowMeans(V$vst[pass, ]))
ann <- ann[!grepl("_PAR_Y$", ens) & !is.na(sym) & nzchar(sym)]
n_dupsym <- ann[, .N, by = sym][N > 1, .N]
ann <- ann[order(sym, -mean_vst)][!duplicated(sym)]
ts("STEP 3 symbols: ", nrow(ann), " unique symbols (", n_dupsym, " symbols had > 1 Ensembl row; highest-mean row kept)")
X <- V$vst[ann$ens, ]; rownames(X) <- ann$sym

## ---- STEP 7 (loaded early): regulon ----
if (!file.exists(F_REG)) stop("missing ", F_REG, " -- run 03_build_aracne_regulon.R first")
reg <- readRDS(F_REG)
ts("STEP 7 regulon ARACNe-BRCA V2: ", length(reg), " regulators | 61-set present: ", sum(TF61 %in% names(reg)),
   " | absent: ", paste(setdiff(TF61, names(reg)), collapse = ", "))

## ---- STEPS 5-6, 8 per contrast (filtered genes) ----
rows <- list(); MRS <- list(); LED <- list()
for (cn in names(LABEL)) {
  ts("=== ", LABEL[[cn]], " ===")
  E1 <- X[, grp == "Basal"]; E2 <- X[, grp == cn]
  ts("  genes: ", nrow(E1), " (in ARACNe-BRCA interactome: ",
     sum(rownames(E1) %in% unique(c(names(reg), unlist(lapply(reg, function(r) names(r$tfmode)))))),
     ") | n Basal = ", ncol(E1), " | n ", cn, " = ", ncol(E2))
  sig <- ck(file.path(STEP, sprintf("sig_%s_%s.rds", cn, F_TAG)), function() {           # STEP 5
    st <- viper::rowTtest(E1, E2)
    setNames((qnorm(st$p.value / 2, lower.tail = FALSE) * sign(st$statistic))[, 1], rownames(E1)) })
  nul <- ck(file.path(STEP, sprintf("null_%s_%s.rds", cn, F_TAG)), function() {          # STEP 6
    set.seed(42); viper::ttestNull(E1, E2, per = 1000, repos = TRUE, verbose = FALSE) })
  ## Non-finite z: with no expression filter, a gene that is constant in both groups
  ## (e.g. zero counts -> VST floor in every sample) has an undefined t-test (NaN). This happens in the
  ## signature (LumB 67 genes, Her2 188) and, far more often, in the permuted null (a rare gene becomes
  ## constant whenever a permutation puts all its non-zero samples in one group: LumA 1,731 genes, LumB
  ## 4,271). aREA ranks genes with base::rank(), which puts NaN LAST -- i.e. treats them as the most
  ## up-regulated -- distorting the enrichment and the null. A constant gene carries no evidence of a
  ## difference, so its z is set to 0 (applied in memory; the cached sig/null stay raw).
  n_bad_sig <- sum(!is.finite(sig)); n_bad_null <- sum(rowSums(!is.finite(nul)) > 0)
  sig[!is.finite(sig)] <- 0; nul[!is.finite(nul)] <- 0
  ts("  STEP 5 signature: ", length(sig), " genes, z ", sprintf("%.1f to %.1f", min(sig), max(sig)),
     " | STEP 6 null: ", nrow(nul), " x ", ncol(nul),
     " | non-finite z set to 0: signature ", n_bad_sig, " genes, null ", n_bad_null, " genes (>= 1 permutation)")
  mrs <- ck(file.path(STEP, sprintf("msviper_%s_%s.rds", cn, MS_TAG)), function()        # STEP 8
    viper::msviper(sig, reg, nul, minsize = MINSIZE, verbose = FALSE))
  MRS[[cn]] <- mrs
  r <- data.table(TF = names(mrs$es$nes), NES = as.numeric(mrs$es$nes), p = as.numeric(mrs$es$p.value),
                  tested_size = as.integer(mrs$es$size[names(mrs$es$nes)]))
  r[, FDR_network := p.adjust(p, "BH")]
  r <- r[TF %in% TF61]; r[, FDR_61 := p.adjust(p, "BH")]
  r[, `:=`(contrast = cn, reliable_size = tested_size >= MIN_RELIABLE)]
  rows[[cn]] <- r
  ts(sprintf("  STEP 8 tested %d/61 (%d with >= %d targets) | p<0.05: %d | FDR_61<0.05: %d | FDR_network<0.05: %d",
             nrow(r), sum(r$reliable_size), MIN_RELIABLE, sum(r$p < 0.05), sum(r$FDR_61 < 0.05), sum(r$FDR_network < 0.05)))
  led <- ck(file.path(STEP, sprintf("ledge_%s_%s.rds", cn, MS_TAG)), function() viper::ledge(mrs))
  LED[[cn]] <- led
  lg <- led$ledge[names(led$ledge) %in% TF61 & lengths(led$ledge) > 0]
  if (length(lg)) fwrite(rbindlist(lapply(names(lg), function(tf) data.table(TF = tf, rank = seq_along(lg[[tf]]), gene = lg[[tf]]))),
                         file.path(OUT, sprintf("ledge_61TF_%s.csv", cn)))
}
R <- rbindlist(rows)
setcolorder(R, c("TF", "contrast", "NES", "p", "FDR_61", "FDR_network", "tested_size", "reliable_size"))
fwrite(R, file.path(OUT, "VIPER_VST_msviper_61TF_4contrasts_ARACNe_BRCA.csv"))
cat("61-set TFs never tested:", paste(setdiff(TF61, R$TF), collapse = ", "), "\n")
cat("tested TFs with < ", MIN_RELIABLE, " targets in any contrast: ",
    paste(sort(unique(R[reliable_size == FALSE, TF])), collapse = ", "), "\n", sep = "")

## ---- STEP 9a: viper's msviper plot (vignette 5.4), 61 TFs, top N_PLOT by p ----
for (cn in names(MRS)) {
  mrs <- MRS[[cn]]; r <- rows[[cn]][order(p)]
  sel <- head(r$TF, N_PLOT); sel <- sel[order(mrs$es$nes[sel], decreasing = TRUE)]
  f <- file.path(FIGP, sprintf("msviper_VST_61TF_%s.png", cn))
  png(f, width = 7, height = 1.4 + 0.3 * length(sel), units = "in", res = 300)
  par(oma = c(0, 0, 2.2, 0))
  ok <- tryCatch({ plot(mrs, mrs = sel, cex = 0.7); TRUE }, error = function(e) { message(cn, ": ", conditionMessage(e)); FALSE })
  if (ok) mtext(sprintf("%s  |  ARACNe-BRCA  |  DESeq2 VST  |  61 TNBC-specific TF genes, top %d by p", LABEL[[cn]], length(sel)),
                side = 3, outer = TRUE, line = 0.6, cex = 0.75, font = 2)
  dev.off(); if (!ok) unlink(f) else ts("wrote ", f)
}

## ---- STEP 9b: summary dot plot (NES fill, -log10 p size, ring p < 0.05, * FDR_61 < 0.05, † < 25 targets) ----
d <- copy(R)
d[, contrast := factor(LABEL[contrast], levels = LABEL)]
small <- d[, .(small = any(!reliable_size)), by = TF]
d <- merge(d, small, by = "TF")
d[, TF_lab := paste0(TF, fifelse(small, " †", ""))]
ord <- d[, .(m = mean(NES)), by = TF_lab][order(m), TF_lab]
d[, TF_lab := factor(TF_lab, levels = ord)]
lim <- max(abs(d$NES))
p <- ggplot(d, aes(contrast, TF_lab)) +
  geom_point(aes(size = -log10(p), fill = NES, colour = p < 0.05), shape = 21, stroke = 0.6) +
  geom_text(data = d[FDR_61 < 0.05], aes(label = "*"), nudge_x = 0.28, size = 4) +
  scale_fill_gradient2(low = "#3B4CC0", mid = "white", high = "#B40426", limits = c(-lim, lim)) +
  scale_colour_manual(values = c(`TRUE` = "black", `FALSE` = "grey75"), name = "p < 0.05") +
  scale_size_continuous(range = c(0.8, 5), name = "-log10 p") +
  labs(x = NULL, y = "TF gene", title = "TF activity, Basal vs comparator (msVIPER, ARACNe-BRCA, DESeq2 VST)",
       caption = paste0("TCGA-BRCA female & white; ring = p < 0.05; * = FDR < 0.05 within the 61-gene set;\n",
                        "† = < 25 regulon targets tested in at least one contrast (low confidence; Alvarez et al. 2016)")) +
  theme_bw(base_size = 10) +
  theme(axis.text = element_text(colour = "black"), axis.text.x = element_text(angle = 40, hjust = 1),
        panel.grid.minor = element_blank(), plot.caption = element_text(hjust = 0, size = 8),
        plot.title.position = "plot", plot.caption.position = "plot")
f <- file.path(FIG, "VIPER_VST_61TF_summary_ARACNe_BRCA.png")
ggsave(f, p, width = 6.8, height = 2.2 + 0.2 * uniqueN(d$TF), dpi = 300, bg = "white"); ts("wrote ", f)

## ---- STEP 9c: VST vs TPM (sensitivity) ----
if (file.exists(F_TPM)) {
  tp <- fread(F_TPM)[, .(TF, contrast, NES_TPM = NES, p_TPM = p)]
  C <- merge(R[, .(TF, contrast, NES_VST = NES, p_VST = p)], tp, by = c("TF", "contrast"))
  S <- C[, .(n = .N, spearman_NES = round(cor(NES_VST, NES_TPM, method = "spearman"), 3),
             sig_both = sum(p_VST < 0.05 & p_TPM < 0.05), sig_VST_only = sum(p_VST < 0.05 & p_TPM >= 0.05),
             sig_TPM_only = sum(p_VST >= 0.05 & p_TPM < 0.05),
             sign_flip_sig = sum(sign(NES_VST) != sign(NES_TPM) & pmin(p_VST, p_TPM) < 0.05)), by = contrast]
  fwrite(C, file.path(OUT, "VIPER_VST_vs_TPM_61TF.csv")); fwrite(S, file.path(OUT, "VIPER_VST_vs_TPM_summary.csv"))
  cat("\nVST vs TPM (sensitivity):\n"); print(S, row.names = FALSE)
  C[, contrast := factor(LABEL[contrast], levels = LABEL)]
  g <- ggplot(C, aes(NES_TPM, NES_VST)) + geom_abline(linetype = 2, colour = "grey60") +
    geom_point(aes(colour = pmin(p_VST, p_TPM) < 0.05), size = 1.6) +
    scale_colour_manual(values = c(`TRUE` = "#B40426", `FALSE` = "grey55"), name = "p < 0.05 in either") +
    facet_wrap(~ contrast, nrow = 1) + coord_equal() +
    labs(x = "NES, log2(TPM + 1)", y = "NES, DESeq2 VST", title = "msVIPER NES: DESeq2 VST vs log2(TPM + 1), 61 TF genes") +
    theme_bw(base_size = 10) + theme(legend.position = "bottom")
  f <- file.path(FIG, "VIPER_VST_vs_TPM_NES.png"); ggsave(f, g, width = 11, height = 3.8, dpi = 300, bg = "white"); ts("wrote ", f)
} else cat("TPM results not found (", F_TPM, ") -- comparison skipped\n")

## ---- STEP 10: leading-edge, shadow and synergy analyses (viper vignette 5.5, 6.1, 6.2) ----
## 10a leading edge  viper::ledge + summary(): the targets that drive each TF's enrichment (GSEA-style
##                   leading edge). Reported for the 61 TFs, ordered by p.
## 10b shadow        viper::shadow: is a TF's enrichment explained by targets it shares with another, more
##                   strongly enriched regulator ("A -> B": part of B's apparent activity is due to A)?
##                   Regulators examined = 61-set TFs with p < 0.05 PLUS the genome-wide top SH_TOP
##                   regulators (so a 61-set TF can be shadowed by any strong master regulator, e.g. FOXA1,
##                   ESR1, E2F1). Package defaults: >= 10 shared targets, shadow p 0.01, 1,000 permutations.
##                   Reported: the shadow pairs, and each 61-set TF's p before vs after the shadow correction.
## 10c synergy       viper::msviperCombinatorial + msviperSynergy: do PAIRS of 61-set TFs (p < 0.05 in that
##                   contrast) have a co-regulon (shared targets) enriched beyond the union of their regulons?
##                   level = 2 (pairs only); co-regulons must reach minsize; synergy p from 1,000 permutations.
## All use the null model stored in the msviper object (same 1,000 sample permutations as step 6).
hr10 <- function(s) cat("\n-------------------- ", s, " --------------------\n", sep = "")
LE <- list(); SHP <- list(); SHR <- list(); SYN <- list()
for (cn in names(MRS)) {
  mrs <- MRS[[cn]]; led <- LED[[cn]]; r <- rows[[cn]][order(p)]; nom <- r[p < 0.05, TF]

  ## 10a leading edge
  hr10(sprintf("10a %s: leading edge, 61 TFs (viper summary)", LABEL[[cn]]))
  s <- as.data.table(summary(led, mrs = r$TF))
  s[, n_ledge := lengths(led$ledge[Regulon])][, contrast := cn]
  LE[[cn]] <- s
  print(head(as.data.frame(s), 15), row.names = FALSE)

  ## 10b shadow
  hr10(sprintf("10b %s: shadow analysis", LABEL[[cn]]))
  if (!length(nom)) { cat("no 61-set TF with p < 0.05 -- shadow skipped\n") } else {
    reg_sh <- union(nom, names(sort(mrs$es$p.value))[seq_len(min(SH_TOP, length(mrs$es$p.value)))])
    msh <- ck(file.path(STEP, sprintf("shadow_%s_%s.rds", cn, MS_TAG)), function()
      viper::shadow(mrs, regulators = reg_sh, targets = 10, shadow = 0.01, per = 1000, seed = 1, verbose = FALSE))
    pr <- if (!is.null(msh$shadow) && length(msh$shadow)) {
      x <- as.data.table(as.data.frame(msh$shadow, stringsAsFactors = FALSE)); setnames(x, 1:2, c("shadowing", "shadowed")); x[, 1:2]
    } else data.table(shadowing = character(), shadowed = character())
    pr[, `:=`(contrast = cn, shadowed_in61 = shadowed %in% TF61, shadowing_in61 = shadowing %in% TF61)]
    SHP[[cn]] <- pr
    af <- data.table(TF = nom, NES_before = mrs$es$nes[nom], p_before = mrs$es$p.value[nom],
                     NES_after_shadow = msh$es$nes[nom], p_after_shadow = msh$es$p.value[nom])
    af[, shadowed_by := vapply(TF, function(t) paste(pr[shadowed == t, shadowing], collapse = ", "), "")]
    af[, contrast := cn]; SHR[[cn]] <- af
    cat(sprintf("regulators examined: %d (%d 61-set TFs with p < 0.05 + genome-wide top %d) | shadow pairs: %d, involving a 61-set TF as shadowed: %d\n",
                length(reg_sh), length(nom), SH_TOP, nrow(pr), sum(pr$shadowed_in61)))
    print(as.data.frame(af[, .(TF, NES_before = round(NES_before, 2), p_before = signif(p_before, 2),
                               NES_after_shadow = round(NES_after_shadow, 2), p_after_shadow = signif(p_after_shadow, 2), shadowed_by)]), row.names = FALSE)
  }

  ## 10c synergy
  hr10(sprintf("10c %s: synergy analysis (pairs of 61-set TFs with p < 0.05)", LABEL[[cn]]))
  if (length(nom) < 2) { cat("fewer than 2 61-set TFs with p < 0.05 -- synergy skipped\n"); next }
  ## msviperCombinatorial STOPS with "No significant co-regulons to evaluate" when no pair of the selected
  ## TFs shares enough targets to form a significant co-regulon (seen for Basal vs LumA, minsize 25). That is
  ## a result ("no synergy testable"), not a failure, so it is caught and recorded instead of aborting.
  syn <- tryCatch(ck(file.path(STEP, sprintf("synergy_%s_%s.rds", cn, MS_TAG)), function() {
    m <- viper::msviperCombinatorial(mrs, regulators = nom, level = 2, verbose = FALSE)
    viper::msviperSynergy(m, per = 1000, seed = 1, verbose = FALSE) }),
    error = function(e) { cat("synergy not testable:", conditionMessage(e), "\n"); NULL })
  if (is.null(syn)) { SYN[[cn]] <- data.table(Regulon = NA_character_, contrast = cn, note = "no significant co-regulon"); next }
  combos <- grep("--", names(syn$es$nes), value = TRUE, fixed = TRUE)
  cat(sprintf("TFs: %d -> possible pairs: %d | co-regulons tested (>= minsize shared targets, significant): %d\n",
              length(nom), choose(length(nom), 2), length(combos)))
  if (!length(combos)) { cat("no co-regulon tested\n"); next }
  ## summary.msviper() treats a length-1 `mrs` as "top N" and calls round() on it, so a single co-regulon
  ## name crashes it (seen: Basal vs NAT, minsize 5). Pass the TFs as well, then keep only the co-regulon rows.
  st <- as.data.table(summary(syn, mrs = c(nom, combos)))[Regulon %in% combos]; st[, contrast := cn]
  SYN[[cn]] <- st
  print(as.data.frame(st), row.names = FALSE)
  sel <- c(nom[order(mrs$es$nes[nom], decreasing = TRUE)], head(combos[order(syn$es$p.value[combos])], 15))
  f <- file.path(FIG, sprintf("msviper_VST_synergy_%s.png", cn))
  png(f, width = 7.5, height = 1.4 + 0.3 * length(sel), units = "in", res = 300); par(oma = c(0, 0, 2.2, 0))
  ok <- tryCatch({ plot(syn, mrs = sel, cex = 0.7); TRUE }, error = function(e) { message("synergy plot ", cn, ": ", conditionMessage(e)); FALSE })
  if (ok) mtext(sprintf("%s  |  synergy: 61-set TFs (p < 0.05) and their co-regulons (A--B)", LABEL[[cn]]),
                side = 3, outer = TRUE, line = 0.6, cex = 0.75, font = 2)
  dev.off(); if (!ok) unlink(f) else ts("wrote ", f)
}
LEt <- rbindlist(LE, fill = TRUE);   fwrite(LEt, file.path(OUT, "VIPER_VST_61TF_leading_edge_summary.csv"))
SHPt <- rbindlist(SHP, fill = TRUE); fwrite(SHPt, file.path(OUT, "VIPER_VST_shadow_pairs.csv"))
SHRt <- rbindlist(SHR, fill = TRUE); fwrite(SHRt, file.path(OUT, "VIPER_VST_61TF_shadow_before_after.csv"))
SYNt <- rbindlist(SYN, fill = TRUE); fwrite(SYNt, file.path(OUT, "VIPER_VST_61TF_synergy.csv"))

## one row per TF x contrast: msviper + leading edge + shadow + synergy partners
ALL <- merge(R, LEt[, .(TF = Regulon, contrast, n_ledge, Ledge)], by = c("TF", "contrast"), all.x = TRUE)
if (nrow(SHRt)) ALL <- merge(ALL, SHRt[, .(TF, contrast, p_after_shadow, shadowed_by)], by = c("TF", "contrast"), all.x = TRUE)
if (nrow(SYNt)) {
  scol <- grep("^Synergy$", names(SYNt), value = TRUE)
  pairs <- SYNt[, .(a = sub("--.*$", "", Regulon), b = sub("^.*--", "", Regulon), contrast,
                    syn_p = if (length(scol)) get(scol) else NA_real_)]
  part <- rbind(pairs[, .(TF = a, partner = b, contrast, syn_p)], pairs[, .(TF = b, partner = a, contrast, syn_p)])
  part <- part[, .(synergy_partners = paste(sprintf("%s (p=%s)", partner, signif(syn_p, 2)), collapse = "; ")), by = .(TF, contrast)]
  ALL <- merge(ALL, part, by = c("TF", "contrast"), all.x = TRUE)
}
fwrite(ALL, file.path(OUT, "VIPER_VST_61TF_msviper_ledge_shadow_synergy.csv"))
ts("wrote leading-edge / shadow / synergy tables to ", OUT)

cat("\nnominal p < 0.05 (VST):\n")
print(R[p < 0.05][order(TF, contrast), .(TF, contrast, NES = round(NES, 2), p = signif(p, 2), FDR_61 = signif(FDR_61, 2),
                                          tested_size, reliable_size)], row.names = FALSE)
ts("done")
