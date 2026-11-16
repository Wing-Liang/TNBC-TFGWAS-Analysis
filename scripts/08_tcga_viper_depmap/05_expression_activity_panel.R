# ==============================================================================
# 05_expression_activity_panel.R
# TF gene expression versus inferred TF activity.
# Produces: Fig 4A-B
# Inputs  : results/tcga/TCGA_BRCA_DESeq2_BasalRef_4contrasts_allgenes.rds
# Run     : Rscript scripts/08_tcga_viper_depmap/05_expression_activity_panel.R   (from the repository root)
# ==============================================================================
## 05_expression_activity_panel.R — Fig 4A-B: TF GENE RNA change (TCGA DESeq2 log2FC of the TF-encoding gene) and TF
## activity (msVIPER NES, inferred from the TF's regulon targets). Labels  # "TF gene RNA", not "mRNA".
## for the 61 TNBC-specific TF genes, side by side, plus their concordance.
## Inputs (no new estimation):
##   results/tf_activity/TCGA_DESeq2_61TF.csv                                      (01; Basal as reference)
##   results/tf_activity/viper_vst/f10n39_min5/VIPER_VST_msviper_61TF_4contrasts_ARACNe_BRCA.csv
##                                     (03; DESeq2 VST, ARACNe-BRCA, minsize 5 --: min5 only;
##                                      TFs with < 25 tested targets in >= 1 contrast are marked † on the TF label)
## Both are Basal vs comparator (LumA, LumB, Her2, NAT) on the same 836 TCGA-BRCA samples (female, white).
## Significance: DESeq2 FDR (BH over all genes) < 0.05; VIPER FDR_61 (BH within the tested 61-gene set) < 0.05.
## "Supported by both" = up in Basal by DESeq2 in all four contrasts (the manuscript's 13 genes) AND higher
## activity in Basal (NES > 0, FDR_61 < 0.05) in at least one contrast.
## dot-plot column labels -> "RNA (log2FC)" and "Inferred / TF activity (NES)".
## Run with R 4.5.1 (see environment/README.md).
## ---- paths (set the environment variables; see README) ----
TNBC_ROOT <- Sys.getenv("TNBC_ROOT", ".")
suppressMessages({library(data.table); library(ggplot2); library(patchwork)})
setwd(TNBC_ROOT)
OUT <- "results/tf_activity/"; FIG <- "results/tf_activity/figures/combined"
dir.create(FIG, showWarnings = FALSE, recursive = TRUE)
CN <- c("LumA", "LumB", "Her2", "NAT"); LABEL <- setNames(paste("Basal vs", CN), CN)

## ---- inputs, long format ----
## DESeq2 FDR is re-adjusted WITHIN the tested TF genes (BH over the raw
## DESeq2 p-values of the 61 TF genes, per contrast), so it uses the same denominator idea as VIPER's
## FDR_61 (BH within the tested TFs). The genome-wide padj (BH over ~28,000 genes) is not used for
## the ring / classification. Raw p-values come from the genome-wide DESeq2 cache.
de <- fread(file.path(OUT, "TCGA_DESeq2_61TF.csv"))
RES <- readRDS("results/tcga/TCGA_BRCA_DESeq2_BasalRef_4contrasts_allgenes.rds")
RES[, contrast := sub("^Basal vs ", "", Contrast)]
PV <- RES[Ensembl %in% de$Ensembl & contrast %in% CN, .(Ensembl, contrast, pvalue)]
PV[, FDR_DE_within := p.adjust(pvalue, "BH"), by = contrast]      # BH within the 61 TF genes, per contrast
PV <- merge(PV, de[, .(Ensembl, TF = toupper(Gene))], by = "Ensembl")
DE <- rbindlist(lapply(CN, function(cn) { cc <- LABEL[[cn]]
  d <- data.table(TF = toupper(de$Gene), Ensembl = de$Ensembl, contrast = cn,
                  value = de[[paste0("log2FoldChange__", cc)]])
  merge(d, PV[contrast == cn, .(Ensembl, q = FDR_DE_within)], by = "Ensembl") }))
DE[, up_all4 := TF %in% DE[, .(ok = all(q < 0.05 & value > 0)), by = TF][ok == TRUE, TF]]  # up in Basal, all 4, within-TF FDR
DE[, Ensembl := NULL]
vp <- fread(file.path(OUT, "viper_vst/f10n39_min5/VIPER_VST_msviper_61TF_4contrasts_ARACNe_BRCA.csv"))
VP <- vp[, .(TF, contrast, value = NES, q = FDR_61, p)]
SMALL <- vp[, .(small = any(!reliable_size)), by = TF][small == TRUE, TF]   # < 25 tested targets in >= 1 contrast
dag <- function(x) ifelse(x %in% SMALL, paste0(x, " †"), x)
stopifnot(uniqueN(DE$TF) == 61)

## ---- concordance table (one row per TF x contrast) ----
C <- merge(DE[, .(TF, contrast, log2FC = value, FDR_DE = q, up_all4)],
           VP[, .(TF, contrast, NES = value, FDR_VIPER = q, p_VIPER = p)], by = c("TF", "contrast"), all.x = TRUE)
cls <- function(v, q) fifelse(is.na(v), "not tested", fifelse(!is.na(q) & q < 0.05, fifelse(v > 0, "up", "down"), "ns"))
C[, `:=`(DE_call = cls(log2FC, FDR_DE), VIPER_call = cls(NES, FDR_VIPER))]
C[, agreement := fifelse(VIPER_call == "not tested", "VIPER not testable",
                 fifelse(DE_call %in% c("up", "down") & DE_call == VIPER_call, "both, same direction",
                 fifelse(DE_call %in% c("up", "down") & VIPER_call %in% c("up", "down"), "both, opposite direction",
                 fifelse(DE_call %in% c("up", "down"), "RNA only",
                 fifelse(VIPER_call %in% c("up", "down"), "activity only", "neither")))))]
both <- C[up_all4 == TRUE & VIPER_call == "up", unique(TF)]
fwrite(C[order(TF, contrast)], file.path(OUT, "TCGA_DE_vs_VIPER_61TF_concordance.csv"))
cat("TFs supported by both (DE up in all 4 AND VIPER activity up, FDR_61 < 0.05, in >= 1 contrast):\n  ",
    paste(sort(both), collapse = ", "), "\n\n")
cat("agreement per contrast:\n"); print(dcast(C, agreement ~ contrast, fun.aggregate = length))
rho <- C[!is.na(NES), .(n = .N, rho = cor(log2FC, NES, method = "spearman"),
                         p = cor.test(log2FC, NES, method = "spearman", exact = FALSE)$p.value), by = contrast]
cat("\nSpearman correlation, log2FC vs NES:\n"); print(rho[, .(contrast, n, rho = round(rho, 2), p = signif(p, 2))])

## ---- Panel A: dot plot, two readouts per contrast ----
L <- rbind(DE[, .(TF, contrast, readout = "RNA\n(log2FC)", value, q)],
           VP[, .(TF, contrast, readout = "Inferred\nTF activity\n(NES)", value, q)])
## circle SIZE = effect size (|log2FC| for TF gene RNA, |NES| for activity); star removed.
## log2FC and NES are on different scales, so colour and size are scaled WITHIN each readout; a few extreme
## log2FCs (HOXB13 vs NAT, 5.8) would otherwise wash out everything else.
## Colour: signed value / 95th percentile of |value| in that readout, clamped to [-1, 1].
## Size:   |value| / 95th percentile of |value| in that readout, capped at 1.
## Ring:   FDR < 0.05. FDR NA (DESeq2 independent filtering of very low counts: CDX4, NKX2-4 in some
##         contrasts) is drawn as not significant.
L[, q_plot := fifelse(is.na(q), 1, q)]
L[, scaled := pmax(-1, pmin(1, value / quantile(abs(value), 0.95, na.rm = TRUE))), by = readout]
L[, eff_size := pmin(1, abs(value) / quantile(abs(value), 0.95, na.rm = TRUE)), by = readout]
L[, contrast := factor(LABEL[contrast], levels = LABEL)]
L[, readout := factor(readout, levels = c("RNA\n(log2FC)", "Inferred\nTF activity\n(NES)"))]
ord <- merge(data.table(TF = unique(DE$TF)), VP[, .(mNES = mean(value)), by = TF], all.x = TRUE)
ord <- merge(ord, DE[, .(mFC = mean(value)), by = TF], by = "TF")
ord <- ord[order(!is.na(mNES), fifelse(is.na(mNES), mFC, mNES))]              # untested at bottom
L[, TF_lab := factor(dag(TF), levels = dag(ord$TF))]
miss <- CJ(TF = setdiff(unique(DE$TF), unique(VP$TF)), contrast = CN)[, `:=`(readout = "Inferred\nTF activity\n(NES)")]
miss[, `:=`(contrast = factor(LABEL[contrast], levels = LABEL), readout = factor(readout, levels = levels(L$readout)),
            TF_lab = factor(dag(TF), levels = levels(L$TF_lab)))]
pA <- ggplot(L[!is.na(value)], aes(readout, TF_lab)) +
  geom_point(aes(fill = scaled, size = eff_size, colour = q_plot < 0.05), shape = 21, stroke = 0.6) +
  geom_point(data = miss, shape = 4, size = 1.2, colour = "grey60") +
  scale_fill_gradient2(low = "#3B4CC0", mid = "white", high = "#B40426", limits = c(-1, 1),
                       name = "Higher in Basal (red)\nLower in Basal (blue)") +
  ## legend shows only the black-border key
  scale_colour_manual(values = c(`TRUE` = "black", `FALSE` = "grey75"), breaks = "TRUE", labels = "FDR < 0.05", name = NULL) +
  scale_size_continuous(range = c(0.5, 4.5), limits = c(0, 1), breaks = c(0.25, 0.5, 1),
                        name = "Effect size |log2FC| or |NES|") +
  ## size keys = solid grey circles; all legend titles to the left of their keys, one row
  guides(fill  = guide_colourbar(order = 1, title.position = "left", title.vjust = 0.5, barwidth = 8, barheight = 0.6,
                                theme = theme(legend.title = element_text(size = 8, margin = margin(r = 16)))),
         size  = guide_legend(order = 2, title.position = "left", title.vjust = 0.5, nrow = 1,
                              override.aes = list(fill = "grey60", colour = "grey60")),
         colour = guide_legend(order = 3, override.aes = list(size = 3, fill = "white"))) +
  facet_wrap(~ contrast, nrow = 1) +
  labs(x = NULL, y = NULL, title = "TF gene RNA (DESeq2) and TF activity (msVIPER), TCGA-BRCA",
       caption = paste0("61 TF genes of the 57 TNBC-specific motifs; Basal vs comparator; female, white (Basal 116, LumA 446, LumB 131, Her2 39, NAT 104).\n",
                        "TF gene RNA: DESeq2 log2FC of the TF-encoding gene. TF activity: msVIPER NES from the TF's regulon targets (ARACNe-BRCA, DESeq2 VST, minsize 5).\n",
                        "FDR (both measures) is Benjamini-Hochberg within the tested TF genes, per contrast; † = < 25 regulon targets tested in >= 1 contrast (low confidence).\n",
                        "Colour (signed effect) and size (|effect|) are scaled separately per measure (each divided by the 95th percentile of its own |effect|, capped at 1).\n",
                        "Black border = FDR < 0.05; x = not testable by VIPER.")) +
  theme_bw(base_size = 9) +
  theme(axis.text = element_text(colour = "black"), panel.grid.minor = element_blank(), strip.text = element_text(face = "bold"),
        plot.caption = element_text(hjust = 0, size = 7), plot.title.position = "plot", plot.caption.position = "plot",
        legend.position = "bottom", legend.box = "horizontal", legend.box.just = "center",
        legend.title = element_text(size = 8, margin = margin(r = 8)),
        legend.text = element_text(size = 7.5), legend.key.height = unit(0.4, "lines"))
fA <- file.path(FIG, "Fig5_combined_DE_VIPER_dotplot.png")
ggsave(fA, pA, width = 8.5, height = 2.2 + 0.19 * 61, dpi = 300, bg = "white"); cat("wrote", fA, "\n")

## ---- Panel B: log2FC vs NES per contrast ----
S <- C[!is.na(NES)]
S[, contrast := factor(LABEL[contrast], levels = LABEL)]
S[, agreement := factor(agreement, levels = c("both, same direction", "RNA only", "activity only", "both, opposite direction", "neither"))]
KEY <- unique(c(both, "FOXA1", "RXRA", "PPARA", "ERG", "EBF1", "TCFL5", "MEF2B"))
## every red point (both, same direction) is labelled in its own panel, in addition to KEY
rl <- rho[, .(contrast = factor(LABEL[contrast], levels = LABEL), lab = sprintf("rho = %.2f, p = %.1g", rho, p))]
pB <- ggplot(S, aes(log2FC, NES)) +
  geom_hline(yintercept = 0, colour = "grey70") + geom_vline(xintercept = 0, colour = "grey70") +
  geom_point(aes(colour = agreement), size = 1.8) +
  ggrepel::geom_text_repel(data = S[TF %in% KEY | agreement == "both, same direction"], aes(label = TF), size = 2.5, max.overlaps = Inf,
                           min.segment.length = 0, segment.size = 0.2, seed = 42) +
  geom_text(data = rl, aes(x = -Inf, y = Inf, label = lab), hjust = -0.05, vjust = 1.4, size = 2.8) +
  ## red = gene expression and TF activity, yellow = gene expression only; legend labels match Fig. 4B legend
  scale_colour_manual(values = c(`both, same direction` = "#B40426", `RNA only` = "#F2C200", `activity only` = "#3B4CC0",
                                 `both, opposite direction` = "#6A3D9A", neither = "grey70"), drop = TRUE,
                      labels = c(`both, same direction` = "Gene expression and TF activity", `RNA only` = "Gene expression only",
                                 `activity only` = "TF activity only", `both, opposite direction` = "Opposite directions",
                                 neither = "Neither"),
                      name = "FDR < 0.05") +   # BH detail moved to the figure legend text
  facet_wrap(~ contrast, ncol = 1, scales = "free_x") +   # one column x 4 rows
  labs(x = "TF gene RNA log2 fold change (DESeq2)", y = "TF activity NES (msVIPER)",
       title = "TF gene RNA change vs TF activity
(Spearman correlation)") +
  guides(colour = guide_legend(nrow = 1, title.position = "top")) +   # legend in one row
  theme_bw(base_size = 9) + theme(legend.position = "bottom", legend.text = element_text(size = 7.5), legend.key.width = unit(0.3, "lines"),
                                  panel.grid.minor = element_blank(), strip.text = element_text(face = "bold"))
fB <- file.path(FIG, "Fig5_combined_DE_VIPER_scatter.png")
ggsave(fB, pB, width = 4.8, height = 13, dpi = 300, bg = "white"); cat("wrote", fB, "\n")
