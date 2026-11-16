# ==============================================================================
# 07_depmap_tf_genes.R
# CRISPR gene effects of the 61 TF genes in TNBC cell lines (DepMap 25Q3).
# Produces: Fig 4C, Table S11
# Inputs  : results/tfgwas/tfgwas_significant_TNBC_unique.csv
# Run     : Rscript scripts/08_tcga_viper_depmap/07_depmap_tf_genes.R   (from the repository root)
# ==============================================================================
## 07_depmap_tf_genes.R — Fig 4C: CRISPR gene dependency (DepMap Public 25Q3, Chronos) of the 61 TF
## genes in TNBC cell lines.
##
## Reproduces Table S11 from the DepMap Public 25Q3 release:
##   gene effect : $DEPMAP_DIR/CRISPRGeneEffect_25Q3.csv (column names "SYMBOL (Entrez)" -> SYMBOL)
##   cell lines  : Model.csv, OncotreeLineage == "Breast" and ModelSubtypeFeatures in
##                 {TNBC, luminal TNBC, basal_A TNBC, basal_B TNBC, basal TNBC} -> 34 models, of which
##                 25 have a 25Q3 CRISPR screen (the manuscript's 25)
##   genes       : the 57 TNBC-specific motif names split on "::" -> 61 symbols (all 61 are in 25Q3)
##   dependent   : gene effect < -0.5 (per line); genes labelled "gene (n dependent lines)"
## Figures: all 61 genes, and the dependent-only subset (>= 1 dependent line; the 13 named in the
## manuscript legend). Run with R 4.5.1.
## ---- paths (set the environment variables; see README) ----
TNBC_ROOT <- Sys.getenv("TNBC_ROOT", ".")
DEPMAP_DIR <- Sys.getenv("DEPMAP_DIR", "data/DepMap")
suppressMessages({library(data.table); library(ggplot2)})
setwd(TNBC_ROOT)
OUT <- "results/tf_activity/"; FIG <- "results/tf_activity/figures/depmap"
dir.create(OUT, showWarnings = FALSE, recursive = TRUE); dir.create(FIG, showWarnings = FALSE, recursive = TRUE)
DEP <- -0.5

motifs <- unique(fread("results/tfgwas/tfgwas_significant_TNBC_unique.csv")$motif.name)
TF61   <- unique(unlist(strsplit(motifs, "::")))
stopifnot(length(TF61) == 61)

ge <- fread(file.path(DEPMAP_DIR, "CRISPRGeneEffect_25Q3.csv"))
setnames(ge, sub(" [(].*[)]$", "", names(ge)))
setnames(ge, 1, "ModelID")
model <- fread(file.path(DEPMAP_DIR, "Model.csv"))
tn <- model[OncotreeLineage == "Breast" &
            ModelSubtypeFeatures %in% c("TNBC", "luminal TNBC", "basal_B TNBC", "basal_A TNBC", "basal TNBC")]
tn <- tn[ModelID %in% ge$ModelID]
cat(sprintf("TNBC cell lines with a 25Q3 CRISPR screen: %d\n", nrow(tn)))
stopifnot(all(TF61 %in% names(ge)))

L <- melt(ge[ModelID %in% tn$ModelID, c("ModelID", TF61), with = FALSE], id.vars = "ModelID",
          variable.name = "gene", value.name = "gene_effect", variable.factor = FALSE)
L <- merge(L, tn[, .(ModelID, CellLineName, ModelSubtypeFeatures, PrimaryOrMetastasis)], by = "ModelID")
S <- L[, .(n_dependent = sum(gene_effect < DEP, na.rm = TRUE), n_lines = sum(!is.na(gene_effect)),
           median_effect = median(gene_effect, na.rm = TRUE), min_effect = min(gene_effect, na.rm = TRUE),
           dependent_lines = paste(sort(CellLineName[gene_effect < DEP & !is.na(gene_effect)]), collapse = ";")),
       by = gene][order(-n_dependent, median_effect)]
fwrite(S, file.path(OUT, "DepMap_25Q3_61TF_summary.csv"))
fwrite(dcast(L, ModelID + CellLineName + ModelSubtypeFeatures + PrimaryOrMetastasis ~ gene, value.var = "gene_effect"),
       file.path(OUT, "DepMap_25Q3_61TF_gene_effect_wide.csv"))

cat("\ndependent genes (effect <", DEP, "in >= 1 line):\n"); print(S[n_dependent > 0, .(gene, n_dependent, median_effect = round(median_effect, 2))])
REC <- c(MZF1 = 7, JUNB = 6, TFAP4 = 5, FOSL1 = 4, ZBTB7B = 2, FOSL2 = 2, NFE2L2 = 1, RUNX3 = 1,
         MYF6 = 1, FEV = 1, PPARA = 1, MAF = 1, JUND = 1)          # counts reported in the manuscript
now <- setNames(S$n_dependent, S$gene)[names(REC)]
cat("reported counts (13 genes) reproduced:", isTRUE(all(now == REC)) && sum(S$n_dependent > 0) == 13, "\n")

## ---- plot (box + jitter by cell line, ordered by median effect) ----
L <- merge(L, S[, .(gene, n_dependent, median_effect)], by = "gene")
## gene on the x-axis (vertical boxes), legend at the bottom in 9 columns.
## Genes ordered left -> right from the most negative median gene effect (strongest dependency).
L[, x_label := paste0(gene, " (", n_dependent, ")")]
L[, x_label := reorder(factor(x_label), median_effect)]
cols <- setNames(grDevices::colorRampPalette(RColorBrewer::brewer.pal(12, "Paired"))(uniqueN(L$CellLineName)),
                 sort(unique(L$CellLineName)))
dep_plot <- function(d, f, w) {
  p <- ggplot(d, aes(x = x_label, y = gene_effect)) +
    geom_boxplot(outlier.shape = NA, fill = "#d3ebeb", color = "navyblue") +
    geom_jitter(aes(color = CellLineName), width = 0.15, height = 0, size = 1.4) +
    geom_hline(yintercept = DEP, linetype = "dotted", color = "red") +
    scale_color_manual(values = cols) +
    labs(x = "Gene (No. of dependent cell lines)", y = "Gene effect (Chronos, DepMap 25Q3)", color = "TNBC cell line") +
    theme_classic() +
    theme(axis.text = element_text(color = "black"), axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "bottom", legend.title = element_text(size = 10), legend.text = element_text(size = 8),
          legend.key.size = unit(0.4, "lines")) +
    guides(color = guide_legend(ncol = 9, title.position = "top", override.aes = list(size = 2)))
  ggsave(f, p, width = w, height = 6.5, dpi = 300, bg = "white")
  cat("wrote", f, "\n")
}
dep_plot(L, file.path(FIG, "DepMap_25Q3_61TF_all.png"), max(12, 1.5 + 0.28 * 61))
dep_plot(L[n_dependent > 0], file.path(FIG, "DepMap_25Q3_61TF_dependent.png"), max(11, 1.5 + 0.5 * uniqueN(L[n_dependent > 0]$gene)))
