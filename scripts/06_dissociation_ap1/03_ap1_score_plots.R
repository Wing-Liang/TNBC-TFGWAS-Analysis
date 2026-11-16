# ==============================================================================
# 03_ap1_score_plots.R
# Violin and scatter plots of the AP-1 score.
# Produces: Fig S9
# Run     : Rscript scripts/06_dissociation_ap1/03_ap1_score_plots.R   (from the repository root)
# ==============================================================================
## 03_ap1_score_plots.R — step 2 of the AP-1 gene-set analysis: same figures as for the full
## dissociation signature, for AP1_score1 (FOS, FOSB, JUN, JUNB, JUND, ATF3):
##   {main,replication}_VlnPlot_AP1_by_celltype.png                  violins, shared y axis across the two sets, Fig. 2 colours
##   {main,replication}_FeatureScatter_AP1_nFeature_byCelltype.png   per-cell-type scatter vs nFeature_RNA, 4 per row,
##                                                                   x and y axes on every panel, Spearman rho in each title
##   spearman_AP1_vs_nFeature_by_celltype.csv
## ---- paths (set the environment variables; see README) ----
TNBC_ROOT <- Sys.getenv("TNBC_ROOT", ".")
suppressPackageStartupMessages({library(data.table); library(ggplot2)})
D <- file.path(TNBC_ROOT, "results/motif_sets/geneset_scores")
fig2 <- c("LASP-1" = "#E64B35CC", "Fb" = "#4DBBD5CC", "LHS" = "#00A087CC", "LASP-2" = "#B09C85CC",
          "B-Myo" = "#F39B7FCC", "Endo" = "#8491B4CC", "T" = "#91D1C2CC", "Pc" = "#f9e79f",
          "Myl" = "#7E6148CC", "B" = "#3C5488CC")
COL <- list(
  main = fig2[c("LASP-1", "LASP-2", "LHS", "B-Myo", "Fb", "Endo", "Pc", "T", "B", "Myl")],
  replication = c("LASP" = fig2[["LASP-1"]], "LHS" = fig2[["LHS"]], "BM" = fig2[["B-Myo"]],
                  "Fibroblasts" = fig2[["Fb"]], "Endo-1" = fig2[["Endo"]], "Endo-2" = fig2[["Endo"]],
                  "Adi-1" = "#BDBDBD", "Adi-2" = "#8C8C8C", "T-cells" = fig2[["T"]], "Macrophages" = fig2[["Myl"]]))
XLAB <- c("BM" = "B-Myo", "Fibroblasts" = "Fb", "T-cells" = "T", "Macrophages" = "Myl")
short <- function(x) ifelse(x %in% names(XLAB), XLAB[x], x)
TITLE <- c(main = "Main analysis set (GSE162798)", replication = "Replication set (GSE244585)")
YLAB <- "AP-1 gene-set score\n(FOS, FOSB, JUN, JUNB, JUND, ATF3)"

md <- rbindlist(lapply(names(COL), function(t) fread(file.path(D, sprintf("%s_AP1_score_per_cell.csv", t)))[, set := t]))
YL <- range(md$AP1_score1); YL <- c(floor(YL[1] * 2) / 2, ceiling(YL[2] * 2) / 2); BR <- seq(YL[1], YL[2], 0.5)
out <- list()
for (t in names(COL)) {
  d <- md[set == t]; stopifnot(all(d[, unique(celltype)] %in% names(COL[[t]])))
  ## violins (shared y axis)
  dv <- copy(d)[, celltype := factor(celltype, levels = names(COL[[t]]))]
  p <- ggplot(dv, aes(celltype, AP1_score1, fill = celltype)) +
    geom_violin(scale = "width", colour = "grey20", linewidth = 0.3) +
    geom_boxplot(width = 0.12, outlier.shape = NA, fill = "white", colour = "grey20", linewidth = 0.3) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50", linewidth = 0.3) +
    scale_fill_manual(values = COL[[t]], guide = "none") +
    scale_x_discrete(labels = short) +
    scale_y_continuous(limits = YL, breaks = BR, labels = function(x) sprintf("%.1f", x)) +
    labs(x = NULL, y = YLAB, title = TITLE[[t]]) + theme_classic(base_size = 12) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, colour = "black"), axis.text.y = element_text(colour = "black"),
          plot.title = element_text(hjust = 0.5, face = "bold", size = 12))
  ggsave(file.path(D, sprintf("%s_VlnPlot_AP1_by_celltype.png", t)), p, width = 6.5, height = 4.2, dpi = 300, bg = "white")
  ## per-cell-type scatter vs number of detected genes
  st <- d[, {ct <- cor.test(AP1_score1, nFeature_RNA, method = "spearman", exact = FALSE)
             .(n = .N, rho = unname(ct$estimate), p = ct$p.value)}, by = celltype]
  st[, lab := sprintf("%s (n = %s)\nρ = %.2f, p = %s", short(celltype), format(n, big.mark = ",", trim = TRUE), rho,
                      formatC(p, format = "g", digits = 2))]
  lv <- st[match(names(COL[[t]]), celltype), lab]
  ds <- merge(d, st[, .(celltype, lab)], by = "celltype")[, lab := factor(lab, levels = lv)]
  q <- ggplot(ds[sample(.N)], aes(AP1_score1, nFeature_RNA, colour = lab)) + geom_point(size = 0.35, alpha = 0.8) +
    facet_wrap(~ lab, ncol = 4, axes = "all", axis.labels = "all") +
    scale_colour_manual(values = setNames(COL[[t]], lv), guide = "none") +
    labs(x = "AP-1 gene-set score", y = "Number of detected genes (nFeature_RNA)", title = TITLE[[t]]) +
    theme_classic(base_size = 11) +
    theme(strip.background = element_blank(), strip.text = element_text(size = 9, face = "bold"),
          plot.title = element_text(hjust = 0.5, face = "bold"), axis.text = element_text(colour = "black"))
  ggsave(file.path(D, sprintf("%s_FeatureScatter_AP1_nFeature_byCelltype.png", t)), q, width = 10, height = 8, dpi = 300, bg = "white")
  out[[t]] <- st[, .(set = t, celltype, n, rho = round(rho, 3), p = signif(p, 2))]
  cat("wrote", t, "violin + scatter\n")
}
res <- rbindlist(out); fwrite(res, file.path(D, "spearman_AP1_vs_nFeature_by_celltype.csv")); print(res)
