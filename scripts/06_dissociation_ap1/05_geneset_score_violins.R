# ==============================================================================
# 05_geneset_score_violins.R
# Violins of dissociation-stress scores, both cohorts.
# Produces: Fig S9
# Run     : Rscript scripts/06_dissociation_ap1/05_geneset_score_violins.R   (from the repository root)
# ==============================================================================
## 05_geneset_score_violins.R — redraw the dissociation-score violins for the two sets.
## Design: (1) the two figures share the same y axis (same limits, breaks and title) so they can be read side
## by side; (2) colours follow the Fig. 2 RNA UMAP cell-type palette (04_annotate_cell_types.R, `colors`, applied
## in cluster order LS-1, Fb, LHr, LS-2, Bs, Endo, T, Pc, Myl, B).
## Replication cell types take the colour of the matching main-set lineage (LASP <- LASP-1, BM <- B-Myo,
## Fibroblasts <- Fb, Endo-1/Endo-2 <- Endo, T-cells <- T, Macrophages <- Myl); adipocytes, absent from the main set,
## get two neutral greys. Cell types are ordered epithelial -> stromal -> immune in both panels.
## Data: the per-cell scores saved by 01_geneset_scores.R (no objects reloaded).
## PNG only.
## ---- paths (set the environment variables; see README) ----
TNBC_ROOT <- Sys.getenv("TNBC_ROOT", ".")
suppressPackageStartupMessages({library(data.table); library(ggplot2)})
D <- file.path(TNBC_ROOT, "results/motif_sets/geneset_scores")
A <- file.path(D, "_archive_violin_v1"); dir.create(A, showWarnings = FALSE)

fig2 <- c("LASP-1" = "#E64B35CC", "Fb" = "#4DBBD5CC", "LHS" = "#00A087CC", "LASP-2" = "#B09C85CC",
          "B-Myo" = "#F39B7FCC", "Endo" = "#8491B4CC", "T" = "#91D1C2CC", "Pc" = "#f9e79f",
          "Myl" = "#7E6148CC", "B" = "#3C5488CC")
COL <- list(
  main = fig2[c("LASP-1", "LASP-2", "LHS", "B-Myo", "Fb", "Endo", "Pc", "T", "B", "Myl")],
  replication = c("LASP" = fig2[["LASP-1"]], "LHS" = fig2[["LHS"]], "BM" = fig2[["B-Myo"]],
                  "Fibroblasts" = fig2[["Fb"]], "Endo-1" = fig2[["Endo"]], "Endo-2" = fig2[["Endo"]],
                  "Adi-1" = "#BDBDBD", "Adi-2" = "#8C8C8C", "T-cells" = fig2[["T"]], "Macrophages" = fig2[["Myl"]]))
## x-axis labels: replication clusters shown with the main-set abbreviations, so both panels have labels of the same
## length and their plotting areas (and y axes) line up exactly.
XLAB <- c("BM" = "B-Myo", "Fibroblasts" = "Fb", "T-cells" = "T", "Macrophages" = "Myl")
TITLE <- c(main = "Main analysis set (GSE162798)", replication = "Replication set (GSE244585)")

md <- rbindlist(lapply(names(COL), function(t)
  fread(file.path(D, sprintf("%s_geneset_scores_per_cell.csv", t)))[, set := t]))
stopifnot(all(md[set == "main", unique(celltype)] %in% names(COL$main)),
          all(md[set == "replication", unique(celltype)] %in% names(COL$replication)))
YL <- range(md$genes_dissoc1); YL <- c(floor(YL[1] * 2) / 2, ceiling(YL[2] * 2) / 2)   # shared limits, 0.5 grid
BR <- seq(YL[1], YL[2], by = 0.5)

for (t in names(COL)) {
  d <- md[set == t][, celltype := factor(celltype, levels = names(COL[[t]]))]
  f <- file.path(D, sprintf("%s_VlnPlot_dissoc_by_celltype.png", t))
  if (file.exists(f) && !file.exists(file.path(A, basename(f)))) file.copy(f, file.path(A, basename(f)), copy.date = TRUE)
  p <- ggplot(d, aes(celltype, genes_dissoc1, fill = celltype)) +
    geom_violin(scale = "width", colour = "grey20", linewidth = 0.3, trim = TRUE) +
    geom_boxplot(width = 0.12, outlier.shape = NA, fill = "white", colour = "grey20", linewidth = 0.3) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50", linewidth = 0.3) +
    scale_fill_manual(values = COL[[t]], guide = "none") +
    scale_x_discrete(labels = function(x) ifelse(x %in% names(XLAB), XLAB[x], x)) +
    scale_y_continuous(limits = YL, breaks = BR, labels = function(x) sprintf("%.1f", x)) +
    labs(x = NULL, y = "Dissociation signature score", title = TITLE[[t]]) +
    theme_classic(base_size = 12) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, colour = "black"), axis.text.y = element_text(colour = "black"),
          plot.title = element_text(hjust = 0.5, face = "bold", size = 12))
  ggsave(f, p, width = 6.5, height = 4.2, dpi = 300, bg = "white")
  cat("wrote", f, "| y limits", YL, "\n")
}
