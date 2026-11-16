# ==============================================================================
# 07_module_score_umap_plot.R
# Draw the TNBC-specific module-score UMAP.
# Produces: Fig 3D
# Inputs  : results/motif_sets/umap_embedding_celltype.csv
#           results/motif_sets/AddModuleScore_per_cell_all_sets.csv
# Run     : Rscript scripts/04_celltype_enrichment/07_module_score_umap_plot.R   (from the repository root)
# ==============================================================================
## 07_module_score_umap_plot.R — draw the TNBC-specific integrated UMAP recoloured to the Fig 3A UpSet TNBC colour
## (#E64B35), from the exported embedding (step 1) + saved per-cell module scores. No Seurat/Signac needed.
## Same hull, cutoffs and floor-at-0 as 05_motif_sets_upset_heatmap.R umap_plot(). Needs fmesher.
## ---- paths (set the environment variables; see README) ----
TNBC_ROOT <- Sys.getenv("TNBC_ROOT", ".")
suppressMessages({library(data.table); library(ggplot2); library(fmesher)})
setwd(TNBC_ROOT)
OUTF <- "results/motif_sets/figures"; HI_COL <- "#E64B35"
epi <- c("LASP-1", "LASP-2", "LHS", "B-Myo")

umap_df <- fread("results/motif_sets/umap_embedding_celltype.csv")
sc_dt   <- fread("results/motif_sets/AddModuleScore_per_cell_all_sets.csv")
sc <- setNames(sc_dt$TNBC_specific, sc_dt$cell)
stopifnot(all(umap_df$cell %in% names(sc)))
lab_df <- umap_df[, .(UMAP_1 = median(UMAP_1), UMAP_2 = median(UMAP_2)), by = celltype]

## hull around epithelial clusters (same call as 05_motif_sets_upset_heatmap.R l.300)
emb  <- as.matrix(umap_df[celltype %in% epi, .(UMAP_1, UMAP_2)])
hull <- fmesher::fm_nonconvex_hull_inla(emb, convex = -0.0021)
hull_seg <- data.frame(x = hull$loc[hull$idx[, 1], 1], y = hull$loc[hull$idx[, 1], 2],
                       xend = hull$loc[hull$idx[, 2], 1], yend = hull$loc[hull$idx[, 2], 2])
blk <- theme(text = element_text(colour = "black"), axis.text = element_text(colour = "black"),
             axis.title = element_text(colour = "black"), legend.text = element_text(colour = "black"))

## umap_plot: identical to 05_motif_sets_upset_heatmap.R except cols = c("grey90", HI_COL)
umap_plot <- function(sc, cells = NULL, title = "", hull = TRUE, legend = TRUE,
                      cols = c("grey90", HI_COL), floor0 = TRUE, min.cutoff = 0.05, max.cutoff = 0.95, pt.size = 0.45) {
  d <- copy(umap_df); d[, score := sc[cell]]
  if (!is.null(cells)) d <- d[cell %in% cells]
  q  <- quantile(d$score, c(min.cutoff, max.cutoff))
  lo <- if (floor0) 0 else q[[1]]
  d[, score := pmin(pmax(score, lo), q[[2]])]
  setorder(d, score)
  p <- ggplot(d, aes(UMAP_1, UMAP_2)) + geom_point(aes(colour = score), size = pt.size, stroke = 0) +
    scale_colour_gradientn(colours = cols, limits = c(lo, q[[2]]), name = "Module\nscore") +
    geom_text(data = lab_df[celltype %in% unique(d$celltype)], aes(label = celltype), size = 3.5, fontface = "bold") +
    labs(x = "UMAP 1", y = "UMAP 2", title = title) + theme_classic(base_size = 11) + blk +
    theme(plot.title = element_text(hjust = 0.5), legend.position = "right")
  if (hull) p <- p + geom_segment(data = hull_seg, aes(x = x, y = y, xend = xend, yend = yend),
                                  colour = "black", linetype = "dotted", linewidth = 0.5)
  if (!legend) p <- p + theme(legend.position = "none")
  p
}
f <- file.path(OUTF, "UMAP_TNBC_specific_AddModuleScore_all_TNBCcol.jpg")
ggsave(f, umap_plot(sc, title = ""), width = 6, height = 5, dpi = 300, bg = "white")
cat("wrote", f, "| high colour", HI_COL, "\n")
