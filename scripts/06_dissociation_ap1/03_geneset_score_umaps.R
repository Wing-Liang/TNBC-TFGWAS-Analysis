# ==============================================================================
# 03_geneset_score_umaps.R
# UMAPs of dissociation-stress and AP-1 scores, both cohorts.
# Produces: Fig S9
# Inputs  : results/motif_sets/umap_embedding_celltype.csv
#           replication/data/cellxgene/embeddings_replication.csv
# Run     : Rscript scripts/06_dissociation_ap1/03_geneset_score_umaps.R   (from the repository root)
# ==============================================================================
## 03_geneset_score_umaps.R — UMAPs of the AddModuleScore gene-set scores:
##   dissociation signature (genes_dissoc1; 32 genes, O'Flanagan et al. 2019) and its AP-1 subset
##   (AP1_score1; FOS, FOSB, JUN, JUNB, JUND, ATF3), for the main analysis set and the replication set.
## Embeddings (RNA UMAP, as in Fig. 2 / Fig. 3D):
##   main        results/motif_sets/umap_embedding_celltype.csv (umap.rna of the 16,451 final cells)
##   replication replication/data/cellxgene/embeddings_replication.csv (CELLxGENE umap.rna; 10,700 final-QC cells)
## Colour scale as the Fig. 3D module-score UMAPs (05_motif_sets_upset_heatmap.R umap_plot): grey90 -> high colour, floored
## at 0 (cells with a score <= 0 are grey), capped at the 95th percentile of that score in that set; cells drawn in
## increasing score order so high-scoring cells sit on top. Scores are per-set (AddModuleScore is relative to control
## genes within each dataset), so each panel has its own colour scale.
## Outputs (PNG only), results/motif_sets/geneset_scores/:
##   {main,replication}_UMAP_dissoc_AP1.png   two panels per set (dissociation | AP-1)
##   UMAP_dissoc_AP1_main_vs_replication.png  2 x 2 (rows: main, replication; columns: dissociation, AP-1)
## ---- paths (set the environment variables; see README) ----
TNBC_ROOT <- Sys.getenv("TNBC_ROOT", ".")
suppressPackageStartupMessages({library(data.table); library(ggplot2); library(patchwork)})
setwd(TNBC_ROOT); D <- "results/motif_sets/geneset_scores"
HI <- "#B2182B"
emb <- list(
  main = fread("results/motif_sets/umap_embedding_celltype.csv")[, .(cell, UMAP_1, UMAP_2)],
  replication = fread("replication/data/cellxgene/embeddings_replication.csv")[, .(cell, UMAP_1 = umap.rna_1, UMAP_2 = umap.rna_2)])
TITLE <- c(main = "Main analysis set (GSE162798)", replication = "Replication set (GSE244585)")
SC <- c(genes_dissoc1 = "Dissociation signature score", AP1_score1 = "AP-1 gene-set score")

umap_panel <- function(d, sc, title) {
  q95 <- quantile(d[[sc]], 0.95)
  d <- copy(d)[, v := pmin(pmax(get(sc), 0), q95)][order(v)]
  ggplot(d, aes(UMAP_1, UMAP_2, colour = v)) + geom_point(size = 0.3) +
    scale_colour_gradient(low = "grey90", high = HI, limits = c(0, q95), name = "Score",
                          labels = function(x) sprintf("%.2f", x)) +
    labs(title = title, x = "UMAP_1", y = "UMAP_2") + coord_equal() +
    theme_classic(base_size = 11) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 11), axis.text = element_text(colour = "black"),
          legend.key.height = unit(0.8, "lines"))
}

P <- list()
for (t in names(emb)) {
  s <- fread(file.path(D, sprintf("%s_AP1_score_per_cell.csv", t)))      # cell, celltype, nFeature_RNA, AP1_score1, genes_dissoc1
  d <- merge(emb[[t]], s, by = "cell")
  stopifnot(nrow(d) == nrow(s), !anyNA(d$genes_dissoc1))
  cat(sprintf("%-12s %d cells matched to the embedding\n", t, nrow(d)))
  P[[t]] <- lapply(names(SC), function(sc) umap_panel(d, sc, SC[[sc]]))
  g <- wrap_plots(P[[t]], ncol = 2) + plot_annotation(title = TITLE[[t]],
         theme = theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 13)))
  ggsave(file.path(D, sprintf("%s_UMAP_dissoc_AP1.png", t)), g, width = 11, height = 5, dpi = 300, bg = "white")
}
row_lab <- function(t) wrap_elements(grid::textGrob(TITLE[[t]], rot = 90, gp = grid::gpar(fontface = "bold", fontsize = 11)))
g2 <- (row_lab("main") + P$main[[1]] + P$main[[2]] + row_lab("replication") + P$replication[[1]] + P$replication[[2]]) +
  plot_layout(ncol = 3, widths = c(0.08, 1, 1))
ggsave(file.path(D, "UMAP_dissoc_AP1_main_vs_replication.png"), g2, width = 11.5, height = 9.5, dpi = 300, bg = "white")
cat("wrote UMAP figures to", D, "\n")
