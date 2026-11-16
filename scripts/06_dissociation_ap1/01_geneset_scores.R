## 01_geneset_scores.R
## Per-cell gene-set scores in both cohorts (Seurat "gene set expression across cells"):
##   (i)  cell cycle -- S and G2/M phase genes (Tirosh et al. 2015; Seurat cc.genes.updated.2019),
##        CellCycleScoring -> S.Score, G2M.Score, Phase;
##   (ii) dissociation stress -- the 32-gene core signature induced by tissue dissociation
##        (O'Flanagan et al. 2019, Genome Biol 20:210),
##        AddModuleScore(ctrl = 20, name = "genes_dissoc") -> genes_dissoc1.
## Each score is plotted against the number of detected genes per cell and summarised by cell type,
## so the dissociation score can be compared with AP-1 motif enrichment by cell type.
##   discovery    results/singlecell/atac_motif.rds, post-QC cells (16,451)
##   replication  replication/rna/replication_RNA_annotated.rds (10,700 cells, CELLxGENE cell types)
##
## Outputs  results/motif_sets/geneset_scores/:
##   {main,replication}_geneset_scores_per_cell.csv     per-cell S.Score, G2M.Score, Phase, genes_dissoc1, nFeature_RNA
##   {main,replication}_geneset_scores_by_celltype.csv  mean/median scores and phase proportions per cell type
##   {main,replication}_FeatureScatter_SScore_nFeature.png, _FeatureScatter_dissoc_nFeature.png,
##   _VlnPlot_dissoc_by_celltype.png
## Run from the repository root:  Rscript scripts/06_dissociation_ap1/01_geneset_scores.R
## ---- paths (set the environment variables; see README) ----
TNBC_ROOT <- Sys.getenv("TNBC_ROOT", ".")
suppressPackageStartupMessages({library(Seurat); library(Signac); library(data.table); library(ggplot2)})
setwd(TNBC_ROOT)
OUT <- "results/motif_sets/geneset_scores"; dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

## ---- gene sets ----
## Cell-cycle markers (Tirosh et al. 2015), split into S-phase and G2/M-phase markers.
s.genes   <- cc.genes.updated.2019$s.genes
g2m.genes <- cc.genes.updated.2019$g2m.genes
## Genes upregulated during dissociation of tissue into single cells.
genes.dissoc <- unique(c("ATF3", "BTG2", "CEBPB", "CEBPD", "CXCL3", "CXCL2", "CXCL1", "DNAJA1", "DNAJB1", "DUSP1",
                         "EGR1", "FOS", "FOSB", "HSP90AA1", "HSP90AB1", "HSPA1A", "HSPA1B", "HSPA1A", "HSPA1B",
                         "HSPA8", "HSPB1", "HSPE1", "HSPH1", "ID3", "IER2", "JUN", "JUNB", "JUND", "MT1X", "NFKBIA",
                         "NR4A1", "PPP1R15A", "SOCS3", "ZFP36"))

## ---- the datasets ----
load_main <- function() {
  fb  <- fread("results/singlecell/final_cell_barcodes.csv", header = TRUE)
  obj <- subset(readRDS("results/singlecell/atac_motif.rds"), cells = fb$x)
  ids <- c("LASP-1", "Fb", "LHS", "LASP-2", "B-Myo", "Endo", "T", "Pc", "Myl", "B"); names(ids) <- levels(obj)
  obj <- RenameIdents(obj, ids); obj$celltype <- as.character(Idents(obj)); obj
}
load_replication <- function() {
  obj <- readRDS("replication/rna/replication_RNA_annotated.rds"); obj$celltype <- as.character(obj$celltype); obj
}

score_one <- function(obj, tag) {
  DefaultAssay(obj) <- "RNA"
  obj <- NormalizeData(obj, verbose = FALSE)                   # scores are computed on the log-normalised data layer
  if (!"nFeature_RNA" %in% colnames(obj@meta.data))           # number of genes detected per cell
    obj$nFeature_RNA <- Matrix::colSums(GetAssayData(obj, assay = "RNA", layer = "counts") > 0)
  cat(sprintf("\n[%s] %d cells | S genes present %d/%d | G2M genes present %d/%d | dissociation genes present %d/%d\n",
              tag, ncol(obj), sum(s.genes %in% rownames(obj)), length(s.genes),
              sum(g2m.genes %in% rownames(obj)), length(g2m.genes),
              sum(genes.dissoc %in% rownames(obj)), length(genes.dissoc)))

  ## 1. Cell-cycle scoring (adds S.Score, G2M.Score and Phase). set.ident = FALSE so the cell-type identities are kept.
  obj <- CellCycleScoring(obj, s.features = s.genes, g2m.features = g2m.genes, set.ident = FALSE)

  ## 2. Dissociation score: average expression of the dissociation genes per cell, minus that of 20 control genes per
  ##    gene drawn from the same expression bins (this is how AddModuleScore controls for technical effects).
  set.seed(1)
  obj <- AddModuleScore(obj, features = list(genes.dissoc), ctrl = 20, name = "genes_dissoc")

  ## How do we know the names of these scores? Look at the meta-data columns.
  cat(sprintf("[%s] meta-data columns:\n  %s\n", tag, paste(colnames(obj@meta.data), collapse = ", ")))
  stopifnot(all(c("S.Score", "G2M.Score", "Phase", "genes_dissoc1") %in% colnames(obj@meta.data)))

  ## 3. Correlation between the number of genes and the scores
  for (sc in c("S.Score", "G2M.Score", "genes_dissoc1")) {
    ct <- cor.test(obj@meta.data[[sc]], obj$nFeature_RNA, method = "spearman", exact = FALSE)
    cat(sprintf("[%s] Spearman rho(%s, nFeature_RNA) = %.3f, p = %.2g\n", tag, sc, ct$estimate, ct$p.value))
  }
  png(file.path(OUT, sprintf("%s_FeatureScatter_SScore_nFeature.png", tag)), width = 6, height = 5, units = "in", res = 300)
  print(FeatureScatter(obj, feature1 = "S.Score", feature2 = "nFeature_RNA", group.by = "celltype", pt.size = 0.3))
  dev.off()
  png(file.path(OUT, sprintf("%s_FeatureScatter_dissoc_nFeature.png", tag)), width = 6, height = 5, units = "in", res = 300)
  print(FeatureScatter(obj, feature1 = "genes_dissoc1", feature2 = "nFeature_RNA", group.by = "celltype", pt.size = 0.3))
  dev.off()

  ## 4. Dissociation score by cell type (does the stress signature follow the cell types with strong AP-1 enrichment?)
  png(file.path(OUT, sprintf("%s_VlnPlot_dissoc_by_celltype.png", tag)), width = 8, height = 4.5, units = "in", res = 300)
  print(VlnPlot(obj, features = "genes_dissoc1", group.by = "celltype", pt.size = 0) + NoLegend() +
          ggtitle(sprintf("Dissociation signature score (%s)", tag)))
  dev.off()

  md <- data.table(cell = colnames(obj), celltype = obj$celltype, nFeature_RNA = obj$nFeature_RNA,
                   S.Score = obj$S.Score, G2M.Score = obj$G2M.Score, Phase = obj$Phase, genes_dissoc1 = obj$genes_dissoc1)
  fwrite(md, file.path(OUT, sprintf("%s_geneset_scores_per_cell.csv", tag)))
  byct <- md[, .(n = .N, dissoc_mean = mean(genes_dissoc1), dissoc_median = median(genes_dissoc1),
                 S_mean = mean(S.Score), G2M_mean = mean(G2M.Score),
                 pct_G1 = 100 * mean(Phase == "G1"), pct_S = 100 * mean(Phase == "S"), pct_G2M = 100 * mean(Phase == "G2M")),
             by = celltype][order(-dissoc_mean)]
  fwrite(byct, file.path(OUT, sprintf("%s_geneset_scores_by_celltype.csv", tag)))
  cat(sprintf("[%s] scores by cell type (highest dissociation score first):\n", tag)); print(byct, digits = 3)
  invisible(NULL)
}

score_one(load_main(), "main"); gc(verbose = FALSE)
score_one(load_replication(), "replication")
cat("\nwrote", OUT, "\n")
