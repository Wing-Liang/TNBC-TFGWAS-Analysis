# ==============================================================================
# 02_dissociation_ap1_scores.R
# Score cells for the AP-1 subset of the dissociation-stress signature.
# Produces: Fig S9
# Inputs  : results/singlecell/final_cell_barcodes.csv
#           results/singlecell/atac_motif.rds
#           replication/rna/replication_RNA_annotated.rds
# Run     : Rscript scripts/06_dissociation_ap1/02_dissociation_ap1_scores.R   (from the repository root)
# ==============================================================================
## 02_dissociation_ap1_scores.R — AP-1 subset of the dissociation signature, scored across cells (tutorial §8.8).
## Same method as 01_geneset_scores.R, restricted to the AP-1 genes of genes.dissoc:
##   FOS, FOSB, JUN, JUNB, JUND (Fos- and Jun-related) and ATF3 (ATF/AP-1 bZIP partner).
## AddModuleScore(features = list(genes.AP1), ctrl = 20, name = "AP1_score") -> meta-data column "AP1_score1".
## Step 1 of 2 (scoring, signac_env). Step 2 = 03_ap1_score_plots.R (R 4.5.1; violins + per-cell-type scatters).
## Outputs, results/motif_sets/geneset_scores/:
##   {main,replication}_AP1_score_per_cell.csv     cell, celltype, nFeature_RNA, AP1_score1, genes_dissoc1 (from the full run)
##   {main,replication}_AP1_score_by_celltype.csv
## ---- paths (set the environment variables; see README) ----
TNBC_ROOT <- Sys.getenv("TNBC_ROOT", ".")
suppressPackageStartupMessages({library(Seurat); library(Signac); library(data.table)})
setwd(TNBC_ROOT)
OUT <- "results/motif_sets/geneset_scores"
genes.AP1 <- c("FOS", "FOSB", "JUN", "JUNB", "JUND", "ATF3")

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
  DefaultAssay(obj) <- "RNA"; obj <- NormalizeData(obj, verbose = FALSE)
  if (!"nFeature_RNA" %in% colnames(obj@meta.data))
    obj$nFeature_RNA <- Matrix::colSums(GetAssayData(obj, assay = "RNA", layer = "counts") > 0)
  stopifnot(all(genes.AP1 %in% rownames(obj)))
  set.seed(1)
  obj <- AddModuleScore(obj, features = list(genes.AP1), ctrl = 20, name = "AP1_score")
  stopifnot("AP1_score1" %in% colnames(obj@meta.data))
  md <- data.table(cell = colnames(obj), celltype = obj$celltype, nFeature_RNA = obj$nFeature_RNA, AP1_score1 = obj$AP1_score1)
  full <- fread(file.path(OUT, sprintf("%s_geneset_scores_per_cell.csv", tag)))[, .(cell, genes_dissoc1)]
  md <- merge(md, full, by = "cell", all.x = TRUE)
  fwrite(md, file.path(OUT, sprintf("%s_AP1_score_per_cell.csv", tag)))
  r1 <- cor.test(md$AP1_score1, md$nFeature_RNA, method = "spearman", exact = FALSE)
  r2 <- cor.test(md$AP1_score1, md$genes_dissoc1, method = "spearman", exact = FALSE)
  cat(sprintf("[%s] %d cells | Spearman rho(AP1_score1, nFeature_RNA) = %.3f (p = %.2g) | rho(AP1_score1, genes_dissoc1) = %.3f\n",
              tag, nrow(md), r1$estimate, r1$p.value, r2$estimate))
  byct <- md[, .(n = .N, AP1_mean = mean(AP1_score1), AP1_median = median(AP1_score1),
                 dissoc_mean = mean(genes_dissoc1)), by = celltype][order(-AP1_mean)]
  fwrite(byct, file.path(OUT, sprintf("%s_AP1_score_by_celltype.csv", tag))); print(byct, digits = 3)
}
score_one(load_main(), "main"); invisible(gc(verbose = FALSE))
score_one(load_replication(), "replication")
