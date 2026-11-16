## 12_motif_annotation_chromvar.R
## JASPAR 2022 motif annotation of the discovery peaks, chromVAR motif activities, differential
## motif activity between cell types, and motif enrichment in cell-type-unique DA peaks.
##
##  1. JASPAR 2022 human motifs (all versions); the 727 tested motifs are the curated list in
##     reference/st4_motif_ids.txt (one version per base motif, heterodimers kept; Table S4).
##  2. AddMotifs (hg19) and RunChromVAR on the annotated ATAC object.
##  3. Differential chromVAR activity (FindAllMarkers, LR test, latent nCount_ATAC; scores are
##     z-scores, so mean.fxn = rowMeans): all cell types, epithelial (LS-1, LS-2, LHr, Bs) and
##     non-epithelial cell types.
##  4. Peak x motif matrix for all peaks.
##  5. Motif enrichment (FindMotifs) in the DA peaks unique to each cell type.
##
## Inputs   results/singlecell/atac_ATAC_RNAlabeled.rds          (05_atac_umap_figures.R)
##          results/singlecell/final_cell_barcodes.csv
##          results/singlecell/peak_annotation.csv               (08_markers_and_da_peaks.R)
##          results/singlecell/combined_peaks_gene_all.csv       (08_markers_and_da_peaks.R)
##          reference/st4_motif_ids.txt
## Outputs  results/singlecell/jaspar2022_motif_list.csv          all human JASPAR 2022 motifs
##          results/singlecell/jaspar2022_motifs_727.csv          the 727 tested motifs
##          results/singlecell/atac_motif.rds                     ATAC object + motifs + chromVAR
##          results/singlecell/Motif_DiffAct_general_cluster.csv
##          results/singlecell/Motif_DiffAct_LS_LHr_Bs.csv
##          results/singlecell/Motif_DiffAct_others.csv
##          results/singlecell/motif_matrix_allPeaks.csv
##          results/singlecell/enriched_motifs_list_uniqDAP_ALL.csv
## Run from the repository root:  Rscript scripts/01_singlecell_annotation/12_motif_annotation_chromvar.R

TNBC_ROOT <- Sys.getenv("TNBC_ROOT", ".")

library(Signac)
library(Seurat)
library(GenomicRanges)
library(BSgenome.Hsapiens.UCSC.hg19)
library(TFBSTools)
library(JASPAR2022)
library(chromVAR)
library(data.table)
library(dplyr)
library(tidyr)

setwd(TNBC_ROOT)
OUT <- "results/singlecell"

## ---- 1. JASPAR 2022 motifs ----
pfm <- getMatrixSet(x = JASPAR2022, opts = list(all_versions = TRUE, species = "Homo sapiens"))
motif_data <- data.frame(
  motif      = sapply(pfm@listData, function(x) x@ID),
  motif.name = sapply(pfm@listData, function(x) x@name),
  family     = sapply(pfm@listData, function(x) if (is.null(x@tags$family)) NA else paste(x@tags$family, collapse = ";"))
)
fwrite(motif_data, file.path(OUT, "jaspar2022_motif_list.csv"), row.names = FALSE)

tested <- readLines("reference/st4_motif_ids.txt"); tested <- tested[nzchar(tested)]
motif_id_name <- motif_data %>%
  filter(motif %in% tested) %>%
  transmute(motif, motif.name,
            base_motif_ids = sub("\\..*$", "", motif),
            motif_version  = as.integer(sub("^.*\\.", "", motif))) %>%
  arrange(motif)
stopifnot(nrow(motif_id_name) == 727)
fwrite(motif_id_name, file.path(OUT, "jaspar2022_motifs_727.csv"), row.names = FALSE)

## ---- 2. motif annotation and chromVAR ----
atac_subset_qc <- readRDS(file.path(OUT, "atac_ATAC_RNAlabeled.rds"))
DefaultAssay(atac_subset_qc) <- "ATAC"
atac_subset_qc <- AddMotifs(object = atac_subset_qc, genome = BSgenome.Hsapiens.UCSC.hg19,
                            pfm = pfm[names(pfm) %in% motif_id_name$motif])
atac_subset_qc <- RunChromVAR(object = atac_subset_qc, genome = BSgenome.Hsapiens.UCSC.hg19)
saveRDS(atac_subset_qc, file.path(OUT, "atac_motif.rds"))

## ---- 3. differential motif activity between cell types ----
Final_cell_barcode <- fread(file.path(OUT, "final_cell_barcodes.csv"), header = TRUE)
atac_subset_qc <- subset(atac_subset_qc, cells = Final_cell_barcode$x)
DefaultAssay(atac_subset_qc) <- "chromvar"

da_motifs <- function(obj) {
  res <- FindAllMarkers(object = obj, only.pos = TRUE, mean.fxn = rowMeans, test.use = "LR",
                        fc.name = "avg_diff", latent.vars = "nCount_ATAC")
  res <- res[res$p_val_adj < 0.05 & res$avg_diff > 0, ] %>% rename("motif" = "gene")
  left_join(motif_id_name, res, by = "motif")
}
write.csv(da_motifs(atac_subset_qc), file.path(OUT, "Motif_DiffAct_general_cluster.csv"), row.names = FALSE)
write.csv(da_motifs(subset(atac_subset_qc, idents = c("LS-1", "LS-2", "LHr", "Bs"))),
          file.path(OUT, "Motif_DiffAct_LS_LHr_Bs.csv"), row.names = FALSE)
write.csv(da_motifs(subset(atac_subset_qc, idents = c("Fb", "Endo", "T", "Pc", "Myl", "B"))),
          file.path(OUT, "Motif_DiffAct_others.csv"), row.names = FALSE)

## ---- 4. peak x motif matrix ----
DefaultAssay(atac_subset_qc) <- "ATAC"
motif_matrix <- GetMotifData(object = atac_subset_qc[["ATAC"]])
motif_matrix_df <- as.data.frame(as.matrix(motif_matrix))
motif_matrix_df$peak <- rownames(motif_matrix_df)
write.csv(motif_matrix_df, file.path(OUT, "motif_matrix_allPeaks.csv"), row.names = FALSE)

## ---- 5. motif enrichment in cell-type-unique DA peaks ----
celltype_levels <- c("LHr", "LS-1", "LS-2", "Bs", "T", "B", "Myl", "Endo", "Pc", "Fb")
Idents(atac_subset_qc) <- factor(Idents(atac_subset_qc), levels = celltype_levels)

## peaks that are differentially accessible in exactly one cell type
peak_anno <- fread(file.path(OUT, "peak_annotation.csv"), header = TRUE) %>% unique() %>%
  mutate(cluster_name = factor(cluster_name, levels = celltype_levels)) %>% arrange(cluster_name)
upset_matrix <- peak_anno %>%
  dplyr::select(peaks, cluster_name) %>% mutate(value = 1) %>%
  pivot_wider(names_from = cluster_name, values_from = value, values_fill = 0) %>%
  as.data.frame()
rownames(upset_matrix) <- upset_matrix$peaks
upset_matrix <- upset_matrix[, -1]
unique_peaks <- function(ct) rownames(upset_matrix)[rowSums(upset_matrix) == 1 & upset_matrix[[ct]] == 1]

combined_peaks_gene_all <- fread(file.path(OUT, "combined_peaks_gene_all.csv"), header = TRUE)
names(combined_peaks_gene_all)[1] <- "peaks"

enriched_motifs_list <- NULL
for (ct in unique(combined_peaks_gene_all$cluster)) {
  features <- combined_peaks_gene_all$peaks[combined_peaks_gene_all$cluster == ct]
  features <- features[features %in% unique_peaks(ct)]
  enriched <- FindMotifs(object = atac_subset_qc, features = features)
  enriched$cluster <- ct
  enriched_motifs_list <- rbind(enriched_motifs_list, enriched)
  cat("Processed cluster:", ct, "\n")
}
## keep -log10(p) finite for plotting
enriched_motifs_list$pvalue <- ifelse(enriched_motifs_list$pvalue == 0, 1e-300,
                                      ifelse(enriched_motifs_list$pvalue == 1, 0.99999999999999999999999,
                                             enriched_motifs_list$pvalue))
write.csv(enriched_motifs_list, file.path(OUT, "enriched_motifs_list_uniqDAP_ALL.csv"), row.names = FALSE)
