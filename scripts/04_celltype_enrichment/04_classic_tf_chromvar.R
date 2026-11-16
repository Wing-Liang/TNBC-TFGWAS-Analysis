# ==============================================================================
# 04_classic_tf_chromvar.R
# chromVAR scores of classic mammary transcription factors across epithelial cell types.
# Produces: Fig S2
# Inputs  : results/singlecell/atac_motif.rds
#           results/singlecell/final_cell_barcodes.csv
#           results/singlecell/jaspar2022_motifs_727.csv
#           results/tfgwas/tfgwas_results_TNBC_ALL.csv
#           results/tfgwas/tfgwas_results_LumA_ALL.csv
#           results/tfgwas/tfgwas_results_LumB_ALL.csv
#           results/tfgwas/tfgwas_results_HER2_ALL.csv
#           ...
# Outputs : results/tfgwas/Heatmap_classic_motif.jpg
#           results/tfgwas/Heatmap_classic_gene_expr.jpg
#           results/tfgwas/heatmap_motifscore_GWASmotif.pdf
#           results/tfgwas/heatmap_motifscore_GWASmotif_LHr_LS_Bs.pdf
#           results/tfgwas/heatmap_motifscore_GWASmotif_uniqMotif_TNBC.pdf
#           results/singlecell/GWASMotif01.all_TNBC_uniqueMotifs_Addmodulescore.jpg
# Run     : Rscript scripts/04_celltype_enrichment/04_classic_tf_chromvar.R   (from the repository root)
# ==============================================================================
## ---- paths (set the environment variables; see README) ----
TNBC_ROOT <- Sys.getenv("TNBC_ROOT", ".")

library(Signac)
library(Seurat)
library(data.table)

library(dplyr)
library(tidyr)

# server
setwd(TNBC_ROOT)

atac_subset_qc<-readRDS("results/singlecell/atac_motif.rds")
Final_cell_barcode<-fread('results/singlecell/final_cell_barcodes.csv',header=TRUE)
#ATAC
atac_subset_qc <- subset(atac_subset_qc , cells = Final_cell_barcode$x)
DefaultAssay(atac_subset_qc)<-'ATAC'

# Extract Motif Accessibility
motif_accessibility <- GetAssayData(atac_subset_qc, assay = "chromvar", layer = "data")

# Extract gene activity scores
RNAgene_expression <- GetAssayData(atac_subset_qc, assay = "RNA", layer = "data")

# Get motif positions
motif.obj <- SeuratObject::GetAssayData(
  object = atac_subset_qc[['ATAC']], layer = "motifs"
)
motif_positions<-Signac::GetMotifData(object = motif.obj)
dim(motif_positions)

# 2. Extract peak counts from the ATAC assay
peak_counts <- GetAssayData(atac_subset_qc, assay = "ATAC", layer = "counts")

# Convert peak names to GRanges if not already
peaks <- granges(atac_subset_qc)

library(SummarizedExperiment)
fragment_counts <- SummarizedExperiment(
  assays = list(counts = peak_counts),
  rowRanges = peaks,
  colData = atac_subset_qc@meta.data
)

# matched peaks with similar GC content and accessibility
library(chromVAR)
library(motifmatchr)
library(BSgenome.Hsapiens.UCSC.hg19)
fragment_counts <- addGCBias(
  object = fragment_counts,
  genome = BSgenome.Hsapiens.UCSC.hg19
)

# Get motif positions from the Seurat object
motif.obj <- SeuratObject::GetAssayData(
  object = atac_subset_qc[['ATAC']],
  layer = "motifs"
)
motif_positions <- Signac::GetMotifData(object = motif.obj)

# Ensure motif_positions aligns with peaks
nrow(fragment_counts) == nrow(motif_positions)

# Compute deviations
dev <- computeDeviations(
  object = fragment_counts,
  annotations = motif_positions)

# Extract deviation scores
motif_scores <- deviations(dev)
# Positive deviation: The motif is more accessible than expected (potentially indicating TF activity).
# Negative deviation: The motif is less accessible than expected.
# Near zero: Accessibility matches the background expectation.

# Add motif deviations to Seurat object
# Create a new assay for motif deviations
motif_assay <- CreateAssayObject(counts =motif_scores)
atac_subset_qc[["MotifDeviations"]] <- motif_assay
dim(atac_subset_qc[["MotifDeviations"]])
# Get unique cell types
cell_types <- unique(atac_subset_qc@meta.data$celltype_RNA)

# Average scores across cell types instead of metacells
DefaultAssay(atac_subset_qc)<-'RNA'
getCellTypeScores <- function(matrix, seurat_obj, cell_types) {
  scores <- lapply(cell_types, function(ct) {
    cells <- rownames(seurat_obj@meta.data)[seurat_obj@meta.data$celltype_RNA == ct]
    if(length(cells) > 1) {
      rowMeans(as.matrix(matrix[, cells, drop = FALSE]))
    } else {
      matrix[, cells, drop = FALSE]
    }
  })
  do.call(cbind, scores)
}

# Calculate scores per cell type
motif_scores_ct <- getCellTypeScores(motif_scores, atac_subset_qc, cell_types)
colnames(motif_scores_ct) <- cell_types
colnames(motif_scores_ct) <-c("Fb","LHS","Endo", "LASP-1","B-Myo","Pc" ,"LASP-2","T","Myl", "B")

####################------------- heat map ----------##################
motif_id_name<-fread('results/singlecell/jaspar2022_motifs_727.csv', header = TRUE)

motif_scores_ct <- subset(motif_scores_ct,rownames(motif_scores_ct) %in%  motif_id_name$motif)
motif_scores_ct <- as.data.frame(motif_scores_ct)
motif_scores_ct$motif <- row.names(motif_scores_ct)

motif_list <- c('ELF5','GRHL1','GRHL2','TFCP2','TBX2','EOMES','TBR1',
                'TBX21','MEIS2','TBX20',"SOX9",'EBF1','PBX2','NHLH1','EBF3','FERD3L',
                'FOXC2','FOXB1','FOXD2','FOXE1','FOXC1',
                'GATA2','GATA4','GATA6','GATA5','GATA3',
                'FOXP1','FOXA2','ESR1','FOXA3','HNF1A',
                'TP53','TP73','TP63','SRF','KLF17','SMAD5','ALX3','ZNF263',
                'PRRX1','MIXL1','LHX9','GBX2','SMAD3','HAND2','SHOX')
motif_list <- data.frame(motif.name=motif_list)
motif_id_list <- left_join(motif_list, motif_id_name, by='motif.name')
motif_score_selected <- left_join(motif_id_list, motif_scores_ct, by='motif')
motif_score_selected$motif.name <- make.unique(motif_score_selected$motif.name)
row.names(motif_score_selected)<-motif_score_selected$motif.name
motif_score_selected <- select(motif_score_selected, `LASP-1`,`LASP-2`,LHS,`B-Myo`)
head(motif_score_selected)

motif_df<-as.data.frame(motif_score_selected)
motif_df$motif<-rownames(motif_df)
motif_df <- motif_df[rownames(motif_df)!='MEIS2.1',]

## filter the breast cancer subtypes
TNBC<-fread("results/tfgwas/tfgwas_results_TNBC_ALL.csv", header = TRUE)
TNBC_selected <- TNBC %>% filter(motif.name %in% motif_list & motif %in% motif_df$motif)
TNBC_selected$TNBC<- -log10(TNBC_selected$TF_p_value)

LumA<-fread("results/tfgwas/tfgwas_results_LumA_ALL.csv", header = TRUE)
LumA_selected <- LumA %>% filter(motif.name %in% motif_list & motif %in% motif_df$motif)
LumA_selected$LumA <- -log10(LumA_selected$TF_p_value)

LumB<-fread("results/tfgwas/tfgwas_results_LumB_ALL.csv", header = TRUE)
LumB_selected <- LumB %>% filter(motif.name %in% motif_list & motif %in% motif_df$motif)
LumB_selected$LumB <- -log10(LumB_selected$TF_p_value)

HER2<-fread("results/tfgwas/tfgwas_results_HER2_ALL.csv", header = TRUE)
HER2_selected <- HER2 %>% filter(motif.name %in% motif_list & motif %in% motif_df$motif)
HER2_selected$HER2 <- -log10(HER2_selected$TF_p_value)

heatmap_df <- left_join(select(TNBC_selected, motif.name, TNBC),
                        select(LumA_selected, motif.name, LumA),
                        by='motif.name')
heatmap_df <- left_join(heatmap_df,
                        select(LumB_selected, motif.name, LumB),
                        by='motif.name')
heatmap_df <- left_join(heatmap_df,
                        select(HER2_selected, motif.name, HER2),
                        by='motif.name')
heatmap_df<- heatmap_df[match(heatmap_rownames, heatmap_df$motif.name), ]

# for the motif score heatmap
motif_df$motif<-NULL
motif_df$motif.name<-NULL
motif_matrix <- as.matrix(motif_df)

scaled_mat <- t(scale(t(motif_matrix)))

library(ComplexHeatmap)
library(dendextend)
set.seed(2025)
heatmap_obj <-Heatmap(scaled_mat,
        name = "Motif score",          # Title of legend
        row_names_gp = gpar(fontsize = 7),
        row_names_side = "left",
        column_names_centered = TRUE,
        width = unit(50, "mm"),
        col = circlize::colorRamp2(c(-1.5,0, 1.5), c("darkblue","grey95","red1")),
        rect_gp = gpar(col = "grey99"), # Set cell border color (instead of border_color)
        cluster_rows = FALSE,  # Preserve the order from motif_list
        cluster_columns = FALSE,  # Cluster columns (or set to FALSE if you want to preserve column order)
        show_row_names = TRUE,  # Show motif names
        show_column_names = TRUE,  # Show column names
        # row_split = 3,                 # Ensure 3 clusters without labels
        column_order = c('LASP-1','LASP-2','LHS','B-Myo'),
        row_title = NULL               # Remove cluster titles/labels (1, 2, 3)
)
jpeg('results/tfgwas/Heatmap_classic_motif.jpg',width=7, height = 8.5, units = 'in', res = 150)
heatmap_obj
dev.off()

# gene expression
library(ggplot2)
DefaultAssay(atac_subset_qc)<-'RNA'
p<-DotPlot(object = atac_subset_qc, features = rev(motif_list), idents=c('LASP-1','LASP-2','LHS','B-Myo')) +
  coord_flip()
jpeg('results/tfgwas/Heatmap_classic_gene_expr.jpg',width=5, height = 8, units = 'in', res = 150)
p
dev.off()

# gene expression data
dotplot_data <- p$data
expr_data <- dotplot_data %>% dplyr::select(features.plot, id, avg.exp.scaled)
names(expr_data)<-c('gene','cluster','avg.exp.scaled')
expr_data$avg.exp.scaled<-ifelse(is.na(expr_data$avg.exp.scaled),0,expr_data$avg.exp.scaled)
library(tidyr)
# Reshape the data to wide format
reshaped_data <- pivot_wider(expr_data,
                             names_from = cluster,
                             values_from = avg.exp.scaled,
                             id_cols = gene)
row.names(reshaped_data)<-reshaped_data$gene
reshaped_data$gene<-NULL
expr_data_matrix <- as.matrix(reshaped_data)

scaled_mat_expr <- t(scale(t(expr_data_matrix)))

library(ComplexHeatmap)
library(dendextend)
set.seed(2025)
heatmap_obj <-Heatmap(scaled_mat_expr,
                      name = "Average gene expression",          # Title of legend
                      row_names_gp = gpar(fontsize = 7),
                      row_names_side = "left",
                      column_names_centered = TRUE,
                      width = unit(50, "mm"),
                      col = circlize::colorRamp2(c(-1.5,0, 1.5), c("brown4","grey95","orange")),
                      rect_gp = gpar(col = "grey99"), # Set cell border color (instead of border_color)
                      cluster_rows = FALSE,  # Preserve the order from motif_list
                      cluster_columns = FALSE,  # Cluster columns (or set to FALSE if you want to preserve column order)
                      show_row_names = TRUE,  # Show motif names
                      show_column_names = TRUE,  # Show column names
                      # row_split = 3,                 # Ensure 3 clusters without labels
                      column_order = c("LS-1", "LS-2", "LHr", "Bs"),
                      row_title = NULL               # Remove cluster titles/labels (1, 2, 3)
)
jpeg('results/tfgwas/Heatmap_classic_gene_expr.jpg',width=5, height = 8, units = 'in', res = 150)
heatmap_obj
dev.off()

# Extract reordered row names
row_order_indices <- row_order(heatmap_obj)
reordered_row_indices <- unlist(row_order_indices)
heatmap_rownames <- rownames(scaled_mat)[reordered_row_indices]

tnbc_heatmap<-Heatmap(heatmap_df$TNBC,
        name = "TNBC", width = unit(5, "mm"),
        cluster_rows = FALSE,
          col = circlize::colorRamp2(c(0, 10), c("grey95", "#ffbf00")))
luma_heatmap<-Heatmap(heatmap_df$LumA,
         name = "LumA", width = unit(5, "mm"),
         cluster_rows = FALSE,
         col = circlize::colorRamp2(c(0, 20), c("grey95","#ff9f00" )))
lumb_heatmap<-Heatmap(heatmap_df$LumB,
                      name = "LumB", width = unit(5, "mm"),
                      cluster_rows = FALSE,
                      col = circlize::colorRamp2(c(0, 10), c("grey95", "#ff8f00" )))
her2_heatmap<-Heatmap(heatmap_df$HER2,
                      name = "HER2", width = unit(5, "mm"),
                      cluster_rows = FALSE,
                      col = circlize::colorRamp2(c(0, 5), c("grey95","#ff8c00")))
combined_heatmap <-  tnbc_heatmap + luma_heatmap + lumb_heatmap + her2_heatmap

heatmap_obj+combined_heatmap

pdf("results/tfgwas/heatmap_motifscore_GWASmotif.pdf", width = 6, height = 6)
draw(heatmap_obj+combined_heatmap)
dev.off()

################ - top 10 unique DAmotif (larger avg_diff if duplicates) across cell types -######
subset.differential.activity<- fread('results/singlecell/Motif_DiffAct_LS_LHr_Bs.csv',header = TRUE)
subset.differential.activity <- subset.differential.activity%>% dplyr::filter(avg_diff>0.5)
subset.differential.activity  <-subset.differential.activity %>%
  group_by(motif.name) %>%
  filter(p_val == min(p_val )) %>%
  ungroup()
table(subset.differential.activity$cluster)
top10_DAMotif <-subset.differential.activity %>%
  group_by(cluster) %>%
  arrange(p_val, -abs(pct.1 -pct.2))%>%
  slice_head(n = 10)  %>%
  ungroup()
table(top10_DAMotif$cluster)

motif_id_name<-fread('results/singlecell/jaspar2022_motifs_727.csv', header = TRUE)

motif_list <- top10_DAMotif$motif
motif_id_list <- motif_id_name[motif_id_name$motif %in% motif_list,]
motif_score_selected <- subset(motif_scores_ct,rownames(motif_scores_ct) %in% motif_id_list$motif)
motif_score_selected <- motif_score_selected[,c(2,4,5,7)]
head(motif_score_selected)

motif_df<-as.data.frame(motif_score_selected)
motif_df$motif<-rownames(motif_df)

motif_df <- left_join(motif_df, dplyr::select(motif_id_list, motif, motif.name), by='motif')
motif_df$motif.name <- make.unique(motif_df$motif.name)
rownames(motif_df) <- motif_df$motif.name

## filter the breast cancer subtypes
TNBC<-fread("results/tfgwas/tfgwas_results_TNBC_ALL.csv", header = TRUE)
TNBC$motif.name <- NULL
TNBC <- unique(TNBC)
TNBC_selected <- TNBC %>% dplyr::filter(motif %in% motif_list & motif.name.orig %in% motif_df$motif.name)
TNBC_selected$TNBC<- -log10(TNBC_selected$TF_p_value)
TNBC_selected <- TNBC_selected %>% rename('motif.name.orig'='motif.name')

LumA<-fread("results/tfgwas/tfgwas_results_LumA_ALL.csv", header = TRUE)
LumA$motif.name <- NULL
LumA <- unique(LumA)
LumA_selected <- LumA %>% filter(motif %in% motif_list & motif.name.orig %in% motif_df$motif.name)
LumA_selected$LumA <- -log10(LumA_selected$TF_p_value)
LumA_selected <- LumA_selected %>% rename('motif.name.orig'='motif.name')

LumB<-fread("results/tfgwas/tfgwas_results_LumB_ALL.csv", header = TRUE)
LumB$motif.name <- NULL
LumB <- unique(LumB)
LumB_selected <- LumB %>% filter(motif %in% motif_list & motif.name.orig %in% motif_df$motif.name)
LumB_selected$LumB <- -log10(LumB_selected$TF_p_value)
LumB_selected <- LumB_selected %>% rename('motif.name.orig'='motif.name')

HER2<-fread("results/tfgwas/tfgwas_results_HER2_ALL.csv", header = TRUE)
HER2$motif.name <- NULL
HER2 <- unique(HER2)
HER2_selected <- HER2 %>% filter(motif %in% motif_list & motif.name.orig %in% motif_df$motif.name)
HER2_selected$HER2 <- -log10(HER2_selected$TF_p_value)
HER2_selected <- HER2_selected %>% rename('motif.name.orig'='motif.name')

heatmap_df <- left_join(dplyr::select(TNBC_selected, motif.name, TNBC),
                        dplyr::select(LumA_selected, motif.name, LumA),
                        by='motif.name')
heatmap_df <- left_join(heatmap_df,
                        dplyr::select(LumB_selected, motif.name, LumB),
                        by='motif.name')
heatmap_df <- left_join(heatmap_df,
                        dplyr::select(HER2_selected, motif.name, HER2),
                        by='motif.name')

# for the motif score heatmap
motif_df$motif<-NULL
motif_df$motif.name<-NULL
motif_matrix <- as.matrix(motif_df)

scaled_mat <- t(scale(t(motif_matrix)))

library(ComplexHeatmap)
library(dendextend)
heatmap_obj <-Heatmap(scaled_mat,
                      name = "Motif score",          # Title of legend
                      row_names_gp = gpar(fontsize = 7),
                      row_names_side = "left",
                      column_names_centered = TRUE,
                      width = unit(50, "mm"),
                      col = circlize::colorRamp2(c(-1.5,0, 1.5), c("#483d8b","grey95","red1")),
                      km = 4,                        # K-means clustering with 10 clusters
                      rect_gp = gpar(col = "grey99"), # Set cell border color (instead of border_color)
                      show_row_dend = FALSE,         # Hide row dendrogram to avoid km labels
                      show_column_dend = FALSE,
                      # row_split = 10,                 # Ensure 3 clusters without labels
                      row_title = NULL               # Remove cluster titles/labels (1, 2, 3)
)

# Extract reordered row names
row_order_indices <- row_order(heatmap_obj)
reordered_row_indices <- unlist(row_order_indices)
heatmap_rownames <- rownames(scaled_mat)[reordered_row_indices]

heatmap_df<- heatmap_df[match(heatmap_rownames, heatmap_df$motif.name), ]

tnbc_heatmap<-Heatmap(heatmap_df$TNBC,
                      name = "TNBC", width = unit(5, "mm"),
                      cluster_rows = FALSE,
                      col = circlize::colorRamp2(c(0, 10), c("grey95", "#ffbf00")))
luma_heatmap<-Heatmap(heatmap_df$LumA,
                      name = "LumA", width = unit(5, "mm"),
                      cluster_rows = FALSE,
                      col = circlize::colorRamp2(c(0, 20), c("grey95","#ff9f00" )))
lumb_heatmap<-Heatmap(heatmap_df$LumB,
                      name = "LumB", width = unit(5, "mm"),
                      cluster_rows = FALSE,
                      col = circlize::colorRamp2(c(0, 10), c("grey95", "#ff8f00" )))
her2_heatmap<-Heatmap(heatmap_df$HER2,
                      name = "HER2", width = unit(5, "mm"),
                      cluster_rows = FALSE,
                      col = circlize::colorRamp2(c(0, 5), c("grey95","#ff8c00")))
combined_heatmap <-  tnbc_heatmap + luma_heatmap + lumb_heatmap + her2_heatmap

heatmap_obj+combined_heatmap

pdf("results/tfgwas/heatmap_motifscore_GWASmotif_LHr_LS_Bs.pdf", width = 6, height = 6)
draw(heatmap_obj+combined_heatmap)
dev.off()

# box plot
library(tidyverse)

# Convert matrix or data.frame from wide to long format
motif_scores_ct <- subset(motif_scores_ct,rownames(motif_scores_ct) %in%  TNBC_unique$motif)
motif_scores_ct <- as.data.frame(motif_scores_ct)  %>% select(`LS-1`, Fb,LHr ,`LS-2`,Bs,Endo,T, Pc,Myl, B)
motif_scores_long <- motif_scores_ct %>%
  rownames_to_column(var = "Motif") %>%
  pivot_longer(-Motif, names_to = "CellType", values_to = "Score")

# Create the box plot
ggplot(motif_scores_long, aes(x = CellType, y = Score)) +
  geom_boxplot(fill = "skyblue", alpha = 0.7, outlier.size = 0.5) +
  labs(title = "Motif Scores by Cell Type",
       x = "Cell Type",
       y = "Motif Score") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  coord_cartesian(ylim = c(-0.1, 0.1))

################ - top 10 unique DAmotif (larger avg_diff if duplicates) across cell types -######
subset.differential.activity<- fread('results/singlecell/Motif_DiffAct_others.csv',header = TRUE)
subset.differential.activity <- subset.differential.activity%>% dplyr::filter(avg_diff>0.5)

subset.differential.activity  <-subset.differential.activity %>%
  group_by(motif.name) %>%
  filter(p_val == min(p_val )) %>%
  ungroup()
table(subset.differential.activity$cluster)
top10_DAMotif <-subset.differential.activity %>%
  group_by(cluster) %>%
  arrange(p_val, -abs(pct.1 -pct.2))%>%
  slice_head(n = 10)  %>%
  ungroup()
table(top10_DAMotif$cluster)

motif_id_name<-fread('results/singlecell/jaspar2022_motifs_727.csv', header = TRUE)

motif_list <- top10_DAMotif$motif
motif_id_list <- motif_id_name[motif_id_name$motif %in% motif_list,]
motif_score_selected <- subset(motif_scores_ct,rownames(motif_scores_ct) %in% motif_id_list$motif)
motif_score_selected <- motif_score_selected[,c(1,3,6,8:10)]
head(motif_score_selected)

motif_df<-as.data.frame(motif_score_selected)
motif_df$motif<-rownames(motif_df)

motif_df <- left_join(motif_df, dplyr::select(motif_id_list, motif, motif.name), by='motif')
motif_df$motif.name <- make.unique(motif_df$motif.name)
rownames(motif_df) <- motif_df$motif.name

## filter the breast cancer subtypes
TNBC<-fread("results/tfgwas/tfgwas_results_TNBC_ALL.csv", header = TRUE)
TNBC$motif.name <- NULL
TNBC <- unique(TNBC)
TNBC_selected <- TNBC %>% dplyr::filter(motif %in% motif_list & motif.name.orig %in% motif_df$motif.name)
TNBC_selected$TNBC<- -log10(TNBC_selected$TF_p_value)
TNBC_selected <- TNBC_selected %>% rename('motif.name.orig'='motif.name')

LumA<-fread("results/tfgwas/tfgwas_results_LumA_ALL.csv", header = TRUE)
LumA$motif.name <- NULL
LumA <- unique(LumA)
LumA_selected <- LumA %>% filter(motif %in% motif_list & motif.name.orig %in% motif_df$motif.name)
LumA_selected$LumA <- -log10(LumA_selected$TF_p_value)
LumA_selected <- LumA_selected %>% rename('motif.name.orig'='motif.name')

LumB<-fread("results/tfgwas/tfgwas_results_LumB_ALL.csv", header = TRUE)
LumB$motif.name <- NULL
LumB <- unique(LumB)
LumB_selected <- LumB %>% filter(motif %in% motif_list & motif.name.orig %in% motif_df$motif.name)
LumB_selected$LumB <- -log10(LumB_selected$TF_p_value)
LumB_selected <- LumB_selected %>% rename('motif.name.orig'='motif.name')

HER2<-fread("results/tfgwas/tfgwas_results_HER2_ALL.csv", header = TRUE)
HER2$motif.name <- NULL
HER2 <- unique(HER2)
HER2_selected <- HER2 %>% filter(motif %in% motif_list & motif.name.orig %in% motif_df$motif.name)
HER2_selected$HER2 <- -log10(HER2_selected$TF_p_value)
HER2_selected <- HER2_selected %>% rename('motif.name.orig'='motif.name')

heatmap_df <- left_join(dplyr::select(TNBC_selected, motif, TNBC),
                        dplyr::select(LumA_selected, motif, LumA),
                        by='motif')
heatmap_df <- left_join(heatmap_df,
                        dplyr::select(LumB_selected, motif, LumB),
                        by='motif')
heatmap_df <- left_join(heatmap_df,
                        dplyr::select(HER2_selected, motif, HER2),
                        by='motif')
heatmap_df <- left_join(dplyr::select(motif_df, motif, motif.name),heatmap_df, by='motif')

# for the motif score heatmap
motif_df$motif<-NULL
motif_df$motif.name<-NULL
motif_matrix <- as.matrix(motif_df)

scaled_mat <- t(scale(t(motif_matrix)))
row_dist <- dist(scaled_mat)
row_clust <- hclust(row_dist)
col_clust = hclust(dist(scaled_mat))

library(ComplexHeatmap)
library(dendextend)
heatmap_obj <-Heatmap(scaled_mat,
                      name = "Motif score",          # Title of legend
                      row_names_gp = gpar(fontsize = 7),
                      row_names_side = "left",
                      column_names_centered = TRUE,
                      width = unit(50, "mm"),
                      col = circlize::colorRamp2(c(-1.5,0, 1.5), c("#483d8b","grey95","red1")),
                      km =6,                        # K-means clustering with 10 clusters
                      rect_gp = gpar(col = "grey99"), # Set cell border color (instead of border_color)
                      show_row_dend = FALSE,         # Hide row dendrogram to avoid km labels
                      show_column_dend = FALSE,
                      # row_split = 10,                 # Ensure 3 clusters without labels
                      row_title = NULL               # Remove cluster titles/labels (1, 2, 3)
)
heatmap_obj

# Extract reordered row names
row_order_indices <- row_order(heatmap_obj)
reordered_row_indices <- unlist(row_order_indices)
heatmap_rownames <- rownames(scaled_mat)[reordered_row_indices]

heatmap_df<- heatmap_df[match(heatmap_rownames, heatmap_df$motif.name), ]

tnbc_heatmap<-Heatmap(heatmap_df$TNBC,
                      name = "TNBC", width = unit(5, "mm"),
                      cluster_rows = FALSE,
                      col = circlize::colorRamp2(c(0, 10), c("grey95", "#ffbf00")))
luma_heatmap<-Heatmap(heatmap_df$LumA,
                      name = "LumA", width = unit(5, "mm"),
                      cluster_rows = FALSE,
                      col = circlize::colorRamp2(c(0, 20), c("grey95","#ff9f00" )))
lumb_heatmap<-Heatmap(heatmap_df$LumB,
                      name = "LumB", width = unit(5, "mm"),
                      cluster_rows = FALSE,
                      col = circlize::colorRamp2(c(0, 10), c("grey95", "#ff8f00" )))
her2_heatmap<-Heatmap(heatmap_df$HER2,
                      name = "HER2", width = unit(5, "mm"),
                      cluster_rows = FALSE,
                      col = circlize::colorRamp2(c(0, 5), c("grey95","#ff8c00")))
combined_heatmap <-  tnbc_heatmap + luma_heatmap + lumb_heatmap + her2_heatmap

heatmap_obj+combined_heatmap

pdf("results/tfgwas/heatmap_motifscore_GWAS
    motif_othercells.pdf", width = 6, height = 7)
draw(heatmap_obj+combined_heatmap)
dev.off()

########## ------ overlap TF between LHr, LS, BS and other cells
subset.differential.activity<- fread('results/singlecell/Motif_DiffAct_LS_LHr_Bs.csv',header = TRUE)
subset.differential.activity <- subset.differential.activity%>% dplyr::filter(avg_diff>0.5)

subset.differential.activity  <-subset.differential.activity %>%
  group_by(motif.name) %>%
  filter(p_val == min(p_val )) %>%
  ungroup()
table(subset.differential.activity$cluster)
top10_DAMotif <-subset.differential.activity %>%
  group_by(cluster) %>%
  arrange(p_val, -abs(pct.1 -pct.2))%>%
  slice_head(n = 10)  %>%
  ungroup()
table(top10_DAMotif$cluster)

subset.differential.activity<- fread('results/singlecell/Motif_DiffAct_others.csv',header = TRUE)
subset.differential.activity <- subset.differential.activity%>% dplyr::filter(avg_diff>0.5)
subset.differential.activity  <-subset.differential.activity %>%
  group_by(motif.name) %>%
  filter(p_val == min(p_val )) %>%
  ungroup()
table(subset.differential.activity$cluster)
top10_DAMotif_other <-subset.differential.activity %>%
  group_by(cluster) %>%
  arrange(p_val, -abs(pct.1 -pct.2))%>%
  slice_head(n = 10)  %>%
  ungroup()
table(top10_DAMotif_other$cluster)

top10_DAMotif[top10_DAMotif$motif.name %in% top10_DAMotif_other$motif.name,]

###############################----------- selected from significant motif in TNBC-----#####
# enriched in the uniPeaks and Unique motifs in TNBC
enriched_motifs_list_TNBC_unique<-fread("results/singlecell/enriched_motifs_list_uniqDAP_TNBCUniqSignificant.csv",header = TRUE)
# filter from the enriched motif (from unique peaks) across cell type (the signitificant motif from GWAS)
# cutoff p val < 0.05
enriched_motifs_select <- enriched_motifs_list_TNBC_unique %>% filter(pvalue <0.05)
enriched_motifs_select <- enriched_motifs_select %>%
  group_by(motif.name.orig) %>%
  filter(fold.enrichment == max(fold.enrichment)) %>%
  ungroup()

motif_id_name<-fread('results/singlecell/jaspar2022_motifs_727.csv', header = TRUE)

motif_list <- enriched_motifs_select$motif
motif_id_list <- motif_id_name[motif_id_name$motif %in% motif_list,]
motif_score_selected <- subset(motif_scores_ct,rownames(motif_scores_ct) %in% motif_id_list$motif)
head(motif_score_selected)

motif_df<-as.data.frame(motif_score_selected)
motif_df$motif<-rownames(motif_df)

motif_df <- left_join(motif_df, dplyr::select(motif_id_list, motif, motif.name), by='motif')
motif_df$motif.name <- make.unique(motif_df$motif.name)
rownames(motif_df) <- motif_df$motif.name

## filter the breast cancer subtypes
TNBC<-fread("results/tfgwas/tfgwas_results_TNBC_ALL.csv", header = TRUE)
TNBC$motif.name <- NULL
TNBC <- unique(TNBC)
TNBC_selected <- TNBC %>% dplyr::filter(motif %in% motif_list & motif.name.orig %in% motif_df$motif.name)
TNBC_selected$TNBC<- -log10(TNBC_selected$TF_p_value)
TNBC_selected <- TNBC_selected %>% rename('motif.name.orig'='motif.name')

LumA<-fread("results/tfgwas/tfgwas_results_LumA_ALL.csv", header = TRUE)
LumA$motif.name <- NULL
LumA <- unique(LumA)
LumA_selected <- LumA %>% filter(motif %in% motif_list & motif.name.orig %in% motif_df$motif.name)
LumA_selected$LumA <- -log10(LumA_selected$TF_p_value)
LumA_selected <- LumA_selected %>% rename('motif.name.orig'='motif.name')

LumB<-fread("results/tfgwas/tfgwas_results_LumB_ALL.csv", header = TRUE)
LumB$motif.name <- NULL
LumB <- unique(LumB)
LumB_selected <- LumB %>% filter(motif %in% motif_list & motif.name.orig %in% motif_df$motif.name)
LumB_selected$LumB <- -log10(LumB_selected$TF_p_value)
LumB_selected <- LumB_selected %>% rename('motif.name.orig'='motif.name')

HER2<-fread("results/tfgwas/tfgwas_results_HER2_ALL.csv", header = TRUE)
HER2$motif.name <- NULL
HER2 <- unique(HER2)
HER2_selected <- HER2 %>% filter(motif %in% motif_list & motif.name.orig %in% motif_df$motif.name)
HER2_selected$HER2 <- -log10(HER2_selected$TF_p_value)
HER2_selected <- HER2_selected %>% rename('motif.name.orig'='motif.name')

heatmap_df <- left_join(dplyr::select(TNBC_selected, motif, TNBC),
                        dplyr::select(LumA_selected, motif, LumA),
                        by='motif')
heatmap_df <- left_join(heatmap_df,
                        dplyr::select(LumB_selected, motif, LumB),
                        by='motif')
heatmap_df <- left_join(heatmap_df,
                        dplyr::select(HER2_selected, motif, HER2),
                        by='motif')
heatmap_df <- left_join(dplyr::select(motif_df, motif, motif.name),heatmap_df, by='motif')

# for the motif score heatmap
motif_df$motif<-NULL
motif_df$motif.name<-NULL
motif_matrix <- as.matrix(motif_df)

scaled_mat <- t(scale(t(motif_matrix)))
row_dist <- dist(scaled_mat)
row_clust <- hclust(row_dist)
col_clust = hclust(dist(scaled_mat))

library(ComplexHeatmap)
library(dendextend)
heatmap_obj <-Heatmap(scaled_mat,
                      name = "Motif score",          # Title of legend
                      row_names_gp = gpar(fontsize = 7),
                      row_names_side = "left",
                      column_names_centered = TRUE,
                      width = unit(50, "mm"),
                      col = circlize::colorRamp2(c(-1.8,0, 1.8), c("#483d8b","grey95","red1")),
                      km =9,                        # K-means clustering with 10 clusters
                      rect_gp = gpar(col = "grey99"), # Set cell border color (instead of border_color)
                      show_row_dend = FALSE,         # Hide row dendrogram to avoid km labels
                      show_column_dend = FALSE,
                      column_split = 4,                 # Ensure 3 clusters without labels
                      row_title = NULL ,
                      column_title = NULL
)
heatmap_obj

pdf("results/tfgwas/heatmap_motifscore_GWASmotif_uniqMotif_TNBC.pdf", width = 4, height = 7)
draw(heatmap_obj)
dev.off()

DefaultAssay(atac_subset_qc) <- "chromvar"

TNBC_unique<-list(unique(enriched_motifs_select$motif))
atac_subset_qc <- AddModuleScore(
  object = atac_subset_qc,
  features = TNBC_unique,
  ctrl = 5,
  name = 'TNBC_unique_motif'
)

jpeg('results/singlecell/GWASMotif01.all_TNBC_uniqueMotifs_Addmodulescore.jpg',width=6, height = 5, units = 'in', res = 150)
FeaturePlot(
  object = atac_subset_qc,
  features = "TNBC_unique_motif1",
  min.cutoff = 'q5',
  max.cutoff = 'q95',
  pt.size = 0.1,
  label = TRUE,
  cols = c("#eeeeee", "purple4", "darkblue"))
dev.off()

####################################################################
##################       Gene expression  ##########################
gene_list <- c('ELF5','GRHL1','GRHL2','TFCP2','TBX2','EOMES','TBR1',
                'TBX21','MEIS2','TBX20','EBF1','PBX2','NHLH1','EBF3','FERD3L',
                'FOXC2','FOXB1','FOXD2','FOXE1','FOXC1',
                'GATA2','GATA4','GATA6','GATA5','GATA3',
                'FOXP1','FOXA2','ESR1','FOXA3','HNF1A',
                'TP53','TP73','TP63','SRF','KLF17','SMAD5','ALX3','ZNF263',
                'PRRX1','MIXL1','LHX9','GBX2','SMAD3','HAND2','SHOX')
DefaultAssay(atac_subset_qc)<-'RNA'
DotPlot(object = atac_subset_qc , features = gene_list) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
