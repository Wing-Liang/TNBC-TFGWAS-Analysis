# ==============================================================================
# 02_motif_enrichment_by_subtype.R
# Motif enrichment in cell-type differentially accessible peaks and module scores, per subtype.
# Produces: Figs S3, S4, S7, S8
# Inputs  : results/singlecell/final_cell_barcodes.csv
#           results/singlecell/atac_motif.rds
#           results/singlecell/peak_annotation.csv
#           results/singlecell/combined_peaks_gene_all.csv
#           results/singlecell/Enriched_motif_Endo.csv
#           results/singlecell/jaspar2022_motifs_727.csv
#           results/tfgwas/tfgwas_results_TNBC_ALL.csv
#           ...
# Outputs : results/singlecell/Enriched_motif_Endo.csv
#           results/singlecell/enriched_motifs_list_uniqDAP_TNBCsignificant.csv
#           results/singlecell/Motif_enricheded_heatmap_uniDAP_TNBC_signMotif.jpg
#           results/singlecell/GWASMotif02_FeaturePlot_TNBC_signicant_all_AddModuleScore.jpg
#           results/singlecell/GWASMotif03_FeaturePlot_TNBC_all_colmeans.jpg
#           results/singlecell/enriched_motifs_list_uniqDAP_LumAsignificant.csv
#           results/singlecell/Motif_enricheded_heatmap_uniDAP_LumA_signMotif.jpg
#           ...
# Run     : Rscript scripts/04_celltype_enrichment/02_motif_enrichment_by_subtype.R   (from the repository root)
# ==============================================================================
## ---- paths (set the environment variables; see README) ----
TNBC_ROOT <- Sys.getenv("TNBC_ROOT", ".")
# Load required packages
library(lme4)      # For generalized mixed models
library(dplyr)     # For data manipulation
library(ggplot2)   # For visualization (optional)

# Load required packages
library(Seurat)        # For scATAC-seq object handling
library(Signac)        # For scATAC-seq utilities
library(dplyr)         # For data manipulation
library(GenomicRanges) # For genomic overlaps
library(motifmatchr)   # For motif matching
library(TFBSTools)     # For motif handling
library(data.table)

setwd(TNBC_ROOT)
#############################################
#         heatmap of motif enrichment       #
##############################################
Final_cell_barcode<-fread('results/singlecell/final_cell_barcodes.csv',header=TRUE)

#ATAC
atac_subset_qc<-readRDS('results/singlecell/atac_motif.rds')
atac_subset_qc <- subset(atac_subset_qc , cells = Final_cell_barcode$x)
DefaultAssay(atac_subset_qc)<-'ATAC'

# re-order the levels sequence
# cell-type display order
fig5D_levels <- c('LHr','LS-1','LS-2','Bs','T','B','Myl','Endo','Pc','Fb')
# Reorder the Idents based on sorted levels
Idents(atac_subset_qc) <- factor(Idents(atac_subset_qc), levels = fig5D_levels)
atac_subset_qc$celltype_ATAC <- Idents(atac_subset_qc)

# enriched motif
# DAP unique in the cell
combined_peaks_gene_all <- fread("results/singlecell/peak_annotation.csv", header = TRUE)
combined_peaks_gene_all <- combined_peaks_gene_all %>% unique()%>%
  mutate(cluster_name = factor(cluster_name, levels = fig5D_levels)) %>%
  arrange(cluster_name)

# Convert data into a binary matrix (presence/absence of peaks in each cluster)
library(tidyr)
binary_matrix <- combined_peaks_gene_all %>%
  dplyr::select(peaks, cluster_name) %>%
  mutate(value = 1) %>%
  pivot_wider(names_from = cluster_name, values_from = value, values_fill = 0)

# Remove peaks column for upset plot input
upset_matrix <- as.data.frame(binary_matrix)
rownames(upset_matrix) <- upset_matrix$peaks
upset_matrix <- upset_matrix[ , -1]
upset_matrix$SumEqualsOne <- rowSums(upset_matrix) == 1
Endo<-upset_matrix %>% dplyr::filter(SumEqualsOne==TRUE & Endo ==1) %>% dplyr::select(Endo)
Endo[,1]<-rownames(Endo)
names(Endo) <- 'peaks'
nrow(Endo)

#top5 DAPs
combined_peaks_gene_all<-fread("results/singlecell/combined_peaks_gene_all.csv", header = TRUE)

combined_peaks_gene_all<-combined_peaks_gene_all%>%rename('gene'='peaks')

# Get the unique cluster names
cluster_names <- unique(combined_peaks_gene_all$cluster)[2:10]

# the first cluster - Endo
feature1<-combined_peaks_gene_all$peaks[combined_peaks_gene_all$cluster==unique(combined_peaks_gene_all$cluster)[1]] # Endo
feature1 <- feature1[feature1 %in% Endo$peaks]
# test enrichment
DefaultAssay(atac_subset_qc)<-'ATAC'
enriched.motifs1 <- FindMotifs(
  object = atac_subset_qc,
  features = feature1)

fwrite(enriched.motifs1, 'results/singlecell/Enriched_motif_Endo.csv', row.names = FALSE)
enriched.motifs1<-fread('results/singlecell/Enriched_motif_Endo.csv', header = TRUE)
#################################################################################################
#                              TNBC significant GWAS motif                           #
#################################################################################################
motif_id_list_final<-fread('results/singlecell/jaspar2022_motifs_727.csv', header = TRUE)
TNBC<-fread("results/tfgwas/tfgwas_results_TNBC_ALL.csv", header = TRUE)
TNBC<- left_join(select(motif_id_list_final,motif, motif.name), TNBC, by=c('motif','motif.name')) %>% dplyr::filter(TF_p_value <= 0.05/727)%>%unique()
length(unique(TNBC$motif))
enriched.motifs1$motif.name<-NULL  # not correct, have make.unique names
enriched_motifs_list_TNBC <-left_join(TNBC, enriched.motifs1, by='motif') %>% unique()

enriched_motifs_list_TNBC$cluster <- unique(combined_peaks_gene_all$cluster)[1]

# Loop through each cluster name
for (cluster in cluster_names) {
  # select the unique peaks
  cluster_peak<-upset_matrix %>% dplyr::filter(SumEqualsOne == TRUE & .data[[cluster]] == 1) %>%
    dplyr::select(cluster)
  cluster_peak[,1]<-rownames(cluster_peak)
  names(cluster_peak) <- 'peaks'

  print(cluster)
  nrow(cluster_peak)

  feature1 <-combined_peaks_gene_all$peaks[combined_peaks_gene_all$cluster == cluster]
  feature1 <- feature1[feature1 %in% cluster_peak$peaks]

  # Find motifs for the current cluster
  enriched_motifs <- FindMotifs(
    object = atac_subset_qc,
    features = feature1)

  # Filter the enriched motifs to include only those in significant GWAS
  enriched_motifs$motif.name<-NULL
  enriched_motifs_filtered <- left_join(TNBC, enriched_motifs, by='motif') %>% unique()
  enriched_motifs_filtered$cluster <- cluster

  # Store the filtered results in the list with the cluster name as the key
  enriched_motifs_list_TNBC <- rbind(enriched_motifs_list_TNBC,enriched_motifs_filtered)

  # Print a message to indicate progress
  cat("Processed cluster:", cluster, "\n")
}

enriched_motifs_list_TNBC$pvalue <- ifelse(enriched_motifs_list_TNBC$pvalue==0,1e-300,
                                           ifelse(enriched_motifs_list_TNBC$pvalue==1,0.99999999999999999999999,enriched_motifs_list_TNBC$pvalue))
# Save the enrichment results for the TNBC-significant motifs
write.csv(enriched_motifs_list_TNBC, file = "results/singlecell/enriched_motifs_list_uniqDAP_TNBCsignificant.csv", row.names = FALSE)
enriched_motifs_list_TNBC<-fread("results/singlecell/enriched_motifs_list_uniqDAP_TNBCsignificant.csv", header=TRUE)

# create a matrix
# Get all unique motif names across all clusters
all_motifs <- unique(enriched_motifs_list_TNBC$motif)

# Get all cluster names
all_clusters <- unique(enriched_motifs_list_TNBC$cluster)

# Initialize an empty matrix
motif_matrix <- matrix(
  0,
  ncol = length(all_clusters),
  nrow = length(all_motifs),
  dimnames = list( all_motifs,all_clusters)
)

# Fill the matrix with p-values for enriched motifs in each cluster
for (cluster1 in all_clusters) {
  # Get the enriched motifs and their p-values for the current cluster
  enriched_data <- enriched_motifs_list_TNBC[enriched_motifs_list_TNBC$cluster==cluster1]

  # Add p-values to the matrix
  for (i in 1:nrow(enriched_data)) {
    motif <- enriched_data$motif[i]
    neglog10pvalue <- -log10(enriched_data$pvalue[i])  # Replace "pvalue" with the actual column name for p-values
    motif_matrix[motif,cluster1] <-neglog10pvalue
  }
}

# Convert the matrix to a data frame for better readability
motif_df<- as.data.frame(motif_matrix)
motif_df$motif<-rownames(motif_df)
motif_df <- left_join(motif_df, select(motif_id_list_final, motif, motif.name), by='motif')
motif_df$motif.name <- make.unique(motif_df$motif.name)
rownames(motif_df) <- motif_df$motif.name
motif_df$motif<-NULL
motif_df$motif.name<-NULL
motif_matrix <- as.matrix(motif_df)

motif_df<-dplyr::select(motif_df,'LHr','LS-1','LS-2','Bs','T','B','Myl','Endo','Pc','Fb')
library(pheatmap)
library(ComplexHeatmap)
library(circlize)

motif_label <- motif_df[apply(motif_df, 1, function(x) any(x >= 300)), ]

motif_df<-motif_df%>%arrange(desc(LHr),desc(`LS-1`), desc(Endo))
TF_list<-rownames(motif_df)
cell_list<-colnames(motif_df)
motif_heatmap<-SingleRasterMap(
  data = t(motif_df),
  raster = TRUE,
  cell.order = cell_list,
  feature.order = rev(TF_list),
  colors = colorRampPalette(c('lightyellow',"purple4","darkblue"))(300),   # Custom color mapping
  disp.min = 0,      # Minimum display value
  disp.max = 300,    # Maximum display value
  limits = c(0, 300) # Set explicit limits
)
num_cols <- ncol(motif_df)
line_positions <- seq(1.5, num_cols - 0.5, by = 1)  # Positions between columns
motif_heatmap<-motif_heatmap +
  geom_vline(xintercept = line_positions, color = "purple4", size = 0.25)+
  scale_y_discrete(
    breaks = rownames(motif_label),                 # Only label significant genes
    labels = rownames(motif_label)) +
  labs(fill = "-log10(pval)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1,size = 12),
        axis.text.y = element_text(angle = 0,size = 5)  # Adjust y-axis font size
  )
motif_heatmap
ggsave("results/singlecell/Motif_enricheded_heatmap_uniDAP_TNBC_signMotif.jpg",
       motif_heatmap, width = 5.5, height = 10, dpi = 200)

#----------------- Set chromvar as default assay---------------------------
DefaultAssay(atac_subset_qc) <- "chromvar"

TNBC_significant_motif_ls<-list(unique(TNBC$motif))
atac_subset_qc <- AddModuleScore(
  object = atac_subset_qc,
  features = TNBC_significant_motif_ls,
  ctrl = 5,
  name = 'TNBC_sinificant_gmotif_AddModuleScore'
)

summary(atac_subset_qc$TNBC_sinificant_gmotif_AddModuleScore1)

jpeg('results/singlecell/GWASMotif02_FeaturePlot_TNBC_signicant_all_AddModuleScore.jpg',width=10, height = 5, units = 'in', res = 150)
a<-FeaturePlot(
  object = subset(atac_subset_qc, idents=c('LHr','LS-1','LS-2','Bs')),
  features = "TNBC_sinificant_gmotif_AddModuleScore1",
  min.cutoff = 'q5',
  max.cutoff = 'q95',
  pt.size = 0.1,
  label = TRUE
  )
b<-FeaturePlot(
  object = subset(atac_subset_qc, idents=c('T','B','Myl','Endo','Pc','Fb')),
  features = "TNBC_sinificant_gmotif_AddModuleScore1",
  min.cutoff = 'q5',
  max.cutoff = 'q95',
  pt.size = 0.1,
  label = TRUE
)
a&NoLegend()|b
dev.off()

# colMean
# Compute average chromVAR score across motifs
motif_accessibility <- GetAssayData(atac_subset_qc, assay = "chromvar", layer = "data")

atac_subset_qc$TNBC_avg_score_all <- Matrix::colMeans(motif_accessibility[unique(TNBC$motif), ])
length(atac_subset_qc$TNBC_avg_score_all)

jpeg('results/singlecell/GWASMotif03_FeaturePlot_TNBC_all_colmeans.jpg',width=10, height = 5, units = 'in', res = 150)
a<-FeaturePlot(
  object = subset(atac_subset_qc, idents=c('LHr','LS-1','LS-2','Bs')),
  features = "TNBC_avg_score_all",
  min.cutoff = 'q5',
  max.cutoff = 'q95',
  pt.size = 0.1,
  label = TRUE,
  cols = c("#eeeeee", "purple4", "darkblue"))
b<-FeaturePlot(
  object = subset(atac_subset_qc, idents=c('T','B','Myl','Endo','Pc','Fb')),
  features = "TNBC_avg_score_all",
  min.cutoff = 'q5',
  max.cutoff = 'q95',
  pt.size = 0.1,
  label = TRUE,
  cols = c("#eeeeee", "purple4", "darkblue")
)
a&NoLegend()|b
dev.off()
#################################################################################################
#                              LumA significant GWAS motif                           #
#################################################################################################
motif_id_list_final<-fread('results/singlecell/jaspar2022_motifs_727.csv', header = TRUE)
LumA<-fread("results/tfgwas/tfgwas_results_LumA_ALL.csv", header = TRUE)
LumA<- left_join(select(motif_id_list_final,motif, motif.name), LumA, by=c('motif','motif.name')) %>% dplyr::filter(TF_p_value <= 0.05/727)
length(unique(LumA$motif))
enriched.motifs1$motif.name<-NULL
enriched_motifs_list_LumA <-left_join(LumA, enriched.motifs1, by='motif') %>% unique()

enriched_motifs_list_LumA$cluster <- unique(combined_peaks_gene_all$cluster)[1]

DefaultAssay(atac_subset_qc)<-'ATAC'
# Loop through each cluster name
for (cluster in cluster_names) {
  # select the unique peaks
  cluster_peak<-upset_matrix %>% dplyr::filter(SumEqualsOne == TRUE & .data[[cluster]] == 1) %>%
    dplyr::select(cluster)
  cluster_peak[,1]<-rownames(cluster_peak)
  names(cluster_peak) <- 'peaks'

  print(cluster)
  nrow(cluster_peak)

  feature1 <-combined_peaks_gene_all$peaks[combined_peaks_gene_all$cluster == cluster]
  feature1 <- feature1[feature1 %in% cluster_peak$peaks]

  # Find motifs for the current cluster
  enriched_motifs <- FindMotifs(
    object = atac_subset_qc,
    features = feature1)

  # Filter the enriched motifs to include only those in significant GWAS
  enriched_motifs$motif.name<-NULL
  enriched_motifs_filtered <- left_join(LumA, enriched_motifs, by='motif') %>% unique()
  enriched_motifs_filtered$cluster <- cluster

  # Store the filtered results in the list with the cluster name as the key
  enriched_motifs_list_LumA <- rbind(enriched_motifs_list_LumA,enriched_motifs_filtered)

  # Print a message to indicate progress
  cat("Processed cluster:", cluster, "\n")
}

enriched_motifs_list_LumA$pvalue <- ifelse(enriched_motifs_list_LumA$pvalue==0,1e-300,
                                           ifelse(enriched_motifs_list_LumA$pvalue==1,0.99999999999999999999999,enriched_motifs_list_LumA$pvalue))

# Save the motif_cor_top5DAP to a file
write.csv(enriched_motifs_list_LumA, file = "results/singlecell/enriched_motifs_list_uniqDAP_LumAsignificant.csv", row.names = FALSE)

# create a matrix
# Get all unique motif names across all clusters
all_motifs <- unique(enriched_motifs_list_LumA$motif)

# Get all cluster names
all_clusters <- unique(enriched_motifs_list_LumA$cluster)

# Initialize an empty matrix
motif_matrix <- matrix(
  0,
  ncol = length(all_clusters),
  nrow = length(all_motifs),
  dimnames = list( all_motifs,all_clusters)
)

# Fill the matrix with p-values for enriched motifs in each cluster
for (cluster1 in all_clusters) {
  # Get the enriched motifs and their p-values for the current cluster
  enriched_data <- enriched_motifs_list_LumA[enriched_motifs_list_LumA$cluster==cluster1]

  # Add p-values to the matrix
  for (i in 1:nrow(enriched_data)) {
    motif <- enriched_data$motif[i]
    neglog10pvalue <- -log10(enriched_data$pvalue[i])  # Replace "pvalue" with the actual column name for p-values
    motif_matrix[motif,cluster1] <-neglog10pvalue
  }
}

# Convert the matrix to a data frame for better readability
motif_df<- as.data.frame(motif_matrix)
motif_df$motif<-rownames(motif_df)
motif_df <- left_join(motif_df, select(motif_id_list_final, motif, motif.name), by='motif')
motif_df$motif.name <- make.unique(motif_df$motif.name)
rownames(motif_df) <- motif_df$motif.name
motif_df$motif<-NULL
motif_df$motif.name<-NULL
motif_matrix <- as.matrix(motif_df)

motif_df<-dplyr::select(motif_df,'LHr','LS-1','LS-2','Bs','T','B','Myl','Endo','Pc','Fb')
library(pheatmap)
library(ComplexHeatmap)
library(circlize)

motif_label <- motif_df[apply(motif_df, 1, function(x) any(x >= 300)), ]

motif_df<-motif_df%>%arrange(desc(LHr),desc(`LS-1`),desc(Endo))
TF_list<-rownames(motif_df)
cell_list<-colnames(motif_df)
motif_heatmap<-SingleRasterMap(
  data = t(motif_df),
  raster = TRUE,
  cell.order = cell_list,
  feature.order = rev(TF_list),
  colors = colorRampPalette(c('lightyellow',"purple4","darkblue"))(300),   # Custom color mapping
  disp.min = 0,      # Minimum display value
  disp.max = 300,    # Maximum display value
  limits = c(0, 300) # Set explicit limits
)
num_cols <- ncol(motif_df)
line_positions <- seq(1.5, num_cols - 0.5, by = 1)  # Positions between columns
motif_heatmap<-motif_heatmap +
  geom_vline(xintercept = line_positions, color = "purple4", size = 0.25)+
  scale_y_discrete(
    breaks = rownames(motif_label),                 # Only label significant genes
    labels = rownames(motif_label)) +
  labs(fill = "-log10(pval)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1,size = 12),
        axis.text.y = element_text(angle = 0,size = 3)  # Adjust y-axis font size
  )
motif_heatmap
ggsave("results/singlecell/Motif_enricheded_heatmap_uniDAP_LumA_signMotif.jpg",
       motif_heatmap, width = 5.5, height = 11, dpi = 200)

#----------------- Set chromvar as default assay---------------------------
DefaultAssay(atac_subset_qc) <- "chromvar"

LumA_significant_motif_ls<-list(unique(LumA$motif))
atac_subset_qc <- AddModuleScore(
  object = atac_subset_qc,
  features = LumA_significant_motif_ls,
  ctrl = 5,
  name = 'LumA_sinificant_gmotif_AddModuleScore'
)

jpeg('results/singlecell/GWASMotif04_FeaturePlot_LumA_signicant_all_AddModuleScore.jpg',width=10, height = 5, units = 'in', res = 150)
a<-FeaturePlot(
  object = subset(atac_subset_qc, idents=c('LHr','LS-1','LS-2','Bs')),
  features = "LumA_sinificant_gmotif_AddModuleScore1",
  min.cutoff = 'q5',
  max.cutoff = 'q95',
  pt.size = 0.1,
  label = TRUE
  )
b<-FeaturePlot(
  object = subset(atac_subset_qc, idents=c('T','B','Myl','Endo','Pc','Fb')),
  features = "LumA_sinificant_gmotif_AddModuleScore1",
  min.cutoff = 'q5',
  max.cutoff = 'q95',
  pt.size = 0.1,
  label = TRUE
)
a&NoLegend()|b
dev.off()

# colMean
# Compute average chromVAR score across motifs
motif_accessibility <- GetAssayData(atac_subset_qc, assay = "chromvar", layer = "data")

atac_subset_qc$LumA_avg_score_all <- Matrix::colMeans(motif_accessibility[unique(LumA$motif), ])
length(atac_subset_qc$LumA_avg_score_all)

jpeg('results/singlecell/GWASMotif05_FeaturePlot_LumA_all_colmeans.jpg',width=10, height = 5, units = 'in', res = 150)
a<-FeaturePlot(
  object = subset(atac_subset_qc, idents=c('LHr','LS-1','LS-2','Bs')),
  features = "LumA_avg_score_all",
  min.cutoff = 'q5',
  max.cutoff = 'q95',
  pt.size = 0.1,
  label = TRUE,
  cols = c("#eeeeee", "purple4", "darkblue"))
b<-FeaturePlot(
  object = subset(atac_subset_qc, idents=c('T','B','Myl','Endo','Pc','Fb')),
  features = "LumA_avg_score_all",
  min.cutoff = 'q5',
  max.cutoff = 'q95',
  pt.size = 0.1,
  label = TRUE,
  cols = c("#eeeeee", "purple4", "darkblue")
)
a&NoLegend()|b
dev.off()
#################################################################################################
#                              LumB significant GWAS motif                           #
#################################################################################################
motif_id_list_final<-fread('results/singlecell/jaspar2022_motifs_727.csv', header = TRUE)
LumB<-fread("results/tfgwas/tfgwas_results_LumB_ALL.csv", header = TRUE)
LumB<- left_join(motif_id_list_final, LumB, by=c('motif','motif.name')) %>% dplyr::filter(TF_p_value <= 0.05/727)
length(unique(LumB$motif))
enriched.motifs1$motif.name<-NULL
enriched_motifs_list_LumB <-left_join(LumB, enriched.motifs1, by='motif') %>% unique()

enriched_motifs_list_LumB$cluster <- unique(combined_peaks_gene_all$cluster)[1]
DefaultAssay(atac_subset_qc)<-'ATAC'
# Loop through each cluster name
for (cluster in cluster_names) {
  # select the unique peaks
  cluster_peak<-upset_matrix %>% dplyr::filter(SumEqualsOne == TRUE & .data[[cluster]] == 1) %>%
    dplyr::select(cluster)
  cluster_peak[,1]<-rownames(cluster_peak)
  names(cluster_peak) <- 'peaks'

  print(cluster)
  nrow(cluster_peak)

  feature1 <-combined_peaks_gene_all$peaks[combined_peaks_gene_all$cluster == cluster]
  feature1 <- feature1[feature1 %in% cluster_peak$peaks]

  # Find motifs for the current cluster
  enriched_motifs <- FindMotifs(
    object = atac_subset_qc,
    features = feature1)

  # Filter the enriched motifs to include only those in significant GWAS
  enriched_motifs$motif.name<-NULL
  enriched_motifs_filtered <- left_join(LumB, enriched_motifs, by='motif') %>% unique()
  enriched_motifs_filtered$cluster <- cluster

  # Store the filtered results in the list with the cluster name as the key
  enriched_motifs_list_LumB <- rbind(enriched_motifs_list_LumB,enriched_motifs_filtered)

  # Print a message to indicate progress
  cat("Processed cluster:", cluster, "\n")
}

enriched_motifs_list_LumB$pvalue <- ifelse(enriched_motifs_list_LumB$pvalue==0,1e-300,
                                           ifelse(enriched_motifs_list_LumB$pvalue==1,0.99999999999999999999999,enriched_motifs_list_LumB$pvalue))

# Save the motif_cor_top5DAP to a file
write.csv(enriched_motifs_list_LumB, file = "results/singlecell/enriched_motifs_list_uniqDAP_LumBsignificant.csv", row.names = FALSE)

# create a matrix
# Get all unique motif names across all clusters
all_motifs <- unique(enriched_motifs_list_LumB$motif)

# Get all cluster names
all_clusters <- unique(enriched_motifs_list_LumB$cluster)

# Initialize an empty matrix
motif_matrix <- matrix(
  0,
  ncol = length(all_clusters),
  nrow = length(all_motifs),
  dimnames = list( all_motifs,all_clusters)
)

# Fill the matrix with p-values for enriched motifs in each cluster
for (cluster1 in all_clusters) {
  # Get the enriched motifs and their p-values for the current cluster
  enriched_data <- enriched_motifs_list_LumB[enriched_motifs_list_LumB$cluster==cluster1]

  # Add p-values to the matrix
  for (i in 1:nrow(enriched_data)) {
    motif <- enriched_data$motif[i]
    neglog10pvalue <- -log10(enriched_data$pvalue[i])  # Replace "pvalue" with the actual column name for p-values
    motif_matrix[motif,cluster1] <-neglog10pvalue
  }
}

# Convert the matrix to a data frame for better readability
motif_df<- as.data.frame(motif_matrix)
motif_df$motif<-rownames(motif_df)
motif_df <- left_join(motif_df, select(motif_id_list_final, motif, motif.name), by='motif')
motif_df$motif.name <- make.unique(motif_df$motif.name)
rownames(motif_df) <- motif_df$motif.name
motif_df$motif<-NULL
motif_df$motif.name<-NULL
motif_matrix <- as.matrix(motif_df)

motif_df<-dplyr::select(motif_df,'LHr','LS-1','LS-2','Bs','T','B','Myl','Endo','Pc','Fb')
library(pheatmap)
library(ComplexHeatmap)
library(circlize)

motif_label <- motif_df[apply(motif_df, 1, function(x) any(x >= 300)), ]

motif_df<-motif_df%>%arrange(desc(LHr),desc(`LS-1`))
TF_list<-rownames(motif_df)
cell_list<-colnames(motif_df)
motif_heatmap<-SingleRasterMap(
  data = t(motif_df),
  raster = TRUE,
  cell.order = cell_list,
  feature.order = rev(TF_list),
  colors = colorRampPalette(c('lightyellow',"purple4","darkblue"))(300),   # Custom color mapping
  disp.min = 0,      # Minimum display value
  disp.max = 300,    # Maximum display value
  limits = c(0, 300) # Set explicit limits
)
num_cols <- ncol(motif_df)
line_positions <- seq(1.5, num_cols - 0.5, by = 1)  # Positions between columns
motif_heatmap<-motif_heatmap +
  geom_vline(xintercept = line_positions, color = "purple4", size = 0.25)+
  scale_y_discrete(
    breaks = rownames(motif_label),                 # Only label significant genes
    labels = rownames(motif_label)) +
  labs(fill = "-log10(pval)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1,size = 12),
        axis.text.y = element_text(angle = 0,size = 5)  # Adjust y-axis font size
  )
motif_heatmap
ggsave("results/singlecell/Motif_enricheded_heatmap_uniDAP_LumB_signMotif.jpg",
       motif_heatmap, width = 5.5, height = 10, dpi = 200)

#----------------- Set chromvar as default assay---------------------------
DefaultAssay(atac_subset_qc) <- "chromvar"

LumB_significant_motif_ls<-list(unique(LumB$motif))
atac_subset_qc <- AddModuleScore(
  object = atac_subset_qc,
  features = LumB_significant_motif_ls,
  ctrl = 5,
  name = 'LumB_sinificant_gmotif_AddModuleScore'
)

jpeg('results/singlecell/GWASMotif06_FeaturePlot_LumB_signicant_all_AddModuleScore.jpg',width=10, height = 5, units = 'in', res = 150)
a<-FeaturePlot(
  object = subset(atac_subset_qc, idents=c('LHr','LS-1','LS-2','Bs')),
  features = "LumB_sinificant_gmotif_AddModuleScore1",
  min.cutoff = 'q5',
  max.cutoff = 'q95',
  pt.size = 0.1,
  label = TRUE
  )
b<-FeaturePlot(
  object = subset(atac_subset_qc, idents=c('T','B','Myl','Endo','Pc','Fb')),
  features = "LumB_sinificant_gmotif_AddModuleScore1",
  min.cutoff = 'q5',
  max.cutoff = 'q95',
  pt.size = 0.1,
  label = TRUE
)
a&NoLegend()|b
dev.off()

# colMean
# Compute average chromVAR score across motifs
motif_accessibility <- GetAssayData(atac_subset_qc, assay = "chromvar", layer = "data")

atac_subset_qc$LumB_avg_score_all <- Matrix::colMeans(motif_accessibility[unique(LumB$motif), ])
length(atac_subset_qc$LumB_avg_score_all)

jpeg('results/singlecell/GWASMotif07_FeaturePlot_LumB_all_colmeans.jpg',width=10, height = 5, units = 'in', res = 150)
a<-FeaturePlot(
  object = subset(atac_subset_qc, idents=c('LHr','LS-1','LS-2','Bs')),
  features = "LumB_avg_score_all",
  min.cutoff = 'q5',
  max.cutoff = 'q95',
  pt.size = 0.1,
  label = TRUE,
  cols = c("#eeeeee", "purple4", "darkblue"))
b<-FeaturePlot(
  object = subset(atac_subset_qc, idents=c('T','B','Myl','Endo','Pc','Fb')),
  features = "LumB_avg_score_all",
  min.cutoff = 'q5',
  max.cutoff = 'q95',
  pt.size = 0.1,
  label = TRUE,
  cols = c("#eeeeee", "purple4", "darkblue")
)
a&NoLegend()|b
dev.off()

#################################################################################################
#                              LumB-HER2Neg significant GWAS motif                           #
#################################################################################################
motif_id_list_final<-fread('results/singlecell/jaspar2022_motifs_727.csv', header = TRUE)
LumB_HER2Neg<-fread("results/tfgwas/tfgwas_results_LumB_HER2Neg_ALL.csv", header = TRUE)
LumB_HER2Neg<- left_join(motif_id_list_final, LumB_HER2Neg, by=c('motif','motif.name')) %>% dplyr::filter(TF_p_value <= 0.05/727)
length(unique(LumB_HER2Neg$motif))
enriched.motifs1$motif.name<-NULL
enriched_motifs_list_LumB_HER2Neg <-left_join(LumB_HER2Neg, enriched.motifs1, by='motif') %>% unique()

enriched_motifs_list_LumB_HER2Neg$cluster <- unique(combined_peaks_gene_all$cluster)[1]
DefaultAssay(atac_subset_qc)<-'ATAC'
# Loop through each cluster name
for (cluster in cluster_names) {
  # select the unique peaks
  cluster_peak<-upset_matrix %>% dplyr::filter(SumEqualsOne == TRUE & .data[[cluster]] == 1) %>%
    dplyr::select(cluster)
  cluster_peak[,1]<-rownames(cluster_peak)
  names(cluster_peak) <- 'peaks'

  print(cluster)
  nrow(cluster_peak)

  feature1 <-combined_peaks_gene_all$peaks[combined_peaks_gene_all$cluster == cluster]
  feature1 <- feature1[feature1 %in% cluster_peak$peaks]

  # Find motifs for the current cluster
  enriched_motifs <- FindMotifs(
    object = atac_subset_qc,
    features = feature1)

  # Filter the enriched motifs to include only those in significant GWAS
  enriched_motifs$motif.name<-NULL
  enriched_motifs_filtered <- left_join(LumB_HER2Neg, enriched_motifs, by='motif') %>% unique()
  enriched_motifs_filtered$cluster <- cluster

  # Store the filtered results in the list with the cluster name as the key
  enriched_motifs_list_LumB_HER2Neg <- rbind(enriched_motifs_list_LumB_HER2Neg,enriched_motifs_filtered)

  # Print a message to indicate progress
  cat("Processed cluster:", cluster, "\n")
}

enriched_motifs_list_LumB_HER2Neg$pvalue <- ifelse(enriched_motifs_list_LumB_HER2Neg$pvalue==0,1e-300,
                                           ifelse(enriched_motifs_list_LumB_HER2Neg$pvalue==1,0.99999999999999999999999,enriched_motifs_list_LumB_HER2Neg$pvalue))

# Save the motif_cor_top5DAP to a file
write.csv(enriched_motifs_list_LumB_HER2Neg, file = "results/singlecell/enriched_motifs_list_uniqDAP_LumB_HER2Negsignificant.csv", row.names = FALSE)

# create a matrix
# Get all unique motif names across all clusters
all_motifs <- unique(enriched_motifs_list_LumB_HER2Neg$motif)

# Get all cluster names
all_clusters <- unique(enriched_motifs_list_LumB_HER2Neg$cluster)

# Initialize an empty matrix
motif_matrix <- matrix(
  0,
  ncol = length(all_clusters),
  nrow = length(all_motifs),
  dimnames = list( all_motifs,all_clusters)
)

# Fill the matrix with p-values for enriched motifs in each cluster
for (cluster1 in all_clusters) {
  # Get the enriched motifs and their p-values for the current cluster
  enriched_data <- enriched_motifs_list_LumB_HER2Neg[enriched_motifs_list_LumB_HER2Neg$cluster==cluster1]

  # Add p-values to the matrix
  for (i in 1:nrow(enriched_data)) {
    motif <- enriched_data$motif[i]
    neglog10pvalue <- -log10(enriched_data$pvalue[i])  # Replace "pvalue" with the actual column name for p-values
    motif_matrix[motif,cluster1] <-neglog10pvalue
  }
}

# Convert the matrix to a data frame for better readability
motif_df<- as.data.frame(motif_matrix)
motif_df$motif<-rownames(motif_df)
motif_df <- left_join(motif_df, select(motif_id_list_final, motif, motif.name), by='motif')
motif_df$motif.name <- make.unique(motif_df$motif.name)
rownames(motif_df) <- motif_df$motif.name
motif_df$motif<-NULL
motif_df$motif.name<-NULL
motif_matrix <- as.matrix(motif_df)

motif_df<-dplyr::select(motif_df,'LHr','LS-1','LS-2','Bs','T','B','Myl','Endo','Pc','Fb')
library(pheatmap)
library(ComplexHeatmap)
library(circlize)

motif_label <- motif_df[apply(motif_df, 1, function(x) any(x >= 300)), ]

motif_df<-motif_df%>%arrange(desc(LHr),desc(`LS-1`))
TF_list<-rownames(motif_df)
cell_list<-colnames(motif_df)
motif_heatmap<-SingleRasterMap(
  data = t(motif_df),
  raster = TRUE,
  cell.order = cell_list,
  feature.order = rev(TF_list),
  colors = colorRampPalette(c('lightyellow',"purple4","darkblue"))(300),   # Custom color mapping
  disp.min = 0,      # Minimum display value
  disp.max = 300,    # Maximum display value
  limits = c(0, 300) # Set explicit limits
)
num_cols <- ncol(motif_df)
line_positions <- seq(1.5, num_cols - 0.5, by = 1)  # Positions between columns
motif_heatmap<-motif_heatmap +
  geom_vline(xintercept = line_positions, color = "purple4", size = 0.25)+
  scale_y_discrete(
    breaks = rownames(motif_label),                 # Only label significant genes
    labels = rownames(motif_label)) +
  labs(fill = "-log10(pval)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1,size = 12),
        axis.text.y = element_text(angle = 0,size = 5)  # Adjust y-axis font size
  )
motif_heatmap
ggsave("results/singlecell/Motif_enricheded_heatmap_uniDAP_LumB_HER2Neg_signMotif.jpg",
       motif_heatmap, width = 5.5, height = 10, dpi = 200)

#----------------- Set chromvar as default assay---------------------------
DefaultAssay(atac_subset_qc) <- "chromvar"

LumB_HER2Neg_significant_motif_ls<-list(unique(LumB_HER2Neg$motif))
atac_subset_qc <- AddModuleScore(
  object = atac_subset_qc,
  features = LumB_HER2Neg_significant_motif_ls,
  ctrl = 5,
  name = 'LumB_HER2Neg_sinificant_gmotif_AddModuleScore'
)

jpeg('results/singlecell/GWASMotif08_FeaturePlot_LumB_HER2Neg_signicant_all_AddModuleScore.jpg',width=10, height = 5, units = 'in', res = 150)
a<-FeaturePlot(
  object = subset(atac_subset_qc, idents=c('LHr','LS-1','LS-2','Bs')),
  features = "LumB_HER2Neg_sinificant_gmotif_AddModuleScore1",
  min.cutoff = 'q5',
  max.cutoff = 'q95',
  pt.size = 0.1,
  label = TRUE
)
b<-FeaturePlot(
  object = subset(atac_subset_qc, idents=c('T','B','Myl','Endo','Pc','Fb')),
  features = "LumB_HER2Neg_sinificant_gmotif_AddModuleScore1",
  min.cutoff = 'q5',
  max.cutoff = 'q95',
  pt.size = 0.1,
  label = TRUE
)
a&NoLegend()|b
dev.off()

#################################################################################################
  #                            HER2 significant GWAS motif                           #
#################################################################################################
motif_id_list_final<-fread('results/singlecell/jaspar2022_motifs_727.csv', header = TRUE)
HER2<-fread("results/tfgwas/tfgwas_results_HER2_ALL.csv", header = TRUE)
HER2<- left_join(select(motif_id_list_final,motif,motif.name), HER2, by=c('motif','motif.name')) %>% dplyr::filter(TF_p_value <= 0.05/727)
length(unique(HER2$motif))
enriched.motifs1$motif.name<-NULL
enriched_motifs_list_HER2 <-left_join(HER2, enriched.motifs1, by='motif') %>% unique()

enriched_motifs_list_HER2$cluster <- unique(combined_peaks_gene_all$cluster)[1]

DefaultAssay(atac_subset_qc)<-'ATAC'
# Loop through each cluster name
for (cluster in cluster_names) {
  # select the unique peaks
  cluster_peak<-upset_matrix %>% dplyr::filter(SumEqualsOne == TRUE & .data[[cluster]] == 1) %>%
    dplyr::select(cluster)
  cluster_peak[,1]<-rownames(cluster_peak)
  names(cluster_peak) <- 'peaks'

  print(cluster)
  nrow(cluster_peak)

  feature1 <-combined_peaks_gene_all$peaks[combined_peaks_gene_all$cluster == cluster]
  feature1 <- feature1[feature1 %in% cluster_peak$peaks]

  # Find motifs for the current cluster
  enriched_motifs <- FindMotifs(
    object = atac_subset_qc,
    features = feature1)

  # Filter the enriched motifs to include only those in significant GWAS
  enriched_motifs$motif.name <-NULL
  enriched_motifs_filtered <- left_join(HER2, enriched_motifs, by='motif') %>% unique()
  enriched_motifs_filtered$cluster <- cluster

  # Store the filtered results in the list with the cluster name as the key
  enriched_motifs_list_HER2 <- rbind(enriched_motifs_list_HER2,enriched_motifs_filtered)

  # Print a message to indicate progress
  cat("Processed cluster:", cluster, "\n")
}

enriched_motifs_list_HER2$pvalue <- ifelse(enriched_motifs_list_HER2$pvalue==0,1e-300,
                                           ifelse(enriched_motifs_list_HER2$pvalue==1,0.99999999999999999999999,enriched_motifs_list_HER2$pvalue))

# Save the motif_cor_top5DAP to a file
write.csv(enriched_motifs_list_HER2, file = "results/singlecell/enriched_motifs_list_uniqDAP_HER2significant.csv", row.names = FALSE)

# create a matrix
# Get all unique motif names across all clusters
all_motifs <- unique(enriched_motifs_list_HER2$motif)

# Get all cluster names
all_clusters <- unique(enriched_motifs_list_HER2$cluster)

# Initialize an empty matrix
motif_matrix <- matrix(
  0,
  ncol = length(all_clusters),
  nrow = length(all_motifs),
  dimnames = list( all_motifs,all_clusters)
)

# Fill the matrix with p-values for enriched motifs in each cluster
for (cluster1 in all_clusters) {
  # Get the enriched motifs and their p-values for the current cluster
  enriched_data <- enriched_motifs_list_HER2[enriched_motifs_list_HER2$cluster==cluster1]

  # Add p-values to the matrix
  for (i in 1:nrow(enriched_data)) {
    motif <- enriched_data$motif[i]
    neglog10pvalue <- -log10(enriched_data$pvalue[i])  # Replace "pvalue" with the actual column name for p-values
    motif_matrix[motif,cluster1] <-neglog10pvalue
  }
}

# Convert the matrix to a data frame for better readability
motif_df<- as.data.frame(motif_matrix)
motif_df$motif<-rownames(motif_df)
motif_df <- left_join(motif_df, select(motif_id_list_final, motif, motif.name), by='motif')
motif_df$motif.name <- make.unique(motif_df$motif.name)
rownames(motif_df) <- motif_df$motif.name
motif_df$motif<-NULL
motif_df$motif.name<-NULL
motif_matrix <- as.matrix(motif_df)

motif_df<-dplyr::select(motif_df,'LHr','LS-1','LS-2','Bs','T','B','Myl','Endo','Pc','Fb')
library(pheatmap)
library(ComplexHeatmap)
library(circlize)

motif_label <- motif_df[apply(motif_df, 1, function(x) any(x >= 300)), ]

motif_df<-motif_df%>%arrange(desc(LHr),desc(`LS-1`))
TF_list<-rownames(motif_df)
cell_list<-colnames(motif_df)
motif_heatmap<-SingleRasterMap(
  data = t(motif_df),
  raster = TRUE,
  cell.order = cell_list,
  feature.order = rev(TF_list),
  colors = colorRampPalette(c('lightyellow',"purple4","darkblue"))(300),   # Custom color mapping
  disp.min = 0,      # Minimum display value
  disp.max = 300,    # Maximum display value
  limits = c(0, 300) # Set explicit limits
)
num_cols <- ncol(motif_df)
line_positions <- seq(1.5, num_cols - 0.5, by = 1)  # Positions between columns
motif_heatmap<-motif_heatmap +
  geom_vline(xintercept = line_positions, color = "purple4", size = 0.25)+
  scale_y_discrete(
    breaks = rownames(motif_label),                 # Only label significant genes
    labels = rownames(motif_label)) +
  labs(fill = "-log10(pval)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1,size = 12),
        axis.text.y = element_text(angle = 0,size = 5)  # Adjust y-axis font size
  )
motif_heatmap
ggsave("results/singlecell/Motif_enricheded_heatmap_uniDAP_HER2_signMotif.jpg",
       motif_heatmap, width = 5.5, height = 5, dpi = 200)

#----------------- Set chromvar as default assay---------------------------
DefaultAssay(atac_subset_qc) <- "chromvar"

HER2_significant_motif_ls<-list(unique(HER2$motif))
atac_subset_qc <- AddModuleScore(
  object = atac_subset_qc,
  features = HER2_significant_motif_ls,
  ctrl = 5,
  name = 'HER2_sinificant_gmotif_AddModuleScore'
)

jpeg('results/singlecell/GWASMotif08_FeaturePlot_HER2_signicant_all_AddModuleScore.jpg',width=10, height = 5, units = 'in', res = 150)
a<-FeaturePlot(
  object = subset(atac_subset_qc, idents=c('LHr','LS-1','LS-2','Bs')),
  features = "HER2_sinificant_gmotif_AddModuleScore1",
  min.cutoff = 'q5',
  max.cutoff = 'q95',
  pt.size = 0.1,
  label = TRUE
  )
b<-FeaturePlot(
  object = subset(atac_subset_qc, idents=c('T','B','Myl','Endo','Pc','Fb')),
  features = "HER2_sinificant_gmotif_AddModuleScore1",
  min.cutoff = 'q5',
  max.cutoff = 'q95',
  pt.size = 0.1,
  label = TRUE
)
a&NoLegend()|b
dev.off()

# colMean
# Compute average chromVAR score across motifs
motif_accessibility <- GetAssayData(atac_subset_qc, assay = "chromvar", layer = "data")

atac_subset_qc$HER2_avg_score_all <- Matrix::colMeans(motif_accessibility[unique(HER2$motif), ])
length(atac_subset_qc$HER2_avg_score_all)

jpeg('results/singlecell/GWASMotif09_FeaturePlot_HER2_all_colmeans.jpg',width=10, height = 5, units = 'in', res = 150)
a<-FeaturePlot(
  object = subset(atac_subset_qc, idents=c('LHr','LS-1','LS-2','Bs')),
  features = "HER2_avg_score_all",
  min.cutoff = 'q5',
  max.cutoff = 'q95',
  pt.size = 0.1,
  label = TRUE,
  cols = c("#eeeeee", "purple4", "darkblue"))
b<-FeaturePlot(
  object = subset(atac_subset_qc, idents=c('T','B','Myl','Endo','Pc','Fb')),
  features = "HER2_avg_score_all",
  min.cutoff = 'q5',
  max.cutoff = 'q95',
  pt.size = 0.1,
  label = TRUE,
  cols = c("#eeeeee", "purple4", "darkblue")
)
a&NoLegend()|b
dev.off()
