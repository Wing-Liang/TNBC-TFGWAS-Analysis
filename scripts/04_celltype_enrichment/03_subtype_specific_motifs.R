# ==============================================================================
# 03_subtype_specific_motifs.R
# Subtype-specific motifs: module scores, Kruskal-Wallis / Dunn tests and enrichment heatmaps.
# Produces: Fig 3E-F, Table S5
# Inputs  : results/singlecell/final_cell_barcodes.csv
#           results/singlecell/atac_motif.rds
#           results/singlecell/peak_annotation.csv
#           results/singlecell/combined_peaks_gene_all.csv
#           results/singlecell/jaspar2022_motifs_727.csv
#           results/tfgwas/tfgwas_results_TNBC_ALL.csv
#           results/tfgwas/tfgwas_results_LumA_ALL.csv
#           ...
# Outputs : results/tfgwas/tfgwas_significant_TNBC_unique.csv
#           results/singlecell/enriched_motifs_list_uniqDAP_TNBCUniqSignificant.csv
#           results/singlecell/Motif_enricheded_heatmap_uniDAP_TNBC_unique_signMotif_52motif.jpg
#           results/tfgwas/heatmap_enrichedmotif_uniqMotif_TNBC_orginal_52motifs.pdf
#           results/tfgwas/heatmap_enrichedmotif_uniqMotif_TNBC_57motifs.pdf
#           results/tfgwas/GSEAplot_TNBC_uniquemotifs.jpg
#           results/tfgwas/GWASMotif01.TNBC_motif_count.jpg
#           ...
# Run     : Rscript scripts/04_celltype_enrichment/03_subtype_specific_motifs.R   (from the repository root)
# ==============================================================================
## ---- paths (set the environment variables; see README) ----
TNBC_ROOT <- Sys.getenv("TNBC_ROOT", ".")

# Load required packages
library(lme4)      # For generalized mixed models
library(dplyr)     # For data manipulation
library(ggplot2)   # For visualization (optional)
library(Seurat)        # For scATAC-seq object handling
library(Signac)        # For scATAC-seq utilities
library(dplyr)         # For data manipulation
library(GenomicRanges) # For genomic overlaps
library(motifmatchr)   # For motif matching
library(TFBSTools)     # For motif handling
library(data.table)
library(ggpubr)
library(fmesher)
library(INLA)

setwd(TNBC_ROOT)
#############################################
#         heatmap of motif enrichment       #
##############################################
Final_cell_barcode<-fread('results/singlecell/final_cell_barcodes.csv',header=TRUE)

#ATAC
atac_subset_qc<-readRDS('results/singlecell/atac_motif.rds')
atac_subset_qc <- subset(atac_subset_qc , cells = Final_cell_barcode$x)
DefaultAssay(atac_subset_qc)<-'RNA'
new.cluster.ids <- c('LASP-1','Fb','LHS','LASP-2','B-Myo',
                     'Endo','T','Pc','Myl','B')
names(new.cluster.ids) <- levels(atac_subset_qc)
atac_subset_qc <- RenameIdents(atac_subset_qc, new.cluster.ids)
atac_subset_qc$celltype_RNA <- Idents(atac_subset_qc)

DefaultAssay(atac_subset_qc)<-'ATAC'
names(new.cluster.ids) <- levels(atac_subset_qc)
atac_subset_qc <- RenameIdents(atac_subset_qc, new.cluster.ids)
atac_subset_qc$celltype_ATAC <- Idents(atac_subset_qc)

# re-order the levels sequence
# cell-type display order
fig5D_levels <- c('LHS','LASP-1','LASP-2','B-Myo','T','B','Myl','Endo','Pc','Fb')
# Reorder the Idents based on sorted levels
Idents(atac_subset_qc) <- factor(Idents(atac_subset_qc), levels = fig5D_levels)
atac_subset_qc$celltype_ATAC <- Idents(atac_subset_qc)

# enriched motif
# DAP unique in the cell
combined_peaks_gene_all <- fread("results/singlecell/peak_annotation.csv", header = TRUE)
combined_peaks_gene_all$cluster_name <- gsub('LS-1','LASP-1',combined_peaks_gene_all$cluster_name )
combined_peaks_gene_all$cluster_name <- gsub('LASP-2','LASP-2',combined_peaks_gene_all$cluster_name )
combined_peaks_gene_all$cluster_name <- gsub('LHr','LHS',combined_peaks_gene_all$cluster_name )
combined_peaks_gene_all$cluster_name <- gsub('Bs','B-Myo',combined_peaks_gene_all$cluster_name )
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
names(Endo) <-'peaks'
nrow(Endo)

#top5 DAPs
combined_peaks_gene_all<-fread("results/singlecell/combined_peaks_gene_all.csv", header = TRUE)
combined_peaks_gene_all$cluster <- gsub('LS-1','LASP-1',combined_peaks_gene_all$cluster )
combined_peaks_gene_all$cluster <- gsub('LASP-2','LASP-2',combined_peaks_gene_all$cluster )
combined_peaks_gene_all$cluster <- gsub('LHr','LHS',combined_peaks_gene_all$cluster )
combined_peaks_gene_all$cluster <- gsub('Bs','B-Myo',combined_peaks_gene_all$cluster )

names(combined_peaks_gene_all)[1]<-'peaks'

# Get the unique cluster names
cluster_names <- unique(combined_peaks_gene_all$cluster)[2:10]
# "LASP-1" "Pc"     "Fb"     "B-Myo"  "LHS"    "T"      "B"      "LASP-2" "Myl"

# the first cluster - Endo
feature1<-combined_peaks_gene_all$peaks[combined_peaks_gene_all$cluster==unique(combined_peaks_gene_all$cluster)[1]] # Endo
feature1 <- feature1[feature1 %in% Endo$peaks]

# test enrichment
DefaultAssay(atac_subset_qc)<-'ATAC'
enriched.motifs1 <- FindMotifs(
  object = atac_subset_qc,
  features = feature1)

#################################################################################################
#                          Unique Motif in TNBC significant GWAS motif                         #
#################################################################################################
motif_id_list_final<-fread('results/singlecell/jaspar2022_motifs_727.csv', header = TRUE)
TNBC<-fread("results/tfgwas/tfgwas_results_TNBC_ALL.csv", header = TRUE)
TNBC <- TNBC %>% select(motif, motif.name, TF_p_value)
TNBC<- left_join(motif_id_list_final, TNBC, by=c('motif','motif.name')) %>% dplyr::filter(TF_p_value < 0.05/727)
length(unique(TNBC$motif))

LumA<-fread("results/tfgwas/tfgwas_results_LumA_ALL.csv", header = TRUE)
LumB<-fread("results/tfgwas/tfgwas_results_LumB_ALL.csv", header = TRUE)
LumB_HER2Neg<-fread("results/tfgwas/tfgwas_results_LumB_HER2Neg_ALL.csv", header = TRUE)
HER2<-fread("results/tfgwas/tfgwas_results_HER2_ALL.csv", header = TRUE)

all<-rbind(LumA, LumB)
all<-rbind(all, LumB_HER2Neg)
all<-rbind(all, HER2)
all<-all %>% dplyr::filter(TF_p_value < 0.05/727)%>% unique()
#                               TF_p_value_LumB>=5e-4, TF_p_value_LumB_HER2Neg>=5e-4, TF_p_value_HER2>=5e-4)
other_sig_motif<-unique(all$motif)
length(other_sig_motif)
#475 sig

TNBC_unique <- TNBC %>% dplyr::filter(!motif %in% other_sig_motif)  %>% unique()
length(unique(TNBC_unique$motif))
fwrite(TNBC_unique,"results/tfgwas/tfgwas_significant_TNBC_unique.csv", row.names = FALSE)

TNBC_unique<-fread("results/tfgwas/tfgwas_significant_TNBC_unique.csv",header = TRUE)
# # TNBC_unique_common<-TNBC_unique[!motif %in% motif_not_in_motifbreakR$motif,]
# TNBC_unique[motif %in% motif_not_in_motifbreakR$motif,motif.name]

enriched.motifs1$motif.name<-NULL
enriched_motifs_list_TNBC_unique <-left_join(TNBC_unique, enriched.motifs1,
                                             by='motif') %>% unique()

enriched_motifs_list_TNBC_unique$cluster <- unique(combined_peaks_gene_all$cluster)[1]

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
  enriched_motifs_filtered <- left_join(TNBC_unique, enriched_motifs, by='motif') %>% unique()
  enriched_motifs_filtered$cluster <- cluster

  # Store the filtered results in the list with the cluster name as the key
  enriched_motifs_list_TNBC_unique <- rbind(enriched_motifs_list_TNBC_unique,enriched_motifs_filtered)

  # Print a message to indicate progress
  cat("Processed cluster:", cluster, "\n")
}

enriched_motifs_list_TNBC_unique$pvalue <- ifelse(enriched_motifs_list_TNBC_unique$pvalue==0,1e-300,
                                                  ifelse(enriched_motifs_list_TNBC_unique$pvalue==1,0.99999999999999999999999,enriched_motifs_list_TNBC_unique$pvalue))

# Save the motif_cor_top5DAP to a file
write.csv(enriched_motifs_list_TNBC_unique, file = "results/singlecell/enriched_motifs_list_uniqDAP_TNBCUniqSignificant.csv", row.names = FALSE)

enriched_motifs_list_TNBC_unique<-fread("results/singlecell/enriched_motifs_list_uniqDAP_TNBCUniqSignificant.csv",header = TRUE)
# create a matrix
# Get all unique motif names across all clusters
all_motifs <- unique(enriched_motifs_list_TNBC_unique$motif.name)

# Get all cluster names
all_clusters <- unique(enriched_motifs_list_TNBC_unique$cluster)

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
  enriched_data <- enriched_motifs_list_TNBC_unique[enriched_motifs_list_TNBC_unique$cluster==cluster1]

  # Add p-values to the matrix
  for (i in 1:nrow(enriched_data)) {
    motif <- enriched_data$motif.name[i]
    neglog10pvalue <- -log10(enriched_data$pvalue[i])  # Replace "pvalue" with the actual column name for p-values
    motif_matrix[motif,cluster1] <-neglog10pvalue
  }
}

# Convert the matrix to a data frame for better readability
motif_df<- as.data.frame(motif_matrix)

motif_df<-dplyr::select(motif_df,'LS-1','LASP-2','LHr','Bs','Fb','Endo','T','B','Myl','Pc')
library(pheatmap)
library(ComplexHeatmap)
library(circlize)

motif_label_Fb <- motif_df[motif_df$Fb>= -log10(0.05),]  #17
rownames(motif_label_Fb)

motif_df<-motif_df%>%arrange(desc(LHr),desc(`LS-1`),desc(Fb))
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
  #   breaks = rownames(motif_label),                 # Only label significant genes
  labs(fill = "-log10(pval)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1,size = 12),
        axis.text.y = element_text(angle = 0,size = 5)  # Adjust y-axis font size
  )
motif_heatmap
ggsave("results/singlecell/Motif_enricheded_heatmap_uniDAP_TNBC_unique_signMotif_52motif.jpg",
       motif_heatmap, width = 4, height = 5, dpi = 200)

# ---------------------heatmap ---------------#
motif_id_name<-fread('results/singlecell/jaspar2022_motifs_727.csv', header = TRUE)

# for the motif score heatmap
motif_matrix <- as.matrix(motif_df)

scaled_mat <- t(scale(t(motif_matrix)))
row_dist <- dist(scaled_mat)
row_clust <- hclust(row_dist)
col_clust = hclust(dist(scaled_mat))

library(ComplexHeatmap)
library(dendextend)
set.seed(2025)
heatmap_obj_orig <-Heatmap(motif_matrix,
                           name = "-log10 pval",          # Title of legend
                           row_names_gp = gpar(fontsize = 7),
                           row_names_side = "left",
                           column_names_centered = TRUE,
                           width = unit(50, "mm"),
                           col = circlize::colorRamp2(c(0,80), c("#483d8b","red1")),
                           rect_gp = gpar(col = "grey99"), # Set cell border color (instead of border_color)
                           show_row_dend = FALSE,         # Hide row dendrogram to avoid km labels
                           show_column_dend = FALSE,
                           # column_split = 4,                 # Ensure 3 clusters without labels
                           row_title = NULL ,
                           column_title = NULL
)
heatmap_obj_orig

pdf("results/tfgwas/heatmap_enrichedmotif_uniqMotif_TNBC_orginal_52motifs.pdf", width = 5, height =7)
draw(heatmap_obj_orig)
dev.off()

# heatmap with specific motif sequence
motif_df<- as.data.frame(motif_matrix)

colSums(motif_matrix >= -log10(0.05 / 727))
# Endo LS-1   Pc   Fb   Bs  LHr    T    B LASP-2  Myl
colSums(motif_matrix >= 80)
# Endo LS-1   Pc   Fb   Bs  LHr    T    B LASP-2  Myl

motif_df <- motif_df[order(-motif_df$Fb),]

rownames(motif_df)
motif_list <- c("NFIC","MYOD1","FIGLA","ELF2","TBX4","TBX5","TBX6","ZBTB6","CDX1","CDX4",
                "TBX1","MGA","TBX15","EOMES","EBF1","NKX2-4","HOXD12::ELK1","FOSL1::JUNB","FOSL2::JUND","FOSB::JUNB",
                "FOSL2::JUNB","MAF::NFE2","RUNX3",
                "BHLHA15","TFAP4","NEUROD1","CEBPD", "CEBPA","MYF6", "CEBPG","HLF","NHLH2","TEF","ATF2","GLI2",
                "TCF21","ERG","FLI1","ELK4" ,"ERF::HOXB13","ELK1::HOXA1",
                "NFE2L2","HOXA6","HOXA4" ,"MZF1","HOXB6",
                "IRF9","TCFL5","MSANTD3","FEV" ,"ELK1::HOXB13","ZBTB7B",  "E2F2","MEF2B","PPARA::RXRA","E2F4","ZNF85")%>%unique()

motif_df<- as.data.frame(motif_matrix)
head(motif_df)
# Reorder motif_df rows based on motif_list
motif_df <- motif_df[match(motif_list, rownames(motif_df)), ]
motif_df <- motif_df %>% select('LS-1','LASP-2','LHr','Bs','Fb','Endo','T','B','Myl','Pc')
motif_matrix <- as.matrix(motif_df)

library(ComplexHeatmap)
library(circlize)

# Create the heatmap
pdf("results/tfgwas/heatmap_enrichedmotif_uniqMotif_TNBC_57motifs.pdf", width = 5, height =7)
Heatmap(motif_matrix,
        col =  circlize::colorRamp2(c(0, 150), c("#483d8b", "red1")),
        row_names_gp = gpar(fontsize = 8),
        cluster_rows = FALSE,  # Preserve the order from motif_list
        cluster_columns = FALSE,  # Cluster columns (or set to FALSE if you want to preserve column order)
        show_row_names = TRUE,  # Show motif names
        show_column_names = TRUE,  # Show column names
        row_names_side = "left",  # Place row names on the left
        heatmap_legend_param = list(title = "-log10 pval"),
        rect_gp = gpar(col = "grey99"),)  # Customize legend title

dev.off()

############################ GSEA bubble plots ###########################################
# Load libraries
library(tidyverse)
library(ggplot2)
library(paletteer)

dat_TNBC<-fread("results/singlecell/enriched_motifs_list_uniqDAP_TNBCUniqSignificant.csv",header = TRUE)
dat_TNBC$cluster <- gsub('LS-1','LASP-1',dat_TNBC$cluster )
dat_TNBC$cluster <- gsub('LASP-2','LASP-2',dat_TNBC$cluster )
dat_TNBC$cluster <- gsub('LHr','LHS',dat_TNBC$cluster )
dat_TNBC$cluster <- gsub('Bs','B-Myo',dat_TNBC$cluster )

dat_TNBC$cluster <- factor(dat_TNBC$cluster, levels = c('LASP-1','LASP-2','LHS','B-Myo','Fb','Endo','T','B','Myl','Pc'))
dat_TNBC$avg.exp.scaled <- ifelse(is.na(dat_TNBC$avg.exp.scaled), 0, dat_TNBC$avg.exp.scaled)

dat_TNBC <- dat_TNBC %>%
  mutate(family = if_else(
    family == 'Thyroid hormone receptor-related factors (NR1);RXR-related receptors (NR2)',
    'Thyroid hormone receptor/RXR-related factor',
    family
  ))
# Thyroid hormone receptor-related factors (NR1);RXR-related receptors (NR2)

jpeg('results/tfgwas/GSEAplot_TNBC_uniquemotifs.jpg', width = 15, height = 18, units = 'in', res = 150)
ggplot1 <- ggplot(dat_TNBC, aes(x = fold.enrichment, y = reorder(motif.name, fold.enrichment)) )+
  geom_point(aes(size = avg.exp.scaled, color = neglog10pval)) +
  coord_cartesian(clip = "off") +
  scale_color_gradient(low = "#F88586", high = "#FC0D0F", name = "Motif enrichment\n-log10(P value)", limits = c(-log10(0.05), 350)) +
  scale_x_continuous(breaks = 1) +  # <--- Only show 1 on x-axis
  theme_bw() +
  scale_size_continuous(name = "Average gene expression\n(normalized)",
                        range = c(0.2, 6)) +
  labs(x = "Fold enrichment\n(Over-representation, enriched by unique peaks in each cell type)", y = "TF motif Name", title = "") +
  facet_grid(family ~ cluster, scales = "free_y",space = "free_y",
             labeller = label_wrap_gen(width = 28)) + # Wrap family names
  theme(axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        axis.text.y = element_text(size = 10, color = 'black'),
        legend.position = "bottom",
        strip.text.y = element_text(size = 10, angle = 0), # Horizontal family names
        strip.text.x = element_text(size = 10) , # Style cluster labels (horizontal)
        panel.spacing.y = unit(0.2, "lines"))
ggplot1
dev.off()

#   geom_point(aes(size = pct.expression, color = neglog10pval)) +
#   # scale_color_paletteer_c(palette = "grDevices::topo.colors",
#   #                         name="-log10 pval(>0.05)",limits = c(-log10(0.05), 300))+
#   # scale_color_gradientn(colors = c("lightblue", "purple3", "red2"), name = "neglog10pval") +  # Continuous color: blue-purple-red
#   theme(axis.title.x = element_text(size = 5),  # Smaller font size for x-label
# ggplot2

dat_TNBC_motif <- unique(select(dat_TNBC, motif, motif.name, family))
dat_TNBC_motif_group <- dat_TNBC_motif %>%
  group_by(family) %>%
  mutate(count = n()) %>%
  select(family, count) %>%
  ungroup() %>%
  unique()

library(stringr)
jpeg('results/tfgwas/GWASMotif01.TNBC_motif_count.jpg',
     width=6, height = 5, units = 'in', res = 250)
ggplot(dat_TNBC_motif_group, aes(x = count, y = reorder(family, count), fill = count)) +
  geom_bar(stat = "identity", fill = "#FF4500") +
  labs(x = "No. of TF motifs", y = "TF motif family", title = "") +
  scale_y_discrete(labels = function(x) str_wrap(x, width = 40)) +
  theme_classic() +
  theme(legend.position = "right") +
  scale_x_continuous(breaks = seq(0, 5, 1), limits = c(0, 6))
dev.off()

dat_TNBC_rsrp<-dat_TNBC%>% filter(motif.name %in% c('TFAP4','FOSL1::JUNB'))
ggplot1 <- ggplot(dat_TNBC_rsrp, aes(x = fold.enrichment, y = reorder(motif.name, fold.enrichment)) )+
  geom_point(aes(size = avg.exp.scaled, color = neglog10pval)) +
  coord_cartesian(clip = "off") +
  scale_color_gradient(low = "#F88586", high = "#FC0D0F", name = "Motif enrichment\n-log10(P value)", limits = c(-log10(0.05), 300)) +
  theme_bw() +
  scale_size_continuous(name = "Average gene expression\n(normalized)") +
  labs(x = "Fold enrichment\n(Over-representation, enriched by unique peaks in each cell type)", y = "TF motif Name", title = "") +
  facet_grid(family ~ cluster, scales = "free_y",space = "free_y",
             labeller = label_wrap_gen(width = 45)) + # Wrap family names
  theme(axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        axis.text.y = element_text(size = 10, color = 'black'),
        legend.position = "bottom",
        strip.text.y = element_text(size = 10, angle = 0), # Horizontal family names
        strip.text.x = element_text(size = 12) , # Style cluster labels (horizontal)
        panel.spacing.y = unit(0.2, "lines"))
ggplot1

##########################################################################################
#                                 feature plot                                           #
##########################################################################################
DefaultAssay(atac_subset_qc) <- "chromvar"

TNBC_unique<-fread("results/tfgwas/tfgwas_significant_TNBC_unique.csv",header = TRUE)
TNBC_unique<-list(unique(TNBC_unique$motif))
set.seed(4242025)
atac_subset_qc <- AddModuleScore(
  object = atac_subset_qc,
  features = TNBC_unique,
  ctrl = 15,
  name = 'TNBC_unique_motif')

p<-FeaturePlot(
  object = atac_subset_qc,
  features = "TNBC_unique_motif1",
  min.cutoff = 'q5',
  max.cutoff = 'q95',
  pt.size = 0.1,
  label = TRUE,
  cols = c( "grey90","orange")
  )+
  xlab("UMAP 1") +
  ylab("UMAP 2")+
  theme(legend.position = "right")

# Extract UMAP coordinates and cluster assignments
umap_data <- as.data.frame(atac_subset_qc@reductions$umap.rna@cell.embeddings)
umap_data$cluster <- atac_subset_qc@meta.data$celltype_RNA

# Filter for the clusters of interest
clusters_to_encircle <- c("LHS", "LASP-1", "LASP-2", "B-Myo")
umap_subset <- umap_data[umap_data$cluster %in% clusters_to_encircle, ]

points_matrix <- as.matrix(umap_subset[, c("rnaUMAP_1", "rnaUMAP_2")])
hull <- inla.nonconvex.hull(points_matrix, convex = -0.0021)
# Adjust convex for tightness (more negative for less smoothing)
segments <- data.frame(
  x = hull$loc[hull$idx[, 1], 1],
  y = hull$loc[hull$idx[, 1], 2],
  xend = hull$loc[hull$idx[, 2], 1],
  yend = hull$loc[hull$idx[, 2], 2]
)

p <- p + geom_segment(
  data = segments,
  aes(x = x, y = y, xend = xend, yend = yend),
  color = "black",  # Outline color
  linetype = "dotted",  # Dotted line
  size = 0.5  # Line thickness
)

jpeg('results/tfgwas/GWASMotif01.all_TNBC_uniqueMotifs_Addmodulescore57motifs.jpg',
     width=6, height = 5, units = 'in', res = 150)
print(p)
dev.off()

jpeg('results/tfgwas/GWASMotif01_TNBC_uniqueMotifs_Addmodulescore_57motifs.jpg',width=10, height = 5, units = 'in', res = 150)
a<-FeaturePlot(
  object = subset(atac_subset_qc, idents=c('LHS','LASP-1','LASP-2','B-Myo')),
  features = "TNBC_unique_motif1",
  min.cutoff = 'q5',
  max.cutoff = 'q95',
  pt.size = 0.1,
  label = TRUE
  # ,cols = c("orange", "darkblue")
  )
b<-FeaturePlot(
  object = subset(atac_subset_qc, idents=c('T','B','Myl','Endo','Pc','Fb')),
  features = "TNBC_unique_motif1",
  min.cutoff = 'q5',
  max.cutoff = 'q95',
  pt.size = 0.1,
  label = TRUE
  # ,cols = c("orange", "darkblue")
)
a &NoLegend()|b
dev.off()

VlnPlot(atac_subset_qc, features = "TNBC_unique_motif1", pt.size = 0.1) +
  ggtitle("Violin Plot of TNBC_unique_motif1")

df <- atac_subset_qc@meta.data

# Convert celltype_RNA to a factor with specified order
celltype_order <- c("LASP-1", "LASP-2", "LHS", "B-Myo", "Fb", "Endo", "T", "B", "Myl", "Pc")
df$celltype_RNA <- factor(df$celltype_RNA, levels = celltype_order)

# Perform statistical comparison (e.g., Kruskal-Wallis test for non-parametric data)
kw<-kruskal.test(TNBC_unique_motif1 ~ celltype_RNA, data = atac_subset_qc@meta.data)

df_main <- df %>% filter(celltype_RNA %in% c("LHS", "LASP-1", "LASP-2", "B-Myo"))
# Step 2: Subset for other subtypes (Fb, Endo, T, B, Myl, Pc)
df_other <- df %>% filter(celltype_RNA %in% c("Fb", "Endo", "T", "B", "Myl", "Pc"))

library(dunn.test)
dunn_test_main <- dunn.test(
  x = df_main$TNBC_unique_motif1,
  g = df_main$celltype_RNA,
  method = "bh" )
dunn_df_main <- data.frame(
  Comparison = dunn_test_main$comparisons,
  Z = dunn_test_main$Z,
  P_unadj = dunn_test_main$P,
  P_adj = dunn_test_main$P.adjusted)
my_comparisons1 <- list(
  c("LASP-1", "LASP-2"), c("LASP-1", "LHS"), c("LASP-1", "B-Myo"),
  c("LASP-2", "LHS"), c("LASP-2", "B-Myo"), c("LHS", "B-Myo")
)
##  3. Normalise Dunn comparison strings (both directions)
dunn_df_main <- dunn_df_main %>%
  tidyr::separate(Comparison, into = c("g1","g2"), sep = " - ", remove = FALSE) %>%
  mutate(
    key_fwd = paste(g1, g2, sep = "-"),
    key_rev = paste(g2, g1, sep = "-"))

##  4. Map to the exact order of my_comparisons
target_keys <- sapply(my_comparisons1, paste, collapse = "-")

dunn_df_main <- dunn_df_main %>%
  mutate(
    match_key = case_when(
      key_fwd %in% target_keys ~ key_fwd,
      key_rev %in% target_keys ~ key_rev,
      TRUE                     ~ NA_character_
    )
  ) %>%
  filter(!is.na(match_key)) %>%                     # keep only the 15 requested pairs
  arrange(match(match_key, target_keys))            # order exactly as target_keys
# Add significance stars
dunn_df_main <- dunn_df_main %>%
  mutate(signif = case_when(
    is.na(P_adj)                  ~ "NA",
    P_adj == 0                    ~ "***",
    P_adj <= 0.001                ~ "***",
    P_adj <= 0.01                 ~ "**",
    P_adj <= 0.05                 ~ "*",
    P_adj <= 0.1                  ~ ".",
    TRUE                          ~ "ns"
  ))

# Plot
p1 <- ggboxplot(df_main, x = "celltype_RNA", y = "TNBC_unique_motif1",
                color = "celltype_RNA", palette = "jco", legend = "none") +
  stat_compare_means(method = "kruskal.test", label.y = 2.5, size = 3) +
  geom_signif(comparisons = my_comparisons1,
              annotations = dunn_df_main$signif,
              y_position = seq(1.2, by = 0.18, length.out = 6),  # adjust spacing
              textsize = 3,
              vjust = 0.5) +
  labs(x = "Cell type", y = "AddModuleScore by TNBC unique TF motifs", title = "") +
  theme(axis.title.y = element_text(size = 8))
jpeg('results/tfgwas/Boxplot_TNBC_uniqueMotifs_Addmodulescore_part1.jpg',width=3.5, height = 5, units = 'in', res = 150)
print(p1)
dev.off()

write.csv(dunn_df_main,'results/tfgwas/dunn_df_main_TNBC.csv',row.names = FALSE)

# other cell type
dunn_test_other <- dunn.test(
  x = df_other$TNBC_unique_motif1,
  g = df_other$celltype_RNA,
  method = "bh" )
dunn_df_other <- data.frame(
  Comparison = dunn_test_other$comparisons,
  Z = dunn_test_other$Z,
  P_unadj = dunn_test_other$P,
  P_adj = dunn_test_other$P.adjusted)
##  2. Desired comparison order (my_comparisons2)
my_comparisons2 <- list(
  c("Fb","Endo"), c("Fb","T"),   c("Fb","B"),   c("Fb","Myl"), c("Fb","Pc"),
  c("Endo","T"),  c("Endo","B"), c("Endo","Myl"),c("Endo","Pc"),
  c("T","B"),     c("T","Myl"),  c("T","Pc"),
  c("B","Myl"),   c("B","Pc"),   c("Myl","Pc")
)
##  3. Normalise Dunn comparison strings (both directions)
dunn_df_other <- dunn_df_other %>%
  tidyr::separate(Comparison, into = c("g1","g2"), sep = " - ", remove = FALSE) %>%
  mutate(
    key_fwd = paste(g1, g2, sep = "-"),
    key_rev = paste(g2, g1, sep = "-"))

##  4. Map to the exact order of my_comparisons2
target_keys <- sapply(my_comparisons2, paste, collapse = "-")

dunn_df_other <- dunn_df_other %>%
  mutate(
    match_key = case_when(
      key_fwd %in% target_keys ~ key_fwd,
      key_rev %in% target_keys ~ key_rev,
      TRUE                     ~ NA_character_
    )
  ) %>%
  filter(!is.na(match_key)) %>%                     # keep only the 15 requested pairs
  arrange(match(match_key, target_keys))            # order exactly as target_keys

##  5. (Optional) Re-create signif column – safe for NA
dunn_df_other <- dunn_df_other %>%
  mutate(
    signif = case_when(
      is.na(P_adj)               ~ "NA",
      P_adj == 0                 ~ "***",
      P_adj <= 0.001             ~ "***",
      P_adj <= 0.01              ~ "**",
      P_adj <= 0.05              ~ "*",
      P_adj <= 0.1               ~ ".",
      TRUE                       ~ "ns") )
##  6. Plot
jpeg('results/tfgwas/Boxplot_TNBC_uniqueMotifs_Addmodulescore_part2.jpg',
     width = 4.5, height = 6, units = 'in', res = 150)
p2 <- ggboxplot(df_other,
                x = "celltype_RNA",
                y = "TNBC_unique_motif1",
                color = "celltype_RNA",
                palette = "YlGnBu",
                legend = "none") +

  ## overall Kruskal-Wallis (small font)
  stat_compare_means(method = "kruskal.test",
                     label.y = 6,
                     size = 3) +
  ## 15 pairwise Dunn stars (BH-adjusted)
  geom_signif(comparisons = my_comparisons2,
              annotations = dunn_df_other$signif,
              y_position = seq(1.4, by = 0.3, length.out = length(my_comparisons2)),
              textsize = 3,
              vjust = 0.5) +
  labs(x = "Cell type",
       y = "AddModuleScore by TNBC unique TF motifs",
       title = "") +
  theme(axis.title.y = element_text(size = 8))
print(p2)
dev.off()
write.csv(dunn_df_other,'results/tfgwas/dunn_df_other_TNBC.csv',row.names = FALSE)

#################################################################################################
#                          Unique Motif in LumA significant GWAS motif                         #
#################################################################################################
motif_id_list_final<-fread('results/singlecell/jaspar2022_motifs_727.csv', header = TRUE)
LumA<-fread("results/tfgwas/tfgwas_results_LumA_ALL.csv", header = TRUE)
motif_id_list_final$motif.name.orig<-motif_id_list_final$motif.name
LumA<- left_join(select(motif_id_list_final,motif,motif.name), LumA, by=c('motif','motif.name')) %>% dplyr::filter(TF_p_value < 0.05/727)
length(unique(LumA$motif))

TNBC<-fread("results/tfgwas/tfgwas_results_TNBC_ALL.csv", header = TRUE)
LumB<-fread("results/tfgwas/tfgwas_results_LumB_ALL.csv", header = TRUE)
LumB_HER2Neg<-fread("results/tfgwas/tfgwas_results_LumB_HER2Neg_ALL.csv", header = TRUE)
HER2<-fread("results/tfgwas/tfgwas_results_HER2_ALL.csv", header = TRUE)

all<-rbind(LumB, TNBC)
all<-rbind(all, LumB_HER2Neg)
all<-rbind(all, HER2)
all<-all %>% dplyr::filter(TF_p_value < 0.05/727)%>% unique()
#                               TF_p_value_LumB>=5e-4, TF_p_value_LumB_HER2Neg>=5e-4, TF_p_value_HER2>=5e-4)
other_sig_motif<-unique(all$motif)
length(other_sig_motif)
#
LumA_unique <- LumA %>% dplyr::filter(!motif %in% other_sig_motif)  %>% unique()
length(unique(LumA_unique$motif))
fwrite(LumA_unique,"results/tfgwas/tfgwas_significant_LumA_unique.csv", row.names = FALSE)

# # LumA_unique_common<-LumA_unique[!motif %in% motif_not_in_motifbreakR$motif,]
# LumA_unique[motif %in% motif_not_in_motifbreakR$motif,motif.name]

enriched.motifs1$motif.name<-NULL
enriched_motifs_list_LumA_unique <-left_join(LumA_unique, enriched.motifs1, by='motif') %>% unique()

enriched_motifs_list_LumA_unique$cluster <- unique(combined_peaks_gene_all$cluster)[1]

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
  enriched_motifs_filtered <- left_join(LumA_unique, enriched_motifs, by='motif') %>% unique()
  enriched_motifs_filtered$cluster <- cluster

  # Store the filtered results in the list with the cluster name as the key
  enriched_motifs_list_LumA_unique <- rbind(enriched_motifs_list_LumA_unique,enriched_motifs_filtered)

  # Print a message to indicate progress
  cat("Processed cluster:", cluster, "\n")
}

enriched_motifs_list_LumA_unique$pvalue <- ifelse(enriched_motifs_list_LumA_unique$pvalue==0,1e-300,
                                                  ifelse(enriched_motifs_list_LumA_unique$pvalue==1,0.99999999999999999999999,enriched_motifs_list_LumA_unique$pvalue))

# Save
write.csv(enriched_motifs_list_LumA_unique, file = "results/singlecell/enriched_motifs_list_uniqDAP_LumAUniqSignificant.csv", row.names = FALSE)

enriched_motifs_list_LumA_unique<-fread("results/singlecell/enriched_motifs_list_uniqDAP_LumAUniqSignificant.csv",header = TRUE)

# create a matrix
# Get all unique motif names across all clusters
all_motifsID <- unique(enriched_motifs_list_LumA_unique$motif)
all_motifs <- motif_id_list_final[motif %in% all_motifsID]
all_motifs$motif.name <- make.unique(all_motifs$motif.name)
all_motifs <- all_motifs$motif.name

# Get all cluster names
all_clusters <- unique(enriched_motifs_list_LumA_unique$cluster)

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
  enriched_data <- enriched_motifs_list_LumA_unique[enriched_motifs_list_LumA_unique$cluster==cluster1]

  # Add p-values to the matrix
  for (i in 1:nrow(enriched_data)) {
    motif <- enriched_data$motif.name[i]
    neglog10pvalue <- -log10(enriched_data$pvalue[i])  # Replace "pvalue" with the actual column name for p-values
    motif_matrix[motif,cluster1] <-neglog10pvalue
  }
}

# Convert the matrix to a data frame for better readability
motif_df<- as.data.frame(motif_matrix)

motif_df<-dplyr::select(motif_df,'LS-1','LASP-2','LHr','Bs','Fb','Endo','T','B','Myl','Pc')

motif_df<-motif_df%>%arrange(desc(LHr),desc(`LS-1`),desc(Bs))
library(pheatmap)
library(ComplexHeatmap)
library(circlize)

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
  #   breaks = rownames(motif_label),                 # Only label significant genes
  labs(fill = "-log10(pval)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1,size = 12),
        axis.text.y = element_text(angle = 0,size = 5)  # Adjust y-axis font size
  )
motif_heatmap
ggsave("results/singlecell/Motif_enricheded_heatmap_uniDAP_LumA_unique_signMotif_96motifs.jpg",
       motif_heatmap, width = 4, height = 8, dpi = 150)

# ---------------------heatmap ---------------#
motif_id_name<-fread('results/singlecell/jaspar2022_motifs_727.csv', header = TRUE)

# for the motif score heatmap
motif_matrix <- as.matrix(motif_df)

library(ComplexHeatmap)
library(dendextend)
set.seed(2025)
heatmap_obj_orig <-Heatmap(motif_matrix,
                           name = "-log10 pval",          # Title of legend
                           row_names_gp = gpar(fontsize = 7),
                           row_names_side = "left",
                           column_names_centered = TRUE,
                           width = unit(80, "mm"),
                           col = circlize::colorRamp2(c(0,80), c("#483d8b","red1")),
                           rect_gp = gpar(col = "grey99"), # Set cell border color (instead of border_color)
                           show_row_dend = FALSE,         # Hide row dendrogram to avoid km labels
                           show_column_dend = FALSE,
                           # column_split = 5,                 # Ensure 3 clusters without labels
                           row_title = NULL ,
                           column_title = NULL
)
heatmap_obj_orig
pdf("results/tfgwas/heatmap_enrichedmotif_uniqMotif_LumA_orginal.pdf", width = 5, height =10)
draw(heatmap_obj_orig)
dev.off()

# heatmap with specific motif sequence
motif_df<- as.data.frame(motif_matrix)

colSums(motif_matrix >= -log10(0.05 / 727))
# LHr LS-1 LASP-2   Bs    T    B  Myl Endo   Pc   Fb
colSums(motif_matrix >= 80)
# LHr LS-1 LASP-2   Bs    T    B  Myl Endo   Pc   Fb
colSums(motif_matrix >= 50)
# LHr LS-1 LASP-2   Bs    T    B  Myl Endo   Pc   Fb

motif_list <- c("SOX2","NFIX","TFAP2A","NFIC","NFIB","TGIF1","SCRT2","PKNOX1","ZNF528",
                "ONECUT1","ZNF24","REST","ONECUT3","ZFP14","ZNF320","PBX2","DMRTA1","POU4F1",
                "MEIS2","ZNF574","MSC","ZNF701","SIX2","ZFP42","AR" ,"NR4A1","TCF7L2",
                "ZNF684","ZNF667","LEF1","TFAP4::ETV1","ZKSCAN1","SOHLH2","BHLHE22","NEUROD2","FOXD1","FOXO3",
                "FOS::JUND","FOSL1::JUN","FOSL2","FOS","BNC2","JUND","NFE2","LMX1B","LMX1A",
                "ALX3","RAX","SHOX","VAX2","KLF17","CEBPG","NKX6-2","RAX2","MEOX2",
                "NFKB1","LHX5","VAX1","LBX2","LHX2","PHOX2B","HOXB4","HOXD4","PPARG",
                "NKX6-3","ESR2","GBX2","NEUROG2","EGR2","EGR4","OLIG2","NFATC4","NRL","MAX::MYC","ZNF610","MAX",
                "HSF4","ZNF317","VEZF1","HSF1","RARA::RXRG","BARX2","HSF2","BSX",
                "GLIS2","LHX1", "DLX1","MSX1","TFEC","ZNF669","EN2",
                "SMAD3","SRF","GLIS3","ZNF76","HOXC8", "NR2F1" ,"DLX5", "CTCF","PAX5","RFX7","THAP11","NFYA","YY1")%>%unique()

motif_df<- as.data.frame(motif_matrix)
head(motif_df)
# Reorder motif_df rows based on motif_list
motif_df <- motif_df[match(motif_list, rownames(motif_df)), ]
motif_df <- motif_df %>% select('LS-1','LASP-2','LHr','Bs','Fb','Endo','T','B','Myl','Pc')
motif_matrix <- as.matrix(motif_df)

library(ComplexHeatmap)
library(circlize)

# Create the heatmap
pdf("results/tfgwas/heatmap_enrichedmotif_uniqMotif_LumA_104motifs.pdf", width = 6, height =11)
Heatmap(motif_matrix,
        col =  circlize::colorRamp2(c(0, 150), c("#483d8b", "red1")),
        row_names_gp = gpar(fontsize = 7.5),
        cluster_rows = FALSE,  # Preserve the order from motif_list
        cluster_columns = FALSE,  # Cluster columns (or set to FALSE if you want to preserve column order)
        show_row_names = TRUE,  # Show motif names
        show_column_names = TRUE,  # Show column names
        row_names_side = "left",  # Place row names on the left
        heatmap_legend_param = list(title = "-log10 pval"),
        rect_gp = gpar(col = "grey99"),)  # Customize legend title

dev.off()

############################ GSEA bubble plots ###########################################
# Load libraries
library(tidyverse)
library(ggplot2)
library(paletteer)
#
#

dat_LumA<-fread("results/singlecell/enriched_motifs_list_uniqDAP_LumAUniqSignificant.csv",header = TRUE)
dat_LumA$cluster <- gsub('LS-1','LASP-1',dat_LumA$cluster )
dat_LumA$cluster <- gsub('LS-2','LASP-2',dat_LumA$cluster )
dat_LumA$cluster <- gsub('LHr','LHS',dat_LumA$cluster )
dat_LumA$cluster <- gsub('Bs','B-Myo',dat_LumA$cluster )

dat_LumA$cluster <- factor(dat_LumA$cluster, levels = c('LASP-1','LASP-2','LHS','B-Myo','Fb','Endo','T','B','Myl','Pc'))
dat_LumA$avg.exp.scaled <- ifelse(is.na(dat_LumA$avg.exp.scaled), 0, dat_LumA$avg.exp.scaled)

dat_LumA1 <- dat_LumA[dat_LumA$family %in% unique(dat_LumA$family)[1:18],]
dat_LumA2 <- dat_LumA[dat_LumA$family %in% unique(dat_LumA$family)[19:40],]

jpeg('results/tfgwas/GSEAplot_LumA_uniquemotifs_part1.jpg',width = 15, height = 18, units = 'in', res = 500)
ggplot1 <- ggplot(dat_LumA1, aes(x = fold.enrichment, y = reorder(motif.name, fold.enrichment)) )+
  geom_point(aes(size = avg.exp.scaled, color = neglog10pval)) +
  coord_cartesian(clip = "off") +
  scale_x_continuous(breaks = 1) +  # <--- Only show 1 on x-axis
  scale_color_gradient(low = "#F88586", high = "#FC0D0F", name = "Motif enrichment\n-log10(P value)", limits = c(-log10(0.05), 300)) +
  theme_bw() +
  scale_size_continuous(name = "Average gene expression\n(normalized)") +
  labs(x = "Fold enrichment\n(Over-representation, enriched by unique peaks in each cell type)", y = "TF motif Name", title = "") +
  facet_grid(family ~ cluster, scales = "free_y",space = "free_y",
             labeller = label_wrap_gen(width = 45)) + # Wrap family names
  theme(axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        axis.text.y = element_text(size = 10, color = 'black'),
        legend.position = "bottom",
        strip.text.y = element_text(size = 10, angle = 0), # Horizontal family names
        strip.text.x = element_text(size = 12) , # Style cluster labels (horizontal)
        panel.spacing.y = unit(0.2, "lines"))
ggplot1
dev.off()

jpeg('results/tfgwas/GSEAplot_LumA_uniquemotifs_part2.jpg',width = 15, height = 15, units = 'in', res = 500)
ggplot1 <- ggplot(dat_LumA2, aes(x = fold.enrichment, y = reorder(motif.name, fold.enrichment)) )+
  geom_point(aes(size = avg.exp.scaled, color = neglog10pval)) +
  coord_cartesian(clip = "off") +
  scale_color_gradient(low = "#F88586", high = "#FC0D0F", name = "Motif enrichment\n-log10(P value)", limits = c(-log10(0.05), 300)) +
  theme_bw() +
  scale_x_continuous(breaks = 1) +  # <--- Only show 1 on x-axis
  scale_size_continuous(name = "Average gene expression\n(normalized)") +
  labs(x = "Fold enrichment\n(Over-representation, enriched by unique peaks in each cell type)", y = "TF motif Name", title = "") +
  facet_grid(family ~ cluster, scales = "free_y",space = "free_y",
             labeller = label_wrap_gen(width = 45)) + # Wrap family names
  theme(axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        axis.text.y = element_text(size = 10, color = 'black'),
        legend.position = "bottom",
        strip.text.y = element_text(size = 10, angle = 0), # Horizontal family names
        strip.text.x = element_text(size = 12) , # Style cluster labels (horizontal)
        panel.spacing.y = unit(0.2, "lines"))
ggplot1
dev.off()

dat_LumA_motif <- unique(select(dat_LumA, motif, motif.name, family))
dat_LumA_motif_group <- dat_LumA_motif %>%
  group_by(family) %>%
  mutate(count = n()) %>%
  select(family, count) %>%
  ungroup() %>%
  unique()
dat_LumA_motif_group<- dat_LumA_motif_group[order(-dat_LumA_motif_group$count),]
dat_LumA1_motif_group<-dat_LumA_motif_group[c(1:20),]
library(stringr)
jpeg('results/tfgwas/GWASMotif01.LumA_motif_count1.jpg',
     width=7, height = 10, units = 'in', res = 250)
ggplot(dat_LumA_motif_group, aes(x = count, y = reorder(family, count), fill = count)) +
  geom_bar(stat = "identity", fill = "#FF4500") +
  labs(x = "No. of TF motifs", y = "TF motif family", title = "") +
  scale_y_discrete(labels = function(x) str_wrap(x, width = 40)) +
  theme_classic() +
  theme(legend.position = "right") +
  scale_x_continuous(breaks = seq(0, 18, 1), limits = c(0, 18))
dev.off()

jpeg('results/tfgwas/GWASMotif01.LumA_motif_count2.jpg',
     width=7, height = 6, units = 'in', res = 250)
ggplot(dat_LumA1_motif_group, aes(x = count, y = reorder(family, count), fill = count)) +
  geom_bar(stat = "identity") +
  scale_fill_gradient(low = "#FFE4B5", high = "#FF4500", name = "Count") +
  labs(x = "Count", y = "TF motif family (top 20)", title = "") +
  scale_y_discrete(labels = function(x) str_wrap(x, width = 40)) +
  theme_classic() +
  theme(legend.position = "right") +
  scale_x_continuous(breaks = seq(0, 18, 1), limits = c(0, 18))
dev.off()
##########################################################################################
#                                 feature plot                                           #
##########################################################################################
DefaultAssay(atac_subset_qc) <- "chromvar"
LumA_unique<-fread("results/tfgwas/tfgwas_significant_LumA_unique.csv", header = TRUE)
LumA_unique<-list(unique(LumA_unique$motif))
set.seed(4242025)
atac_subset_qc <- AddModuleScore(
  object = atac_subset_qc,
  features = LumA_unique,
  ctrl = 15,
  name = 'LumA_unique_motif'
)

jpeg('results/tfgwas/GWASMotif01.all_LumA_uniqueMotifs_Addmodulescore.jpg',width=6, height = 5, units = 'in', res = 150)
p<-FeaturePlot(
  object = atac_subset_qc,
  features = "LumA_unique_motif1",
  min.cutoff = 'q5',
  max.cutoff = 'q95',
  pt.size = 0.1,
  label = TRUE,
  ,cols = c( "grey90","orange"))+
  xlab("UMAP 1") +
  ylab("UMAP 2")+
  theme(legend.position = "right")

# Extract UMAP coordinates and cluster assignments
umap_data <- as.data.frame(atac_subset_qc@reductions$umap.rna@cell.embeddings)
umap_data$cluster <- atac_subset_qc@meta.data$celltype_RNA

# Filter for the clusters of interest
clusters_to_encircle <- c("LHS", "LASP-1", "LASP-2", "B-Myo")
umap_subset <- umap_data[umap_data$cluster %in% clusters_to_encircle, ]

points_matrix <- as.matrix(umap_subset[, c("rnaUMAP_1", "rnaUMAP_2")])
hull <- inla.nonconvex.hull(points_matrix, convex = -0.0021)
# Adjust convex for tightness (more negative for less smoothing)
segments <- data.frame(
  x = hull$loc[hull$idx[, 1], 1],
  y = hull$loc[hull$idx[, 1], 2],
  xend = hull$loc[hull$idx[, 2], 1],
  yend = hull$loc[hull$idx[, 2], 2]
)

p <- p + geom_segment(
  data = segments,
  aes(x = x, y = y, xend = xend, yend = yend),
  color = "black",  # Outline color
  linetype = "dotted",  # Dotted line
  size = 0.5  # Line thickness
)
print(p)
dev.off()

jpeg('results/tfgwas/GWASMotif01_LumA_uniqueMotifs_Addmodulescore.jpg',width=10, height = 5, units = 'in', res = 150)
a<-FeaturePlot(
  object = subset(atac_subset_qc, idents=c('LHr','LS-1','LASP-2','Bs')),
  features = "LumA_unique_motif1",
  min.cutoff = 'q5',
  max.cutoff = 'q95',
  pt.size = 0.1,
  label = TRUE)
b<-FeaturePlot(
  object = subset(atac_subset_qc, idents=c('T','B','Myl','Endo','Pc','Fb')),
  features = "LumA_unique_motif1",
  min.cutoff = 'q5',
  max.cutoff = 'q95',
  pt.size = 0.1,
  label = TRUE
)
a &NoLegend()|b
dev.off()

df <- atac_subset_qc@meta.data
celltype_order <- c("LASP-1", "LASP-2", "LHS", "B-Myo", "Fb", "Endo", "T", "B", "Myl", "Pc")

# Convert celltype_RNA to a factor with specified order
levels(df$celltype_RNA)
df$celltype_RNA <- factor(df$celltype_RNA, levels = celltype_order)

# Perform statistical comparison (e.g., Kruskal-Wallis test for non-parametric data)
kw<-kruskal.test(LumA_unique_motif1 ~ celltype_RNA, data = atac_subset_qc@meta.data)

df_main <- df %>% filter(celltype_RNA %in% c("LHS", "LASP-1", "LASP-2", "B-Myo"))
# Step 2: Subset for other subtypes (Fb, Endo, T, B, Myl, Pc)
df_other <- df %>% filter(celltype_RNA %in% c("Fb", "Endo", "T", "B", "Myl", "Pc"))

library(dunn.test)
dunn_test_main <- dunn.test(
  x = df_main$LumA_unique_motif1,
  g = df_main$celltype_RNA,
  method = "bh" )
dunn_df_main <- data.frame(
  Comparison = dunn_test_main$comparisons,
  Z = dunn_test_main$Z,
  P_unadj = dunn_test_main$P,
  P_adj = dunn_test_main$P.adjusted)
my_comparisons1 <- list(
  c("LASP-1", "LASP-2"), c("LASP-1", "LHS"), c("LASP-1", "B-Myo"),
  c("LASP-2", "LHS"), c("LASP-2", "B-Myo"), c("LHS", "B-Myo")
)
##  3. Normalise Dunn comparison strings (both directions)
dunn_df_main <- dunn_df_main %>%
  tidyr::separate(Comparison, into = c("g1","g2"), sep = " - ", remove = FALSE) %>%
  mutate(
    key_fwd = paste(g1, g2, sep = "-"),
    key_rev = paste(g2, g1, sep = "-"))

##  4. Map to the exact order of my_comparisons
target_keys <- sapply(my_comparisons1, paste, collapse = "-")

dunn_df_main <- dunn_df_main %>%
  mutate(
    match_key = case_when(
      key_fwd %in% target_keys ~ key_fwd,
      key_rev %in% target_keys ~ key_rev,
      TRUE                     ~ NA_character_
    )
  ) %>%
  filter(!is.na(match_key)) %>%                     # keep only the 15 requested pairs
  arrange(match(match_key, target_keys))            # order exactly as target_keys
# Add significance stars
dunn_df_main <- dunn_df_main %>%
  mutate(signif = case_when(
    is.na(P_adj)                  ~ "NA",
    P_adj == 0                    ~ "***",
    P_adj <= 0.001                ~ "***",
    P_adj <= 0.01                 ~ "**",
    P_adj <= 0.05                 ~ "*",
    TRUE                          ~ "ns"
  ))

# Plot
p1 <- ggboxplot(df_main, x = "celltype_RNA", y = "LumA_unique_motif1",
                color = "celltype_RNA", palette = "jco", legend = "none") +
  stat_compare_means(method = "kruskal.test", label.y = 2, size = 3) +
  geom_signif(comparisons = my_comparisons1,
              annotations = dunn_df_main$signif,
              y_position = seq(0.6, by = 0.18, length.out = 6),  # adjust spacing
              textsize = 3,
              vjust = 0.5) +
  labs(x = "Cell type", y = "AddModuleScore by LumA unique TF motifs", title = "") +
  theme(axis.title.y = element_text(size = 8))
jpeg('results/tfgwas/Boxplot_LumA_uniqueMotifs_Addmodulescore_part1.jpg',width=3.5, height = 5, units = 'in', res = 150)
print(p1)
dev.off()

write.csv(dunn_df_main,'results/tfgwas/dunn_df_main_LumA.csv',row.names = FALSE)

# other cell type
dunn_test_other <- dunn.test(
  x = df_other$LumA_unique_motif1,
  g = df_other$celltype_RNA,
  method = "bh" )
dunn_df_other <- data.frame(
  Comparison = dunn_test_other$comparisons,
  Z = dunn_test_other$Z,
  P_unadj = dunn_test_other$P,
  P_adj = dunn_test_other$P.adjusted)
##  2. Desired comparison order (my_comparisons2)
my_comparisons2 <- list(
  c("Fb","Endo"), c("Fb","T"),   c("Fb","B"),   c("Fb","Myl"), c("Fb","Pc"),
  c("Endo","T"),  c("Endo","B"), c("Endo","Myl"),c("Endo","Pc"),
  c("T","B"),     c("T","Myl"),  c("T","Pc"),
  c("B","Myl"),   c("B","Pc"),   c("Myl","Pc")
)
##  3. Normalise Dunn comparison strings (both directions)
dunn_df_other <- dunn_df_other %>%
  tidyr::separate(Comparison, into = c("g1","g2"), sep = " - ", remove = FALSE) %>%
  mutate(
    key_fwd = paste(g1, g2, sep = "-"),
    key_rev = paste(g2, g1, sep = "-"))

##  4. Map to the exact order of my_comparisons2
target_keys <- sapply(my_comparisons2, paste, collapse = "-")

dunn_df_other <- dunn_df_other %>%
  mutate(
    match_key = case_when(
      key_fwd %in% target_keys ~ key_fwd,
      key_rev %in% target_keys ~ key_rev,
      TRUE                     ~ NA_character_
    )
  ) %>%
  filter(!is.na(match_key)) %>%                     # keep only the 15 requested pairs
  arrange(match(match_key, target_keys))            # order exactly as target_keys

##  5. (Optional) Re-create signif column – safe for NA
dunn_df_other <- dunn_df_other %>%
  mutate(
    signif = case_when(
      is.na(P_adj)               ~ "NA",
      P_adj == 0                 ~ "***",
      P_adj <= 0.001             ~ "***",
      P_adj <= 0.01              ~ "**",
      P_adj <= 0.05              ~ "*",
      TRUE                       ~ "ns") )
##  6. Plot
jpeg('results/tfgwas/Boxplot_LumA_uniqueMotifs_Addmodulescore_part2.jpg',
     width = 4.5, height = 6, units = 'in', res = 150)
p2 <- ggboxplot(df_other,
                x = "celltype_RNA",
                y = "LumA_unique_motif1",
                color = "celltype_RNA",
                palette = "YlGnBu",
                legend = "none") +

  ## overall Kruskal-Wallis (small font)
  stat_compare_means(method = "kruskal.test",
                     label.y = 6,
                     size = 3) +
  ## 15 pairwise Dunn stars (BH-adjusted)
  geom_signif(comparisons = my_comparisons2,
              annotations = dunn_df_other$signif,
              y_position = seq(1.4, by = 0.3, length.out = length(my_comparisons2)),
              textsize = 3,
              vjust = 0.5) +
  labs(x = "Cell type",
       y = "AddModuleScore by LumA unique TF motifs",
       title = "") +
  theme(axis.title.y = element_text(size = 8))
print(p2)
dev.off()
write.csv(dunn_df_other,'results/tfgwas/dunn_df_other_LumA.csv',row.names = FALSE)

#################################################################################################
#                          Unique Motif in LumB significant GWAS motif                         #
#################################################################################################
motif_id_list_final<-fread('results/singlecell/jaspar2022_motifs_727.csv', header = TRUE)
LumB<-fread("results/tfgwas/tfgwas_results_LumB_ALL.csv", header = TRUE)
LumB<- left_join(motif_id_list_final, LumB, by=c('motif','motif.name')) %>% dplyr::filter(TF_p_value < 0.05/727)
length(unique(LumB$motif))

LumA<-fread("results/tfgwas/tfgwas_results_LumA_ALL.csv", header = TRUE)
TNBC<-fread("results/tfgwas/tfgwas_results_TNBC_ALL.csv", header = TRUE)
LumB_HER2Neg<-fread("results/tfgwas/tfgwas_results_LumB_HER2Neg_ALL.csv", header = TRUE)
HER2<-fread("results/tfgwas/tfgwas_results_HER2_ALL.csv", header = TRUE)

all<-rbind(TNBC, LumA)
all<-rbind(all, LumB_HER2Neg)
all<-rbind(all, HER2)
all<-all %>% dplyr::filter(TF_p_value < 0.05/727)%>% unique()
#                               TF_p_value_TNBC>=5e-4, TF_p_value_LumB_HER2Neg>=5e-4, TF_p_value_HER2>=5e-4)
other_sig_motif<-unique(all$motif)
length(other_sig_motif)
# 5 new criteria

LumB_unique <- LumB %>% dplyr::filter(!motif %in% other_sig_motif)  %>% unique()
# #13
length(unique(LumB_unique$motif))
fwrite(LumB_unique,"results/tfgwas/tfgwas_significant_LumB_unique.csv", row.names = FALSE)
#

enriched.motifs1$motif.name<-NULL
enriched_motifs_list_LumB_unique <-left_join(LumB_unique, enriched.motifs1, by='motif') %>% unique()

enriched_motifs_list_LumB_unique$cluster <- unique(combined_peaks_gene_all$cluster)[1]

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
  enriched_motifs_filtered <- left_join(LumB_unique, enriched_motifs, by='motif') %>% unique()
  enriched_motifs_filtered$cluster <- cluster

  # Store the filtered results in the list with the cluster name as the key
  enriched_motifs_list_LumB_unique <- rbind(enriched_motifs_list_LumB_unique,enriched_motifs_filtered)

  # Print a message to indicate progress
  cat("Processed cluster:", cluster, "\n")
}

enriched_motifs_list_LumB_unique$pvalue <- ifelse(enriched_motifs_list_LumB_unique$pvalue==0,1e-300,
                                                  ifelse(enriched_motifs_list_LumB_unique$pvalue==1,0.99999999999999999999999,enriched_motifs_list_LumB_unique$pvalue))

# Save
write.csv(enriched_motifs_list_LumB_unique, file = "results/singlecell/enriched_motifs_list_uniqDAP_LumBUniqSignificant.csv", row.names = FALSE)

enriched_motifs_list_LumB_unique<-fread("results/singlecell/enriched_motifs_list_uniqDAP_LumBUniqSignificant.csv",header = TRUE)

# create a matrix
# Get all unique motif names across all clusters
all_motifsID <- unique(enriched_motifs_list_LumB_unique$motif)
all_motifs <- motif_id_list_final[motif %in% all_motifsID]
all_motifs$motif.name <- make.unique(all_motifs$motif.name)
all_motifs <- all_motifs$motif.name

# Get all cluster names
all_clusters <- unique(enriched_motifs_list_LumB_unique$cluster)

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
  enriched_data <- enriched_motifs_list_LumB_unique[enriched_motifs_list_LumB_unique$cluster==cluster1]

  # Add p-values to the matrix
  for (i in 1:nrow(enriched_data)) {
    motif <- enriched_data$motif.name[i]
    neglog10pvalue <- -log10(enriched_data$pvalue[i])  # Replace "pvalue" with the actual column name for p-values
    motif_matrix[motif,cluster1] <-neglog10pvalue
  }
}

# Convert the matrix to a data frame for better readability
motif_df<- as.data.frame(motif_matrix)

motif_df<-dplyr::select(motif_df,'LS-1','LASP-2','LHr','Bs','Fb','Endo','T','B','Myl','Pc')

motif_df<-motif_df%>%arrange(desc(LHr),desc(`LS-1`),desc(Bs))
library(pheatmap)
library(ComplexHeatmap)
library(circlize)

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
  #   breaks = rownames(motif_label),                 # Only label significant genes
  labs(fill = "-log10(pval)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1,size = 12),
        axis.text.y = element_text(angle = 0,size = 5)  # Adjust y-axis font size
  )
motif_heatmap
ggsave("results/singlecell/Motif_enricheded_heatmap_uniDAP_LumB_unique_signMotif.jpg",
       motif_heatmap, width = 4, height = 5, dpi = 150)

# ---------------------heatmap ---------------#
motif_id_name<-fread('results/singlecell/jaspar2022_motifs_727.csv', header = TRUE)

# for the motif score heatmap
motif_matrix <- as.matrix(motif_df)

library(ComplexHeatmap)
library(dendextend)
set.seed(2025)
heatmap_obj_orig <-Heatmap(motif_matrix,
                           name = "-log10 pval",          # Title of legend
                           row_names_gp = gpar(fontsize = 7),
                           row_names_side = "left",
                           column_names_centered = TRUE,
                           width = unit(80, "mm"),
                           col = circlize::colorRamp2(c(0,80), c("#483d8b","red1")),
                           rect_gp = gpar(col = "grey99"), # Set cell border color (instead of border_color)
                           show_row_dend = FALSE,         # Hide row dendrogram to avoid km labels
                           show_column_dend = FALSE,
                           # column_split = 4,                 # Ensure 3 clusters without labels
                           row_title = NULL ,
                           column_title = NULL
)
heatmap_obj_orig
pdf("results/tfgwas/heatmap_enrichedmotif_uniqMotif_LumB_orginal.pdf", width = 5, height =7)
draw(heatmap_obj_orig)
dev.off()

# ---------heatmap with specific motif sequence----------------
motif_df<- as.data.frame(motif_matrix)
motif_df <- motif_df %>% select('LS-1','LASP-2','LHr','Bs','Fb','Endo','T','B','Myl','Pc')

colSums(motif_matrix >= -log10(0.05 / 727))

colSums(motif_matrix >= 80)

colSums(motif_matrix >= 50)

motif_df <- motif_df[order(-motif_df$`LS-1`),]
rownames(motif_df)
motif_list <- c("SNAI2","TFAP2C" ,"HES6", "ARGFX", "PRRX2" , "LHX6",
                "JUN::JUNB","NR1I3","FOSB::JUNB","NR1H4::RXRA","OVOL1","CREM", "OTX1")%>% unique()

# create the new data frame
motif_df<- as.data.frame(motif_matrix)
head(motif_df)
# Reorder motif_df rows based on motif_list
motif_df <- motif_df[match(motif_list, rownames(motif_df)), ]
motif_df <- motif_df %>% select('LS-1','LASP-2','LHr','Bs','Fb','Endo','T','B','Myl','Pc')
motif_matrix <- as.matrix(motif_df)

library(ComplexHeatmap)
library(circlize)

# Create the heatmap
pdf("results/tfgwas/heatmap_enrichedmotif_uniqMotif_LumB_13motifs.pdf", width = 5, height =2.5)
Heatmap(motif_matrix,
        col =  circlize::colorRamp2(c(0, 80), c("#483d8b", "red1")),
        row_names_gp = gpar(fontsize = 7.5),
        cluster_rows = FALSE,  # Preserve the order from motif_list
        cluster_columns = FALSE,  # Cluster columns (or set to FALSE if you want to preserve column order)
        show_row_names = TRUE,  # Show motif names
        show_column_names = TRUE,  # Show column names
        row_names_side = "left",  # Place row names on the left
        heatmap_legend_param = list(title = "-log10 pval"),
        rect_gp = gpar(col = "grey99"),)  # Customize legend title

dev.off()

############################ GSEA bubble plots ###########################################
# Load libraries
library(tidyverse)
library(ggplot2)
library(paletteer)
library(tidyverse)
library(ggplot2)
library(paletteer)
dat_LumB<-fread("results/singlecell/enriched_motifs_list_uniqDAP_LumBUniqSignificant.csv",header = TRUE)
dat_LumB$cluster <- gsub('LS-1','LASP-1',dat_LumB$cluster )
dat_LumB$cluster <- gsub('LS-2','LASP-2',dat_LumB$cluster )
dat_LumB$cluster <- gsub('LHr','LHS',dat_LumB$cluster )
dat_LumB$cluster <- gsub('Bs','B-Myo',dat_LumB$cluster )
#
#
dat_LumB$cluster <- factor(dat_LumB$cluster, levels = c('LASP-1','LASP-2','LHS','B-Myo','Fb','Endo','T','B','Myl','Pc'))
dat_LumB$avg.exp.scaled <- ifelse(is.na(dat_LumB$avg.exp.scaled), 0, dat_LumB$avg.exp.scaled)
#

jpeg('results/tfgwas/GSEAplot_LumB_uniquemotifs.jpg',width = 15, height = 8, units = 'in', res = 500)
ggplot1 <- ggplot(dat_LumB, aes(x = fold.enrichment, y = reorder(motif.name, fold.enrichment)) )+
  geom_point(aes(size = avg.exp.scaled, color = neglog10pval)) +
  coord_cartesian(clip = "off") +
  scale_color_gradient(low = "#F88586", high = "#FC0D0F", name = "Motif enrichment\n-log10(P value)", limits = c(-log10(0.05), 300)) +
  theme_bw() +
  scale_x_continuous(breaks = 1) +
  scale_size_continuous(name = "Average gene expression\n(normalized)") +
  labs(x = "Fold enrichment\n(Over-representation, enriched by unique peaks in each cell type)", y = "TF motif Name", title = "") +
  facet_grid(family ~ cluster, scales = "free_y",space = "free_y",
             labeller = label_wrap_gen(width = 45)) + # Wrap family names
  theme(axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        axis.text.y = element_text(size = 10, color = 'black'),
        legend.position = "bottom",
        strip.text.y = element_text(size = 10, angle = 0), # Horizontal family names
        strip.text.x = element_text(size = 12) , # Style cluster labels (horizontal)
        panel.spacing.y = unit(0.2, "lines"))
ggplot1
dev.off()

dat_LumB_motif <- unique(select(dat_LumB, motif, motif.name, family))
dat_LumB_motif_group <- dat_LumB_motif %>%
  group_by(family) %>%
  mutate(count = n()) %>%
  select(family, count) %>%
  ungroup() %>%
  unique()

library(stringr)
jpeg('results/tfgwas/GWASMotif01.LumB_motif_count.jpg',
     width=6, height = 5, units = 'in', res = 250)
ggplot(dat_LumB_motif_group, aes(x = count, y = reorder(family, count), fill = count)) +
  geom_bar(stat = "identity", fill = "#FF4500") +
  labs(x = "No. of TF motifs", y = "TF motif family", title = "") +
  scale_y_discrete(labels = function(x) str_wrap(x, width = 40)) +
  theme_classic() +
  theme(legend.position = "right") +
  scale_x_continuous(breaks = seq(0, 3, 1), limits = c(0, 4))
dev.off()
##########################################################################################
#                                 feature plot                                           #
##########################################################################################
DefaultAssay(atac_subset_qc) <- "chromvar"

LumB_unique<-fread("results/tfgwas/tfgwas_significant_LumB_unique.csv", header=TRUE)
LumB_unique<-list(unique(LumB_unique$motif))
set.seed(4242025)
atac_subset_qc <- AddModuleScore(
  object = atac_subset_qc,
  features = LumB_unique,
  ctrl = 15,
  name = 'LumB_unique_motif'
)

jpeg('results/tfgwas/GWASMotif01.all_LumB_uniqueMotifs_Addmodulescore.jpg',width=6, height = 5, units = 'in', res = 150)
p<-FeaturePlot(
  object = atac_subset_qc,
  features = "LumB_unique_motif1",
  min.cutoff = 'q5',
  max.cutoff = 'q95',
  pt.size = 0.1,
  label = TRUE,
  cols = c( "grey90","orange")
)+
  xlab("UMAP 1") +
  ylab("UMAP 2")+
  theme(legend.position = "right")

# Extract UMAP coordinates and cluster assignments
umap_data <- as.data.frame(atac_subset_qc@reductions$umap.rna@cell.embeddings)
umap_data$cluster <- atac_subset_qc@meta.data$celltype_RNA

# Filter for the clusters of interest
clusters_to_encircle <- c("LHS", "LASP-1", "LASP-2", "B-Myo")
umap_subset <- umap_data[umap_data$cluster %in% clusters_to_encircle, ]

points_matrix <- as.matrix(umap_subset[, c("rnaUMAP_1", "rnaUMAP_2")])
hull <- inla.nonconvex.hull(points_matrix, convex = -0.0021)
# Adjust convex for tightness (more negative for less smoothing)
segments <- data.frame(
  x = hull$loc[hull$idx[, 1], 1],
  y = hull$loc[hull$idx[, 1], 2],
  xend = hull$loc[hull$idx[, 2], 1],
  yend = hull$loc[hull$idx[, 2], 2]
)

p <- p + geom_segment(
  data = segments,
  aes(x = x, y = y, xend = xend, yend = yend),
  color = "black",  # Outline color
  linetype = "dotted",  # Dotted line
  size = 0.5  # Line thickness
)
print(p)
dev.off()

jpeg('results/tfgwas/GWASMotif01_LumB_uniqueMotifs_Addmodulescore.jpg',width=11, height = 5, units = 'in', res = 150)
a<-FeaturePlot(
  object = subset(atac_subset_qc, idents=c('LHS','LASP-1','LASP-2','B-Myo')),
  features = "LumB_unique_motif1",
  min.cutoff = 'q5',
  max.cutoff = 'q95',
  pt.size = 0.1,
  label = TRUE)
b<-FeaturePlot(
  object = subset(atac_subset_qc, idents=c('T','B','Myl','Endo','Pc','Fb')),
  features = "LumB_unique_motif1",
  min.cutoff = 'q5',
  max.cutoff = 'q95',
  pt.size = 0.1,
  label = TRUE
)
a &NoLegend()|b
dev.off()

df <- atac_subset_qc@meta.data
celltype_order <- c("LASP-1", "LASP-2", "LHS", "B-Myo", "Fb", "Endo", "T", "B", "Myl", "Pc")

# Convert celltype_RNA to a factor with specified order
levels(df$celltype_RNA)
df$celltype_RNA <- factor(df$celltype_RNA, levels = celltype_order)

# Perform statistical comparison (e.g., Kruskal-Wallis test for non-parametric data)
kw<-kruskal.test(LumB_unique_motif1 ~ celltype_RNA, data = atac_subset_qc@meta.data)

df_main <- df %>% filter(celltype_RNA %in% c("LHS", "LASP-1", "LASP-2", "B-Myo"))
# Step 2: Subset for other subtypes (Fb, Endo, T, B, Myl, Pc)
df_other <- df %>% filter(celltype_RNA %in% c("Fb", "Endo", "T", "B", "Myl", "Pc"))

library(dunn.test)
dunn_test_main <- dunn.test(
  x = df_main$LumB_unique_motif1,
  g = df_main$celltype_RNA,
  method = "bh" )
dunn_df_main <- data.frame(
  Comparison = dunn_test_main$comparisons,
  Z = dunn_test_main$Z,
  P_unadj = dunn_test_main$P,
  P_adj = dunn_test_main$P.adjusted)
my_comparisons1 <- list(
  c("LASP-1", "LASP-2"), c("LASP-1", "LHS"), c("LASP-1", "B-Myo"),
  c("LASP-2", "LHS"), c("LASP-2", "B-Myo"), c("LHS", "B-Myo")
)
##  3. Normalise Dunn comparison strings (both directions)
dunn_df_main <- dunn_df_main %>%
  tidyr::separate(Comparison, into = c("g1","g2"), sep = " - ", remove = FALSE) %>%
  mutate(
    key_fwd = paste(g1, g2, sep = "-"),
    key_rev = paste(g2, g1, sep = "-"))

##  4. Map to the exact order of my_comparisons
target_keys <- sapply(my_comparisons1, paste, collapse = "-")

dunn_df_main <- dunn_df_main %>%
  mutate(
    match_key = case_when(
      key_fwd %in% target_keys ~ key_fwd,
      key_rev %in% target_keys ~ key_rev,
      TRUE                     ~ NA_character_
    )
  ) %>%
  filter(!is.na(match_key)) %>%                     # keep only the 15 requested pairs
  arrange(match(match_key, target_keys))            # order exactly as target_keys
# Add significance stars
dunn_df_main <- dunn_df_main %>%
  mutate(signif = case_when(
    is.na(P_adj)                  ~ "NA",
    P_adj == 0                    ~ "***",
    P_adj <= 0.001                ~ "***",
    P_adj <= 0.01                 ~ "**",
    P_adj <= 0.05                 ~ "*",
    TRUE                          ~ "ns"
  ))

# Plot
p1 <- ggboxplot(df_main, x = "celltype_RNA", y = "LumB_unique_motif1",
                color = "celltype_RNA", palette = "jco", legend = "none") +
  stat_compare_means(method = "kruskal.test", label.y = 3.5, size = 3) +
  geom_signif(comparisons = my_comparisons1,
              annotations = dunn_df_main$signif,
              y_position = seq(1.6, by = 0.19, length.out = 6),  # adjust spacing
              textsize = 3,
              vjust = 0.5) +
  labs(x = "Cell type", y = "AddModuleScore by LumB unique TF motifs", title = "") +
  theme(axis.title.y = element_text(size = 8))
jpeg('results/tfgwas/Boxplot_LumB_uniqueMotifs_Addmodulescore_part1.jpg',width=3.5, height = 5, units = 'in', res = 150)
print(p1)
dev.off()

write.csv(dunn_df_main,'results/tfgwas/dunn_df_main_LumB.csv',row.names = FALSE)

# other cell type
dunn_test_other <- dunn.test(
  x = df_other$LumB_unique_motif1,
  g = df_other$celltype_RNA,
  method = "bh" )
dunn_df_other <- data.frame(
  Comparison = dunn_test_other$comparisons,
  Z = dunn_test_other$Z,
  P_unadj = dunn_test_other$P,
  P_adj = dunn_test_other$P.adjusted)
##  2. Desired comparison order (my_comparisons2)
my_comparisons2 <- list(
  c("Fb","Endo"), c("Fb","T"),   c("Fb","B"),   c("Fb","Myl"), c("Fb","Pc"),
  c("Endo","T"),  c("Endo","B"), c("Endo","Myl"),c("Endo","Pc"),
  c("T","B"),     c("T","Myl"),  c("T","Pc"),
  c("B","Myl"),   c("B","Pc"),   c("Myl","Pc")
)
##  3. Normalise Dunn comparison strings (both directions)
dunn_df_other <- dunn_df_other %>%
  tidyr::separate(Comparison, into = c("g1","g2"), sep = " - ", remove = FALSE) %>%
  mutate(
    key_fwd = paste(g1, g2, sep = "-"),
    key_rev = paste(g2, g1, sep = "-"))

##  4. Map to the exact order of my_comparisons2
target_keys <- sapply(my_comparisons2, paste, collapse = "-")

dunn_df_other <- dunn_df_other %>%
  mutate(
    match_key = case_when(
      key_fwd %in% target_keys ~ key_fwd,
      key_rev %in% target_keys ~ key_rev,
      TRUE                     ~ NA_character_
    )
  ) %>%
  filter(!is.na(match_key)) %>%                     # keep only the 15 requested pairs
  arrange(match(match_key, target_keys))            # order exactly as target_keys

##  5. (Optional) Re-create signif column – safe for NA
dunn_df_other <- dunn_df_other %>%
  mutate(
    signif = case_when(
      is.na(P_adj)               ~ "NA",
      P_adj == 0                 ~ "***",
      P_adj <= 0.001             ~ "***",
      P_adj <= 0.01              ~ "**",
      P_adj <= 0.05              ~ "*",
      TRUE                       ~ "ns") )
##  6. Plot
jpeg('results/tfgwas/Boxplot_LumB_uniqueMotifs_Addmodulescore_part2.jpg',
     width = 4.5, height = 6, units = 'in', res = 150)
p2 <- ggboxplot(df_other,
                x = "celltype_RNA",
                y = "LumB_unique_motif1",
                color = "celltype_RNA",
                palette = "YlGnBu",
                legend = "none") +

  ## overall Kruskal-Wallis (small font)
  stat_compare_means(method = "kruskal.test",
                     label.y = 7,
                     size = 3) +
  ## 15 pairwise Dunn stars (BH-adjusted)
  geom_signif(comparisons = my_comparisons2,
              annotations = dunn_df_other$signif,
              y_position = seq(2, by = 0.3, length.out = length(my_comparisons2)),
              textsize = 3,
              vjust = 0.5) +
  labs(x = "Cell type",
       y = "AddModuleScore by LumB unique TF motifs",
       title = "") +
  theme(axis.title.y = element_text(size = 8))
print(p2)
dev.off()
write.csv(dunn_df_other,'results/tfgwas/dunn_df_other_LumB.csv',row.names = FALSE)

#################################################################################################
#                          Unique Motif in LumB_HER2Neg significant GWAS motif                         #
#################################################################################################
motif_id_list_final<-fread('results/singlecell/jaspar2022_motifs_727.csv', header = TRUE)
LumB_HER2Neg<-fread("results/tfgwas/tfgwas_results_LumB_HER2Neg_ALL.csv", header = TRUE)
LumB_HER2Neg<- left_join(select(motif_id_list_final,motif,motif.name), LumB_HER2Neg, by=c('motif','motif.name')) %>%
  dplyr::filter(TF_p_value < 0.05/727)
length(unique(LumB_HER2Neg$motif))
#433 previous Bonferroni corrected

TNBC<-fread("results/tfgwas/tfgwas_results_TNBC_ALL.csv", header = TRUE)
LumB<-fread("results/tfgwas/tfgwas_results_LumB_ALL.csv", header = TRUE)
LumA<-fread("results/tfgwas/tfgwas_results_LumA_ALL.csv", header = TRUE)
HER2<-fread("results/tfgwas/tfgwas_results_HER2_ALL.csv", header = TRUE)

all<-rbind(TNBC, LumA)
all<-rbind(all, LumB)
all<-rbind(all, HER2)
all<-all %>% dplyr::filter(TF_p_value < 0.05/727)%>% unique()
#                               TF_p_value_LumA>=5e-4, TF_p_value_LumB>=5e-4, TF_p_value_HER2>=5e-4)
other_sig_motif<-unique(all$motif)
length(other_sig_motif)
# 7 new criteria

LumB_HER2Neg_unique <- LumB_HER2Neg %>% dplyr::filter(!motif %in% other_sig_motif)  %>% unique()
length(unique(LumB_HER2Neg_unique$motif))
fwrite(LumB_HER2Neg_unique,"results/tfgwas/tfgwas_significant_LumB_HER2Neg_unique.csv", row.names = FALSE)

enriched.motifs1$motif.name<-NULL
enriched_motifs_list_LumB_HER2Neg_unique <-left_join(LumB_HER2Neg_unique, enriched.motifs1, by='motif') %>% unique()

enriched_motifs_list_LumB_HER2Neg_unique$cluster <- unique(combined_peaks_gene_all$cluster)[1]

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
  enriched_motifs_filtered <- left_join(LumB_HER2Neg_unique, enriched_motifs, by='motif') %>% unique()
  enriched_motifs_filtered$cluster <- cluster

  # Store the filtered results in the list with the cluster name as the key
  enriched_motifs_list_LumB_HER2Neg_unique <- rbind(enriched_motifs_list_LumB_HER2Neg_unique,enriched_motifs_filtered)

  # Print a message to indicate progress
  cat("Processed cluster:", cluster, "\n")
}

enriched_motifs_list_LumB_HER2Neg_unique$pvalue <- ifelse(enriched_motifs_list_LumB_HER2Neg_unique$pvalue==0,1e-300,
                                                  ifelse(enriched_motifs_list_LumB_HER2Neg_unique$pvalue==1,0.99999999999999999999999,enriched_motifs_list_LumB_HER2Neg_unique$pvalue))

# Save
write.csv(enriched_motifs_list_LumB_HER2Neg_unique, file = "results/singlecell/enriched_motifs_list_uniqDAP_LumB_HER2NegUniqSignificant.csv", row.names = FALSE)
write.csv(enriched_motifs_list_LumB_HER2Neg_unique, file = "results/singlecell/enriched_motifs_list_uniqDAP_LumB_HER2Neg_7motifs.csv", row.names = FALSE)

enriched_motifs_list_LumB_HER2Neg_unique<-fread("results/singlecell/enriched_motifs_list_uniqDAP_LumB_HER2NegUniqSignificant.csv",header = TRUE)

# create a matrix
# Get all unique motif names across all clusters
all_motifsID <- unique(enriched_motifs_list_LumB_HER2Neg_unique$motif)
all_motifs <- motif_id_list_final[motif %in% all_motifsID]
all_motifs$motif.name <- make.unique(all_motifs$motif.name)
all_motifs <- all_motifs$motif.name

# Get all cluster names
all_clusters <- unique(enriched_motifs_list_LumB_HER2Neg_unique$cluster)

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
  enriched_data <- enriched_motifs_list_LumB_HER2Neg_unique[enriched_motifs_list_LumB_HER2Neg_unique$cluster==cluster1]

  # Add p-values to the matrix
  for (i in 1:nrow(enriched_data)) {
    motif <- enriched_data$motif.name[i]
    neglog10pvalue <- -log10(enriched_data$pvalue[i])  # Replace "pvalue" with the actual column name for p-values
    motif_matrix[motif,cluster1] <-neglog10pvalue
  }
}

# Convert the matrix to a data frame for better readability
motif_df<- as.data.frame(motif_matrix)

motif_df<-dplyr::select(motif_df,'LS-1','LASP-2','LHr','Bs','Fb','Endo','T','B','Myl','Pc')

motif_df<-motif_df%>%arrange(desc(LHr),desc(`LS-1`),desc(Bs))
library(pheatmap)
library(ComplexHeatmap)
library(circlize)

TF_list<-rownames(motif_df)
cell_list<-colnames(motif_df)
motif_heatmap<-SingleRasterMap(
  data = t(motif_df),
  raster = TRUE,
  cell.order = cell_list,
  feature.order = rev(TF_list),
  colors = colorRampPalette(c('lightyellow',"purple4","darkblue"))(300),   # Custom color mapping
  disp.min = 0,      # Minimum display value
  disp.max = 100,    # Maximum display value
  limits = c(0, 100) # Set explicit limits
)
num_cols <- ncol(motif_df)
line_positions <- seq(1.5, num_cols - 0.5, by = 1)  # Positions between columns
motif_heatmap<-motif_heatmap +
  geom_vline(xintercept = line_positions, color = "purple4", size = 0.25)+
  #   breaks = rownames(motif_label),                 # Only label significant genes
  labs(fill = "-log10(pval)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1,size = 12),
        axis.text.y = element_text(angle = 0,size = 5)  # Adjust y-axis font size
  )
motif_heatmap
ggsave("results/singlecell/Motif_enricheded_heatmap_uniDAP_LumB_HER2Neg_unique_signMotif.jpg",
       motif_heatmap, width = 4, height = 5, dpi = 150)

# ---------------------heatmap ---------------#
motif_id_name<-fread('results/singlecell/jaspar2022_motifs_727.csv', header = TRUE)

# for the motif score heatmap
motif_matrix <- as.matrix(motif_df)

library(ComplexHeatmap)
library(dendextend)
set.seed(2025)
heatmap_obj_orig <-Heatmap(motif_matrix,
                           name = "-log10 pval",          # Title of legend
                           row_names_gp = gpar(fontsize = 7),
                           row_names_side = "left",
                           column_names_centered = TRUE,
                           width = unit(80, "mm"),
                           col = circlize::colorRamp2(c(0,60), c("#483d8b","red1")),
                           rect_gp = gpar(col = "grey99"), # Set cell border color (instead of border_color)
                           show_row_dend = FALSE,         # Hide row dendrogram to avoid km labels
                           show_column_dend = FALSE,
                           # column_split = 4,                 # Ensure 3 clusters without labels
                           row_title = NULL ,
                           column_title = NULL
)
heatmap_obj_orig
pdf("results/tfgwas/heatmap_enrichedmotif_uniqMotif_LumB_HER2Neg_orginal.pdf", width = 5, height =7)
draw(heatmap_obj_orig)
dev.off()

# ---------heatmap with specific motif sequence----------------
motif_df<- as.data.frame(motif_matrix)

colSums(motif_matrix >= -log10(0.05 / 727))
# LHr LS-1 LASP-2   Bs    T    B  Myl Endo   Pc   Fb
colSums(motif_matrix >= 10)
# LHr LS-1 LASP-2   Bs    T    B  Myl Endo   Pc   Fb

motif_df <- motif_df[order(-motif_df$`LHr`),]
rownames(motif_df)
motif_list <- c("ESR1","MTF1","POU6F1","RXRG","FOS","VENTX","NR6A1", "RARB","GSC2" )%>% unique()

# create the new data frame
motif_df<- as.data.frame(motif_matrix)
head(motif_df)
# Reorder motif_df rows based on motif_list
motif_df <- motif_df[match(motif_list, rownames(motif_df)), ]
motif_df <- motif_df %>% select('LS-1','LASP-2','LHr','Bs','Fb','Endo','T','B','Myl','Pc')
motif_matrix <- as.matrix(motif_df)

library(ComplexHeatmap)
library(circlize)

# Create the heatmap
pdf("results/tfgwas/heatmap_enrichedmotif_uniqMotif_LumB_HER2Neg_9motifs.pdf", width = 5, height =2)
Heatmap(motif_matrix,
        col =  circlize::colorRamp2(c(0, 40), c("#483d8b", "red1")),
        row_names_gp = gpar(fontsize = 7.5),
        cluster_rows = FALSE,  # Preserve the order from motif_list
        cluster_columns = FALSE,  # Cluster columns (or set to FALSE if you want to preserve column order)
        show_row_names = TRUE,  # Show motif names
        show_column_names = TRUE,  # Show column names
        row_names_side = "left",  # Place row names on the left
        heatmap_legend_param = list(title = "-log10 pval"),
        rect_gp = gpar(col = "grey99"),)  # Customize legend title

dev.off()

############################ GSEA bubble plots ###########################################
# Load libraries
library(tidyverse)
library(ggplot2)
library(paletteer)
library(tidyverse)
library(ggplot2)
library(paletteer)
dat_LumB_HER2Neg<-fread("results/singlecell/enriched_motifs_list_uniqDAP_LumB_HER2NegUniqSignificant.csv",header = TRUE)
dat_LumB_HER2Neg$cluster <- gsub('LS-1','LASP-1',dat_LumB_HER2Neg$cluster )
dat_LumB_HER2Neg$cluster <- gsub('LS-2','LASP-2',dat_LumB_HER2Neg$cluster )
dat_LumB_HER2Neg$cluster <- gsub('LHr','LHS',dat_LumB_HER2Neg$cluster )
dat_LumB_HER2Neg$cluster <- gsub('Bs','B-Myo',dat_LumB_HER2Neg$cluster )
#
#
dat_LumB_HER2Neg$cluster <- factor(dat_LumB_HER2Neg$cluster, levels = c('LASP-1','LASP-2','LHS','B-Myo','Fb','Endo','T','B','Myl','Pc'))
dat_LumB_HER2Neg$avg.exp.scaled <- ifelse(is.na(dat_LumB_HER2Neg$avg.exp.scaled), 0, dat_LumB_HER2Neg$avg.exp.scaled)

jpeg('results/tfgwas/GSEAplot_LumB_HER2Neg_uniquemotifs.jpg',width = 15, height = 6, units = 'in', res = 500)
ggplot1 <- ggplot(dat_LumB_HER2Neg, aes(x = fold.enrichment, y = reorder(motif.name, fold.enrichment)) )+
  geom_point(aes(size = avg.exp.scaled, color = neglog10pval)) +
  coord_cartesian(clip = "off") +
  scale_color_gradient(low = "#F88586", high = "#FC0D0F", name = "Motif enrichment\n-log10(P value)", limits = c(-log10(0.05), 20)) +
  theme_bw() +
  scale_x_continuous(breaks = 1) +
  scale_size_continuous(name = "Average gene expression\n(normalized)") +
  labs(x = "Fold enrichment\n(Over-representation, enriched by unique peaks in each cell type)", y = "TF motif Name", title = "") +
  facet_grid(family ~ cluster, scales = "free_y",space = "free_y",
             labeller = label_wrap_gen(width = 45)) + # Wrap family names
  theme(axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        axis.text.y = element_text(size = 10, color = 'black'),
        legend.position = "bottom",
        strip.text.y = element_text(size = 10, angle = 0), # Horizontal family names
        strip.text.x = element_text(size = 12) , # Style cluster labels (horizontal)
        panel.spacing.y = unit(0.2, "lines"))
ggplot1
dev.off()

dat_LumB_HER2Neg_motif <- unique(select(dat_LumB_HER2Neg, motif, motif.name, family))
dat_LumB_HER2Neg_motif_group <- dat_LumB_HER2Neg_motif %>%
  group_by(family) %>%
  mutate(count = n()) %>%
  select(family, count) %>%
  ungroup() %>%
  unique()

library(stringr)
jpeg('results/tfgwas/GWASMotif01.LumB_HER2Neg_motif_count.jpg',
     width=6, height = 5, units = 'in', res = 250)
ggplot(dat_LumB_HER2Neg_motif_group, aes(x = count, y = reorder(family, count), fill = count)) +
  geom_bar(stat = "identity", fill = "#FF4500") +
  labs(x = "No. of TF motifs", y = "TF motif family", title = "") +
  scale_y_discrete(labels = function(x) str_wrap(x, width = 40)) +
  theme_classic() +
  theme(legend.position = "right") +
  scale_x_continuous(breaks = seq(0, 3, 1), limits = c(0, 3))
dev.off()
##########################################################################################
#                                 feature plot                                           #
##########################################################################################
DefaultAssay(atac_subset_qc) <- "chromvar"

LumB_HER2Neg_unique<-fread("results/tfgwas/tfgwas_significant_LumB_HER2Neg_unique.csv", header=TRUE)
LumB_HER2Neg_unique<-list(unique(LumB_HER2Neg_unique$motif))
set.seed(4242025)
atac_subset_qc <- AddModuleScore(
  object = atac_subset_qc,
  features = LumB_HER2Neg_unique,
  ctrl = 15,
  name = 'LumB_HER2Neg_unique_motif'
)

jpeg('results/tfgwas/GWASMotif01.all_LumB_HER2Neg_uniqueMotifs_Addmodulescore.jpg',width=6, height = 5, units = 'in', res = 150)
p<-FeaturePlot(
  object = atac_subset_qc,
  features = "LumB_HER2Neg_unique_motif1",
  min.cutoff = 'q5',
  max.cutoff = 'q95',
  pt.size = 0.1,
  label = TRUE,
  cols = c( "grey90","orange")
)+
  xlab("UMAP 1") +
  ylab("UMAP 2")+
  theme(legend.position = "right")

# Extract UMAP coordinates and cluster assignments
umap_data <- as.data.frame(atac_subset_qc@reductions$umap.rna@cell.embeddings)
umap_data$cluster <- atac_subset_qc@meta.data$celltype_RNA

# Filter for the clusters of interest
clusters_to_encircle <- c("LHS", "LASP-1", "LASP-2", "B-Myo")
umap_subset <- umap_data[umap_data$cluster %in% clusters_to_encircle, ]

points_matrix <- as.matrix(umap_subset[, c("rnaUMAP_1", "rnaUMAP_2")])
hull <- inla.nonconvex.hull(points_matrix, convex = -0.0021)
# Adjust convex for tightness (more negative for less smoothing)
segments <- data.frame(
  x = hull$loc[hull$idx[, 1], 1],
  y = hull$loc[hull$idx[, 1], 2],
  xend = hull$loc[hull$idx[, 2], 1],
  yend = hull$loc[hull$idx[, 2], 2]
)

p <- p + geom_segment(
  data = segments,
  aes(x = x, y = y, xend = xend, yend = yend),
  color = "black",  # Outline color
  linetype = "dotted",  # Dotted line
  size = 0.5  # Line thickness
)
print(p)
dev.off()

jpeg('results/tfgwas/GWASMotif01_LumB_HER2Neg_uniqueMotifs_Addmodulescore.jpg',width=10, height = 5, units = 'in', res = 150)
a<-FeaturePlot(
  object = subset(atac_subset_qc, idents=c('LHr','LS-1','LASP-2','Bs')),
  features = "LumB_HER2Neg_unique_motif1",
  min.cutoff = 'q5',
  max.cutoff = 'q95',
  pt.size = 0.1,
  label = TRUE)
b<-FeaturePlot(
  object = subset(atac_subset_qc, idents=c('T','B','Myl','Endo','Pc','Fb')),
  features = "LumB_HER2Neg_unique_motif1",
  min.cutoff = 'q5',
  max.cutoff = 'q95',
  pt.size = 0.1,
  label = TRUE
)
a &NoLegend()|b
dev.off()

df <- atac_subset_qc@meta.data
celltype_order <- c("LASP-1", "LASP-2", "LHS", "B-Myo", "Fb", "Endo", "T", "B", "Myl", "Pc")

# Convert celltype_RNA to a factor with specified order
levels(df$celltype_RNA)
df$celltype_RNA <- factor(df$celltype_RNA, levels = celltype_order)

# Perform statistical comparison (e.g., Kruskal-Wallis test for non-parametric data)
kw<-kruskal.test(LumB_HER2Neg_unique_motif1 ~ celltype_RNA, data = atac_subset_qc@meta.data)

df_main <- df %>% filter(celltype_RNA %in% c("LHS", "LASP-1", "LASP-2", "B-Myo"))
# Step 2: Subset for other subtypes (Fb, Endo, T, B, Myl, Pc)
df_other <- df %>% filter(celltype_RNA %in% c("Fb", "Endo", "T", "B", "Myl", "Pc"))

library(dunn.test)
dunn_test_main <- dunn.test(
  x = df_main$LumB_HER2Neg_unique_motif1,
  g = df_main$celltype_RNA,
  method = "bh" )
dunn_df_main <- data.frame(
  Comparison = dunn_test_main$comparisons,
  Z = dunn_test_main$Z,
  P_unadj = dunn_test_main$P,
  P_adj = dunn_test_main$P.adjusted)
my_comparisons1 <- list(
  c("LASP-1", "LASP-2"), c("LASP-1", "LHS"), c("LASP-1", "B-Myo"),
  c("LASP-2", "LHS"), c("LASP-2", "B-Myo"), c("LHS", "B-Myo")
)
##  3. Normalise Dunn comparison strings (both directions)
dunn_df_main <- dunn_df_main %>%
  tidyr::separate(Comparison, into = c("g1","g2"), sep = " - ", remove = FALSE) %>%
  mutate(
    key_fwd = paste(g1, g2, sep = "-"),
    key_rev = paste(g2, g1, sep = "-"))

##  4. Map to the exact order of my_comparisons
target_keys <- sapply(my_comparisons1, paste, collapse = "-")

dunn_df_main <- dunn_df_main %>%
  mutate(
    match_key = case_when(
      key_fwd %in% target_keys ~ key_fwd,
      key_rev %in% target_keys ~ key_rev,
      TRUE                     ~ NA_character_
    )
  ) %>%
  filter(!is.na(match_key)) %>%                     # keep only the 15 requested pairs
  arrange(match(match_key, target_keys))            # order exactly as target_keys
# Add significance stars
dunn_df_main <- dunn_df_main %>%
  mutate(signif = case_when(
    is.na(P_adj)                  ~ "NA",
    P_adj == 0                    ~ "***",
    P_adj <= 0.001                ~ "***",
    P_adj <= 0.01                 ~ "**",
    P_adj <= 0.05                 ~ "*",
    TRUE                          ~ "ns"
  ))

# Plot
p1 <- ggboxplot(df_main, x = "celltype_RNA", y = "LumB_HER2Neg_unique_motif1",
                color = "celltype_RNA", palette = "jco", legend = "none") +
  stat_compare_means(method = "kruskal.test", label.y = 3.5, size = 3) +
  geom_signif(comparisons = my_comparisons1,
              annotations = dunn_df_main$signif,
              y_position = seq(1.6, by = 0.2, length.out = 6),  # adjust spacing
              textsize = 3,
              vjust = 0.5) +
  labs(x = "Cell type", y = "AddModuleScore by LumB/Her2Neg unique TF motifs", title = "") +
  theme(axis.title.y = element_text(size = 8))
jpeg('results/tfgwas/Boxplot_LumB_HER2Neg_uniqueMotifs_Addmodulescore_part1.jpg',width=3.5, height = 5, units = 'in', res = 150)
print(p1)
dev.off()

write.csv(dunn_df_main,'results/tfgwas/dunn_df_main_LumB_HER2Neg.csv',row.names = FALSE)

# other cell type
dunn_test_other <- dunn.test(
  x = df_other$LumB_HER2Neg_unique_motif1,
  g = df_other$celltype_RNA,
  method = "bh" )
dunn_df_other <- data.frame(
  Comparison = dunn_test_other$comparisons,
  Z = dunn_test_other$Z,
  P_unadj = dunn_test_other$P,
  P_adj = dunn_test_other$P.adjusted)
##  2. Desired comparison order (my_comparisons2)
my_comparisons2 <- list(
  c("Fb","Endo"), c("Fb","T"),   c("Fb","B"),   c("Fb","Myl"), c("Fb","Pc"),
  c("Endo","T"),  c("Endo","B"), c("Endo","Myl"),c("Endo","Pc"),
  c("T","B"),     c("T","Myl"),  c("T","Pc"),
  c("B","Myl"),   c("B","Pc"),   c("Myl","Pc")
)
##  3. Normalise Dunn comparison strings (both directions)
dunn_df_other <- dunn_df_other %>%
  tidyr::separate(Comparison, into = c("g1","g2"), sep = " - ", remove = FALSE) %>%
  mutate(
    key_fwd = paste(g1, g2, sep = "-"),
    key_rev = paste(g2, g1, sep = "-"))

##  4. Map to the exact order of my_comparisons2
target_keys <- sapply(my_comparisons2, paste, collapse = "-")

dunn_df_other <- dunn_df_other %>%
  mutate(
    match_key = case_when(
      key_fwd %in% target_keys ~ key_fwd,
      key_rev %in% target_keys ~ key_rev,
      TRUE                     ~ NA_character_
    )
  ) %>%
  filter(!is.na(match_key)) %>%                     # keep only the 15 requested pairs
  arrange(match(match_key, target_keys))            # order exactly as target_keys

##  5. (Optional) Re-create signif column – safe for NA
dunn_df_other <- dunn_df_other %>%
  mutate(
    signif = case_when(
      is.na(P_adj)               ~ "NA",
      P_adj == 0                 ~ "***",
      P_adj <= 0.001             ~ "***",
      P_adj <= 0.01              ~ "**",
      P_adj <= 0.05              ~ "*",
      TRUE                       ~ "ns") )
##  6. Plot
jpeg('results/tfgwas/Boxplot_LumB_HER2Neg_uniqueMotifs_Addmodulescore_part2.jpg',
     width = 4.5, height = 6, units = 'in', res = 150)
p2 <- ggboxplot(df_other,
                x = "celltype_RNA",
                y = "LumB_HER2Neg_unique_motif1",
                color = "celltype_RNA",
                palette = "YlGnBu",
                legend = "none") +

  ## overall Kruskal-Wallis (small font)
  stat_compare_means(method = "kruskal.test",
                     label.y = 7,
                     size = 3) +
  ## 15 pairwise Dunn stars (BH-adjusted)
  geom_signif(comparisons = my_comparisons2,
              annotations = dunn_df_other$signif,
              y_position = seq(2, by = 0.3, length.out = length(my_comparisons2)),
              textsize = 3,
              vjust = 0.5) +
  labs(x = "Cell type",
       y = "AddModuleScore by LumB/Her2Neg unique TF motifs",
       title = "") +
  theme(axis.title.y = element_text(size = 8))
print(p2)
dev.off()
write.csv(dunn_df_other,'results/tfgwas/dunn_df_other_LumB_HER2Neg.csv',row.names = FALSE)

#################################################################################################
#                          Unique Motif in HER2 significant GWAS motif                         #
#################################################################################################
motif_id_list_final<-fread('results/singlecell/jaspar2022_motifs_727.csv', header = TRUE)
HER2<-fread("results/tfgwas/tfgwas_results_HER2_ALL.csv", header = TRUE)
HER2<- left_join(select(motif_id_list_final,motif,motif.name), HER2, by=c('motif','motif.name')) %>%
  dplyr::filter(TF_p_value < 0.05/727)
length(unique(HER2$motif))
#38 previous Bonferroni corrected

TNBC<-fread("results/tfgwas/tfgwas_results_TNBC_ALL.csv", header = TRUE)
LumB<-fread("results/tfgwas/tfgwas_results_LumB_ALL.csv", header = TRUE)
LumB_HER2Neg<-fread("results/tfgwas/tfgwas_results_LumB_HER2Neg_ALL.csv", header = TRUE)
LumA<-fread("results/tfgwas/tfgwas_results_LumA_ALL.csv", header = TRUE)

all<-rbind(TNBC, LumA)
all<-rbind(all, LumB)
all<-rbind(all, LumB_HER2Neg)
all<-all %>% dplyr::filter(TF_p_value < 0.05/727)%>% unique()
#                               TF_p_value_LumB>=5e-4, TF_p_value_LumB_HER2Neg>=5e-4, TF_p_value_LumA>=5e-4)
other_sig_motif<-unique(all$motif)
length(other_sig_motif)
# # 11 new criteria
# #518
#
HER2_unique <- HER2 %>% dplyr::filter(!motif %in% other_sig_motif)  %>% unique()
length(unique(HER2_unique$motif))
fwrite(HER2_unique,"results/tfgwas/tfgwas_significant_HER2_unique.csv", row.names = FALSE)

enriched.motifs1$motif.name<-NULL
enriched_motifs_list_HER2_unique <-left_join(HER2_unique, enriched.motifs1, by='motif') %>% unique()

enriched_motifs_list_HER2_unique$cluster <- unique(combined_peaks_gene_all$cluster)[1]

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
  enriched_motifs_filtered <- left_join(HER2_unique, enriched_motifs, by='motif') %>% unique()
  enriched_motifs_filtered$cluster <- cluster

  # Store the filtered results in the list with the cluster name as the key
  enriched_motifs_list_HER2_unique <- rbind(enriched_motifs_list_HER2_unique,enriched_motifs_filtered)

  # Print a message to indicate progress
  cat("Processed cluster:", cluster, "\n")
}

enriched_motifs_list_HER2_unique$pvalue <- ifelse(enriched_motifs_list_HER2_unique$pvalue==0,1e-300,
                                                  ifelse(enriched_motifs_list_HER2_unique$pvalue==1,0.99999999999999999999999,enriched_motifs_list_HER2_unique$pvalue))

# Save
write.csv(enriched_motifs_list_HER2_unique, file = "results/singlecell/enriched_motifs_list_uniqDAP_HER2.csv", row.names = FALSE)

enriched_motifs_list_HER2_unique<-fread("results/singlecell/enriched_motifs_list_uniqDAP_HER2.csv",header = TRUE)

# create a matrix
# Get all unique motif names across all clusters
all_motifsID <- unique(enriched_motifs_list_HER2_unique$motif)
all_motifs <- motif_id_list_final[motif %in% all_motifsID]
all_motifs$motif.name <- make.unique(all_motifs$motif.name)
all_motifs <- all_motifs$motif.name

# Get all cluster names
all_clusters <- unique(enriched_motifs_list_HER2_unique$cluster)

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
  enriched_data <- enriched_motifs_list_HER2_unique[enriched_motifs_list_HER2_unique$cluster==cluster1]

  # Add p-values to the matrix
  for (i in 1:nrow(enriched_data)) {
    motif <- enriched_data$motif.name[i]
    neglog10pvalue <- -log10(enriched_data$pvalue[i])  # Replace "pvalue" with the actual column name for p-values
    motif_matrix[motif,cluster1] <-neglog10pvalue
  }
}

# Convert the matrix to a data frame for better readability
motif_df<- as.data.frame(motif_matrix)

motif_df<-dplyr::select(motif_df,'LS-1','LASP-2','LHr','Bs','Fb','Endo','T','B','Myl','Pc')

motif_df<-motif_df%>%arrange(desc(`LASP-2`),desc(Endo))
library(pheatmap)
library(ComplexHeatmap)
library(circlize)

TF_list<-rownames(motif_df)
cell_list<-colnames(motif_df)
motif_heatmap<-SingleRasterMap(
  data = t(motif_df),
  raster = TRUE,
  cell.order = cell_list,
  feature.order = rev(TF_list),
  colors = colorRampPalette(c('lightyellow',"purple4","darkblue"))(300),   # Custom color mapping
  disp.min = 0,      # Minimum display value
  disp.max = 100,    # Maximum display value
  limits = c(0, 100) # Set explicit limits
)
num_cols <- ncol(motif_df)
line_positions <- seq(1.5, num_cols - 0.5, by = 1)  # Positions between columns
motif_heatmap<-motif_heatmap +
  geom_vline(xintercept = line_positions, color = "purple4", size = 0.25)+
  #   breaks = rownames(motif_label),                 # Only label significant genes
  labs(fill = "-log10(pval)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1,size = 12),
        axis.text.y = element_text(angle = 0,size = 5)  # Adjust y-axis font size
  )
motif_heatmap
ggsave("results/singlecell/Motif_enricheded_heatmap_uniDAP_HER2_unique_signMotif.jpg",
       motif_heatmap, width = 4, height = 5, dpi = 150)

# ---------------------heatmap ---------------#
motif_id_name<-fread('results/singlecell/jaspar2022_motifs_727.csv', header = TRUE)

# for the motif score heatmap
motif_matrix <- as.matrix(motif_df)

library(ComplexHeatmap)
library(dendextend)
set.seed(2025)
heatmap_obj_orig <-Heatmap(motif_matrix,
                           name = "-log10 pval",          # Title of legend
                           row_names_gp = gpar(fontsize = 7),
                           row_names_side = "left",
                           column_names_centered = TRUE,
                           width = unit(50, "mm"),
                           col = circlize::colorRamp2(c(0,50), c("#483d8b","red1")),
                           rect_gp = gpar(col = "grey99"), # Set cell border color (instead of border_color)
                           show_row_dend = FALSE,         # Hide row dendrogram to avoid km labels
                           show_column_dend = FALSE,
                           # column_split = 4,                 # Ensure 3 clusters without labels
                           row_title = NULL ,
                           column_title = NULL
)
heatmap_obj_orig
pdf("results/tfgwas/heatmap_enrichedmotif_uniqMotif_HER2_orginal.pdf", width = 5, height =7)
draw(heatmap_obj_orig)
dev.off()

# ---------heatmap with specific motif sequence----------------
motif_df<- as.data.frame(motif_matrix)

colSums(motif_matrix >= -log10(0.05 / 727))
# Endo LS-1   Pc   Fb   Bs  LHr    T    B LASP-2  Myl
colSums(motif_matrix >= 10)
# Endo LS-1   Pc   Fb   Bs  LHr    T    B LASP-2  Myl

motif_df <- motif_df[order(-motif_df$`Endo`),]
rownames(motif_df)
motif_list <- c("PDX1","DRGX","EMX1","HOXA2","HOXB5","HOXA5",
                "ATF7","FOSL2::JUNB","FOSL2::JUN","JUNB","GMEB2","PAX9","PAX1","MYBL2")%>% unique()
motif_df$motif <- rownames(motif_df)
motif_df <- motif_df[!motif_df$motif %in% motif_list,]
motif_df <- motif_df[order(-motif_df$Fb),]
motif_df <- motif_df[order(-motif_df$`LASP-2`),]

# create the new data frame
motif_df<- as.data.frame(motif_matrix)
head(motif_df)
# Reorder motif_df rows based on motif_list
motif_df <- motif_df[match(motif_list, rownames(motif_df)), ]
motif_df <- motif_df %>% select('LS-1','LASP-2','LHr','Bs','Fb','Endo','T','B','Myl','Pc')
motif_matrix <- as.matrix(motif_df)

library(ComplexHeatmap)
library(circlize)

# Create the heatmap
pdf("results/tfgwas/heatmap_enrichedmotif_uniqMotif_HER2_14motifs.pdf", width = 5, height =2.5)
Heatmap(motif_matrix,
        col =  circlize::colorRamp2(c(0, 40), c("#483d8b", "red1")),
        row_names_gp = gpar(fontsize = 7.5),
        cluster_rows = FALSE,  # Preserve the order from motif_list
        cluster_columns = FALSE,  # Cluster columns (or set to FALSE if you want to preserve column order)
        show_row_names = TRUE,  # Show motif names
        show_column_names = TRUE,  # Show column names
        row_names_side = "left",  # Place row names on the left
        heatmap_legend_param = list(title = "-log10 pval"),
        rect_gp = gpar(col = "grey99"),)  # Customize legend title

dev.off()

############################ GSEA bubble plots ###########################################
# Load libraries
library(tidyverse)
library(ggplot2)
library(paletteer)
library(tidyverse)
library(ggplot2)
library(paletteer)
dat_HER2<-fread("results/singlecell/enriched_motifs_list_uniqDAP_HER2.csv",header = TRUE)

dat_HER2$neglog10pval<--log10(dat_HER2$pvalue)
## scaled average RNA expression of each motif's TF gene per cell type (DotPlot avg.exp.scaled);
## for a heterodimer motif (A::B) the lower of the two partners is used
tf_gene_expression <- function(dat) {
  obj <- readRDS("results/singlecell/atac_motif.rds")
  obj <- subset(obj, cells = fread("results/singlecell/final_cell_barcodes.csv", header = TRUE)$x)
  DefaultAssay(obj) <- "RNA"
  ids <- c("LASP-1", "Fb", "LHS", "LASP-2", "B-Myo", "Endo", "T", "Pc", "Myl", "B"); names(ids) <- levels(obj)
  obj <- RenameIdents(obj, ids)
  Idents(obj) <- factor(Idents(obj), levels = c("LASP-1", "LASP-2", "LHS", "B-Myo", "Fb", "Endo", "T", "B", "Myl", "Pc"))
  motifs <- rev(levels(reorder(dat$motif.name, dat$observed)))
  genes  <- unique(unlist(strsplit(motifs, "::")))
  p <- DotPlot(obj, features = rev(genes)) + coord_flip()
  m <- p$data %>% dplyr::select(features.plot, id, avg.exp.scaled) %>%
    pivot_wider(names_from = id, values_from = avg.exp.scaled)
  names(m)[1] <- "feature"
  for (pair in strsplit(motifs[grepl("::", motifs)], "::")) {
    if (all(pair %in% m$feature)) {
      r1 <- m[m$feature == pair[1], -1]; r2 <- m[m$feature == pair[2], -1]
      m <- bind_rows(m, data.frame(feature = paste(pair, collapse = "::"), Map(min, r1, r2), check.names = FALSE))
    }
  }
  names(m)[2:3] <- c("LS-1", "LS-2")
  m %>% rename(motif.name = feature) %>%
    pivot_longer(cols = -motif.name, names_to = "cluster", values_to = "avg.exp.scaled")
}
percent_expr <- tf_gene_expression(dat_HER2)
dat_HER2 <- left_join(dat_HER2, percent_expr, by=c('motif.name','cluster'))
motif_data<-fread('results/singlecell/jaspar2022_motif_list.csv', header=TRUE)
dat_HER2 <- left_join(dat_HER2, motif_data, by=c('motif.name','motif'))

dat_HER2$cluster <- gsub('LS-1','LASP-1',dat_HER2$cluster )
dat_HER2$cluster <- gsub('LS-2','LASP-2',dat_HER2$cluster )
dat_HER2$cluster <- gsub('LHr','LHS',dat_HER2$cluster )
dat_HER2$cluster <- gsub('Bs','B-Myo',dat_HER2$cluster )
dat_HER2$cluster <- factor(dat_HER2$cluster, levels = c('LASP-1','LASP-2','LHS','B-Myo','Fb','Endo','T','B','Myl','Pc'))
dat_HER2$avg.exp.scaled <- ifelse(is.na(dat_HER2$avg.exp.scaled), 0, dat_HER2$avg.exp.scaled)

jpeg('results/tfgwas/GSEAplot_HER2_uniquemotifs.jpg',width = 15, height = 6, units = 'in', res = 500)
ggplot1 <- ggplot(dat_HER2, aes(x = fold.enrichment, y = reorder(motif.name, fold.enrichment)) )+
  geom_point(aes(size = avg.exp.scaled, color = neglog10pval)) +
  coord_cartesian(clip = "off") +
  scale_color_gradient(low = "#F88586", high = "#FC0D0F", name = "Motif enrichment\n-log10(P value)", limits = c(-log10(0.05), 30)) +
  theme_bw() +
  scale_x_continuous(breaks = 1) +
  scale_size_continuous(name = "Average gene expression\n(normalized)") +
  labs(x = "Fold enrichment\n(Over-representation, enriched by unique peaks in each cell type)", y = "TF motif Name", title = "") +
  facet_grid(family ~ cluster, scales = "free_y",space = "free_y",
             labeller = label_wrap_gen(width = 45)) + # Wrap family names
  theme(axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        axis.text.y = element_text(size = 10, color = 'black'),
        legend.position = "bottom",
        strip.text.y = element_text(size = 10, angle = 0), # Horizontal family names
        strip.text.x = element_text(size = 12) , # Style cluster labels (horizontal)
        panel.spacing.y = unit(0.2, "lines"))
ggplot1
dev.off()

dat_HER2_motif <- unique(select(dat_HER2, motif, motif.name, family))
dat_HER2_motif_group <- dat_HER2_motif %>%
  group_by(family) %>%
  mutate(count = n()) %>%
  select(family, count) %>%
  ungroup() %>%
  unique()

library(stringr)
jpeg('results/tfgwas/GWASMotif01.HER2_motif_count.jpg',
     width=6, height = 5, units = 'in', res = 250)
ggplot(dat_HER2_motif_group, aes(x = count, y = reorder(family, count), fill = count)) +
  geom_bar(stat = "identity", fill = "#FF4500") +
  labs(x = "No. of TF motifs", y = "TF motif family", title = "") +
  scale_y_discrete(labels = function(x) str_wrap(x, width = 40)) +
  theme_classic() +
  theme(legend.position = "right") +
  scale_x_continuous(breaks = seq(0, 3, 1), limits = c(0, 3))
dev.off()
##########################################################################################
#                                 feature plot                                           #
##########################################################################################
HER2_unique<-fread("results/tfgwas/tfgwas_significant_HER2_unique.csv",header = TRUE)
HER2_unique<-list(unique(HER2_unique$motif))
DefaultAssay(atac_subset_qc) <- "chromvar"
set.seed(4242025)
atac_subset_qc <- AddModuleScore(
  object = atac_subset_qc,
  features = HER2_unique,
  ctrl = 15,
  name = 'HER2_unique_motif'
)

jpeg('results/tfgwas/GWASMotif01.all_HER2_uniqueMotifs_Addmodulescore.jpg',width=6, height = 5, units = 'in', res = 150)
p<-FeaturePlot(
  object = atac_subset_qc,
  features = "HER2_unique_motif1",
  min.cutoff = 'q5',
  max.cutoff = 'q95',
  pt.size = 0.1,
  label = TRUE,
  cols = c( "grey90","orange")
)+
  xlab("UMAP 1") +
  ylab("UMAP 2")+
  theme(legend.position = "right")

# Extract UMAP coordinates and cluster assignments
umap_data <- as.data.frame(atac_subset_qc@reductions$umap.rna@cell.embeddings)
umap_data$cluster <- atac_subset_qc@meta.data$celltype_RNA

# Filter for the clusters of interest
clusters_to_encircle <- c("LHS", "LASP-1", "LASP-2", "B-Myo")
umap_subset <- umap_data[umap_data$cluster %in% clusters_to_encircle, ]

points_matrix <- as.matrix(umap_subset[, c("rnaUMAP_1", "rnaUMAP_2")])
hull <- inla.nonconvex.hull(points_matrix, convex = -0.0021)
# Adjust convex for tightness (more negative for less smoothing)
segments <- data.frame(
  x = hull$loc[hull$idx[, 1], 1],
  y = hull$loc[hull$idx[, 1], 2],
  xend = hull$loc[hull$idx[, 2], 1],
  yend = hull$loc[hull$idx[, 2], 2]
)

p <- p + geom_segment(
  data = segments,
  aes(x = x, y = y, xend = xend, yend = yend),
  color = "black",  # Outline color
  linetype = "dotted",  # Dotted line
  size = 0.5  # Line thickness
)
print(p)
dev.off()

jpeg('results/tfgwas/GWASMotif01_HER2_uniqueMotifs_Addmodulescore.jpg',width=10, height = 5, units = 'in', res = 150)
a<-FeaturePlot(
  object = subset(atac_subset_qc, idents=c('LHr','LS-1','LASP-2','Bs')),
  features = "HER2_unique_motif1",
  min.cutoff = 'q5',
  max.cutoff = 'q95',
  pt.size = 0.1,
  label = TRUE)
b<-FeaturePlot(
  object = subset(atac_subset_qc, idents=c('T','B','Myl','Endo','Pc','Fb')),
  features = "HER2_unique_motif1",
  min.cutoff = 'q5',
  max.cutoff = 'q95',
  pt.size = 0.1,
  label = TRUE
)
a &NoLegend()|b
dev.off()

df <- atac_subset_qc@meta.data
celltype_order <- c("LASP-1", "LASP-2", "LHS", "B-Myo", "Fb", "Endo", "T", "B", "Myl", "Pc")

# Convert celltype_RNA to a factor with specified order
levels(df$celltype_RNA)
df$celltype_RNA <- factor(df$celltype_RNA, levels = celltype_order)

# Perform statistical comparison (e.g., Kruskal-Wallis test for non-parametric data)
kw<-kruskal.test(HER2_unique_motif1 ~ celltype_RNA, data = atac_subset_qc@meta.data)

df_main <- df %>% filter(celltype_RNA %in% c("LHS", "LASP-1", "LASP-2", "B-Myo"))
# Step 2: Subset for other subtypes (Fb, Endo, T, B, Myl, Pc)
df_other <- df %>% filter(celltype_RNA %in% c("Fb", "Endo", "T", "B", "Myl", "Pc"))

library(dunn.test)
dunn_test_main <- dunn.test(
  x = df_main$HER2_unique_motif1,
  g = df_main$celltype_RNA,
  method = "bh" )
dunn_df_main <- data.frame(
  Comparison = dunn_test_main$comparisons,
  Z = dunn_test_main$Z,
  P_unadj = dunn_test_main$P,
  P_adj = dunn_test_main$P.adjusted)
my_comparisons1 <- list(
  c("LASP-1", "LASP-2"), c("LASP-1", "LHS"), c("LASP-1", "B-Myo"),
  c("LASP-2", "LHS"), c("LASP-2", "B-Myo"), c("LHS", "B-Myo")
)
##  3. Normalise Dunn comparison strings (both directions)
dunn_df_main <- dunn_df_main %>%
  tidyr::separate(Comparison, into = c("g1","g2"), sep = " - ", remove = FALSE) %>%
  mutate(
    key_fwd = paste(g1, g2, sep = "-"),
    key_rev = paste(g2, g1, sep = "-"))

##  4. Map to the exact order of my_comparisons
target_keys <- sapply(my_comparisons1, paste, collapse = "-")

dunn_df_main <- dunn_df_main %>%
  mutate(
    match_key = case_when(
      key_fwd %in% target_keys ~ key_fwd,
      key_rev %in% target_keys ~ key_rev,
      TRUE                     ~ NA_character_
    )
  ) %>%
  filter(!is.na(match_key)) %>%                     # keep only the 15 requested pairs
  arrange(match(match_key, target_keys))            # order exactly as target_keys
# Add significance stars
dunn_df_main <- dunn_df_main %>%
  mutate(signif = case_when(
    is.na(P_adj)                  ~ "NA",
    P_adj == 0                    ~ "***",
    P_adj <= 0.001                ~ "***",
    P_adj <= 0.01                 ~ "**",
    P_adj <= 0.05                 ~ "*",
    TRUE                          ~ "ns"
  ))

# Plot
p1 <- ggboxplot(df_main, x = "celltype_RNA", y = "HER2_unique_motif1",
                color = "celltype_RNA", palette = "jco", legend = "none") +
  stat_compare_means(method = "kruskal.test", label.y = 4.5, size = 3) +
  geom_signif(comparisons = my_comparisons1,
              annotations = dunn_df_main$signif,
              y_position = seq(2.4, by = 0.22, length.out = 6),  # adjust spacing
              textsize = 3,
              vjust = 0.5) +
  labs(x = "Cell type", y = "AddModuleScore by Her2 unique TF motifs", title = "") +
  theme(axis.title.y = element_text(size = 8))
jpeg('results/tfgwas/Boxplot_HER2_uniqueMotifs_Addmodulescore_part1.jpg',width=3.5, height = 5, units = 'in', res = 150)
print(p1)
dev.off()

write.csv(dunn_df_main,'results/tfgwas/dunn_df_main_HER2.csv',row.names = FALSE)

# other cell type
dunn_test_other <- dunn.test(
  x = df_other$HER2_unique_motif1,
  g = df_other$celltype_RNA,
  method = "bh" )
dunn_df_other <- data.frame(
  Comparison = dunn_test_other$comparisons,
  Z = dunn_test_other$Z,
  P_unadj = dunn_test_other$P,
  P_adj = dunn_test_other$P.adjusted)
##  2. Desired comparison order (my_comparisons2)
my_comparisons2 <- list(
  c("Fb","Endo"), c("Fb","T"),   c("Fb","B"),   c("Fb","Myl"), c("Fb","Pc"),
  c("Endo","T"),  c("Endo","B"), c("Endo","Myl"),c("Endo","Pc"),
  c("T","B"),     c("T","Myl"),  c("T","Pc"),
  c("B","Myl"),   c("B","Pc"),   c("Myl","Pc")
)
##  3. Normalise Dunn comparison strings (both directions)
dunn_df_other <- dunn_df_other %>%
  tidyr::separate(Comparison, into = c("g1","g2"), sep = " - ", remove = FALSE) %>%
  mutate(
    key_fwd = paste(g1, g2, sep = "-"),
    key_rev = paste(g2, g1, sep = "-"))

##  4. Map to the exact order of my_comparisons2
target_keys <- sapply(my_comparisons2, paste, collapse = "-")

dunn_df_other <- dunn_df_other %>%
  mutate(
    match_key = case_when(
      key_fwd %in% target_keys ~ key_fwd,
      key_rev %in% target_keys ~ key_rev,
      TRUE                     ~ NA_character_
    )
  ) %>%
  filter(!is.na(match_key)) %>%                     # keep only the 15 requested pairs
  arrange(match(match_key, target_keys))            # order exactly as target_keys

##  5. (Optional) Re-create signif column – safe for NA
dunn_df_other <- dunn_df_other %>%
  mutate(
    signif = case_when(
      is.na(P_adj)               ~ "NA",
      P_adj == 0                 ~ "***",
      P_adj <= 0.001             ~ "***",
      P_adj <= 0.01              ~ "**",
      P_adj <= 0.05              ~ "*",
      TRUE                       ~ "ns") )
##  6. Plot
jpeg('results/tfgwas/Boxplot_HER2_uniqueMotifs_Addmodulescore_part2.jpg',
     width = 4.5, height = 6, units = 'in', res = 150)
p2 <- ggboxplot(df_other,
                x = "celltype_RNA",
                y = "HER2_unique_motif1",
                color = "celltype_RNA",
                palette = "YlGnBu",
                legend = "none") +

  ## overall Kruskal-Wallis (small font)
  stat_compare_means(method = "kruskal.test",
                     label.y = 8,
                     size = 3) +
  ## 15 pairwise Dunn stars (BH-adjusted)
  geom_signif(comparisons = my_comparisons2,
              annotations = dunn_df_other$signif,
              y_position = seq(2.8, by = 0.32, length.out = length(my_comparisons2)),
              textsize = 3,
              vjust = 0.5) +
  labs(x = "Cell type",
       y = "AddModuleScore by Her2 unique TF motifs",
       title = "") +
  theme(axis.title.y = element_text(size = 8))
print(p2)
dev.off()
write.csv(dunn_df_other,'results/tfgwas/dunn_df_other_HER2.csv',row.names = FALSE)

###################################################################################
#                          Overlap                                                #
###################################################################################
LumA<-fread("results/tfgwas/tfgwas_results_LumA_ALL.csv", header = TRUE)
LumB<-fread("results/tfgwas/tfgwas_results_LumB_ALL.csv", header = TRUE)
LumB_HER2Neg<-fread("results/tfgwas/tfgwas_results_LumB_HER2Neg_ALL.csv", header = TRUE)

# Filter significant motifs based on TF_p_value
significant_LumA <- LumA[LumA$TF_p_value < 0.05/727 ,.(motif)]
significant_LumB <- LumB[TF_p_value < 0.05/727, .(motif)]
significant_LumB_HER2Neg <- LumB_HER2Neg[TF_p_value < 0.05/727,.(motif)]

# Find overlapping motifs across all three datasets
common_motifs <- Reduce(intersect, list(
  significant_LumA$motif,
  significant_LumB$motif,
  significant_LumB_HER2Neg$motif
))
common_motifs_ds <- data.frame(motif=common_motifs)
dim(common_motifs_ds)

fwrite(common_motifs_ds,"results/tfgwas/tfgwas_significant_LumOverlaps_unique.csv",row.names = FALSE)
enriched.motifs1$motif.name<-NULL
enriched_motifs_list_common_unique <-left_join(common_motifs_ds, enriched.motifs1, by='motif') %>% unique()

enriched_motifs_list_common_unique$cluster <- unique(combined_peaks_gene_all$cluster)[1]

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
  enriched_motifs_filtered <- left_join(common_motifs_ds, enriched_motifs, by='motif') %>% unique()
  enriched_motifs_filtered$cluster <- cluster

  # Store the filtered results in the list with the cluster name as the key
  enriched_motifs_list_common_unique <- rbind(enriched_motifs_list_common_unique,enriched_motifs_filtered)

  # Print a message to indicate progress
  cat("Processed cluster:", cluster, "\n")
}

enriched_motifs_list_common_unique$pvalue <- ifelse(enriched_motifs_list_common_unique$pvalue==0,1e-300,
                                                          ifelse(enriched_motifs_list_common_unique$pvalue==1,0.99999999999999999999999,enriched_motifs_list_common_unique$pvalue))

# create a matrix
# Get all unique motif names across all clusters
motif_id_list_final<-fread('results/singlecell/jaspar2022_motifs_727.csv', header = TRUE)
enriched_motifs_list_common_unique <- left_join(enriched_motifs_list_common_unique, select(motif_id_list_final, motif, motif.name), by='motif')
# Save
write.csv(enriched_motifs_list_common_unique, file = "results/singlecell/enriched_motifs_list_uniqDAP_LuminalcommonUniqSignificant.csv", row.names = FALSE)

enriched_motifs_list_common_unique<-fread("results/singlecell/enriched_motifs_list_uniqDAP_LuminalcommonUniqSignificant.csv",header = TRUE)

all_motifsID <- unique(enriched_motifs_list_common_unique$motif)
all_motifs <- motif_id_list_final[motif %in% all_motifsID]
all_motifs$motif.name <- make.unique(all_motifs$motif.name)
all_motifs <- all_motifs$motif.name

# Get all cluster names
all_clusters <- unique(enriched_motifs_list_common_unique$cluster)

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
  enriched_data <- enriched_motifs_list_common_unique[enriched_motifs_list_common_unique$cluster==cluster1,]

  # Add p-values to the matrix
  for (i in 1:nrow(enriched_data)) {
    motif <- enriched_data$motif.name[i]
    neglog10pvalue <- -log10(enriched_data$pvalue[i])  # Replace "pvalue" with the actual column name for p-values
    motif_matrix[motif,cluster1] <-neglog10pvalue
  }
}

# ---------heatmap with specific motif sequence----------------
motif_df<- as.data.frame(motif_matrix)

colSums(motif_matrix >= -log10(0.05 / 727))
# Hr LS-1 LASP-2   Bs    T    B  Myl Endo   Pc   Fb
colSums(motif_matrix >= 100)
# LHr LS-1 LASP-2   Bs    T    B  Myl Endo   Pc   Fb

motif_df <- motif_df[order(-motif_df$`LS-1`),]
rownames(motif_df)
motif_list <- c("IKZF1", "TFAP2A","ZEB1","TCF3","EHF","ELF3" ,"TFAP2B",
                "TFAP2C","TCF4","TFAP2E","TCF12","ZNF680","EWSR1-FLI1","ZNF454",
                "FOXO1::ELK1","ELF1","SPI1","FOXJ2::ELF1", "TBX2","ETV5::FOXI1","ETV5::FOXO1",
                "TBX21","ZNF530","ZNF257","ZNF682","SPIB","GABPA","FOXO1::ELF1",
                "ZNF418","ERF::FOXI1","FLI1::FOXI1","ETV1","PLAGL2","RBPJ","ZNF354C",
                "ERF::FOXO1","ZNF354A","PRDM9","ZNF707","ETV2::FOXI1","TBX18","ZNF341",
                "ZIM3","FOXN3","FOXA1","FOXC2","FOXD2","FOXC1","FOXB1","FOXD3",
                "TEAD3" ,"TEAD4" ,"TEAD1","TEAD2","POU1F1","ETV2::HOXB13","POU3F3","POU3F1","BCL6","STAT3",
                "ZKSCAN5","ZNF136","ZBTB26","HES5","HEY1","POU3F4","POU2F1" ,
                "POU2F2","POU2F3","NR5A1","POU4F2","ZSCAN29","ZNF274",
                "POU3F2","POU5F1B","ZNF460",
                "ATOH1","TWIST1","SP5","KLF4","SP4","ZNF148","MAZ",
                "KLF5" ,"PATZ1","KLF1","KLF14","RORB","KLF7","KLF15","SP8",
                "FOS::JUN","CREB1","ZNF263","SP2","KLF2","SATB1","SP3","NR2C2", "JUN",
                "RUNX1" ,
                "LBX1","NKX6-1","KLF6","HOXB3" ,"HNF1B","REL","RELA","HESX1",
                "HNF1A","DLX6","STAT1","STAT1::STAT2","IRF1","POU5F1","HES1","ZBTB14","HES2","HEY2",
                "ESRRB","ESRRA" ,"PPARG::RXRA","NR1H2::RXRA",
                "ZNF675","THRB","THRA","NKX2-2","NR2F1","HNF4A","RFX2","TFDP1","NR1H4","HINFP",
                "TFAP2A.1","TFAP2B.1","NR2C2.1","THRB.1" )%>% unique()

# create the new data frame
motif_df<- as.data.frame(motif_matrix)
head(motif_df)
# Reorder motif_df rows based on motif_list
motif_df <- motif_df[match(motif_list, rownames(motif_df)), ]
motif_df <- motif_df %>% select('LS-1','LASP-2','LHr','Bs','Fb','Endo','T','B','Myl','Pc')
motif_matrix <- as.matrix(motif_df)

library(ComplexHeatmap)
library(circlize)

# Create the heatmap
pdf("results/tfgwas/heatmap_enrichedmotif_uniqMotif_LumOverlap_137motifs.pdf", width = 5, height =13)
Heatmap(motif_matrix,
        col =  circlize::colorRamp2(c(0, 100), c("#483d8b", "red1")),
        row_names_gp = gpar(fontsize = 7),
        cluster_rows = FALSE,  # Preserve the order from motif_list
        cluster_columns = FALSE,  # Cluster columns (or set to FALSE if you want to preserve column order)
        show_row_names = TRUE,  # Show motif names
        show_column_names = TRUE,  # Show column names
        row_names_side = "left",  # Place row names on the left
        heatmap_legend_param = list(title = "-log10 pval"),
        rect_gp = gpar(col = "grey99"),)  # Customize legend title

dev.off()

############################ GSEA bubble plots ###########################################
# Load libraries
library(tidyverse)
library(ggplot2)
library(paletteer)
library(tidyverse)
library(ggplot2)
library(paletteer)
dat_LumOverlap<-fread("results/singlecell/enriched_motifs_list_uniqDAP_LuminalcommonUniqSignificant.csv",header = TRUE)
dat_LumOverlap$neglog10pval<--log10(dat_LumOverlap$pvalue)

motif_data<-fread('results/singlecell/jaspar2022_motif_list.csv', header=TRUE)
dat_LumOverlap <- left_join(dat_LumOverlap, motif_data, by=c('motif.name','motif'))
#
dat_LumOverlap$cluster <- gsub('LS-1','LASP-1',dat_LumOverlap$cluster )
dat_LumOverlap$cluster <- gsub('LS-2','LASP-2',dat_LumOverlap$cluster )
dat_LumOverlap$cluster <- gsub('LHr','LHS',dat_LumOverlap$cluster )
dat_LumOverlap$cluster <- gsub('Bs','B-Myo',dat_LumOverlap$cluster )

dat_LumOverlap$cluster <- factor(dat_LumOverlap$cluster, levels = c('LASP-1','LASP-2','LHS','B-Myo','Fb','Endo','T','B','Myl','Pc'))
dat_LumOverlap$avg.exp.scaled <- ifelse(is.na(dat_LumOverlap$avg.exp.scaled), 0, dat_LumOverlap$avg.exp.scaled)

fwrite(dat_LumOverlap,"results/singlecell/enriched_motifs_list_uniqDAP_HER2UniqSignificant.csv",row.names=FALSE)

jpeg('results/tfgwas/GSEAplot_LuminalOverlaps_uniquemotifs.jpg',width = 18, height = 26, units = 'in', res = 500)
ggplot1 <- ggplot(dat_LumOverlap, aes(x = fold.enrichment, y = reorder(motif.name, fold.enrichment)) )+
  geom_point(aes(size = avg.exp.scaled, color = neglog10pval)) +
  coord_cartesian(clip = "off") +
  scale_color_gradient(low = "#F88586", high = "#FC0D0F", name = "Motif enrichment\n-log10(P value)",
                       limits = c(-log10(0.05), 310)) +
  theme_bw() +
  scale_x_continuous(breaks = 1) +
  scale_size_continuous(name = "Average gene expression\n(normalized)") +
  labs(x = "Fold enrichment\n(Over-representation, enriched by unique peaks in each cell type)", y = "TF motif Name", title = "") +
  facet_grid(family ~ cluster, scales = "free_y",space = "free_y",
             labeller = label_wrap_gen(width = 45)) + # Wrap family names
  theme(axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        axis.text.y = element_text(size = 10, color = 'black'),
        legend.position = "bottom",
        strip.text.y = element_text(size = 10, angle = 0), # Horizontal family names
        strip.text.x = element_text(size = 12) , # Style cluster labels (horizontal)
        panel.spacing.y = unit(0.2, "lines"))
ggplot1
dev.off()

dat_LumOverlap<-fread("results/singlecell/enriched_motifs_list_uniqDAP_LuminalcommonUniqSignificant.csv",header = TRUE)
dat_LumOverlap$neglog10pval<--log10(dat_LumOverlap$pvalue)
motif_data<-fread('results/singlecell/jaspar2022_motif_list.csv', header=TRUE)
dat_LumOverlap <- left_join(dat_LumOverlap, motif_data, by=c('motif.name','motif'))

dat_LumOverlap_motif <- unique(select(dat_LumOverlap, motif, motif.name, family))
dat_LumOverlap_motif_group <- dat_LumOverlap_motif %>%
  group_by(family) %>%
  mutate(count = n()) %>%
  select(family, count) %>%
  ungroup() %>%
  unique()

library(stringr)
jpeg('results/tfgwas/GWASMotif01.LuminalOverlap_motif_count.jpg',
     width=8, height = 10, units = 'in', res = 250)
ggplot(dat_LumOverlap_motif_group, aes(x = count, y = reorder(family, count), fill = count)) +
  geom_bar(stat = "identity", fill = "#FF4500") +
  labs(x = "No. of TF motifs", y = "TF motif family", title = "") +
  scale_y_discrete(labels = function(x) str_wrap(x, width = 40)) +
  theme_classic() +
  theme(legend.position = "right") +
  scale_x_continuous(breaks = seq(0, 22, 1), limits = c(0, 22))
dev.off()

####################### GSEA Bubble plot - motif score####################
dat_LumOverlap<-fread("results/singlecell/enriched_motifs_list_uniqDAP_LuminalcommonUniqSignificant.csv",header = TRUE)
#

dat_LumOverlap$cluster <- factor(dat_LumOverlap$cluster, levels = c('LS-1','LASP-2','LHr','Bs','Fb','Endo','T','B','Myl','Pc'))
dat_LumOverlap$avg.exp.scaled <- ifelse(is.na(dat_LumOverlap$avg.exp.scaled), 0, dat_LumOverlap$avg.exp.scaled)

motif_id_name<-fread('results/singlecell/jaspar2022_motifs_727.csv', header = TRUE)
motif_scores_ct <- subset(motif_scores_ct,rownames(motif_scores_ct) %in%  motif_id_name$motif)
motif_scores_ct <- as.data.frame(motif_scores_ct)
motif_scores_ct$motif <- row.names(motif_scores_ct)

motif_id_list <- left_join(select(dat_LumOverlap,motif,motif.name), motif_id_name, by=c('motif','motif.name')) %>% unique()
motif_score_selected <- left_join(motif_id_list, motif_scores_ct, by='motif') %>% as.data.frame()
# for the motif score heatmap
rownames(motif_score_selected)<-motif_score_selected$motif
motif_score_selected <- motif_score_selected[,c(5:14)]

motif_matrix <- as.matrix(motif_score_selected)
scaled_mat <- t(scale(t(motif_matrix)))
scaled_mat <- as.data.frame(scaled_mat)
scaled_mat$motif <- row.names(scaled_mat)

motif_scores_long <- scaled_mat %>%
  pivot_longer(cols = -motif, names_to = "cluster", values_to = "motif_scores")
dat_LumOverlap <- left_join(dat_LumOverlap, motif_scores_long, by=c('cluster','motif'))

# Convert the cluster column to a factor with the specified order
cluster_order <- c("LS-1", "LASP-2", "LHr", "Bs", "Fb", "Endo", "T", "B", "Myl", "Pc")
dat_LumOverlap$cluster <- factor(dat_LumOverlap$cluster, levels = cluster_order)

jpeg('results/tfgwas/GSEAplot_LumOverlap_motif_score.jpg',width=15, height = 20, units = 'in', res = 140)
ggplot1<-ggplot(dat_LumOverlap, aes(x = fold.enrichment, y = reorder(motif.name, avg.exp.scaled))) +
  geom_point(aes(size = avg.exp.scaled, color = motif_scores)) +
  coord_cartesian(clip = "off") +
  scale_color_gradient2(low="darkblue", high="red",mid = "grey90", name="Motif score",midpoint = 0) +
  theme_bw() +
  scale_size_continuous(name = "Average RNA expression(scaled)") +
  labs(x = "Fold.Enrichment", y = "Motif Name", title = "LumOverlap unique motifs") +
  facet_grid(. ~ cluster, scales = "free_y")+  # Facet by cluster, free y-scales for independent motif lists
  theme(axis.title.x = element_text(size = 8),  # Smaller font size for x-label
        axis.title.y = element_text(size = 8))
ggplot1
dev.off()
##########################################################################################
#                                 feature plot                                           #
##########################################################################################
LumOverlap_unique<-fread("results/tfgwas/tfgwas_significant_LumOverlaps_unique.csv",header = TRUE)
LumOverlap_unique<-list(unique(LumOverlap_unique$motif))
DefaultAssay(atac_subset_qc) <- "chromvar"
set.seed(4242025)
atac_subset_qc <- AddModuleScore(
  object = atac_subset_qc,
  features = LumOverlap_unique,
  ctrl = 15,
  name = 'LumOverlap_unique_motif'
)

library(INLA)
jpeg('results/tfgwas/GWASMotif01.all_LumOverlap_uniqueMotifs_Addmodulescore.jpg',width=6, height = 5, units = 'in', res = 150)
p<-FeaturePlot(
  object = atac_subset_qc,
  features = "LumOverlap_unique_motif1",
  min.cutoff = 'q5',
  max.cutoff = 'q95',
  pt.size = 0.1,
  label = TRUE,
  cols = c( "grey90","orange")
)+
  xlab("UMAP 1") +
  ylab("UMAP 2")+
  theme(legend.position = "right")

# Extract UMAP coordinates and cluster assignments
umap_data <- as.data.frame(atac_subset_qc@reductions$umap.rna@cell.embeddings)
umap_data$cluster <- atac_subset_qc@meta.data$celltype_RNA

# Filter for the clusters of interest
clusters_to_encircle <- c("LHS", "LASP-1", "LASP-2", "B-Myo")
umap_subset <- umap_data[umap_data$cluster %in% clusters_to_encircle, ]

points_matrix <- as.matrix(umap_subset[, c("rnaUMAP_1", "rnaUMAP_2")])
hull <- inla.nonconvex.hull(points_matrix, convex = -0.0021)
# Adjust convex for tightness (more negative for less smoothing)
segments <- data.frame(
  x = hull$loc[hull$idx[, 1], 1],
  y = hull$loc[hull$idx[, 1], 2],
  xend = hull$loc[hull$idx[, 2], 1],
  yend = hull$loc[hull$idx[, 2], 2]
)

p <- p + geom_segment(
  data = segments,
  aes(x = x, y = y, xend = xend, yend = yend),
  color = "black",  # Outline color
  linetype = "dotted",  # Dotted line
  size = 0.5  # Line thickness
)
print(p)
dev.off()

df <- atac_subset_qc@meta.data
celltype_order <- c("LASP-1", "LASP-2", "LHS", "B-Myo", "Fb", "Endo", "T", "B", "Myl", "Pc")

# Convert celltype_RNA to a factor with specified order
levels(df$celltype_RNA)
df$celltype_RNA <- factor(df$celltype_RNA, levels = celltype_order)

# Perform statistical comparison (e.g., Kruskal-Wallis test for non-parametric data)
kw<-kruskal.test(LumOverlap_unique_motif1 ~ celltype_RNA, data = atac_subset_qc@meta.data)

df_main <- df %>% filter(celltype_RNA %in% c("LHS", "LASP-1", "LASP-2", "B-Myo"))
# Step 2: Subset for other subtypes (Fb, Endo, T, B, Myl, Pc)
df_other <- df %>% filter(celltype_RNA %in% c("Fb", "Endo", "T", "B", "Myl", "Pc"))

jpeg('results/tfgwas/Boxplot_LumOverlap_uniqueMotifs_Addmodulescore_part1.jpg',width=6, height = 6, units = 'in', res = 150)
my_comparisons1 <- list( c("LASP-1", "LASP-2"), c("LASP-1", "LHS"), c("LASP-1", "B-Myo"),
                         c("LASP-2","LHS"),c("LASP-2","B-Myo"),c("LHS","B-Myo"))
p1<-ggboxplot(df_main, x = "celltype_RNA", y = "LumOverlap_unique_motif1",
              color = "celltype_RNA", palette = "jco", legend = "none")+
  stat_compare_means(comparisons = my_comparisons1 ,label = "p.signif")+
  stat_compare_means(label.y = 3)  +
  labs(x = "Cell type", y = "AddModuleScore by overlapped TF motifs in luminal subtypes", title = "")
p1
dev.off()

jpeg('results/tfgwas/Boxplot_LumOverlap_uniqueMotifs_Addmodulescore_part2.jpg',width=6, height = 7, units = 'in', res = 150)
my_comparisons2 <- list(c("Fb","Endo"),c("Fb","T"),c("Fb","B"),c("Fb","Myl"),c("Fb","Pc"),
                        c("Endo","T"),c("Endo","B"),c("Endo","Myl"),c("Endo","Pc"),
                        c("T","B"),c("T","Myl"),c("T","Pc"),
                        c("B","Myl"),c("B","Pc"),
                        c("Myl","Pc"))
p2<-ggboxplot(df_other, x = "celltype_RNA", y = "LumOverlap_unique_motif1",
              color = "celltype_RNA", palette = "YlGnBu", legend = "none")+
  stat_compare_means(comparisons = my_comparisons2 ,label = "p.signif")+
  stat_compare_means(label.y = 5)+
labs(x = "Cell type", y = "AddModuleScore by overlapped TF motifs in luminal subtypes", title = "")
# p1+p2
p2
dev.off()
###################################################################################
#                          Overlap- TNBC, Lum                                           #
###################################################################################
TNBC<-fread("results/tfgwas/tfgwas_results_TNBC_ALL.csv", header = TRUE)
LumA<-fread("results/tfgwas/tfgwas_results_LumA_ALL.csv", header = TRUE)
LumB<-fread("results/tfgwas/tfgwas_results_LumB_ALL.csv", header = TRUE)
LumB_HER2Neg<-fread("results/tfgwas/tfgwas_results_LumB_HER2Neg_ALL.csv", header = TRUE)

# Filter significant motifs based on TF_p_value
significant_TNBC <- TNBC[TNBC$TF_p_value < 0.05/727,.(motif)]
significant_LumA <- LumA[LumA$TF_p_value < 0.05/727,.(motif)]
significant_LumB <- LumB[TF_p_value < 0.05/727, .(motif)]
significant_LumB_HER2Neg <- LumB_HER2Neg[TF_p_value < 0.05/727,.(motif)]

# Find overlapping motifs across all three datasets
common_motifs <- Reduce(intersect, list(
  significant_TNBC$motif,
  significant_LumA$motif,
  significant_LumB$motif,
  significant_LumB_HER2Neg$motif
))
common_motifs_ds <- data.frame(motif=common_motifs)
dim(common_motifs_ds)

enriched.motifs1$motif.name<-NULL
enriched_motifs_list_common_unique <-left_join(common_motifs_ds, enriched.motifs1, by='motif') %>% unique()

enriched_motifs_list_common_unique$cluster <- unique(combined_peaks_gene_all$cluster)[1]

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
  enriched_motifs_filtered <- left_join(common_motifs_ds, enriched_motifs, by='motif') %>% unique()
  enriched_motifs_filtered$cluster <- cluster

  # Store the filtered results in the list with the cluster name as the key
  enriched_motifs_list_common_unique <- rbind(enriched_motifs_list_common_unique,enriched_motifs_filtered)

  # Print a message to indicate progress
  cat("Processed cluster:", cluster, "\n")
}

enriched_motifs_list_common_unique$pvalue <- ifelse(enriched_motifs_list_common_unique$pvalue==0,1e-300,
                                                    ifelse(enriched_motifs_list_common_unique$pvalue==1,0.99999999999999999999999,enriched_motifs_list_common_unique$pvalue))

# Save
write.csv(enriched_motifs_list_common_unique, file = "results/singlecell/enriched_motifs_list_uniqDAP_TNBCLuminalcommonUniqSig.csv", row.names = FALSE)

enriched_motifs_list_common_unique<-fread("results/singlecell/enriched_motifs_list_uniqDAP_TNBCLuminalcommonUniqSig.csv",header = TRUE)

# create a matrix
# Get all unique motif names across all clusters
motif_id_list_final<-fread('results/singlecell/jaspar2022_motifs_727.csv', header = TRUE)
enriched_motifs_list_common_unique <- left_join(enriched_motifs_list_common_unique, select(motif_id_list_final, motif, motif.name), by='motif')
all_motifsID <- unique(enriched_motifs_list_common_unique$motif)
all_motifs <- motif_id_list_final[motif %in% all_motifsID]
all_motifs$motif.name <- make.unique(all_motifs$motif.name)
all_motifs <- all_motifs$motif.name

# Get all cluster names
all_clusters <- unique(enriched_motifs_list_common_unique$cluster)

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
  enriched_data <- enriched_motifs_list_common_unique[enriched_motifs_list_common_unique$cluster==cluster1,]

  # Add p-values to the matrix
  for (i in 1:nrow(enriched_data)) {
    motif <- enriched_data$motif.name[i]
    neglog10pvalue <- -log10(enriched_data$pvalue[i])  # Replace "pvalue" with the actual column name for p-values
    motif_matrix[motif,cluster1] <-neglog10pvalue
  }
}

# ---------heatmap with specific motif sequence----------------
motif_df<- as.data.frame(motif_matrix)

colSums(motif_matrix >= -log10(0.05 / 727))
# LHr LS-1 LASP-2   Bs    T    B  Myl Endo   Pc   Fb
colSums(motif_matrix >= 100)
# LHr LS-1 LASP-2   Bs    T    B  Myl Endo   Pc   Fb

motif_df <- motif_df[order(-motif_df$LHr),]
rownames(motif_df)
motif_list <- c("TFAP2A","TFAP2B","EHF","ELF3","IKZF1","ZNF454","ZNF680","ELF1","EWSR1-FLI1",
                "FOXJ2::ELF1","FOXO1::ELK1","SPI1","GABPA","SPIB","ETV1","ZNF257","ZNF418",
                "FLI1::FOXI1","ERF::FOXI1","PLAGL2","ERF::FOXO1","TBX21","TBX2","ZNF354A","PRDM9",
                "RBPJ","ETV2::FOXI1","TBX18","POU3F1","ZNF354C","ZIM3",
                "FOXC1","FOXD3","FOXA1" ,"FOXC2","FOXN3" ,
                "POU1F1","HES5",
                "POU3F3","NR5A1","STAT3","POU2F1","ZNF460","POU3F2",
                "KLF4","RUNX1","KLF5" ,"KLF1","KLF2","TEAD3","TEAD1","KLF6",
                "NKX6-1","SP3","POU2F2","ZKSCAN5","HNF1A","SP8",
                "PATZ1","NR2C2","ZNF263","KLF14","SP5","ESRRB","KLF7","SP2" , "KLF15",
                "STAT1::STAT2","IRF1","POU2F3","POU5F1","HES1","ZBTB14","ATOH1","TFDP1",
                "THRA","NR2F1","NR1H4" ,"RFX2","NKX2-2","TFAP2B.1" )%>% unique()

motif_df$motif <- rownames(motif_df)
motif_df <- motif_df[!motif_df$motif %in% motif_list,]
motif_df <- motif_df[order(-motif_df$`LS-1`),]
motif_df <- motif_df[order(-motif_df$`LASP-2`),]
motif_df <- motif_df[order(-motif_df$`Endo`),]

motif_df <- motif_df[order(-motif_df$B),]

# create the new data frame
motif_df<- as.data.frame(motif_matrix)
head(motif_df)
# Reorder motif_df rows based on motif_list
motif_df <- motif_df[match(motif_list, rownames(motif_df)), ]
motif_df <- motif_df %>% select('LS-1','LASP-2','LHr','Bs','Fb','Endo','T','B','Myl','Pc')
motif_matrix <- as.matrix(motif_df)

library(ComplexHeatmap)
library(circlize)

# Create the heatmap
pdf("results/tfgwas/heatmap_enrichedmotif_uniqMotif_TNBCLumOverlap_137motifs.pdf", width = 5, height =9)
Heatmap(motif_matrix,
        col =  circlize::colorRamp2(c(0, 100), c("#483d8b", "red1")),
        row_names_gp = gpar(fontsize = 7),
        cluster_rows = FALSE,  # Preserve the order from motif_list
        cluster_columns = FALSE,  # Cluster columns (or set to FALSE if you want to preserve column order)
        show_row_names = TRUE,  # Show motif names
        show_column_names = TRUE,  # Show column names
        row_names_side = "left",  # Place row names on the left
        heatmap_legend_param = list(title = "-log10 pval"),
        rect_gp = gpar(col = "grey99"),)  # Customize legend title

dev.off()

##########################################################################################
#                                 feature plot                                           #
##########################################################################################
DefaultAssay(atac_subset_qc) <- "chromvar"

common_motifs_ds<-list(unique(common_motifs_ds$motif))
set.seed(4242025)
atac_subset_qc <- AddModuleScore(
  object = atac_subset_qc,
  features = common_motifs_ds,
  ctrl = 5,
  name = 'common_motifs_motif'
)

jpeg('results/singlecell/GWASMotif01.all_TNBCcommon_Motifs_Addmodulescore.jpg',width=6, height = 5, units = 'in', res = 150)
FeaturePlot(
  object = atac_subset_qc,
  features = "common_motifs_motif1",
  min.cutoff = 'q5',
  max.cutoff = 'q95',
  pt.size = 0.1,
  label = TRUE,
  cols = c("orange", "darkblue")
)
dev.off()

############################################################################################
#                                 Manhattan plot                                          #
############################################################################################
BCAC_all<-fread("results/tfgwas/variant_peak_motif_subtypes_ALL.csv", header = TRUE)

# TNBC
TNBC_unique<-fread("results/tfgwas/tfgwas_significant_TNBC_unique.csv",header = TRUE)

BCAC_TNBC_unique <- BCAC_all%>% select(Chr, Pos, chrpos37,TNBC_pval,TNBC_unique$motif)%>%unique()
BCAC_TNBC_unique$TNBC_pval<-ifelse(BCAC_TNBC_unique$TNBC_pval==0,1e-300,
                                 ifelse(BCAC_TNBC_unique$TNBC_pval==1,0.99999999999999999999999,
                                        BCAC_TNBC_unique$TNBC_pval))
head(BCAC_TNBC_unique)
ma_cols <- grep("^MA", names(BCAC_TNBC_unique), value = TRUE)
all_zero_MA <- rowSums(BCAC_TNBC_unique[, ma_cols] != 0) == 0

# rsid
snp_data<-fread("snp_index.txt", header = TRUE)
snp_data$chrpos37 <- paste0('chr',snp_data$CHR,':',snp_data$POS)
names(snp_data)[1] <- 'rsid'
snp_data  <- snp_data  %>% dplyr::select(chrpos37, rsid, REF, ALT)

# join rsid
BCAC_TNBC_unique <- left_join(BCAC_TNBC_unique, snp_data, by='chrpos37')
table(is.na(BCAC_TNBC_unique$rsid))
BCAC_TNBC_unique$Chr <- as.numeric(gsub('chr','',BCAC_TNBC_unique$Chr))
BCAC_TNBC_unique<-BCAC_TNBC_unique[!is.na(rsid)&!is.na(TNBC_pval),]

# filter the TNBC sig snps
ma_cols <- grep("^MA[0-9]+\\.[0-9]+", names(BCAC_TNBC_unique), value = TRUE)
# Create a logical vector: TRUE if any MA column == 1 for that row
has_ma_1 <- BCAC_TNBC_unique[, rowSums(.SD == 1) > 0, .SDcols = ma_cols]
# Subset the original data to keep only rows where at least one MA column is 1
TNBC_sig_snp <- BCAC_TNBC_unique[has_ma_1]

library(qqman)
png("results/tfgwas/manhattan_TNBC.png", width = 14, height = 8, units = "in", res = 150)
par(mar = c(5, 5, 1, 1))
manhattan_local(BCAC_TNBC_unique,
          chr = "Chr",
          bp = "Pos",
          p = "TNBC_pval",
          snp = "rsid",
          col = c("black", "grey"),
          chrlabs = NULL,
          suggestiveline =FALSE,
          genomewideline = -log10(5e-08),
          highlight = TNBC_sig_snp$rsid,
          logp = TRUE,
          annotateTop = FALSE)
dev.off()

#############################################
# Luminal A
LumA_unique<-fread("results/tfgwas/tfgwas_significant_LumA_unique.csv",header = TRUE)
BCAC_LumA_unique <- BCAC_all%>% select(Chr, Pos, chrpos37,LuminalA_pval,LumA_unique$motif) %>%unique()
BCAC_LumA_unique$LuminalA_pval<-ifelse(BCAC_LumA_unique$LuminalA_pval==0,1e-300,
                                   ifelse(BCAC_LumA_unique$LuminalA_pval==1,0.99999999999999999999999,
                                          BCAC_LumA_unique$LuminalA_pval))
# join rsid
BCAC_LumA_unique <- left_join(BCAC_LumA_unique, snp_data, by='chrpos37')
table(is.na(BCAC_LumA_unique$rsid))
BCAC_LumA_unique$Chr <- as.numeric(gsub('chr','',BCAC_LumA_unique$Chr))
BCAC_LumA_unique<-BCAC_LumA_unique[!is.na(rsid)&!is.na(LuminalA_pval),]

# filter the LumA sig snps
ma_cols <- grep("^MA[0-9]+\\.[0-9]+", names(BCAC_LumA_unique), value = TRUE)
# Create a logical vector: TRUE if any MA column == 1 for that row
has_ma_1 <- BCAC_LumA_unique[, rowSums(.SD == 1) > 0, .SDcols = ma_cols]
# Subset the original data to keep only rows where at least one MA column is 1
LumA_sig_snp <- BCAC_LumA_unique[has_ma_1]

library(qqman)
png("results/tfgwas/manhattan_LumA.png", width = 14, height = 8, units = "in", res = 150)
par(mar = c(5, 5, 1, 1))
manhattan_local(BCAC_LumA_unique,
          chr = "Chr",
          bp = "Pos",
          p = "LuminalA_pval",
          snp = "rsid",
          col = c("black", "grey"),
          chrlabs = NULL,
          suggestiveline =FALSE,
          genomewideline = -log10(5e-08),
          highlight = LumA_sig_snp$rsid,
          logp = TRUE,
          annotateTop = FALSE)
dev.off()

#############################################
# Luminal B
LumB_unique<-fread("results/tfgwas/tfgwas_significant_LumB_unique.csv",header = TRUE)
BCAC_LumB_unique <- BCAC_all%>% select(Chr, Pos, chrpos37,LuminalB_P_meta,LumB_unique$motif) %>%unique()
BCAC_LumB_unique$LuminalB_P_meta<-ifelse(BCAC_LumB_unique$LuminalB_P_meta==0,1e-300,
                                       ifelse(BCAC_LumB_unique$LuminalB_P_meta==1,0.99999999999999999999999,
                                              BCAC_LumB_unique$LuminalB_P_meta))
# join rsid
BCAC_LumB_unique <- left_join(BCAC_LumB_unique, snp_data, by='chrpos37')
table(is.na(BCAC_LumB_unique$rsid))
BCAC_LumB_unique$Chr <- as.numeric(gsub('chr','',BCAC_LumB_unique$Chr))
BCAC_LumB_unique<-BCAC_LumB_unique[!is.na(rsid)&!is.na(LuminalB_P_meta),]

# filter the LumB sig snps
ma_cols <- grep("^MA[0-9]+\\.[0-9]+", names(BCAC_LumB_unique), value = TRUE)
# Create a logical vector: TRUE if any MA column == 1 for that row
has_ma_1 <- BCAC_LumB_unique[, rowSums(.SD == 1) > 0, .SDcols = ma_cols]
# Subset the original data to keep only rows where at least one MA column is 1
LumB_sig_snp <- BCAC_LumB_unique[has_ma_1]

library(qqman)
png("results/tfgwas/manhattan_LumB.png", width = 14, height = 8, units = "in", res = 150)
par(mar = c(4, 4, 1, 1))
manhattan_local(BCAC_LumB_unique,
          chr = "Chr",
          bp = "Pos",
          p = "LuminalB_P_meta",
          snp = "rsid",
          col = c("black", "grey"),
          chrlabs = NULL,
          suggestiveline =FALSE,
          genomewideline = -log10(5e-08),
          highlight = LumB_sig_snp$rsid,
          logp = TRUE,
          annotateTop = FALSE)
dev.off()

#############################################
# Luminal B HER2neg
LumB_HER2Neg_unique<-fread("results/tfgwas/tfgwas_significant_LumB_HER2Neg_unique.csv",header = TRUE)
BCAC_LumB_HER2Neg_unique <- BCAC_all%>% select(Chr, Pos, chrpos37,LuminalB_HER2Neg_pval,LumB_HER2Neg_unique$motif) %>%unique()
BCAC_LumB_HER2Neg_unique$LuminalB_HER2Neg_pval<-ifelse(BCAC_LumB_HER2Neg_unique$LuminalB_HER2Neg_pval==0,1e-300,
                                         ifelse(BCAC_LumB_HER2Neg_unique$LuminalB_HER2Neg_pval==1,0.99999999999999999999999,
                                                BCAC_LumB_HER2Neg_unique$LuminalB_HER2Neg_pval))
# join rsid
BCAC_LumB_HER2Neg_unique <- left_join(BCAC_LumB_HER2Neg_unique, snp_data, by='chrpos37')
table(is.na(BCAC_LumB_HER2Neg_unique$rsid))
BCAC_LumB_HER2Neg_unique$Chr <- as.numeric(gsub('chr','',BCAC_LumB_HER2Neg_unique$Chr))
BCAC_LumB_HER2Neg_unique<-BCAC_LumB_HER2Neg_unique[!is.na(rsid)&!is.na(LuminalB_HER2Neg_pval),]

# filter the LumB sig snps
ma_cols <- grep("^MA[0-9]+\\.[0-9]+", names(BCAC_LumB_HER2Neg_unique), value = TRUE)
# Create a logical vector: TRUE if any MA column == 1 for that row
has_ma_1 <- BCAC_LumB_HER2Neg_unique[, rowSums(.SD == 1) > 0, .SDcols = ma_cols]
# Subset the original data to keep only rows where at least one MA column is 1
LumB_HER2Neg_sig_snp <- BCAC_LumB_HER2Neg_unique[has_ma_1]

library(qqman)
png("results/tfgwas/manhattan_LumB_HER2neg.png", width = 14, height = 8, units = "in", res = 150)
par(mar = c(5, 5, 1, 1))
manhattan_local(BCAC_LumB_HER2Neg_unique,
          chr = "Chr",
          bp = "Pos",
          p = "LuminalB_HER2Neg_pval",
          snp = "rsid",
          col = c("black", "grey"),
          chrlabs = NULL,
          suggestiveline =FALSE,
          genomewideline = -log10(5e-08),
          highlight = LumB_HER2Neg_sig_snp$rsid,
          logp = TRUE,
          annotateTop = FALSE)
dev.off()

#############################################
# HER2
HER2_unique<-fread("results/tfgwas/tfgwas_significant_HER2_unique.csv",header = TRUE)
BCAC_HER2_unique <- BCAC_all%>% select(Chr, Pos, chrpos37,HER2_pval,HER2_unique$motif) %>%unique()
BCAC_HER2_unique$HER2_pval<-ifelse(BCAC_HER2_unique$HER2_pval==0,1e-300,
                                                 ifelse(BCAC_HER2_unique$HER2_pval==1,0.99999999999999999999999,
                                                        BCAC_HER2_unique$HER2_pval))
# join rsid
BCAC_HER2_unique <- left_join(BCAC_HER2_unique, snp_data, by='chrpos37')
table(is.na(BCAC_HER2_unique$rsid))
BCAC_HER2_unique$Chr <- as.numeric(gsub('chr','',BCAC_HER2_unique$Chr))
BCAC_HER2_unique<-BCAC_HER2_unique[!is.na(rsid)&!is.na(HER2_pval),]

# filter the LumB sig snps
ma_cols <- grep("^MA[0-9]+\\.[0-9]+", names(BCAC_HER2_unique), value = TRUE)
# Create a logical vector: TRUE if any MA column == 1 for that row
has_ma_1 <- BCAC_HER2_unique[, rowSums(.SD == 1) > 0, .SDcols = ma_cols]
# Subset the original data to keep only rows where at least one MA column is 1
HER2_sig_snp <- BCAC_HER2_unique[has_ma_1]

library(qqman)
png("results/tfgwas/manhattan_Her2.png", width = 14, height = 8, units = "in", res = 150)
par(mar = c(5, 5, 1, 1))
manhattan_local(BCAC_HER2_unique,
          chr = "Chr",
          bp = "Pos",
          p = "HER2_pval",
          snp = "rsid",
          col = c("black", "grey"),
          chrlabs = NULL,
          suggestiveline =FALSE,
          genomewideline = -log10(5e-08),
          highlight = HER2_sig_snp$rsid,
          logp = TRUE,
          annotateTop = FALSE)
dev.off()

########################################################
# Suppl TABLE
TNBC_unique<-fread("results/tfgwas/tfgwas_significant_TNBC_unique.csv",header = TRUE)
LumA_unique<-fread("results/tfgwas/tfgwas_significant_LumA_unique.csv",header = TRUE)
LumB_unique<-fread("results/tfgwas/tfgwas_significant_LumB_unique.csv",header = TRUE)
LumB_HER2Neg_unique<-fread("results/tfgwas/tfgwas_significant_LumB_HER2Neg_unique.csv",header = TRUE)
HER2_unique<-fread("results/tfgwas/tfgwas_significant_HER2_unique.csv",header = TRUE)

TNBC<-fread("results/tfgwas/tfgwas_results_TNBC_ALL.csv", header = TRUE)
LumA<-fread("results/tfgwas/tfgwas_results_LumA_ALL.csv", header = TRUE)
LumB<-fread("results/tfgwas/tfgwas_results_LumB_ALL.csv", header = TRUE)
LumB_HER2Neg<-fread("results/tfgwas/tfgwas_results_LumB_HER2Neg_ALL.csv", header = TRUE)
HER2<-fread("results/tfgwas/tfgwas_results_HER2_ALL.csv",header = TRUE)

# TNBC
TNBC$Bonf_sig<-ifelse(TNBC$TF_p_value<0.05/727, 'Y','')
TNBC$unique<-ifelse(TNBC$TF %in% TNBC_unique$motif,'Y',"")
TNBC_select<-TNBC %>% select(motif, motif.name, TF_p_value, Bonf_sig, unique)
names(TNBC_select)[3:5]<-paste0("TNBC_",names(TNBC_select)[3:5])
# LumA
LumA$Bonf_sig<-ifelse(LumA$TF_p_value<0.05/727, 'Y','')
LumA$unique<-ifelse(LumA$TF %in% LumA_unique$motif,'Y',"")
LumA_select<-LumA %>% select(motif, motif.name, TF_p_value, Bonf_sig, unique)
names(LumA_select)[3:5]<-paste0("LumA_",names(LumA_select)[3:5])
# LumB
LumB$Bonf_sig<-ifelse(LumB$TF_p_value<0.05/727, 'Y','')
LumB$unique<-ifelse(LumB$TF %in% LumB_unique$motif,'Y',"")
LumB_select<-LumB %>% select(motif, motif.name, TF_p_value, Bonf_sig, unique)
names(LumB_select)[3:5]<-paste0("LumB_",names(LumB_select)[3:5])
# LumB_HER2Neg
LumB_HER2Neg$Bonf_sig<-ifelse(LumB_HER2Neg$TF_p_value<0.05/727, 'Y','')
LumB_HER2Neg$unique<-ifelse(LumB_HER2Neg$TF %in% LumB_HER2Neg_unique$motif,'Y',"")
LumB_HER2Neg_select<-LumB_HER2Neg %>% select(motif, motif.name, TF_p_value, Bonf_sig, unique)
names(LumB_HER2Neg_select)[3:5]<-paste0("LumB_HER2Neg_",names(LumB_HER2Neg_select)[3:5])
# HER2
HER2$Bonf_sig<-ifelse(HER2$TF_p_value<0.05/727, 'Y','')
HER2$unique<-ifelse(HER2$TF %in% HER2_unique$motif,'Y',"")
HER2_select<-HER2 %>% select(motif, motif.name, TF_p_value, Bonf_sig, unique)
names(HER2_select)[3:5]<-paste0("HER2_",names(HER2_select)[3:5])

# combine
All<-left_join(TNBC_select, LumA_select, by=c("motif","motif.name"))
All<-left_join(All, LumB_select, by=c("motif","motif.name"))
All<-left_join(All, LumB_HER2Neg_select, by=c("motif","motif.name"))
All<-left_join(All, HER2_select, by=c("motif","motif.name"))

common_motifs_ds<-fread("results/tfgwas/tfgwas_significant_LumOverlaps_unique.csv",header=TRUE)
common_motifs_ds$Overlaping<-'Y'
All<-left_join(All, common_motifs_ds, by='motif')

write.csv(All,"results/tfgwas/tfgwas_combined_all_subtypes.csv", row.names = FALSE)
