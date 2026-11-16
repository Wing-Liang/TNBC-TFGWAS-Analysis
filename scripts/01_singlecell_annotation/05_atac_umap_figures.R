# ==============================================================================
# 05_atac_umap_figures.R
# scATAC-seq dimensionality reduction, chromVAR motif activity and Fig 2 panels.
# Produces: Fig 2C-E
# Inputs  : results/singlecell/atac_RNA_annotated.rds
#           results/singlecell/final_cell_barcodes.csv
#           results/singlecell/atac_ATAC_RNAlabeled.rds
#           results/singlecell/ATAC_da_results_LS1_2.csv
#           results/singlecell/Markers_ATAC_Subclusters.csv
#           results/singlecell/atac_ATAC_subcluster_annoatated.rds
#           results/singlecell/ATAC_RNAlabelled_filterMT.rds
#           ...
# Outputs : results/singlecell/ATAC01_DepthCor_ATAC.jpg
#           results/singlecell/ATAC02_Dimplot_ATAC_bAnnot_10Clusters.jpg
#           results/singlecell/atac_ATAC_RNAlabeled.rds
#           results/singlecell/ATAC03.0_DimPlot_ATAC_coloredbyOWN.jpg
#           results/singlecell/ATAC03.1_DimPlot_ATAC_coloredbyATAC.jpg
#           results/singlecell/ATAC03.2_DimPlot_ATAC_coloredbyRNA.jpg
#           results/singlecell/ATAC04_DimPlot_ATAC.jpg
#           ...
# Run     : Rscript scripts/01_singlecell_annotation/05_atac_umap_figures.R   (from the repository root)
# ==============================================================================
#             '#edbb99','#b2d8d8',"#B09C85B2","#F39b7fb2","#8491B4B2")
colors <-c( "#E64B35CC","#4DBBD5CC", "#00A087CC", "#B09C85CC","#F39B7FCC", "#8491B4CC", "#91D1C2CC",'#f9e79f', "#7E6148CC",  "#3C5488CC")
####################################################
###                    scATAC-seq                 ###
####################################################
#------------------- Gene expression data processing------------------------------
options(future.globals.maxSize = 8000 * 1024^2) # avoid Error in getGlobalsAndPackages(expr, envir = envir, globals = globals)
atac_subset_qc.rna <- readRDS('results/singlecell/atac_RNA_annotated.rds')
DefaultAssay(atac_subset_qc.rna)<-'RNA'

atac_subset_qc.rna[['RNA']]<-as(atac_subset_qc.rna[['RNA']], Class="Assay5")
Final_cell_barcode<-fread('results/singlecell/final_cell_barcodes.csv',header=TRUE)
atac_subset_qc.rna <- subset(atac_subset_qc.rna , cells = Final_cell_barcode$x)
DefaultAssay(atac_subset_qc.rna)<-'RNA'

# Perform standard analysis of each modality independently RNA analysis
atac_subset_qc.rna <- NormalizeData(atac_subset_qc.rna)
atac_subset_qc.rna <- FindVariableFeatures(atac_subset_qc.rna)
atac_subset_qc.rna <- ScaleData(atac_subset_qc.rna)
atac_subset_qc.rna <- RunPCA(atac_subset_qc.rna)
atac_subset_qc.rna <- RunUMAP(atac_subset_qc.rna, dims = 1:30)

levels(atac_subset_qc.rna)
new.cluster.ids <- c('LS-1','Fb','LHr','LS-2','Bs',
  'Endo','T','Pc','Myl','B')
names(new.cluster.ids) <- levels(atac_subset_qc.rna)
atac_subset_qc.rna <- RenameIdents(atac_subset_qc.rna, new.cluster.ids)

atac_subset_qc.rna$celltype_RNA <- Idents(atac_subset_qc.rna)

# ------------------- DNA accessibility data processing -------------------
# process the DNA accessibility assay the same way we would process a scATAC-seq dataset, by performing latent semantic indexing (LSI).
atac_subset_qc <- readRDS('results/singlecell/atac_RNA_annotated.rds')
DefaultAssay(atac_subset_qc) <- "ATAC"

# We exclude the first dimension as this is typically correlated with sequencing depth
atac_subset_qc <- RunTFIDF(atac_subset_qc)
atac_subset_qc <- FindTopFeatures(atac_subset_qc, min.cutoff = "q5")
atac_subset_qc <- RunSVD(atac_subset_qc)

jpeg('results/singlecell/ATAC01_DepthCor_ATAC.jpg',width=7, height = 4, units = 'in', res = 300)
DepthCor(atac_subset_qc)
dev.off()

atac_subset_qc <- RunUMAP(atac_subset_qc, reduction = "lsi", dims = 2:30,
                          reduction.name = "umap.atac", reduction.key = "atacUMAP_")
atac_subset_qc <- FindNeighbors(object = atac_subset_qc, reduction = 'lsi', dims = 2:30)
atac_subset_qc <- FindClusters(object = atac_subset_qc, algorithm = 3, resolution = 0.15)
# Number of communities: 10

jpeg('results/singlecell/ATAC02_Dimplot_ATAC_bAnnot_10Clusters.jpg',width=6, height = 5, units = 'in', res = 100)
DimPlot(object = atac_subset_qc,reduction = "umap.atac", label = TRUE)+ ggtitle("ATAC")+ scale_color_manual(values = colors)+ NoLegend()
dev.off()

# ------------------- Replace the atac idents with the rna annoated idents-------------------
Idents(atac_subset_qc) <- Idents(atac_subset_qc.rna)
atac_subset_qc$celltype_ATAC <- Idents(atac_subset_qc)

#######################################################
#        Integrating scRNA-seq and scATAC-seq data    #
#######################################################
# Generate gene activities
gene.activities <- GeneActivity(atac_subset_qc,
                                max.width = 5e+08, biotypes =NULL, assay = 'ATAC')
#                                                   'PRKG1', 'MCTP1', 'ZNF385D',
#                                                   'MECOM', 'GLIS3', 'TPT1-AS1',
#                                                   'INPP4B', 'ERBB4'), biotypes =NULL)

dim(gene.activities)

# add gene activities as a new assay
atac_subset_qc[["ACTIVITY"]] <- CreateAssayObject(counts = gene.activities)

# normalize gene activities
DefaultAssay(atac_subset_qc) <- "ACTIVITY"
atac_subset_qc <- NormalizeData(atac_subset_qc)
atac_subset_qc <- ScaleData(atac_subset_qc)

saveRDS(atac_subset_qc, 'results/singlecell/atac_ATAC_RNAlabeled.rds')

atac_subset_qc <-readRDS('results/singlecell/atac_ATAC_RNAlabeled.rds')
Final_cell_barcode<-fread('results/singlecell/final_cell_barcodes.csv',header=TRUE)
atac_subset_qc <- subset(atac_subset_qc , cells = Final_cell_barcode$x)
DefaultAssay(atac_subset_qc)<-'ATAC'

# UMAP plots
jpeg('results/singlecell/ATAC03.0_DimPlot_ATAC_coloredbyOWN.jpg',width=10, height = 4, units = 'in', res = 200)
p12 <- DimPlot(atac_subset_qc, reduction = "umap.rna", group.by ='celltype_RNA', label = TRUE, label.size = 4, repel = TRUE)  +
  ggtitle("RNA")+ scale_color_manual(values = colors)+NoLegend()
p13 <- DimPlot(atac_subset_qc, reduction = "umap.atac", group.by ='celltype_ATAC',  label = TRUE, label.size = 4, repel = TRUE)  +
  ggtitle("ATAC")+ scale_color_manual(values = colors)+NoLegend()
p12+p13
dev.off()

jpeg('results/singlecell/ATAC03.1_DimPlot_ATAC_coloredbyATAC.jpg',width=10, height = 4, units = 'in', res = 200)
p12 <- DimPlot(atac_subset_qc, reduction = "umap.rna", group.by ='celltype_ATAC', label = TRUE, label.size = 4, repel = TRUE)  +
  ggtitle("RNA colored by ATAC")+ scale_color_manual(values = colors)+NoLegend()
p13 <- DimPlot(atac_subset_qc, reduction = "umap.atac", group.by ='celltype_ATAC',  label = TRUE, label.size = 4, repel = TRUE)  +
  ggtitle("ATAC")+ scale_color_manual(values = colors)+NoLegend()
p12+p13
dev.off()

jpeg('results/singlecell/ATAC03.2_DimPlot_ATAC_coloredbyRNA.jpg',width=10, height = 4, units = 'in', res = 200)
p12 <- DimPlot(atac_subset_qc, reduction = "umap.rna", group.by ='celltype_RNA', label = TRUE, label.size = 4, repel = TRUE)  +
  ggtitle("RNA")+ scale_color_manual(values = colors)+NoLegend()
p13 <- DimPlot(atac_subset_qc, reduction = "umap.atac", group.by ='celltype_RNA',  label = TRUE, label.size = 4, repel = TRUE)  +
  ggtitle("ATAC colored by RNA")+ scale_color_manual(values = colors)+NoLegend()
p12+p13
dev.off()

jpeg('results/singlecell/ATAC04_DimPlot_ATAC.jpg',width=6, height = 5, units = 'in', res = 300)
DimPlot(atac_subset_qc, reduction = "umap.atac", group.by ='celltype_ATAC',  label = TRUE, pt.size = 0.5)+
  ggtitle("ATAC")+ scale_color_manual(values = colors)+ NoLegend()
dev.off()

jpeg('results/singlecell/ATAC05_DimPlot_ATAC_Doublet.jpg',width=10, height = 4, units = 'in', res = 300)
p14 <- DimPlot(atac_subset_qc, reduction = "umap.atac", group.by='celltype_ATAC',label = TRUE, label.size = 4, repel = TRUE)  +
  ggtitle("ATAC")+ scale_color_manual(values = colors)+NoLegend()
p15 <- FeaturePlot(atac_subset_qc, features = "scDblFinder_score", reduction = "umap.atac") +
  ggtitle("Doublet Score")+
  theme(plot.title = element_text(hjust = 0.5), legend.position = "right",
        legend.title = element_text(size = 9, face = "bold"),
        legend.text = element_text(size = 8))  +
  labs(color = "DoubletScore")
p14+p15
dev.off()

# ---------LS1 VS LS2-------------------
DefaultAssay(atac_subset_qc)<-'ACTIVITY'
da_peaks_LS1_LS2 <- FindMarkers(
  object = atac_subset_qc,
  ident.1 = 'LS-1',
  ident.2 = 'LS-2',  # Compare against another level
  test.use = 'LR',
  latent.vars = "nCount_ATAC")
da_peaks_LS1_LS2$gene<-row.names(da_peaks_LS1_LS2)
write.csv(da_peaks_LS1_LS2, "results/singlecell/ATAC_da_results_LS1_2.csv", row.names = FALSE)

da_peaks_LS1_LS2<-fread("results/singlecell/ATAC_da_results_LS1_2.csv",header = TRUE)
top20_LS1_2 <- da_peaks_LS1_LS2 %>%
  arrange(p_val_adj) %>%  # Sort by significance
  head(20)  # Take the top 10 genes

# Loading relevant libraries
# Add a significance column for visualization
da_peaks_LS1_LS2 <- da_peaks_LS1_LS2 %>%
  mutate(Significance = case_when(
    p_val_adj < 0.05 & avg_log2FC > 1 ~ "Upregulated",
    p_val_adj < 0.05 & avg_log2FC < -1 ~ "Downregulated",
    TRUE ~ "Not Significant"
  ))
ggplot(da_peaks_LS1_LS2, aes(x = avg_log2FC, y = -log10(p_val_adj), color = Significance)) +
  geom_point(alpha = 0.8) +  # Plot points with transparency
  scale_color_manual(values = c("Upregulated" = "red4", "Downregulated" = "skyblue3", "Not Significant" = "gray")) +
  labs(title = "Volcano Plot of DA Peaks (LS-1 vs LS-2)", x = "log2 Fold Change", y = "-log10 Adjusted P-value") +
  theme_minimal() +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black") +  # FC cutoff
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black")+  # P-value cutoff
  geom_text(data = top20_LS1_2, aes(label = rownames(top20_LS1_2)),
            vjust = -0.5, size = 2.5, color = "black")  # Label top 10 genes

ggsave("results/singlecell/ATAC_volcano_plot_LS1_vs_LS2.png", width = 8, height = 6, dpi = 300)

#######################################################
#                     subset of atac                 #
#######################################################
atac_subset_qc<-readRDS('results/singlecell/atac_ATAC_RNAlabeled.rds')

DefaultAssay(atac_subset_qc)<-'ATAC'
# List all available graphs in the Seurat object
atac_subset_qc@graphs

table(atac_subset_qc$seurat_clusters)

# We exclude the first dimension as this is typically correlated with sequencing depth
atac_subset_qc <- RunTFIDF(atac_subset_qc)
atac_subset_qc <- FindTopFeatures(atac_subset_qc, min.cutoff = "q5")
atac_subset_qc <- RunSVD(atac_subset_qc)
atac_subset_qc <- RunUMAP(atac_subset_qc, reduction = "lsi", dims = 2:30,
                          reduction.name = "umap.atac", reduction.key = "atacUMAP_")
atac_subset_qc <- FindNeighbors(object = atac_subset_qc, reduction = 'lsi', dims = 2:30)
atac_subset_qc <- FindClusters(object = atac_subset_qc, algorithm = 3, resolution = 0.15)
# Number of communities: 10

atac_subset_qc_subcluster <- FindSubCluster(atac_subset_qc, cluster = 0, graph.name = "ATAC_snn", algorithm = 3, resolution =0.2)
Idents(atac_subset_qc_subcluster) <- "sub.cluster"
factor(atac_subset_qc_subcluster$sub.cluster)
atac_subset_qc_subcluster <- FindSubCluster(atac_subset_qc_subcluster, cluster = 1, graph.name = "ATAC_snn", algorithm = 3, resolution =0.1)
Idents(atac_subset_qc_subcluster) <- "sub.cluster"
factor(atac_subset_qc_subcluster$sub.cluster)
atac_subset_qc_subcluster <- FindSubCluster(atac_subset_qc_subcluster, cluster = 2, graph.name = "ATAC_snn", algorithm = 3, resolution =0.08)
Idents(atac_subset_qc_subcluster) <- "sub.cluster"
factor(atac_subset_qc_subcluster$sub.cluster)
atac_subset_qc_subcluster <- FindSubCluster(atac_subset_qc_subcluster, cluster = 3, graph.name = "ATAC_snn", algorithm = 3, resolution =0.08)
Idents(atac_subset_qc_subcluster) <- "sub.cluster"
factor(atac_subset_qc_subcluster$sub.cluster)

jpeg('results/singlecell/ATAC08_subcluster_b4Anno.jpg',width=6, height = 5, units = 'in', res = 100)
DimPlot(object = atac_subset_qc_subcluster,reduction = "umap.atac", group.by = 'sub.cluster',label = TRUE)+ ggtitle("ATAC")+ scale_color_manual(values = colors)+ NoLegend()
dev.off()

gene.activities_subcluster <- GeneActivity(atac_subset_qc_subcluster,
                                max.width = 5e+08, biotypes =NULL, assay = 'ATAC')
#                                                   'PRKG1', 'MCTP1', 'ZNF385D',
#                                                   'MECOM', 'GLIS3', 'TPT1-AS1',
#                                                   'INPP4B', 'ERBB4'), biotypes =NULL)

dim(gene.activities_subcluster)

# add gene activities as a new assay
atac_subset_qc_subcluster[["ACTIVITY"]] <- CreateAssayObject(counts = gene.activities_subcluster)

# normalize gene activities
DefaultAssay(atac_subset_qc_subcluster) <- "ACTIVITY"
atac_subset_qc_subcluster <- NormalizeData(atac_subset_qc_subcluster)
atac_subset_qc_subcluster <- ScaleData(atac_subset_qc_subcluster)

# find markers of the gene activities
levels(Idents(atac_subset_qc_subcluster))
atac_subset_qc.atac_markers_subcluster <- FindAllMarkers(atac_subset_qc_subcluster)

write.csv(atac_subset_qc.atac_markers_subcluster,'results/singlecell/Markers_ATAC_Subclusters.csv',row.names = FALSE)
atac_subset_qc.atac_markers_subcluster <-fread('results/singlecell/Markers_ATAC_Subclusters.csv',header = TRUE)

# Assigning cell type identity to clusters
atac_subset_qc_subcluster <- RenameIdents(atac_subset_qc_subcluster,
                               '0_0'='LS-1.1','0_1'='LS-1.2','0_2'='LS-1.3','1_0' = 'Fb-1','1_1'='Fb-2',
                               '2_0'='LHr-1','2_1'='LHr-2','3_0'='Bs-1','3_1'='Bs-2','4'='LS-2',
                               '5'='Endo','6'='T','7'='Pc','8'='Myl','9'='B')
atac_subset_qc_subcluster$celltype_ATAC_subcluster <- Idents(atac_subset_qc_subcluster)

atac_subset_qc_subcluster

saveRDS(atac_subset_qc_subcluster, 'results/singlecell/atac_ATAC_subcluster_annoatated.rds')

atac_subset_qc_subcluster<-readRDS('results/singlecell/atac_ATAC_subcluster_annoatated.rds')
DefaultAssay(atac_subset_qc_subcluster) <- "ACTIVITY"

jpeg('results/singlecell/ATAC09_DimPlot_ATAC.jpg',width=6, height = 5, units = 'in', res = 300)
DimPlot(atac_subset_qc_subcluster,  reduction = "umap.atac",group.by='celltype_ATAC_subcluster',label = TRUE, pt.size = 0.5)+
  ggtitle("ATAC")+ scale_color_manual(values = colors)+ NoLegend()
dev.off()

# ------------------------- Mitochondrial proportion ------------------------------
atac_subset_qc<-readRDS('results/singlecell/atac_ATAC_RNAlabeled.rds')
DefaultAssay(atac_subset_qc)<-'RNA'
mt_gene <- c("MT-ND1", "MT-ND2", "MT-ND3", "MT-ND4", "MT-ND4L", "MT-ND5", "MT-ND6",
             "MT-CO1", "MT-CO2", "MT-CO3", "MT-ATP6", "MT-ATP8", "MT-CYB")
mt_gene_det <- mt_gene[mt_gene %in% rownames(atac_subset_qc)]
atac_subset_qc[["percent.mt"]] <- PercentageFeatureSet(atac_subset_qc, pattern = "^MT-")

summary(atac_subset_qc$percent.mt)
# Min.  1st Qu.   Median     Mean  3rd Qu.     Max.

table(atac_subset_qc$percent.mt >= 20)

table(atac_subset_qc$celltype_RNA,atac_subset_qc$percent.mt>=20)
# luminal secretory-1 (LS-1)  4277  158
# fibroblast (Fb)             2566    7
# luminal HR (LHr)            2404    7
# luminal secretory-2 (LS-2)  1781    5
# basal epithelial (Bs)       1726   28
# endothelial (Endo)          1498    3
# T-cell (T)                   897    3
# Pericyte (Pc)                794   13
# Myeloid (Myl)                292    2
# B-cell (B)                   216    0

atac_filt <- subset(atac_subset_qc, subset = percent.mt < 20)
atac_filt

saveRDS(atac_filt , 'results/singlecell/ATAC_RNAlabelled_filterMT.rds')

atac_filt<-readRDS('results/singlecell/ATAC_RNAlabelled_filterMT.rds')

###################################################################################
# filter the seurat object by cell barcode
Final_cell_barcode <- colnames(atac_filt)
write.csv(Final_cell_barcode  , 'results/singlecell/final_cell_barcodes.csv',row.names = FALSE)

setMethod("seqinfo", "ChromatinAssay", function(x) x@seqinfo)
setMethod("genome",  "ChromatinAssay", function(x) {
  si <- x@seqinfo
  if (is.null(si)) NULL else GenomeInfoDb::genome(si)
})
library(data.table)
Final_cell_barcode<-fread('results/singlecell/final_cell_barcodes.csv',header=TRUE)
atac_subset_qc<-readRDS('results/singlecell/atac_ATAC_RNAlabeled.rds')
atac_subset_qc <- subset(atac_subset_qc , cells = Final_cell_barcode$x)
DefaultAssay(atac_subset_qc)<-'ATAC'

atac_subset_qc.rna <- readRDS('results/singlecell/atac_RNA_annotated.rds')
DefaultAssay(atac_subset_qc.rna)<-'RNA'
atac_subset_qc.rna[['RNA']]<-as(atac_subset_qc.rna[['RNA']], Class="Assay5")
atac_subset_qc.rna <- subset(atac_subset_qc.rna , cells = Final_cell_barcode$x)

levels(atac_subset_qc.rna)
#                      'Endo','T','Pc','Myl','B')
new.cluster.ids <- c('LASP-1','Fb','LHS','LASP-2','B-Myo',
                     'Endo','T','Pc','Myl','B')
names(new.cluster.ids) <- levels(atac_subset_qc.rna)
atac_subset_qc.rna <- RenameIdents(atac_subset_qc.rna, new.cluster.ids)
atac_subset_qc.rna$celltype_RNA <- Idents(atac_subset_qc.rna)

names(new.cluster.ids) <- levels(atac_subset_qc)
atac_subset_qc <- RenameIdents(atac_subset_qc, new.cluster.ids)
atac_subset_qc$celltype_ATAC <- Idents(atac_subset_qc)

# ------------------ Dot plot ----------------------
# cell-type display order
fig5D_levels <- c('LASP-1','LASP-2','LHS','B-Myo','Fb','Endo','T','B','Myl','Pc')
# Reorder the Idents based on sorted levels
Idents(atac_subset_qc) <- factor(Idents(atac_subset_qc), levels = fig5D_levels)
Idents(atac_subset_qc.rna) <- factor(Idents(atac_subset_qc.rna), levels = fig5D_levels)
atac_subset_qc$celltype_ATAC <- Idents(atac_subset_qc)
atac_subset_qc.rna$celltype_RNA <- Idents(atac_subset_qc.rna)

# UMAP + stacked barplot
library(ggplot2)
library(dplyr)

rna_data <- data.frame(
  rnaUMAP_1 = atac_subset_qc.rna@reductions$umap.rna@cell.embeddings[, 1],
  rnaUMAP_2 = atac_subset_qc.rna@reductions$umap.rna@cell.embeddings[, 2],
  celltype_RNA = atac_subset_qc.rna$celltype_RNA)
centroids_rna <- rna_data %>%
  group_by(celltype_RNA) %>%
  summarize(
    rnaUMAP_1 = mean(rnaUMAP_1),
    rnaUMAP_2 = mean(rnaUMAP_2))

atac_data <- data.frame(
  atacUMAP_1 = atac_subset_qc@reductions$umap.atac@cell.embeddings[, 1],
  atacUMAP_2 = atac_subset_qc@reductions$umap.atac@cell.embeddings[, 2],
  celltype_ATAC = atac_subset_qc$celltype_ATAC)
centroids_atac <- atac_data %>%
  group_by(celltype_ATAC) %>%
  summarize(
    atacUMAP_1 = mean(atacUMAP_1),
    atacUMAP_2 = mean(atacUMAP_2))

# lollipod figure

# Extract metadata from the Seurat object
metadata <- atac_subset_qc.rna@meta.data

# Count the number of cells per cell type
summary_data <- metadata %>%
  group_by(celltype_RNA) %>%
  summarise(num_cells = n())
summary_data$celltype_RNA <- factor(summary_data$celltype_RNA, levels = rev(fig5D_levels))

# Create the lollipop dot chart
library(ggpubr)
ggdot<-ggdotchart(summary_data, x = "celltype_RNA", y = "num_cells",
           add = "segments",                             # Add segments from y = 0 to dots
           add.params = list(color = "lightgray", size = 0.8), # Change segment color and size
           color = colors,                                # Color by groups
           sorting = "none",
           rotate = TRUE,                                # Rotate vertically
           dot.size = 6,                                 # Large dot size
           font.label = list(color = "white", size = 7,
                             vjust = 0.5),               # Adjust label parameters
           ggtheme = theme_pubr(),                      # Use ggpubr theme
           legend = "right"                             # Add legend on the right-hand side
)+
  labs(x = "",
       y='',
      title = "Number of cells")+
  theme(panel.border = element_rect(color = "black", fill = NA, size = 0.7))  # Add border       # Center the title

axis <- ggh4x::guide_axis_truncated(
  trunc_lower = unit(0, "npc"),
  trunc_upper = unit(2.5, "cm")
)

library(patchwork)
p_rna <- ggplot(rna_data, aes(x = rnaUMAP_1, y = rnaUMAP_2, color = celltype_RNA)) +
  geom_point()+
  geom_text(
    data = centroids_rna,  # Use centroids for labeling
    aes(label = celltype_RNA),  # Use celltype_RNA as labels
    color = "black",  # Label color
    size = 5,  # Label size
    vjust = -1.5, hjust = -0.8  # Adjust label position relative to points
  ) +
  theme(axis.line = element_line(arrow = arrow()),
        axis.title = element_text(hjust = 0))+
  theme(
    axis.line = element_line(arrow = arrow(type = 'closed', length = unit(10, 'pt'))),  # Add closed arrows to axes
    panel.grid = element_blank(),  # Remove background grid
    panel.background = element_blank(),  # Remove panel background color
    plot.background = element_blank(),  # Remove plot background color
    axis.text.x = element_blank(),  # Remove x-axis tick labels
    axis.text.y = element_blank(),
    axis.ticks.x = element_blank(),
    axis.ticks.y = element_blank(),
    legend.position = "botoom"  # Remove legend
  ) +
  labs(x = "UMAP 1", y = "UMAP 2")+
  scale_color_manual(values = colors) +  # Use custom colors
  ggtitle("scRNA-seq")

p_atac <- ggplot(atac_data, aes(x = atacUMAP_1, y = atacUMAP_2, color = celltype_ATAC)) +
  geom_point()+
  geom_text(
    data = centroids_atac,  # Use centroids for labeling
    aes(label = celltype_ATAC),  # Use celltype_RNA as labels
    color = "black",  # Label color
    size = 5,  # Label size
    vjust = -1.0, hjust = -0.1  # Adjust label position relative to points
  )  +
  theme(axis.line = element_line(arrow = arrow()),
        axis.title = element_text(hjust = 0))+
  theme(
    axis.line = element_line(arrow = arrow(type = 'closed', length = unit(10, 'pt'))),  # Add closed arrows to axes
    panel.grid = element_blank(),  # Remove background grid
    panel.background = element_blank(),  # Remove panel background color
    plot.background = element_blank(),  # Remove plot background color
    axis.text.x = element_blank(),  # Remove x-axis tick labels
    axis.text.y = element_blank(),
    axis.ticks.x = element_blank(),
    axis.ticks.y = element_blank(),
    legend.position = "none"  # Remove legend
  ) +
  labs(x = "UMAP 1", y = "UMAP 2") +
  scale_color_manual(values = colors) +  # Use custom colors
  ggtitle("scATAC-seq")

jpeg('results/singlecell/ATACfiltMT_UMAP_barchart.jpg',width=13, height = 5.5, units = 'in', res = 200)
(p_rna| p_atac | ggdot) + plot_layout(widths = c(1.6, 1.6, 0.8))  # Set the ratio of widths
dev.off()

jpeg('results/singlecell/ATACfiltMT_scRNA.jpg',width=5.5, height = 5.6, units = 'in', res = 200)
p_rna
dev.off()

jpeg('results/singlecell/ATACfiltMT_scATAC.jpg',width=6, height = 6, units = 'in', res = 200)
p_atac
dev.off()

jpeg('results/singlecell/ATACfiltMT_barchart.jpg',width=3.5, height = 6, units = 'in', res = 200)
ggdot
dev.off()

# ------------------ Dot plot ----------------------
Final_cell_barcode<-fread('results/singlecell/final_cell_barcodes.csv',header=TRUE)

atac_subset_qc<-readRDS('results/singlecell/atac_ATAC_RNAlabeled.rds')
atac_subset_qc <- subset(atac_subset_qc , cells = Final_cell_barcode$x)
DefaultAssay(atac_subset_qc)<-'ATAC'

atac_subset_qc.rna <- readRDS('results/singlecell/atac_RNA_annotated.rds')
DefaultAssay(atac_subset_qc.rna)<-'RNA'
atac_subset_qc.rna[['RNA']]<-as(atac_subset_qc.rna[['RNA']], Class="Assay5")
atac_subset_qc.rna <- subset(atac_subset_qc.rna , cells = Final_cell_barcode$x)

atac_subset_qc.rna_markers<-fread('results/singlecell/Markers_RNA_DEGs_perpaper.csv',header=TRUE)
.rowNamesDF(atac_subset_qc.rna_markers, make.names=TRUE)<-atac_subset_qc.rna_markers$gene

levels(atac_subset_qc.rna)
new.cluster.ids <- c('LASP-1','Fb','LHS','LASP-2','B-Myo',
                     'Endo','T','Pc','Myl','B')
names(new.cluster.ids) <- levels(atac_subset_qc.rna)
atac_subset_qc.rna <- RenameIdents(atac_subset_qc.rna, new.cluster.ids)

atac_subset_qc.rna$celltype_RNA <- Idents(atac_subset_qc.rna)

# ------------------ Dot plot ----------------------
# cell-type display order
fig5D_levels <- c('LASP-1','LASP-2','LHS','B-Myo','Fb','Endo','T','B','Myl','Pc')
# Reorder the Idents based on sorted levels
Idents(atac_subset_qc) <- factor(Idents(atac_subset_qc), levels = fig5D_levels)
Idents(atac_subset_qc.rna) <- factor(Idents(atac_subset_qc.rna), levels = fig5D_levels)
atac_subset_qc$celltype_ATAC <- Idents(atac_subset_qc)
atac_subset_qc.rna$celltype_RNA <- Idents(atac_subset_qc.rna)

# remove pct=0 or not available in activity dataset or duplicated gene
gene_protcoding_lincRNA<-fread("results/singlecell/gene_list_proteincoding_lincRNA.csv",header = TRUE)
atac_subset_qc.rna_markers <- atac_subset_qc.rna_markers[atac_subset_qc.rna_markers$gene %in% gene_protcoding_lincRNA$gene_name,]

atac_subset_qc.rna_markers.filt<-atac_subset_qc.rna_markers[!atac_subset_qc.rna_markers$gene%in%c('TPT1-AS1'),]

# top 5 markers
atac_subset_qc.rna_markers.filt %>%
  group_by(cluster) %>%
  arrange(p_val_adj,desc(avg_log2FC))%>%
  slice_head(n = 5)  %>%
  ungroup() -> top5

cluster<-c('LASP-1','Fb','LHS','LASP-2','B-Myo','Endo','T','Pc','Myl','B')
top5$cluster_name <- rep(cluster, each = 5)
top5 <- top5 %>% arrange(factor(cluster_name, levels = rev(fig5D_levels)))

# DotPlot for RNA Expression
DefaultAssay(atac_subset_qc.rna) <- "RNA"
p16_data <- DotPlot(atac_subset_qc.rna, features = top5$gene, group.by = "celltype_RNA")$data
p16_data$category <- "RNA expression"

library(ggplot2)
p16<-ggplot(p16_data, aes(x = id, y = features.plot)) +
  geom_point(aes(size = pct.exp, color = avg.exp.scaled)) +
  scale_color_gradient2(low = "plum3", mid = "grey90", high = "red2", midpoint = 0) +
  scale_size(range = c(0.6, 6)) +
  facet_grid(. ~ category, scales = "free", space = "free")  +
  theme_minimal() +  # Start with minimal theme
  theme(
    axis.ticks = element_line(color = "black"),
    axis.text = element_text(color = "black"),
    axis.title = element_text(color = "black"),
    strip.text = element_text(size = 14, face = "bold", color = "black"),
    panel.grid = element_blank(),  # Remove grid lines
    panel.spacing = unit(0.5, "lines"),
    # panel.background = element_blank(),  # Remove background
    # panel.border = element_blank(),  # No default panel border
    legend.position = "bottom"
  )   + labs(
    color = "",  # Title for color legend
    size = "Pct"   # Title for size legend
  )

p16 <- p16 +
  # Add border lines for the ATAC panel (all sides)
  annotation_custom(grid::linesGrob(y = c(1, 1), x = c(0, 1), gp = grid::gpar(col = "black", lwd = 1))) +  # Bottom line
  annotation_custom(grid::linesGrob(y = c(0, 1), x = c(1, 1), gp = grid::gpar(col = "black", lwd = 1))) +  # Right line
  annotation_custom(grid::linesGrob(y = c(0, 0), x = c(0, 1), gp = grid::gpar(col = "black", lwd = 1))) +  # Top line
  annotation_custom(grid::linesGrob(y = c(0, 1), x = c(0, 0), gp = grid::gpar(col = "black", lwd = 1)))  # Left line for ATAC panel

jpeg('results/singlecell/ATAC010_Dotplot_top5DEG_RNA.jpg', width = 6, height = 10, units = 'in', res = 200)
p16
dev.off()

# DotPlot for Gene Activity
DefaultAssay(atac_subset_qc) <- "ACTIVITY"
p17_data <- DotPlot(atac_subset_qc, features = top5$gene, group.by = "celltype_ATAC")$data
# The following requested variables were not found: LAMA2, MEG3, PRKG1, MCTP1, ZNF385D, MECOM, GLIS3, TPT1-AS1, INPP4B, ERBB4
p17_data <- p17_data %>% filter(!is.na(features.plot)) # Remove NA values
p17_data$category <- "Gene Activity (ATAC)"

p17<-ggplot(p17_data, aes(x = id, y = features.plot)) +
  geom_point(aes(size = pct.exp, color = avg.exp.scaled)) +
  scale_color_gradient2(low = "plum3", mid = "grey90", high = "red2", midpoint = 0) +
  scale_size(range = c(0.1, 6)) +
  facet_grid(. ~ category, scales = "free", space = "free")  +
  theme_minimal() +  # Start with minimal theme
  theme(
    axis.ticks = element_line(color = "black"),  # Keep x-axis ticks
    axis.text.x = element_text(color = "black"),  # Keep x-axis text
    axis.text.y = element_blank(),  # Remove y-axis text
    axis.title.x = element_text(color = "black"),  # Keep x-axis title
    axis.title.y = element_blank(),  # Remove y-axis title
    axis.ticks.y = element_blank(),  # Remove y-axis ticks
    strip.text = element_text(size = 14, face = "bold", color = "black"),
    panel.grid = element_blank(),  # Remove grid lines
    panel.spacing = unit(0.5, "lines"),
    legend.position = "bottom"
  )+ labs(
    x = "Cell type",
    y = "Gene",
    color = "",
    size = "Pct"
  )

p17 <- p17 +
  # Add border lines for the ATAC panel (all sides)
  annotation_custom(grid::linesGrob(y = c(1, 1), x = c(0, 1), gp = grid::gpar(col = "black", lwd = 1))) +  # Bottom line
  annotation_custom(grid::linesGrob(y = c(0, 1), x = c(1, 1), gp = grid::gpar(col = "black", lwd = 1))) +  # Right line
  annotation_custom(grid::linesGrob(y = c(0, 0), x = c(0, 1), gp = grid::gpar(col = "black", lwd = 1))) +  # Top line
  annotation_custom(grid::linesGrob(y = c(0, 1), x = c(0, 0), gp = grid::gpar(col = "black", lwd = 1)))  # Left line for ATAC panel

library(patchwork)
jpeg('results/singlecell/ATAC011_Dotplot_top5DEG_RNA_ATAC.jpg', width = 11, height = 11, units = 'in', res = 200)
(p16+p17) +
  plot_layout(widths = c(0.93,1))  # Set the ratio of widths
dev.off()

# Combine datasets
combined_data <- rbind(p16_data, p17_data)
combined_data$category <- factor(combined_data$category, levels = c("RNA expression", "Gene Activity (ATAC)"))

# Create Faceted DotPlot
p16_17<-ggplot(combined_data, aes(x = id, y = features.plot)) +
  geom_point(aes(size = pct.exp, color = avg.exp.scaled)) +
  scale_color_gradient2(low = "plum3", mid = "grey90", high = "red2", midpoint = 0) +
  scale_size(range = c(0.6, 6)) +
  facet_grid(. ~ category, scales = "free", space = "free")  +
  theme_minimal() +  # Start with minimal theme
  theme(
    axis.ticks = element_line(color = "black"),
    axis.text = element_text(color = "black"),
    axis.title = element_text(color = "black"),
    strip.text = element_text(size = 14, face = "bold", color = "black"),
    panel.grid = element_blank(),  # Remove grid lines
    panel.spacing = unit(0.5, "lines"),
    panel.background = element_blank(),  # Remove background
    panel.border = element_blank(),  # No default panel border
    legend.position = "bottom"
  )

p16_17 <- p16_17 +
  # Add border lines for the ATAC panel (all sides)
  annotation_custom(grid::linesGrob(y = c(1, 1), x = c(0, 1), gp = grid::gpar(col = "black", lwd = 1))) +  # Bottom line
  annotation_custom(grid::linesGrob(y = c(0, 0), x = c(0, 1), gp = grid::gpar(col = "black", lwd = 1))) +  # Top line
  annotation_custom(grid::linesGrob(y = c(0, 1), x = c(0, 0), gp = grid::gpar(col = "black", lwd = 1)))  # Left line for ATAC panel

jpeg('results/singlecell/ATAC09_Dotplot_top5_RNA_ATAC.jpg', width = 9, height = 10, units = 'in', res = 200)
p16_17
dev.off()

########################### rotated panel #########################
library(ggplot2)
library(dplyr)

combined_data <- rbind(p16_data, p17_data)
combined_data$category <- factor(combined_data$category,
                                 levels = c("RNA expression", "Gene Activity (ATAC)"))  # RNA on top
combined_data$id <- factor(combined_data$id, levels = rev(fig5D_levels))      # LASP-1 on top
combined_data$features.plot <- factor(combined_data$features.plot,
                                      levels = rev(unique(top5$gene)))         # CCL28 leftmost
p16_17 <- ggplot(combined_data, aes(x = features.plot, y = id)) +
  geom_point(aes(size = pct.exp, color = avg.exp.scaled)) +
  scale_color_gradient2(low = "plum3", mid = "grey90", high = "red2", midpoint = 0) +
  scale_size(range = c(0.6, 6)) +
  facet_grid(category ~ ., scales = "free", space = "free", switch = "y") +  # strips on the right
  theme_minimal() +
  theme(
    axis.ticks       = element_line(color = "black"),
    axis.text.x      = element_text(color = "black", angle = 45, hjust = 1, vjust = 1,
                                    face = "italic"),
    axis.text.y      = element_text(color = "black"),
    axis.title       = element_text(color = "black"),
    strip.text.y     = element_text(size = 12, face = "bold", color = "black", angle = 0),
    strip.placement  = "outside",
    panel.grid       = element_blank(),
    panel.spacing    = unit(0.5, "lines"),
    panel.background = element_blank(),
    panel.border     = element_rect(fill = NA, colour = "black", linewidth = 0.6),
    legend.position  = "bottom",
    legend.box       = "horizontal"
  ) +
  labs(x = "Gene", y = "Cell type", color = "Gene expression (Avg.)", size = "Pct.")

jpeg('results/singlecell/ATAC09_Dotplot_top5_RNA_ATAC_rotated.jpg',
     width = 14, height = 7, units = 'in', res = 300)
p16_17
dev.off()

################################ volcano plot ####################
cluster_ids <- c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9)
labels <- c('LS-1', 'Fb', 'LHr', 'LS-2', 'Bs', 'Endo', 'T', 'Pc', 'Myl', 'B')

# Create a named vector for mapping clusters to labels
cluster_mapping <- setNames(labels, cluster_ids)

# Assign labels to the cluster column in the data frame
atac_subset_qc.rna_markers.filt$cluster_label <- cluster_mapping[as.character(atac_subset_qc.rna_markers.filt$cluster)]

all_DEGs<-fread('results/singlecell/Markers_RNA_DEGs_all.csv', header = TRUE)
LS1 <- all_DEGs %>% filter(cluster =='LS-1')
LS1$p_val_adj <-ifelse(LS1$p_val_adj==0,1e-300,LS1$p_val_adj)

top20_LS1 <- LS1 %>%
  arrange(p_val_adj) %>%  # Sort by significance
  head(20)  # Take the top 20 genes

# Add a significance column for visualization
LS1 <- LS1 %>%
  mutate(Significance = case_when(
    p_val_adj < 0.05 & avg_log2FC > 1.5 ~ "Upregulated",
    p_val_adj < 0.05 & avg_log2FC < -1.5 ~ "Downregulated",
    TRUE ~ "Not Significant"
  ))
ggplot(LS1, aes(x = avg_log2FC, y = -log10(p_val_adj), color = Significance)) +
  geom_point(alpha = 0.8) +  # Plot points with transparency
  scale_color_manual(values = c("Upregulated" = "red4", "Downregulated" = "skyblue3", "Not Significant" = "gray")) +
  labs(title = "Volcano Plot of DEG (LS-1 vs others)", x = "log2 Fold Change", y = "-log10 Adjusted P-value") +
  theme_minimal() +
  geom_vline(xintercept = c(-1.5, 1.5), linetype = "dashed", color = "black") +  # FC cutoff
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black")+  # P-value cutoff
  geom_text(data = top20_LS1, aes(label = top20_LS1$gene),
            vjust = -0.5, size = 2.5, color = "black")  # Label top 10 genes

ggsave("results/singlecell/ATAC_volcano_plot_LS1_vs_LS2.png", width = 8, height = 6, dpi = 300)

#------------------ cell-type marker panels ------------------------------
features <- rev(c("ANKRD30A", "ERBB4","SYTL2", "INPP4B","AFF3",
                  "IGF2BP2","ALDH1A3", "COBL", "PIGR", "CCL28",
                  "NRG1","KRT14", "ACTA2", "PTPRT", "MYLK",
                  "PTPRC", "PARP8","FYN", "STAT4","AOAH",
                  "ENAM", "SSR4","MZB1","DERL3",  "HERPUD1",
                  "HLA-DRA","CD74", "HLA-DPA1", "HLA-DRB1","HLA-DPB1",
                  "MCTP1", "ADAMTS9", "ZNF385D","ELTD1",'SELE',
                  'C11orf96',"MT1A",'RGS6','IGFBP5','PRKG1',
                  'DCN',"APOD", 'HPSE2',"CFD","LAMA2"))

# DotPlot for RNA Expression
DefaultAssay(atac_subset_qc.rna) <- "RNA"
p18_data <- DotPlot(atac_subset_qc.rna, features = features, group.by = "celltype_RNA")$data
p18_data$category <- "RNA expression"

# DotPlot for Gene Activity
DefaultAssay(atac_subset_qc) <- "ACTIVITY"
p19_data <- DotPlot(atac_subset_qc, features = features, group.by = "celltype_ATAC")$data
p19_data <- p19_data %>% filter(!is.na(features.plot)) # Remove NA values
p19_data$category <- "Gene Activity (ATAC)"

#   The following requested variables were not found (10 out of 11 shown): LAMA2, HPSE2, PRKG1, RGS6, ZNF385D, MCTP1, PTPRT, NRG1, AFF3, INPP4B

# Combine datasets
combined_data <- rbind(p18_data, p19_data)
combined_data$category <- factor(combined_data$category, levels = c("RNA expression", "Gene Activity (ATAC)"))

# Create Faceted DotPlot
jpeg('results/singlecell/ATAC10_Dotplot_Fig5D.jpg', width = 8, height = 12, units = 'in', res = 200)
p18_19<-ggplot(combined_data, aes(x = id, y = features.plot)) +
  geom_point(aes(size = pct.exp, color = avg.exp.scaled)) +
  scale_color_gradient2(low = "plum3", mid = "grey90", high = "red", midpoint = 0) +
  scale_size(range = c(0.3, 6)) +
  facet_grid(. ~ category, scales = "free", space = "free")  +
  theme_minimal() +  # Start with minimal theme
  theme(
    axis.ticks = element_line(color = "black"),
    axis.text = element_text(color = "black"),
    axis.title = element_text(color = "black"),
    strip.text = element_text(size = 14, face = "bold", color = "black"),
    panel.grid = element_blank(),  # Remove grid lines
    panel.spacing = unit(0.5, "lines"),
    panel.background = element_blank(),  # Remove background
    panel.border = element_blank(),  # No default panel border
    legend.position = "bottom"
  )

p18_19 <- p18_19 +
  # Add border lines for the ATAC panel (all sides)
  annotation_custom(grid::linesGrob(y = c(1, 1), x = c(0, 1), gp = grid::gpar(col = "black", lwd = 1))) +  # Bottom line
  annotation_custom(grid::linesGrob(y = c(0, 0), x = c(0, 1), gp = grid::gpar(col = "black", lwd = 1))) +  # Top line
  annotation_custom(grid::linesGrob(y = c(0, 1), x = c(0, 0), gp = grid::gpar(col = "black", lwd = 1)))  # Left line for ATAC panel

p18_19
dev.off()

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# Heatmap showing Pearson correlation coefficient (cor) of RNA expression and gene activity across different cell types
library(pheatmap)

# Ensure correct assay selection
DefaultAssay(atac_subset_qc.rna) <- "RNA"
rna_variable_genes <- FindVariableFeatures(atac_subset_qc.rna,
                                          nfeatures = 3000)
rna_variable_genes <- VariableFeatures(rna_variable_genes)

DefaultAssay(atac_subset_qc) <- "ACTIVITY"  # Gene activity assay
ga_variable_genes <- FindVariableFeatures(atac_subset_qc,
                                          nfeatures = 3000)
ga_variable_genes <- VariableFeatures(ga_variable_genes)

common_genes <- intersect(rna_variable_genes, ga_variable_genes)[1:850]

# Extract RNA expression per cell type
rna_expr <- AggregateExpression(
  atac_subset_qc.rna
)$RNA

# Aggregate (pseudo-bulk) gene activity per cell type
gene_activity <- AggregateExpression(
  atac_subset_qc
)$ACTIVITY

# Filter RNA and ATAC matrices to include only the top 850 variable genes
rna_expr <- as.matrix(rna_expr[common_genes, ])
gene_activity <- as.matrix(gene_activity[common_genes, ])

fig5D_levels <- c('LHr','LS-1','LS-2','Bs','T','B','Myl','Endo','Pc','Fb')
rna_expr <-rna_expr[,rev(fig5D_levels)]
gene_activity <-gene_activity[,(fig5D_levels)]

# Compute Pearson correlation between RNA expression and ATAC gene activity
cor <- cor(rna_expr, gene_activity,
                  method = "pearson")
cor <- as.matrix(cor)

# Define color palette
color_palette <- colorRampPalette(c( "blue1","#FFFFCC","red3"))

heatmap(cor,
        Rowv = NA,
        Colv = NA,
        col = color_palette(500),
        xlab = "Gene activity (from ATAC)",
        ylab =  "Gene expression (RNA)")
