# ==============================================================================
# 04_annotate_cell_types.R
# Cluster the scRNA-seq data and annotate cell types from canonical markers.
# Produces: Fig 2A
# Inputs  : results/singlecell/atac_postQC.rds
#           results/singlecell/Markers_RNA_DEGs_perpaper.csv
#           results/singlecell/atac_RNA_annotated.rds
#           results/singlecell/final_cell_barcodes.csv
#           results/singlecell/Markers_RNA_DEGs_perpaper_labelled.csv
# Outputs : results/singlecell/RNA01_ElbowPlot_RNA.jpg
#           results/singlecell/RNA01.2_ElbowPlot_lookforPC.jpg
#           results/singlecell/RNA01.3_clustree_testing_resolution.jpg
#           results/singlecell/RNA02_Dimplot_RNA_bAnnot_10Clusters.jpg
#           results/singlecell/Markers_RNA_DEGs_perpaper.csv
#           results/singlecell/Markers_RNA_DEGs_all.csv
#           results/singlecell/Markers_RNA_Compared_cluster3.csv
#           ...
# Run     : Rscript scripts/01_singlecell_annotation/04_annotate_cell_types.R   (from the repository root)
# ==============================================================================
#   "#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd",
#   "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22",
#   "#17becf",
#   "#393b79", "#637939", "#8c6d31", "#843c39", "#7b4173",
#   "#5254a3", "#6b6ecf", "#9c9ede", "#e7969c", "#e7ba52",
#   "#bd9e39"
colors <-c( "#E64B35CC","#4DBBD5CC", "#00A087CC", "#B09C85CC","#F39B7FCC", "#8491B4CC", "#91D1C2CC",'#f9e79f', "#7E6148CC",  "#3C5488CC")

#####################################################################################
# ---------- Finding differentially expressed features (cluster biomarkers) -----  ##
#####################################################################################
#------------------- Gene expression data processing------------------------------
atac_subset_qc <- readRDS('results/singlecell/atac_postQC.rds')
DefaultAssay(atac_subset_qc) <- "RNA"
options(future.globals.maxSize = 8000 * 1024^2) # avoid Error in getGlobalsAndPackages(expr, envir = envir, globals = globals)

atac_subset_qc <- SCTransform(atac_subset_qc)
atac_subset_qc <- RunPCA(atac_subset_qc)

jpeg('results/singlecell/RNA01_ElbowPlot_RNA.jpg',width=8, height = 5, units = 'in', res = 200)
ElbowPlot(atac_subset_qc, ndims = 50)
dev.off()

# Deriving the elbow
# Determine percent of variation associated with each PC
pct <- atac_subset_qc[["pca"]]@stdev/sum(atac_subset_qc[["pca"]]@stdev) * 100
# Calculate cumulative percents for each PC
cumu <- cumsum(pct)
# Determine which PC exhibits cumulative percent greater than 90% and %
# variation associated with the PC as less than 5
co1 <- which(cumu > 90 & pct < 5)[1]
# Determine the difference between variation of PC and subsequent PC, last
# point where change of % of variation is more than 0.1%.
co2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1
pc <- min(co1, co2)
pc
plot_df <- data.frame(pct = pct, cumu = cumu, rank = 1:length(pct))
jpeg('results/singlecell/RNA01.2_ElbowPlot_lookforPC.jpg',width=8, height = 5, units = 'in', res = 200)
ggplot(plot_df, aes(cumu, pct, label = rank, color = rank > pc)) + geom_text() +
  geom_vline(xintercept = 90, color = "grey") + geom_hline(yintercept = min(pct[pct >
                                                                                  5]), color = "grey") + theme_bw()
dev.off()

atac_subset_qc <- FindNeighbors(atac_subset_qc, dims = 1:20)
atac_subset_qc <- RunUMAP(atac_subset_qc, reduction.name = 'umap.rna', reduction.key = 'rnaUMAP_',dims = 1:20)

# Testing resolutions
reso <- c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1)
reso_res <- lapply(1:length(reso), function(x, atac_subset_qc, reso) {
  atac_subset_qc <- FindClusters(atac_subset_qc, resolution = reso[x])
  clust <- setNames(atac_subset_qc$seurat_clusters, Cells(atac_subset_qc))
  return(clust)
}, atac_subset_qc, reso)
k_tab <- do.call(cbind, reso_res)
k_dat <- as.data.frame(k_tab)
head(k_dat, 2)
names(k_dat) <- paste('k',reso,sep = '')

jpeg('results/singlecell/RNA01.3_clustree_testing_resolution.jpg',width=8, height = 5, units = 'in', res = 200)
clustree(k_dat, prefix = "k", node_colour = "sc3_stability")
dev.off()

atac_subset_qc <- FindClusters(atac_subset_qc, resolution = 0.1)
# 10 clusters

# Look at cluster IDs of the first 5 cells
head(Idents(atac_subset_qc), 5)
# ACAGCGGGTAACGTGA-1 ACAGCGGGTCGGTTTG-1 ACAGCGGGTTCCTCAC-1 ACAGCGGGTTGTTCCC-1 ACAGCGGGTGCGTACC-1

jpeg('results/singlecell/RNA02_Dimplot_RNA_bAnnot_10Clusters.jpg',width=6, height = 5, units = 'in', res = 100)
DimPlot(atac_subset_qc, reduction = "umap.rna", label = TRUE)+ scale_color_manual(values = colors)+ NoLegend()
dev.off()

# ------------ find markers for every cluster compared to all remaining cells ------------
# Find the minimum number of cells per cluster to downsample
cell_counts <- table(Idents(atac_subset_qc))
DefaultAssay(atac_subset_qc)<-'RNA'
atac_subset_qc.rna_markers <- FindAllMarkers(
  object = atac_subset_qc,
  min.diff.pct = 0.5  # Ensuring genes are significantly different between groups
)

atac_subset_qc.rna_markers <- atac_subset_qc.rna_markers %>%
  filter(p_val_adj < 0.05 & abs(avg_log2FC) > 1.5)

write.csv(atac_subset_qc.rna_markers,'results/singlecell/Markers_RNA_DEGs_perpaper.csv',row.names = FALSE)

atac_subset_qc.rna_markers <- FindAllMarkers(
  object = atac_subset_qc)
write.csv(atac_subset_qc.rna_markers,'results/singlecell/Markers_RNA_DEGs_all.csv',row.names = FALSE)

atac_subset_qc.rna_markers<-fread('results/singlecell/Markers_RNA_DEGs_perpaper.csv',header=TRUE)
.rowNamesDF(atac_subset_qc.rna_markers, make.names=TRUE)<-atac_subset_qc.rna_markers$gene
head(row.names(atac_subset_qc.rna_markers))

# find all markers distinguishing cluster 3 from clusters 0
cluster3.markers <- FindMarkers(atac_subset_qc, ident.1 = 3, ident.2 = 0)
head(cluster3.markers, n = 10)
# --------- doublet score<=0.95 --------------------
# p_val avg_log2FC pct.1 pct.2 p_val_adj
# PTN            0   4.629776 0.966 0.166         0
# MET            0   4.136893 0.945 0.176         0
# GLIS3          0   4.177809 0.962 0.210         0
# KRT14          0   4.821995 0.767 0.142         0
# ST6GALNAC5     0   3.659692 0.729 0.134         0
# ETS1           0   2.680570 0.798 0.217         0
# TSHZ2          0   2.897592 0.785 0.209         0
# WFDC2          0   2.334933 0.859 0.335         0
# TRABD2B        0   4.637473 0.566 0.049         0
# CAV1           0   3.978085 0.577 0.063         0

write.csv(cluster3.markers,'results/singlecell/Markers_RNA_Compared_cluster3.csv',row.names = FALSE)

jpeg('results/singlecell/RNA04_FeaturePlot_FOXA1_TP63.jpg',width=9, height = 4, units = 'in', res = 150)
FeaturePlot(atac_subset_qc,features = c('FOXA1','TP63'))
dev.off()

jpeg('results/singlecell/RNA05_FeaturePlot_figS5D.jpg',width=18, height = 13, units = 'in', res = 50)
FeaturePlot(atac_subset_qc,features = c('ANKRD30A','ERBB4','PROM1','PTN','KRT14',
                                        'ACTA2','RGS5','PTPRC','CD3E','IL7R','CD79A',
                                        'IGJ','HLA-DRA','VWF','COL1A2','LUM'))
dev.off()

atac_subset_qc.rna_markers %>%
  group_by(cluster) %>%
  dplyr::filter(adj_p_val,des(avg_log2FC)) %>%
  slice_head(n = 5)  %>%
  ungroup() -> top5

jpeg('results/singlecell/RNA06_DoHeatmap_top5gene.jpg',width=10, height = 15, units = 'in', res = 150)
DoHeatmap(atac_subset_qc, features = top5$gene) + NoLegend()
dev.off()

# Assigning cell type identity to clusters
atac_subset_qc <- RenameIdents(atac_subset_qc,
                               '0'='luminal secretory-1 (LS-1)','1' = 'fibroblast (Fb)',
                               '2'='luminal HR (LHr)','3'='luminal secretory-2 (LS-2)','4'='basal epithelial (Bs)',
                               '5'='endothelial (Endo)','6'='T-cell (T)','7'='Pericyte (Pc)','8'='Myeloid (Myl)','9'='B-cell (B)')
atac_subset_qc$celltype_RNA <- Idents(atac_subset_qc)

saveRDS(atac_subset_qc, 'results/singlecell/atac_RNA_annotated.rds')

# restrict to specific cell type
atac_subset_qc<-readRDS('results/singlecell/atac_RNA_annotated.rds')
Final_cell_barcode<-fread('results/singlecell/final_cell_barcodes.csv',header=TRUE)
atac_subset_qc <- subset(atac_subset_qc , cells = Final_cell_barcode$x)

cell_counts <- table(Idents(atac_subset_qc))
# luminal secretory-1 (LS-1)            fibroblast (Fb)           luminal HR (LHr) luminal secretory-2 (LS-2)      basal epithelial (Bs)
# endothelial (Endo)                 T-cell (T)              Pericyte (Pc)              Myeloid (Myl)                 B-cell (B)
atac_subset_qc.rna_markers <- FindAllMarkers(
  object = atac_subset_qc,
  # min.diff.pct = 0.5,  # Ensuring genes are significantly different between groups
)

atac_subset_qc.rna_markers <- atac_subset_qc.rna_markers %>%
  filter(p_val_adj < 0.05 & avg_log2FC > 1 )

write.csv(atac_subset_qc.rna_markers,'results/singlecell/Markers_RNA_DEGs_perpaper_labelled.csv',row.names = FALSE)

# restrict to LS, LHr, Bs
atac_BCcells <- subset(atac_subset_qc, idents = c("luminal secretory-1 (LS-1)",
                                                  "luminal secretory-2 (LS-2)",
                                                  "luminal HR (LHr)","basal epithelial (Bs)"))
atac_BCcells.rna_markers <- FindAllMarkers(
  object = atac_BCcells,
  # min.diff.pct = 0.5,  # Ensuring genes are significantly different between groups
)

atac_BCcells.rna_markers <- atac_BCcells.rna_markers %>%
  filter(p_val_adj < 0.05 & avg_log2FC > 1 )

write.csv(atac_BCcells.rna_markers,'results/singlecell/Markers_RNA_DEGs_LHr_LS_Bs_labelled.csv',row.names = FALSE)

jpeg('results/singlecell/RNA07_DimPlot_RNA_10clusters.jpg',width=6, height = 5, units = 'in', res = 300)
DimPlot(atac_subset_qc, reduction = "umap.rna",  label = TRUE, pt.size = 0.5)+ ggtitle("RNA")+ scale_color_manual(values = colors)+ NoLegend()
dev.off()

atac_subset_qc<-readRDS('results/singlecell/atac_RNA_annotated.rds')

jpeg('results/singlecell/RNA08_DimPlot_RNA_Doublet.jpg',width=12, height = 4, units = 'in', res = 300)
p8 <- DimPlot(atac_subset_qc, reduction = "umap.rna", group.by='celltype_RNA',label = FALSE, label.size = 3.5, repel = TRUE)  + ggtitle("RNA")+ scale_color_manual(values = colors)
p9 <- FeaturePlot(atac_subset_qc, features = "scDblFinder_score", reduction = "umap.rna") +
  ggtitle("Doublet Score")+
  theme(plot.title = element_text(hjust = 0.5), legend.position = "right",
        legend.title = element_text(size = 9, face = "bold"),
        legend.text = element_text(size = 8))  +
  labs(color = "DoubletScore")
p8+p9
dev.off()

jpeg('results/singlecell/RNA09_FeaturePlot_specific_genes.jpg',width=20, height = 20, units = 'in', res = 50)
a<-DimPlot(atac_subset_qc, reduction = "umap.rna",  label = TRUE, ,label.size = 5)+ ggtitle("RNA")+ scale_color_manual(values = colors)+ NoLegend()
b<-FeaturePlot(atac_subset_qc,features = c('FOXA1','TP63','EPCAM','KIT','ITGA6','ESR1','GATA3'))
design<-"##AAA###
         BBBBBBBB
         BBBBBBBB"
wrap_plots(A=a,B=b,ncol = 1,nrow = 3, design = design)
dev.off()

jpeg('results/singlecell/RNA10_FeaturePlot_figS5D_annotated.jpg',width=18, height = 24, units = 'in', res = 50)
a<-DimPlot(atac_subset_qc, reduction = "umap.rna",  label = TRUE, ,label.size = 5)+ ggtitle("RNA")+ scale_color_manual(values = colors)+ NoLegend()
b<-FeaturePlot(atac_subset_qc,features = c('ANKRD30A','ERBB4','PROM1','PTN','KRT14',
                                        'ACTA2','RGS5','PTPRC','CD3E','IL7R','CD79A',
                                        'IGJ','HLA-DRA','VWF','COL1A2','LUM'))
design<-"#AA#
         BBBB
         BBBB
         BBBB"
wrap_plots(A=a,B=b,ncol = 1,nrow = 4, design = design)
dev.off()

#########################################
#   check the mitochondrial proportion  #
#########################################
atac_subset_qc<-readRDS('results/singlecell/atac_RNA_annotated.rds')

# 13 protein-coding genes
mt_gene <- c("MT-ND1", "MT-ND2", "MT-ND3", "MT-ND4", "MT-ND4L", "MT-ND5", "MT-ND6",
             "MT-CO1", "MT-CO2", "MT-CO3", "MT-ATP6", "MT-ATP8", "MT-CYB")
mt_gene_det <- mt_gene[mt_gene %in% rownames(atac_subset_qc)]

atac_subset_qc[["percent.mt"]] <- PercentageFeatureSet(atac_subset_qc, features = mt_gene_det)
summary(atac_subset_qc$percent.mt)
# > summary(atac_subset_qc$percent.mt)
# Min.  1st Qu.   Median     Mean  3rd Qu.     Max.

atac_subset_qc[["percent.mt2"]] <- PercentageFeatureSet(atac_subset_qc, pattern = "^MT-")
summary(atac_subset_qc$percent.mt2)

library(ggplot2)
dat <- data.frame(byPattern = atac_subset_qc$percent.mt, byGene = atac_subset_qc$percent.mt2,
                  stringsAsFactors = FALSE)

cor_val <- cor.test(dat$byPattern, dat$byGene, method = "spearman")

ggplot(dat, aes(x = byPattern, y = byGene)) + geom_point() + geom_smooth() + labs(x = "% of MT, est by genes",
                                                                                  y = "% of MT,est by pattern", subtitle = paste0("rho=", round(cor_val$estimate,
                                                                                                                                                3), "; p-value=", cor_val$p.value[1])) + theme_classic()
# Potential cell debris would show low UMI counts (nCount_RNA) and high percentage of mitochondrial genes (percent.mt)
FeatureScatter(atac_subset_qc, feature1 = "nCount_RNA", feature2 = "percent.mt")+ scale_color_manual(values = colors)

FeatureScatter(atac_subset_qc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")+ scale_color_manual(values = colors)

# In most scRNA cases, we set the percent.mt > 10 (95% of overall population)
# total 16677

table(atac_subset_qc$percent.mt >= 20)

table(atac_subset_qc$celltype_RNA,atac_subset_qc$percent.mt >= 20)
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

FeaturePlot(atac_subset_qc, features = "percent.mt", reduction = "umap.rna") +
  ggtitle("percent.mt")+
  theme(plot.title = element_text(hjust = 0.5), legend.position = "right",
        legend.title = element_text(size = 9, face = "bold"),
        legend.text = element_text(size = 8))  +
  labs(color = "percent.mt")

########### violine plots ################
atac_subset_qc.rna_markers<-fread('results/singlecell/Markers_RNA_DEGs_perpaper_labelled.csv',header = TRUE)
.rowNamesDF(atac_subset_qc.rna_markers, make.names=TRUE)<-atac_subset_qc.rna_markers$gene
unique(atac_subset_qc.rna_markers$cluster)
