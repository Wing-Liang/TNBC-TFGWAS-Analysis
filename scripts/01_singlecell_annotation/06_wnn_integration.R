# ==============================================================================
# 06_wnn_integration.R
# Weighted-nearest-neighbour integration of RNA and ATAC.
# Produces: Fig 2
# Inputs  : results/singlecell/atac_ATAC_RNAlabeled.rds
#           results/singlecell/atac_RNA_annotated.rds
# Outputs : results/singlecell/WNN01_DimPlot_RNA_ATAC_WNN_b4annotated.jpg
#           results/singlecell/WNN_Markers_10clusters.csv
#           results/singlecell/atac_WNN_annotated.rds
#           results/singlecell/WNN02_DimPlot_RNA_ATAC_WNN_annotated.jpg
#           results/singlecell/WNN03_DimPlot_Doublets.jpg
#           results/singlecell/WNN04_DimPlot_RNA_ATAC_WNN_annotated.jpg
#           results/singlecell/umap_combined.csv
#           ...
# Run     : Rscript scripts/01_singlecell_annotation/06_wnn_integration.R   (from the repository root)
# ==============================================================================
colors <- c(
  "#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd",
  "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf",
  "#393b79", "#637939", "#8c6d31", "#843c39", "#7b4173",
  "#5254a3", "#6b6ecf", "#9c9ede", "#e7969c", "#e7ba52",
  "#bd9e39","black"
)
########################################################
###         Weighted Nearest Neighbor Analysis       ###
########################################################
# RNA analysis
options(future.globals.maxSize = 8000 * 1024^3) # avoid Error in getGlobalsAndPackages(expr, envir = envir, globals = globals)
atac_subset_qc <- readRDS('results/singlecell/atac_ATAC_RNAlabeled.rds')

DefaultAssay(atac_subset_qc) <- "RNA"
atac_subset_qc <- SCTransform(atac_subset_qc, verbose = FALSE) %>%
  RunPCA() %>%
  RunUMAP(dims = 1:20, reduction.name = 'umap.rna_wnn', reduction.key = 'rnaUMAP_')

# ATAC analysis
# We exclude the first dimension as this is typically correlated with sequencing depth
DefaultAssay(atac_subset_qc) <- "ATAC"
atac_subset_qc <- RunTFIDF(atac_subset_qc)
atac_subset_qc <- FindTopFeatures(atac_subset_qc, min.cutoff = 'q5')
atac_subset_qc <- RunSVD(atac_subset_qc)
atac_subset_qc <- RunUMAP(atac_subset_qc, reduction = 'lsi', dims = 2:30, reduction.name = "umap.atac_wnn",
                          reduction.key = "atacUMAP_")

# Find MultiModal
atac_subset_qc <- FindMultiModalNeighbors(atac_subset_qc, reduction.list = list("pca",'lsi'),
                                          dims.list = list(1:20, 2:30))
atac_subset_qc <- RunUMAP(atac_subset_qc, nn.name = "weighted.nn", reduction.name = "wnn.umap",
                          reduction.key = "wnnUMAP_", verbose = TRUE)
atac_subset_qc <- FindClusters(atac_subset_qc, graph.name = "wsnn", algorithm = 3,
                               resolution = 0.06)
# Number of communities: 10

jpeg('results/singlecell/WNN01_DimPlot_RNA_ATAC_WNN_b4annotated.jpg',width=10, height = 4, units = 'in', res = 300)
p1 <- DimPlot(atac_subset_qc, reduction = "umap.rna_wnn",label = FALSE, label.size = 3.5, repel = TRUE) + ggtitle("RNA")  + scale_color_manual(values = colors) +NoLegend()
p2 <- DimPlot(atac_subset_qc, reduction = "umap.atac_wnn", label = FALSE, label.size = 3.5, repel = TRUE) + ggtitle("ATAC")+ scale_color_manual(values = colors) +NoLegend()
p3 <- DimPlot(atac_subset_qc, reduction = "wnn.umap", label = FALSE, label.size = 3.5, repel = TRUE) + ggtitle("WNN")+ scale_color_manual(values = colors)
(p1 + p2 + p3) & theme(plot.title = element_text(hjust = 0.5))
dev.off()

# Find markers
DefaultAssay(atac_subset_qc) <- "RNA"
atac_subset_qc.wnn_markers <- FindAllMarkers(atac_subset_qc)

# add annotations
atac_subset_qc <- RenameIdents(atac_subset_qc, '0'='LS-1','1' = 'Fb',
                               '2'='LHr','3'='LS-2','4'='Bs',
                               '5'='Endo','6'='T','7'='Pc','8'='Myl','9'='B')

atac_subset_qc$celltype_WNN <- Idents(atac_subset_qc)

write.csv(atac_subset_qc.wnn_markers,'results/singlecell/WNN_Markers_10clusters.csv',row.names = FALSE)

saveRDS(atac_subset_qc, 'results/singlecell/atac_WNN_annotated.rds')

jpeg('results/singlecell/WNN02_DimPlot_RNA_ATAC_WNN_annotated.jpg',width=13, height = 4, units = 'in', res = 300)
p4 <- DimPlot(atac_subset_qc, reduction = "umap.rna_wnn",group.by='celltype_WNN', label = TRUE, label.size = 3.5, repel = TRUE) + ggtitle("RNA")+ scale_color_manual(values = colors) +NoLegend()
p5 <- DimPlot(atac_subset_qc, reduction = "umap.atac_wnn",group.by='celltype_WNN',label = TRUE, label.size = 3.5, repel = TRUE) + ggtitle("ATAC")+ scale_color_manual(values = colors) +NoLegend()
p6 <- DimPlot(atac_subset_qc, reduction = "wnn.umap", group.by='celltype_WNN',label =TRUE, label.size = 3.5, repel = TRUE) + ggtitle("WNN")+ scale_color_manual(values = colors)+NoLegend()
p4 + p5 + p6
dev.off()

jpeg('results/singlecell/WNN03_DimPlot_Doublets.jpg',width=10, height = 4, units = 'in', res = 300)
p6 <- DimPlot(atac_subset_qc, reduction = "wnn.umap", group.by='celltype_WNN',label = TRUE, label.size = 3.5, repel = TRUE)  + ggtitle("WNN")+ scale_color_manual(values = colors)+NoLegend()
p7 <- FeaturePlot(atac_subset_qc, features = "scDblFinder_score", reduction = "wnn.umap") +
  ggtitle("Doublet Score")+
  theme(plot.title = element_text(hjust = 0.5), legend.position = "right",
        legend.title = element_text(size = 9, face = "bold"),
        legend.text = element_text(size = 8))  +
  labs(color = "DoubletScore")
p6+p7
dev.off()

jpeg('results/singlecell/WNN04_DimPlot_RNA_ATAC_WNN_annotated.jpg',width=13, height = 4, units = 'in', res = 300)
p8 <- DimPlot(atac_subset_qc.rna, reduction = "umap.rna", group.by='celltype_RNA', label = TRUE, label.size = 3.5, repel = TRUE) + ggtitle("RNA")+ scale_color_manual(values = colors) +NoLegend()
p9 <- DimPlot(atac_subset_qc, reduction = "umap.atac",group.by='celltype_ATAC',label = TRUE, label.size = 3.5, repel = TRUE) + ggtitle("ATAC")+ scale_color_manual(values = colors) +NoLegend()
p10 <- DimPlot(atac_subset_qc, reduction = "wnn.umap", group.by='celltype_WNN',label =TRUE, label.size = 3.5, repel = TRUE) + ggtitle("WNN")+ scale_color_manual(values = colors)+NoLegend()
p8 + p9 + p10
dev.off()

##############################################
######         Combination UMAP plots    #####
##############################################

# rna
atac_subset_qc.rna <- readRDS('results/singlecell/atac_RNA_annotated.rds')
DefaultAssay(atac_subset_qc.rna)<-'RNA'
levels(atac_subset_qc.rna)
new.cluster.ids <- c('LS-1','Fb','LHr','LS-2','Bs',
                     'Endo','T','Pc','Myl','B')
names(new.cluster.ids) <- levels(atac_subset_qc.rna)
atac_subset_qc.rna <- RenameIdents(atac_subset_qc.rna, new.cluster.ids)

atac_subset_qc.rna$celltype_RNA <- Idents(atac_subset_qc.rna)

prna <- DimPlot(atac_subset_qc.rna, reduction = "umap.rna", group.by='celltype_RNA',label = TRUE, label.size = 3.5, repel = TRUE)
patac <- DimPlot(atac_subset_qc, reduction = "umap.atac",group.by='celltype_ATAC',label = TRUE, label.size = 3.5, repel = TRUE)
pwnn <- DimPlot(atac_subset_qc, reduction = "wnn.umap", group.by='celltype_WNN',label =TRUE, label.size = 3.5, repel = TRUE)
umap_rna<-prna$data
umap_atac<-patac$data
umap_wnn<-pwnn$data

# Ensure row names (cell barcodes) are preserved as a column for merging
umap_rna$Barcode <- rownames(umap_rna)
umap_atac$Barcode <- rownames(umap_atac)
umap_wnn$Barcode <- rownames(umap_wnn)

# Merge the three data frames by "Barcode"
umap_combined <- merge(umap_rna, umap_atac, by = "Barcode", suffixes = c("_RNA", "_ATAC"))%>%
  merge( umap_wnn, by = "Barcode", suffixes = c("", "_WNN"))

umap_combined$RNA_WNN<-ifelse(umap_combined$celltype_RNA==umap_combined$celltype_WNN,'TRUE','FALSE')
table(umap_combined$RNA_WNN)

umap_combined$RNA_ATAC<-ifelse(umap_combined$celltype_RNA==umap_combined$celltype_ATAC,'TRUE','FALSE')
table(umap_combined$RNA_ATAC)

# View the combined data frame
names(umap_combined)

umap_combined<-select(umap_combined,'Barcode','celltype_RNA','celltype_ATAC','celltype_WNN','RNA_WNN')
# Save to CSV (optional)
write.csv(umap_combined, file = "results/singlecell/umap_combined.csv", row.names = FALSE)

###########
jpeg('results/singlecell/Joint01_DimPlot_coloredRNA.jpg',width=12, height = 4, units = 'in', res = 300)
p20 <- DimPlot(atac_subset_qc, reduction = "umap.rna", group.by='celltype_RNA',label = FALSE, label.size = 3.5, repel = TRUE)  +
  ggtitle("RNA")+ scale_color_manual(values = colors)+NoLegend()
p21 <- DimPlot(atac_subset_qc, reduction = "umap.atac", group.by='celltype_RNA',label = FALSE, label.size = 3.5, repel = TRUE)  +
  ggtitle("ATAC colored by RNA")+ scale_color_manual(values = colors)+NoLegend()
p22 <- DimPlot(atac_subset_qc, reduction = "wnn.umap", group.by='celltype_RNA',label = FALSE, label.size = 3.5, repel = TRUE)  +
  ggtitle("WNN colored by RNA")+ scale_color_manual(values = colors)
p20+p21+p22
dev.off()

jpeg('results/singlecell/Joint02_DimPlot_coloredATAC.jpg',width=10, height = 4, units = 'in', res = 300)
p23 <- DimPlot(atac_subset_qc, reduction = "umap.rna", group.by='celltype_ATAC',label = FALSE, label.size = 3.5, repel = TRUE)  +
  ggtitle("RNA colored by ATAC")+ scale_color_manual(values = colors)+NoLegend()
p24 <- DimPlot(atac_subset_qc, reduction = "umap.atac", group.by='celltype_ATAC',label = FALSE, label.size = 3.5, repel = TRUE)  +
  ggtitle("ATAC")+ scale_color_manual(values = colors)+NoLegend()
p25 <- DimPlot(atac_subset_qc, reduction = "wnn.umap", group.by='celltype_ATAC',label = FALSE, label.size = 3.5, repel = TRUE)  +
  ggtitle("WNN colored by ATAC")+ scale_color_manual(values = colors)
p23+p24+p25
dev.off()

jpeg('results/singlecell/Joint03_DimPlot_coloredWNN.jpg',width=10, height = 4, units = 'in', res = 300)
p23 <- DimPlot(atac_subset_qc, reduction = "umap.rna", group.by='celltype_WNN',label = FALSE, label.size = 3.5, repel = TRUE)  +
  ggtitle("RNA colored by WNN")+ scale_color_manual(values = colors)+NoLegend()
p24 <- DimPlot(atac_subset_qc, reduction = "umap.atac", group.by='celltype_WNN',label = FALSE, label.size = 3.5, repel = TRUE)  +
  ggtitle("ATAC colored by WNN")+ scale_color_manual(values = colors)+NoLegend()
p25 <- DimPlot(atac_subset_qc, reduction = "wnn.umap", group.by='celltype_WNN',label = FALSE, label.size = 3.5, repel = TRUE)  +
  ggtitle("WNN")+ scale_color_manual(values = colors)
p23+p24+p25
dev.off()
