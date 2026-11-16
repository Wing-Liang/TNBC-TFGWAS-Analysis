# ==============================================================================
# 08_markers_and_da_peaks.R
# Differentially accessible peaks, peak annotation and top-5 marker genes per cell type.
# Produces: Fig 2F, Table S1; DA peaks used for Fig 3F
# Inputs  : results/singlecell/atac_ATAC_RNAlabeled.rds
#           results/singlecell/final_cell_barcodes.csv
#           results/singlecell/da_results_all_FindAllMarkers.csv
#           results/singlecell/combined_peaks_gene_all.csv
#           results/singlecell/combined_peaks_gene_proteincoding.csv
#           results/singlecell/peak_annotation.csv
#           results/singlecell/combined_peaks_gene_top5.csv
#           ...
# Outputs : results/singlecell/da_results_all_FindAllMarkers.csv
#           results/singlecell/gene_list_proteincoding_lincRNA.csv
#           results/singlecell/dap_top5/DAP_CovPlot_Bs_LLNLF-65H9.1.jpg
#           results/singlecell/Upset_plot_DApeaks.jpg
#           results/singlecell/Stackedplot_plot_DApeaks.jpg
#           results/singlecell/Table_DApeak.jpg
#           results/singlecell/ATAC_volcano_plot_LS1_vs_LS2.jpg
#           ...
# Run     : Rscript scripts/01_singlecell_annotation/08_markers_and_da_peaks.R   (from the repository root)
# ==============================================================================
## ---- paths (set the environment variables; see README) ----
TNBC_ROOT <- Sys.getenv("TNBC_ROOT", ".")

library(Signac)
library(Seurat)

library(GenomicRanges)
library(ggplot2)
library(patchwork)
library(BSgenome.Hsapiens.UCSC.hg19)
library(EnsDb.Hsapiens.v75)
library(Herper)
library(devtools)
library(reticulate)

library(BiocParallel)

library(data.table)

library(dplyr)
library(tidyr)

library(TFBSTools)

library(patchwork)

library(JASPAR2022)

library(presto)

library(chromVAR)

library(ggseqlogo)
library(data.table)

# Peak annotation
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(org.Hs.eg.db)
library(ggupset)

# server
setwd(TNBC_ROOT)
#######################################################
#                          DAPs                       #
#######################################################
atac_subset_qc<-readRDS('results/singlecell/atac_ATAC_RNAlabeled.rds')
Final_cell_barcode<-fread('results/singlecell/final_cell_barcodes.csv',header=TRUE)
#ATAC
atac_subset_qc <- subset(atac_subset_qc , cells = Final_cell_barcode$x)
DefaultAssay(atac_subset_qc)<-'ATAC'

da_peaks <- FindAllMarkers(
    object = atac_subset_qc,
    test.use = 'LR',
    latent.vars = 'nCount_ATAC'
  )
da_peaks$peaks <- rownames(da_peaks)
write.csv(da_peaks, file = "results/singlecell/da_results_all_FindAllMarkers.csv", row.names = FALSE)

# ------------------- Closest Feature -----------------
# Find the closest gene to each peak (if within 100,000 bp)
da_peaks_FindAllMarkers <- fread("results/singlecell/da_results_all_FindAllMarkers.csv",
                                 header = TRUE)
da_peaks_FindAllMarkers <- da_peaks_FindAllMarkers[da_peaks_FindAllMarkers$p_val_adj < 0.05 &
                                                     abs(da_peaks_FindAllMarkers$avg_log2FC) > 1.5 ]
da_peaks_FindAllMarkers$peaks<-NULL

closest_genes_all <- ClosestFeature(
  object = atac_subset_qc,
  regions = da_peaks_FindAllMarkers$gene)

subset_closest_gene_all <- closest_genes_all[closest_genes_all$distance<100000,]

subset_closest_gene_all <- subset_closest_gene_all %>%
  rename(gene = query_region)

combined_peaks_gene_all <- merge(subset_closest_gene_all,
                                 da_peaks_FindAllMarkers,
                                 by = "gene", all = TRUE)

# Filter out rows where 'gene_id' is NA
combined_peaks_gene_all <- combined_peaks_gene_all %>%
  filter(!is.na(gene_id)) %>%
  unique()

write.csv(combined_peaks_gene_all,
          file = "results/singlecell/combined_peaks_gene_all.csv",
          row.names = FALSE)

######################################
# only keep the protein-coding genes #
######################################
combined_peaks_gene_all <- fread("results/singlecell/combined_peaks_gene_all.csv",
                                 header = TRUE)
combined_peaks_gene_all <- combined_peaks_gene_all[combined_peaks_gene_all$p_val_adj < 0.05 &
                                                     abs(combined_peaks_gene_all$avg_log2FC) > 1.5 ]
combined_peaks_gene_all$peaks<-NULL

table(combined_peaks_gene_all$gene_biotype)
# lincRNA processed_transcript       protein_coding                 rRNA
head(unique(combined_peaks_gene_all$gene_name)[combined_peaks_gene_all$gene_biotype=='lincRNA'])
head(unique(combined_peaks_gene_all$gene_name)[combined_peaks_gene_all$gene_biotype=='processed_transcript'])
head(unique(combined_peaks_gene_all$gene_name)[combined_peaks_gene_all$gene_biotype=='protein_coding'])
head(unique(combined_peaks_gene_all$gene_name)[combined_peaks_gene_all$gene_biotype=='rRNA'])

combined_peaks_gene_all <- combined_peaks_gene_all %>% dplyr::filter(gene_biotype %in% c('protein_coding','lincRNA'))
gene_protcoding_lincRNA<-select(combined_peaks_gene_all,'gene_name','gene_biotype')
gene_protcoding_lincRNA <- unique(gene_protcoding_lincRNA)
write.csv(gene_protcoding_lincRNA, "results/singlecell/gene_list_proteincoding_lincRNA.csv", row.names = FALSE)   # used by 05_atac_umap_figures.R

combined_peaks_gene_all <- combined_peaks_gene_all %>% dplyr::rename(peaks=gene,
                                                                     gene=gene_name,
                                                                     cluster_name=cluster)
write.csv(combined_peaks_gene_all,
          file = "results/singlecell/combined_peaks_gene_proteincoding.csv",
          row.names = FALSE)
######################################
#             peak annotation       #
######################################
combined_peaks_gene_all <- fread("results/singlecell/combined_peaks_gene_proteincoding.csv",
                                 header = TRUE)
# Profile of ChIP peaks binding to TSS regions
peak <- granges(atac_subset_qc)
edb <- EnsDb.Hsapiens.v75
seqlevelsStyle(edb) <- "UCSC"

txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
peakAnno <- annotatePeak(peak , tssRegion=c(-3000, 3000),
                             TxDb=txdb, annoDb="org.Hs.eg.db")
plotAnnoPie(peakAnno)

# Convert annotated peak object to a data frame
peak_anno_df <- as.data.frame(peakAnno)

# View the first few rows
head(peak_anno_df)

# Create a new column `annotation_category` to classify genomic regions
peak_anno_df <- peak_anno_df %>%
  mutate(annotation_category = case_when(
    grepl("Promoter", annotation) ~ "Promoter",
    grepl("5' UTR", annotation) ~ "5' UTR",
    grepl("3' UTR", annotation) ~ "3' UTR",

    # First Exon and First Intron
    grepl("exon", annotation, ignore.case = TRUE) ~ "Exon",
    grepl("intron", annotation, ignore.case = TRUE) ~ "Intron",

    grepl("Downstream \\(<=300bp\\)", annotation) ~ "Downstream (<=300bp)",
    grepl("Distal Intergenic", annotation) ~ "Distal Intergenic",
    TRUE ~ "Other"
  ))

# View the first few rows to verify results
head(peak_anno_df[, c("annotation", "annotation_category")])
table(peak_anno_df$annotation_category)

peak_anno_df <-  peak_anno_df %>%
  mutate(peaks = paste0(seqnames, "-", start, "-", end))

names(peak_anno_df)[1:18]<-paste(names(peak_anno_df)[1:18],'anno',sep = '_')
combined_peaks_gene_all <- left_join(combined_peaks_gene_all, peak_anno_df,by='peaks')

all.equal(combined_peaks_gene_all$gene_id,combined_peaks_gene_all$ENSEMBL_anno)
table(is.na(combined_peaks_gene_all$gene_id))
table(is.na(combined_peaks_gene_all$ENSEMBL_anno))

combined_peaks_gene_all$diff_gene<-ifelse(combined_peaks_gene_all$gene==combined_peaks_gene_all$SYMBOL_anno,'TRUE','FALSE')
write.csv(combined_peaks_gene_all,
          file = "results/singlecell/peak_annotation.csv",
          row.names = FALSE)
# ------------------- Closest Feature - top5  -----------------
# Find the closest gene to each peak (if within 100,000 bp)
combined_peaks_gene_all<-fread("results/singlecell/peak_annotation.csv",header = TRUE)
top5_dap_all <-combined_peaks_gene_all %>%
  group_by(cluster_name) %>%
  arrange(p_val, -abs(pct.1 -pct.2))%>%
  slice_head(n = 5)  %>%
  ungroup()

# Filter out rows where 'gene_id' is NA
combined_peaks_gene_top5 <- top5_dap_all %>%
  filter(!is.na(gene_id)) %>%
  unique()

write.csv(combined_peaks_gene_top5,
          file = "results/singlecell/combined_peaks_gene_top5.csv",
          row.names = FALSE)

##########################################################################################
# --------------------------------- Genomic region plot ---------------------------------
##########################################################################################
Final_cell_barcode<-fread('results/singlecell/final_cell_barcodes.csv',header=TRUE)

combined_peaks_gene_top5<-fread("results/singlecell/combined_peaks_gene_top5.csv", header = TRUE)
atac_subset_qc <- readRDS("results/singlecell/atac_links_general_all.rds")
links<-Links(atac_subset_qc)
atac_subset_qc <- subset(atac_subset_qc , cells = Final_cell_barcode$x)
links<-links[mcols(links)$score >= 0.1]
Links(atac_subset_qc)<-links
DefaultAssay(atac_subset_qc)<-'ATAC'

# Plotting genomic regions
atac_subset_qc <- SortIdents(atac_subset_qc)

# Load required libraries
library(Seurat)
library(Signac)
library(GenomicRanges)

# Define constants ------------------------------------------------------------
OUTPUT_DIR <- "results/singlecell/dap_top5/"

# Extract genes of interest per cluster
GENES_OF_INTEREST_BY_CLUSTER <- split(combined_peaks_gene_top5$gene, combined_peaks_gene_top5$cluster_name)

COMPARISONS <- list(
  list(cluster = "LHr", suffix = "LHr"),
  list(cluster = "LS-1", suffix = "LS1"),
  list(cluster = "LS-2", suffix = "LS2"),
  list(cluster = "Bs", suffix = "Bs"),
  list(cluster = "Fb", suffix = "Fb"),
  list(cluster = "Pc", suffix = "Pc"),
  list(cluster = "Endo", suffix = "Endo"),
  list(cluster = "Myl", suffix = "Myl"),
  list(cluster = "T", suffix = "T"),
  list(cluster = "B", suffix = "B")
)

# Define plotting parameters ---------------------------------------------------
PLOT_WIDTH <- 8  # inches
PLOT_HEIGHT <- 8 # inches
PLOT_RES <- 200   # dpi
EXTEND_UPSTREAM <- 65000
EXTEND_DOWNSTREAM <- 65000

# Core plotting function -------------------------------------------------------
generate_comparison_plots <- function(comparison, combined_data, atac_obj) {
  # Extract genes specific to the cluster

  genes_of_interest <- GENES_OF_INTEREST_BY_CLUSTER[[comparison$cluster]]

  # Validate input data
  if(is.null(genes_of_interest) || length(genes_of_interest) == 0) {
    warning(paste("No genes found for comparison:", comparison$suffix))
    return(NULL)
  }

  # Create output directory
  if(!dir.exists(OUTPUT_DIR)) {
    dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)
  }

  # Generate plots for each gene in the cluster
  plots <- lapply(genes_of_interest, function(gene) {
    tryCatch({
      # Find overlapping regions
      gene_coords <- LookupGeneCoords(atac_obj, gene)
      gene_coords <- resize(gene_coords, width = width(gene_coords) + 2 * 65000, fix = "center")
      peak_ranges <- StringToGRanges(combined_data$peaks)
      highlight_regions <- subsetByOverlaps(peak_ranges, gene_coords)
      highlight_regions$color <- "orange"

      # Generate coverage plot
      CoveragePlot(
        object = atac_obj,
        region = gene,
        extend.upstream = EXTEND_UPSTREAM,
        extend.downstream = EXTEND_DOWNSTREAM,
        annotation = TRUE,
        peaks = TRUE,
        features = gene,
        scale.factor = 1e6,
        region.highlight = highlight_regions,
        expression.assay = "RNA",
        links = TRUE
      )
    }, error = function(e) {
      message(paste("Error processing", gene, ":", e$message))
      NULL
    })
  })

  # Save valid plots
  save_plots(plots, comparison$suffix, genes_of_interest)
}

# Helper function for saving plots ---------------------------------------------
save_plots <- function(plots, comparison_suffix, genes) {
  invisible(lapply(seq_along(plots), function(i) {
    if(!is.null(plots[[i]])) {
      output_file <- file.path(
        OUTPUT_DIR,
        paste0("DAP_CovPlot_links_", comparison_suffix, "_", genes[i], ".jpg")
      )

      jpeg(
        filename = output_file,
        width = PLOT_WIDTH,
        height = PLOT_HEIGHT,
        units = "in",
        res = PLOT_RES
      )
      print(plots[[i]])
      dev.off()
    }
  }))
}

# Execute analysis pipeline ----------------------------------------------------
# Initialize output directory
if(!dir.exists(OUTPUT_DIR)) {
  dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)
}

# Process all comparisons
invisible(lapply(COMPARISONS, function(comp) {
  generate_comparison_plots(
    comparison = comp,
    combined_data = combined_peaks_gene_top5,
    atac_obj = atac_subset_qc
  )
}))

cat("All coverage plots successfully generated for each cluster-specific gene set.\n")

################################################################
 #                          Outstanding genes                  #
################################################################
combined_peaks_gene_top5<-fread("results/singlecell/combined_peaks_gene_top5.csv",header = TRUE)

# find DA peaks overlapping gene of interest
# Bs
gene_coords <- LookupGeneCoords(atac_subset_qc, 'LLNLF-65H9.1')
peak_ranges <- StringToGRanges(combined_peaks_gene_top5[combined_peaks_gene_top5$cluster=='Bs' &
                                                          combined_peaks_gene_top5$gene_name=='LLNLF-65H9.1',]$gene)
gene_coords <- resize(gene_coords, width = width(gene_coords) + 2 * 1.5e5, fix = "center")
highlight_regions <- subsetByOverlaps(peak_ranges, gene_coords)
highlight_regions$color <- "orange"

jpeg("results/singlecell/dap_top5/DAP_CovPlot_Bs_LLNLF-65H9.1.jpg")
CoveragePlot(
  object = atac_subset_qc,
  region = 'LLNLF-65H9.1',
  region.highlight = highlight_regions,
  annotation = TRUE,
  peaks = TRUE,
  links = TRUE)
dev.off()

ExpressionPlot(
  object = atac_subset_qc,
  features = 'LLNLF-65H9.1',
  assay = "RNA"
)

# Check if the gene exists in the RNA assay
"LLNLF-65H9.1" %in% rownames(GetAssayData(atac_subset_qc, assay = "RNA"))

#############################################
#             Upset plot for DAPs          #
#############################################
# Load required libraries
library(UpSetR)
library(dplyr)
library(tidyr)
library(ggplot2)

# Load dataset
combined_peaks_gene_all <- fread("results/singlecell/peak_annotation.csv", header = TRUE)
combined_peaks_gene_all$cluster_name <- gsub("LS", "LASP", combined_peaks_gene_all$cluster_name)
combined_peaks_gene_all$cluster_name <- gsub("LHr", "LHS", combined_peaks_gene_all$cluster_name)
combined_peaks_gene_all$cluster_name <- gsub("Bs", "B-Myo", combined_peaks_gene_all$cluster_name)

table(combined_peaks_gene_all$cluster_name)
fig5D_levels <- c('LASP-1','LASP-2','LHS','B-Myo','Fb','Endo','T','B','Myl','Pc')

combined_peaks_gene_all <- combined_peaks_gene_all %>% unique()%>%
  mutate(cluster_name = factor(cluster_name, levels = fig5D_levels)) %>%
  arrange(cluster_name)

# Convert data into a binary matrix (presence/absence of peaks in each cluster)
binary_matrix <- combined_peaks_gene_all %>%
  dplyr::select(peaks, cluster_name) %>%
  mutate(value = 1) %>%
  pivot_wider(names_from = cluster_name, values_from = value, values_fill = 0)

# Remove peaks column for upset plot input
upset_matrix <- as.data.frame(binary_matrix)
rownames(upset_matrix) <- upset_matrix$peaks
upset_matrix <- upset_matrix[ , -1]

# Define colors
colors <-c( "#E64B35CC","#4DBBD5CC", "#00A087CC", "#B09C85CC","#F39B7FCC", "#8491B4CC", "#91D1C2CC",'#f9e79f', "#7E6148CC",  "#3C5488CC")

# Create Upset Plot
jpeg("results/singlecell/Upset_plot_DApeaks.jpg", width = 11.5, height = 5, units = "in", res = 200)
upset(
  upset_matrix,
  sets = rev(colnames(upset_matrix)),
  order.by = "freq",
  mainbar.y.label = "No. of Peaks",
  sets.bar.color = rev(colors[1:length(colnames(upset_matrix))]),
  matrix.color = "orange2",  # Color the matrix
  main.bar.color = "steelblue", # Color the intersection bars
  keep.order = TRUE
)
dev.off()

cat("Upset Diagram saved as: results/singlecell/upset_plot_DApeaks.jpg\n")

# peak annotation
library(ggplot2)
library(dplyr)
library(data.table)
library(scales)
library(gridExtra)

head(combined_peaks_gene_all)

# Create a summary count by cluster_name and annotation_category_anno
plot_data <- combined_peaks_gene_all %>%
  dplyr::count(cluster_name, annotation_category_anno) %>%
  group_by(cluster_name) %>%
  mutate(Percentage = (n / sum(n)) * 100,  # Convert to percentage
         Label = paste0(n, " (", round(Percentage, 1), "%)"))  # Create labels

# Plot stacked bar chart
# Plot stacked bar chart
jpeg("results/singlecell/Stackedplot_plot_DApeaks.jpg",
     width = 7, height = 5, units = "in", res = 300)

# Stacked histogram plot
stacked_plot <- ggplot(plot_data, aes(x = cluster_name, y = Percentage, fill = annotation_category_anno)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c(
    "3' UTR" = "#2066a8",
    "5' UTR" = "#8ec1da",
    "Distal Intergenic" = "#cde1ec",
    "Downstream (<=300bp)" = "#ededed",
    "Exon" = "#f6d6c2",
    "Intron" = "#f4a582",  # Added a reasonable color (was missing)
    "Promoter" = "#ae282c"
  )) +
  theme_minimal() +
  labs(
    x = "Cell type",
    y = "Percentage",
    fill = "Annotation Category"
  ) +
  scale_y_continuous(labels = scales::percent_format(scale = 1)) +
  theme(
    # Move legend to bottom
    legend.position = "bottom",
    # All text in black
    text = element_text(color = "black"),
    axis.title = element_text(color = "black"),
    axis.text = element_text(color = "black"),
    legend.title = element_text(color = "black"),
    legend.text = element_text(color = "black"),
    plot.title = element_text(color = "black")
  )

# Print plot
print(stacked_plot)
dev.off()

library(grid)  # For grid.draw function
# Create a summary table for display
table_data <- plot_data %>%
  select(cluster_name, annotation_category_anno, Label) %>%
  pivot_wider(names_from = cluster_name, values_from = Label, values_fill = list(Label = "0 (0%)")) %>%
  rename(Annotation = annotation_category_anno)

# Save the table as an independent plot
jpeg("results/singlecell/Table_DApeak.jpg", width = 12, height = 4, units = "in", res = 300)
print(table_data)
dev.off()

#############################################
#                     dot plot              #
#############################################
atac_subset_qc<-readRDS('results/singlecell/atac_ATAC_RNAlabeled.rds')
DefaultAssay(atac_subset_qc)<-'ATAC'

atac_subset_qc.rna_markers<-fread('results/singlecell/Markers_RNA_DEGs_perpaper.csv',header=TRUE)
.rowNamesDF(atac_subset_qc.rna_markers, make.names=TRUE)<-atac_subset_qc.rna_markers$gene
atac_subset_qc.rna <- readRDS('results/singlecell/atac_RNA_annotated.rds')
DefaultAssay(atac_subset_qc.rna)<-'RNA'

atac_subset_qc.rna[['RNA']]<-as(atac_subset_qc.rna[['RNA']], Class="Assay5")

levels(atac_subset_qc.rna)
new.cluster.ids <- c('LS-1','Fb','LHr','LS-2','Bs',
                     'Endo','T','Pc','Myl','B')
names(new.cluster.ids) <- levels(atac_subset_qc.rna)
atac_subset_qc.rna <- RenameIdents(atac_subset_qc.rna, new.cluster.ids)

atac_subset_qc.rna$celltype_RNA <- Idents(atac_subset_qc.rna)

combined_peaks_gene_top5<-fread("results/singlecell/combined_peaks_gene_top5.csv",header = TRUE)

# cell-type display order
fig5D_levels <- c('LHr','LS-1','LS-2','Bs','T','B','Myl','Endo','Pc','Fb')
# Reorder the Idents based on sorted levels
Idents(atac_subset_qc) <- factor(Idents(atac_subset_qc), levels = fig5D_levels)
Idents(atac_subset_qc.rna) <- factor(Idents(atac_subset_qc.rna), levels = fig5D_levels)
atac_subset_qc$celltype_ATAC <- Idents(atac_subset_qc)
atac_subset_qc.rna$celltype_RNA <- Idents(atac_subset_qc.rna)

combined_peaks_gene_top5 <- combined_peaks_gene_top5 %>%
  arrange(factor(cluster, levels = rev(fig5D_levels)))

# DotPlot for RNA Expression

# Define output directory
OUTPUT_DIR <- "results/singlecell/DAP_top5_all"

# Ensure output directory exists
if (!dir.exists(OUTPUT_DIR)) {
  dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)
}

# ------------------------------
# DotPlot for RNA Expression
# ------------------------------
DefaultAssay(atac_subset_qc.rna) <- "RNA"

p1_data <- DotPlot(atac_subset_qc.rna,
                   features = combined_peaks_gene_top5$gene_name,
                   group.by = "celltype_RNA")$data

p1_data$category <- "RNA Expression"

# ------------------------------
# DotPlot for ATAC
# ------------------------------
DefaultAssay(atac_subset_qc) <- "ATAC"

p2_data <- DotPlot(atac_subset_qc,
                   features = combined_peaks_gene_top5$gene,
                   group.by = "celltype_ATAC")$data

# Remove NA values (genes not found in dataset)
p2_data <- p2_data %>% filter(!is.na(features.plot))
p2_data$category <- "ATAC"

# ------------------------------
# Plot 1: RNA Expression DotPlot
# ------------------------------
jpeg(file.path(OUTPUT_DIR, "Dotplot01_top5_RNA.jpg"), width = 4, height = 8, units = 'in', res = 200)

p1 <- ggplot(p1_data, aes(x = id, y = features.plot)) +
  geom_point(aes(size = pct.exp, color = avg.exp.scaled)) +
  scale_color_gradient2(low = "plum3", mid = "grey90", high = "red", midpoint = 0) +
  scale_size(range = c(0.3, 6)) +
  labs(title = "RNA Expression", x = "Cell Type", y = "Genes") +
  theme_minimal() +
  theme(
    axis.ticks = element_line(color = "black"),
    axis.text = element_text(color = "black"),
    axis.title = element_text(color = "black"),
    strip.text = element_text(size = 10, face = "bold", color = "black"),
    panel.grid = element_blank(),
    legend.position = "bottom"
  )

print(p1)
dev.off()

# ------------------------------
# Plot 2: ATAC Gene Activity DotPlot
# ------------------------------
jpeg(file.path(OUTPUT_DIR, "Dotplot02_top5_ATAC.jpg"), width = 6, height = 8, units = 'in', res = 200)

p2 <- ggplot(p2_data, aes(x = id, y = features.plot)) +
  geom_point(aes(size = pct.exp, color = avg.exp.scaled)) +
  scale_color_gradient2(low = "dodgerblue3", mid = "grey90", high = "firebrick3", midpoint = 0) +
  scale_size(range = c(0.3, 6)) +
  labs(title = "ATAC", x = "Cell Type", y = "Peaks") +
  theme_minimal() +
  theme(
    axis.ticks = element_line(color = "black"),
    axis.text = element_text(color = "black"),
    axis.title = element_text(color = "black"),
    strip.text = element_text(size = 10, face = "bold", color = "black"),
    panel.grid = element_blank(),
    legend.position = "bottom"
  )

print(p2)
dev.off()

cat("DotPlots successfully generated: RNA Expression and ATAC.\n")

#############################################
#                     dot plot              #
#############################################
atac_subset_qc<-readRDS('results/singlecell/atac_ATAC_RNAlabeled.rds')
DefaultAssay(atac_subset_qc)<-'ATAC'

atac_subset_qc.rna_markers<-fread('results/singlecell/Markers_RNA_DEGs_perpaper.csv',header=TRUE)
.rowNamesDF(atac_subset_qc.rna_markers, make.names=TRUE)<-atac_subset_qc.rna_markers$gene
atac_subset_qc.rna <- readRDS('results/singlecell/atac_RNA_annotated.rds')
DefaultAssay(atac_subset_qc.rna)<-'RNA'

atac_subset_qc.rna[['RNA']]<-as(atac_subset_qc.rna[['RNA']], Class="Assay5")

levels(atac_subset_qc.rna)
new.cluster.ids <- c('LS-1','Fb','LHr','LS-2','Bs',
                     'Endo','T','Pc','Myl','B')
names(new.cluster.ids) <- levels(atac_subset_qc.rna)
atac_subset_qc.rna <- RenameIdents(atac_subset_qc.rna, new.cluster.ids)

atac_subset_qc.rna$celltype_RNA <- Idents(atac_subset_qc.rna)

atac_subset_qc.rna_markers$cluster <- as.character(atac_subset_qc.rna_markers$cluster)
atac_subset_qc.rna_markers %>%
  group_by(cluster) %>%
  slice_head(n = 5)  %>%
  ungroup() -> top5
cluster<-c('LS-1','Fb','LHr','LS-2','Bs','Endo','T','Pc','Myl','B')
top5$cluster_name <- rep(cluster, each = 5)

# cell-type display order
fig5D_levels <- c('LHr','LS-1','LS-2','Bs','T','B','Myl','Endo','Pc','Fb')
# Reorder the Idents based on sorted levels
top5 <- top5 %>% arrange(factor(cluster_name, levels = rev(fig5D_levels)))

Idents(atac_subset_qc) <- factor(Idents(atac_subset_qc), levels = fig5D_levels)
Idents(atac_subset_qc.rna) <- factor(Idents(atac_subset_qc.rna), levels = fig5D_levels)
atac_subset_qc$celltype_ATAC <- Idents(atac_subset_qc)
atac_subset_qc.rna$celltype_RNA <- Idents(atac_subset_qc.rna)

# DotPlot for RNA Expression
# Define output directory
OUTPUT_DIR <- "results/singlecell/DAP_top5_all"

# Ensure output directory exists
if (!dir.exists(OUTPUT_DIR)) {
  dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)
}

# ------------------------------
# DotPlot for RNA Expression
# ------------------------------
DefaultAssay(atac_subset_qc.rna) <- "RNA"

p1_data <- DotPlot(atac_subset_qc.rna,
                   features = top5$gene,
                   group.by = "celltype_RNA")$data

p1_data$category <- "RNA Expression"

# ------------------------------
# DotPlot for ATAC
# ------------------------------
DefaultAssay(atac_subset_qc) <- "ATAC"
combined_peaks_gene_all<-fread("results/singlecell/combined_peaks_gene_all.csv",header = TRUE)
table(combined_peaks_gene_all$gene_biotype)
# lincRNA processed_transcript       protein_coding                 rRNA
combined_peaks_gene_all <- combined_peaks_gene_all %>% dplyr::filter(gene_biotype=='protein_coding')
combined_peaks_gene_all <- combined_peaks_gene_all %>% dplyr::rename(peaks=gene,
                                                              gene=gene_name,
                                                              cluster_name=cluster)
unique_top_peaks <- combined_peaks_gene_all %>%
  group_by(cluster_name, gene) %>%
  slice_head(n = 1)  %>%
  ungroup()
top5<-left_join(top5, select(unique_top_peaks ,'gene','cluster_name','peaks','distance'),
                by=c('gene','cluster_name')) %>% unique()

p2_data <- DotPlot(atac_subset_qc,
                   features = top5$peaks,
                   group.by = "celltype_ATAC")$data

# Remove NA values (genes not found in dataset)
p2_data <- p2_data %>% filter(!is.na(features.plot))
p2_data$category <- "ATAC"

# ------------------------------
# Plot 1: RNA Expression DotPlot
# ------------------------------
jpeg(file.path(OUTPUT_DIR, "Dotplot01_top5DEGs_RNA.jpg"), width = 6, height = 8, units = 'in', res = 200)

p1 <- ggplot(p1_data, aes(x = id, y = features.plot)) +
  geom_point(aes(size = pct.exp, color = avg.exp.scaled)) +
  scale_color_gradient2(low = "plum3", mid = "grey90", high = "red", midpoint = 0) +
  scale_size(range = c(0.3, 6)) +
  labs(title = "RNA Expression", x = "Cell Type", y = "Genes") +
  theme_minimal() +
  theme(
    axis.ticks = element_line(color = "black"),
    axis.text = element_text(color = "black"),
    axis.title = element_text(color = "black"),
    strip.text = element_text(size = 6, face = "bold", color = "black"),
    panel.grid = element_blank(),
    legend.position = "bottom"
  )

print(p1)
dev.off()

# ------------------------------
# Plot 2: ATAC Gene Activity DotPlot
# ------------------------------
jpeg(file.path(OUTPUT_DIR, "Dotplot02_top5DEGs_ATAC.jpg"), width = 6, height = 8, units = 'in', res = 200)

p2 <- ggplot(p2_data, aes(x = id, y = features.plot)) +
  geom_point(aes(size = pct.exp, color = avg.exp.scaled)) +
  scale_color_gradient2(low = "dodgerblue3", mid = "grey90", high = "firebrick3", midpoint = 0) +
  scale_size(range = c(0.3, 6)) +
  labs(title = "ATAC", x = "Cell Type", y = "Peaks") +
  theme_minimal() +
  theme(
    axis.ticks = element_line(color = "black"),
    axis.text = element_text(color = "black"),
    axis.title = element_text(color = "black"),
    strip.text = element_text(size = 6, face = "bold", color = "black"),
    panel.grid = element_blank(),
    legend.position = "bottom"
  )

print(p2)
dev.off()

cat("DotPlots successfully generated: RNA Expression and ATAC.\n")

################################################################
#                      volcano plot different lv              #
################################################################

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

ggsave("results/singlecell/ATAC_volcano_plot_LS1_vs_LS2.jpg", width = 8, height = 6, dpi = 300)

################################################################
#         heatmap of RNA exp and correlated gene in ATAC       #
################################################################
Final_cell_barcode<-fread('results/singlecell/final_cell_barcodes.csv',header=TRUE)

#ATAC
atac_subset_qc<-readRDS('results/singlecell/atac_ATAC_RNAlabeled.rds')
atac_subset_qc <- subset(atac_subset_qc , cells = Final_cell_barcode$x)
DefaultAssay(atac_subset_qc)<-'ATAC'

#RNA
atac_subset_qc.rna <- readRDS('results/singlecell/atac_RNA_annotated.rds')
DefaultAssay(atac_subset_qc.rna)<-'RNA'
atac_subset_qc.rna[['RNA']]<-as(atac_subset_qc.rna[['RNA']], Class="Assay5")
atac_subset_qc.rna <- subset(atac_subset_qc.rna , cells = Final_cell_barcode$x)

levels(atac_subset_qc.rna)
new.cluster.ids <- c('LS-1','Fb','LHr','LS-2','Bs',
                     'Endo','T','Pc','Myl','B')
names(new.cluster.ids) <- levels(atac_subset_qc.rna)
atac_subset_qc.rna <- RenameIdents(atac_subset_qc.rna, new.cluster.ids)
atac_subset_qc.rna$celltype_RNA <- Idents(atac_subset_qc.rna)

# re-order the levels sequence
# cell-type display order
fig5D_levels <- c('LHr','LS-1','LS-2','Bs','T','B','Myl','Endo','Pc','Fb')
# Reorder the Idents based on sorted levels
Idents(atac_subset_qc) <- factor(Idents(atac_subset_qc), levels = fig5D_levels)
Idents(atac_subset_qc.rna) <- factor(Idents(atac_subset_qc.rna), levels = fig5D_levels)
atac_subset_qc$celltype_ATAC <- Idents(atac_subset_qc)
atac_subset_qc.rna$celltype_RNA <- Idents(atac_subset_qc.rna)

# peak gene correlation: Spearman correlation between each top-5 DA peak (ATAC counts) and
# every gene (normalised RNA), across the post-QC cells. Slow; computed once and cached.
F_PEAKCOR <- "results/singlecell/Peak_RNA_correlation_matrix.csv"
if (!file.exists(F_PEAKCOR)) {
  top5 <- fread("results/singlecell/combined_peaks_gene_top5.csv", header = TRUE)
  peak_accessibility <- as.matrix(GetAssayData(atac_subset_qc, assay = "ATAC", layer = "counts")[top5$peaks, ])
  rna_expression <- as.matrix(GetAssayData(atac_subset_qc, assay = "RNA", layer = "data"))
  correlation_matrix <- matrix(0, nrow = nrow(peak_accessibility), ncol = nrow(rna_expression),
                               dimnames = list(rownames(peak_accessibility), rownames(rna_expression)))
  for (i in seq_len(nrow(peak_accessibility))) {
    for (j in seq_len(nrow(rna_expression))) {
      correlation_matrix[i, j] <- cor(peak_accessibility[i, ], rna_expression[j, ], method = "spearman")
    }
  }
  correlation_df <- as.data.frame(correlation_matrix)
  correlation_df$peak <- rownames(correlation_matrix)
  write.csv(correlation_df, F_PEAKCOR, row.names = FALSE)
}
peak_cor<-fread(F_PEAKCOR,header=TRUE)
# Ensure peak_cor is a data frame (if it's a matrix, convert it)
if (is.matrix(peak_cor)) {
  peak_cor <- as.data.frame(peak_cor)
}

rownames(peak_cor)<-peak_cor$peak

# Initialize a data frame to store results
peak_cor_top5DAP <- data.frame(
  Peak = rownames(peak_cor),  # Peaks (rows)
  Gene = character(nrow(peak_cor)),  # Gene with max correlation for each peak
  Correlation = numeric(nrow(peak_cor)),  # Max correlation value for each peak
  stringsAsFactors = FALSE
)

# Loop through each row (peak) to find the maximum correlation for each peak
for (i in 1:nrow(peak_cor)) {
  # Extract the correlation values for the current peak
  peak_cor_values <-as.vector(peak_cor[i, 1:32738])  # Exclude the "peak" column

  # Find the index of the maximum correlation value
  max_corr_index <- which.max(abs(peak_cor[i, 1:32738]))

  # Extract the gene name and correlation value
  peak_cor_top5DAP$Peak[i] <- peak_cor$peak[i]
  peak_cor_top5DAP$Gene[i] <- colnames(peak_cor)[max_corr_index]
  peak_cor_top5DAP$Correlation[i] <- peak_cor_values[max_corr_index][[1]]
}

# Save the peak_cor_top5DAP to a file
write.csv(peak_cor_top5DAP, file = "results/singlecell/max_correlation_per_peak.csv", row.names = FALSE)

peak_cor_top5DAP<-peak_cor_top5DAP%>% rename(peaks=Peak)

#top5 DAPs
combined_peaks_gene_top5<-fread("results/singlecell/combined_peaks_gene_top5.csv", header = TRUE)

combined_peaks_gene_top5<-left_join(select(combined_peaks_gene_top5,'peaks','gene','cluster_name'),
                                    peak_cor_top5DAP, by='peaks')
fig5D_levels <- c('LHr','LS-1','LS-2','Bs','T','B','Myl','Endo','Pc','Fb')
cluster_data <- combined_peaks_gene_top5 %>% unique()%>%
  mutate(cluster_name = factor(cluster_name, levels = fig5D_levels)) %>%
  arrange(cluster_name)

rna_data <- LayerData(atac_subset_qc.rna, assay = "RNA", layer = "scale.data")
atac_data <- LayerData(atac_subset_qc, assay = "ATAC", layer = "counts")

rna_data_range <- range(rna_data, na.rm = TRUE)
print(rna_data_range)

atac_data_range <- range(atac_data, na.rm = TRUE)
print(atac_data_range)

# Create a heatmap for gene expression (RNA assay)
# Load necessary libraries
library(Seurat)
colors <-c( "#E64B35CC","#4DBBD5CC", "#00A087CC", "#B09C85CC","#F39B7FCC", "#8491B4CC", "#91D1C2CC",'#f9e79f', "#7E6148CC",  "#3C5488CC")

# Identify significant genes (correlation > 0.5)
significant_genes <- unique(cluster_data$Gene[cluster_data$Correlation > 0])
significant_genes
# Create a heatmap for RNA expression (RNA assay)
rna_heatmap <- DoHeatmap(
  object = atac_subset_qc.rna,                  # Seurat object
  features = unique(cluster_data$Gene),         # Plot all genes
  assay = "RNA",                                # Use RNA assay
  group.by = "celltype_RNA",                    # Group cells by metadata (e.g., cell type)
  disp.min = -2,
  disp.max = 3,
  label = FALSE                                # Disable default gene labels
) + ggtitle("RNA expression") + scale_fill_gradientn(colors = c( "#002366", "#195905","yellow"))
# Customize the RNA heatmap theme
rna_heatmap <- rna_heatmap +
  theme(
    legend.position = "right"                   # Move legend to the right
  )
# Add custom y-axis labels for significant genes
rna_heatmap <- rna_heatmap +
  scale_y_discrete(
    breaks = significant_genes,                 # Only label significant genes
    labels = significant_genes
  )

# 6. Create a heatmap for peak accessibility (ATAC assay)
atac_heatmap <- DoHeatmap(
  object = atac_subset_qc,                  # Seurat object
  features = cluster_data$peaks,            # Significant peaks
  assay = "ATAC",                          # Use ATAC assay
  slot = 'counts',
  group.by = "celltype_ATAC",
  disp.min = 0,
  disp.max = 2,
  label = FALSE                            # Label peaks
) + ggtitle("scATAC-seq") + scale_fill_gradientn(colors = c( "#002366", "#195905","yellow"))
atac_heatmap <- atac_heatmap +
  theme(
    axis.text.y = element_blank(),  # Remove y-axis text (gene names)
    legend.position = "right"        # Remove legend
  )

range(cluster_data$Correlation)
# 7. Save the heatmaps
ggsave("results/singlecell/PeaknGene_corr_heatmap_RNA.jpg", rna_heatmap, width = 5.5, height = 6, dpi = 300)
ggsave("results/singlecell/PeaknGene_corr_heatmap_ATAC.jpg", atac_heatmap, width = 4.5, height = 5.5, dpi = 300)

# 8. Display the heatmaps
print(rna_heatmap)
print(atac_heatmap)
