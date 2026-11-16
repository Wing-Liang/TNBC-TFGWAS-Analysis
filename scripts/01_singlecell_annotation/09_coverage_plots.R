# ==============================================================================
# 09_coverage_plots.R
# Coverage plots of marker-gene loci with co-accessibility links.
# Produces: supporting plots
# Inputs  : results/singlecell/atac_links_general_all.rds
#           results/singlecell/final_cell_barcodes.csv
# Outputs : results/singlecell/CovPlot_<gene>.jpg
# Run     : Rscript scripts/01_singlecell_annotation/09_coverage_plots.R   (from the repository root)
# ==============================================================================
## ---- paths (set the environment variables; see README) ----
TNBC_ROOT <- Sys.getenv("TNBC_ROOT", ".")
library(Signac)
library(Seurat)
library(ggplot2)
library(data.table)

setwd(TNBC_ROOT)
## ATAC object with Cicero co-accessibility links (07_coaccessibility_cicero.R)
atac_subset_qc<-readRDS('results/singlecell/atac_links_general_all.rds')
Final_cell_barcode<-fread('results/singlecell/final_cell_barcodes.csv',header=TRUE)
atac_subset_qc <- subset(atac_subset_qc , cells = Final_cell_barcode$x)
DefaultAssay(atac_subset_qc)<-'ATAC'


#######################################################
#          Visualization of genomic regions           #
#######################################################
# List of genes to plot
genes <- c('FOXA1', 'TP63', 'EPCAM', 'KIT', 'ITGA6', 'ESR1', 'GATA3')

# Ensure the output directory exists
dir.create("results/singlecell/", recursive = TRUE, showWarnings = FALSE)

# Loop through each region and generate a CoveragePlot
for (g in genes) {
  # Generate the CoveragePlot
  cov_plot <- CoveragePlot(
    object = atac_subset_qc,
    region = g,
    extend.upstream = 500,
    extend.downstream = 500,
    annotation = TRUE,
    peaks = TRUE,
    features = g
  )

  # Save the plot as a JPEG file
  jpeg(filename = paste0("results/singlecell/CovPlot_", g, ".jpg"),
       width = 8, height = 10, units = "in", res = 200)
  print(cov_plot)
  dev.off()
}

CoveragePlot(
  object = atac_subset_qc,
  region = "chr16-79622735-79639634",
  annotation = TRUE,
  peaks = TRUE,
  features = "MAF",
  links = TRUE
)
