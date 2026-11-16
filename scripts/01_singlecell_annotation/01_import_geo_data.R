# ==============================================================================
# 01_import_geo_data.R
# Import the GSE162798 scRNA-seq and scATAC-seq matrices, fragments and cell metadata.
# Produces: input to all single-cell steps
# Inputs  : data/GSE162798/scATAC/barcodes.tsv.gz
#           data/GSE162798/scATAC/peaks.bed.gz
#           data/GSE162798/scRNA/
#           data/GSE162798/GSM5484484_lib.coassay.std_cellname_dictionary.csv
# Outputs : results/singlecell/clean_ATAC.rds
#           results/singlecell/clean_RNA.rds
# Run     : Rscript scripts/01_singlecell_annotation/01_import_geo_data.R   (from the repository root)
# ==============================================================================
## ---- paths (set the environment variables; see README) ----
TNBC_ROOT <- Sys.getenv("TNBC_ROOT", ".")

setwd(TNBC_ROOT)
library(Signac)
library(Seurat)
library(GenomeInfoDb)
library(ggplot2)
library(patchwork)
library(dplyr)
library(Matrix)
library(GenomicRanges)

set.seed(2023)

### Pre-processing workflow

## Load in ATAC data

barcodes <- read.table("data/GSE162798/scATAC/barcodes.tsv.gz", header = F, sep = "\t")
matrix0 <- readMM("data/GSE162798/scATAC/matrix.mtx.gz")
peaks0 <- read.table("data/GSE162798/scATAC/peaks.bed.gz", header = FALSE, sep = "\t")
peaks1 = paste0(peaks0$V1, ':', peaks0$V2, '-', peaks0$V3)

colnames(matrix0) = barcodes$V1
rownames(matrix0) =peaks1

metadata <- read.csv(
  file = "data/GSE162798/GSM5484482_SNuBarARC-HBCA.scATAC.cell_metainfo.csv",
  header = TRUE,
  row.names = 1
)

chrom_assay <- CreateChromatinAssay(
  counts = matrix0,
  sep = c(":", "-"),
  genome = 'hg19',
  fragments = 'data/GSE162798/fragments/fragments.tsv.gz',
  min.cells = 1
)

HBCA_ATAC <- CreateSeuratObject(
  counts = chrom_assay,
  assay = 'peaks',
  project = 'ATAC',
  meta.data = metadata
)

HBCA_ATAC

## Load in RNA data
rna_matrix <- Read10X(data.dir = 'data/GSE162798/scRNA/')
dim(rna_matrix)
HBCA_RNA <- CreateSeuratObject(counts = rna_matrix, assay = 'RNA', min.cells = 5, min.features = 200, project = 'RNA')

## Find overlap barcodes between RNA and ATAC
dic0 = read.csv('data/GSE162798/GSM5484484_lib.coassay.std_cellname_dictionary.csv')
dic0$RNA = paste0(dic0$RNA_original_cell_name, '-1')
dic0$ATAC = paste0(dic0$ATAC_original_cell_name, '-1')

dic_ATAC = dic0[is.element(dic0$ATAC, colnames(HBCA_ATAC)),]
dic_overlap = dic_ATAC[is.element(dic_ATAC$RNA, colnames(HBCA_RNA)),]

## Subsetting data to common barcodes and sorting
HBCA_ATAC = subset(HBCA_ATAC, cells=dic_overlap$ATAC)
HBCA_RNA = subset(HBCA_RNA, cells=dic_overlap$RNA)
colnames(HBCA_RNA) = colnames(HBCA_ATAC)
HBCA_RNA = RenameCells(HBCA_RNA, new.names = colnames(HBCA_ATAC))

saveRDS(HBCA_ATAC, 'results/singlecell/clean_ATAC.rds')
saveRDS(HBCA_RNA, 'results/singlecell/clean_RNA.rds')
