# ==============================================================================
# 03_qc_doublets.R
# ATAC and RNA quality control and doublet removal (scDblFinder).
# Produces: final cell set
# Inputs  : data/GSE162798/fragments/fragments.tsv.gz
#           results/singlecell/peakcalling_df.bed
#           results/singlecell/atac_preQC.rds
#           results/singlecell/atac_ATAC_RNAlabeled.rds
#           results/singlecell/final_cell_barcodes.csv
# Outputs : results/singlecell/atac_preQC.rds
#           results/singlecell/01_DensityScatter_ATAC.jpg
#           results/singlecell/02_histogram_nucleosome_signal_NS1.jpg
#           results/singlecell/03_histogram_nucleosome_signal_NS1.2.jpg
#           results/singlecell/04_histogram_nucleosome_signal_NS1.1.jpg
#           results/singlecell/05_VlnPlot_before_filter.jpg
#           results/singlecell/atac_postQC.rds
#           ...
# Run     : Rscript scripts/01_singlecell_annotation/03_qc_doublets.R   (from the repository root)
# ==============================================================================
# Analyzing chrom_assay scATAC-seq: 10x multiomic
# https://stuartlab.org/signac/articles/chrom_assay_vignette
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

# read the peakcalling bed file
library(regioneR)

options(future.globals.maxSize=1e9)

library(dplyr)

library(IRanges)

# server
setwd(TNBC_ROOT)

################################################################
###               (1) Pre-processing workflow                ###
################################################################
##  -----------------Load ATAC data create a new seurat object -----------------
# start by creating a Seurat object using the peak/cell matrix and cell metadata
# Create a fragment object
frags <- CreateFragmentObject(path = "data/GSE162798/fragments/fragments.tsv.gz",
                              cells = colnames(atac_subset))

# First call peaks on the dataset
# If you already have a set of peaks you can skip this step
peakcalling <-  as.data.frame(read.table("results/singlecell/peakcalling_df.bed",header = TRUE, sep="\t",stringsAsFactors=FALSE, quote=""))
peakcalling <- toGRanges(peakcalling, format="BED", header=TRUE)
peakcalling #230904

# remove peaks on nonstandard chromosomes and in genomic blacklist regions
peakcalling  <- keepStandardChromosomes(peakcalling , pruning.mode = "coarse")

# Quantify fragments in each peak
macs2_counts <- FeatureMatrix(fragments = frags,
                              features = peakcalling,
                              cells = colnames(atac_subset))

all.equal(colnames(atac_subset), colnames(macs2_counts))
# ---------------- annotation------------------------
# get gene annotations for hg19
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v75)
seqlevels(annotations) <- paste0('chr', seqlevels(annotations))
genome(annotations) <- "hg19"
seqlevelsStyle(annotations) <- 'UCSC'

# remove the old ATAC assay
DefaultAssay(atac_subset)<-'RNA'
atac_subset[['ATAC']]<-NULL
atac_subset

# create a new assay using the MACS2 peak set and add it to the Seurat object
atac_subset[['ATAC']] <- CreateChromatinAssay(
  counts = macs2_counts,
  sep = c(":", "-"),
  fragments = "data/GSE162798/fragments/fragments.tsv.gz",
  annotation = annotations,
  min.cells = 1)
DefaultAssay(atac_subset)<-'ATAC'
atac_subset

# remove the features that correspond to chromosome scaffolds or other sequences instead of the (22+2) standard chromosomes.
peaks.keep <- seqnames(granges(atac_subset)) %in% standardChromosomes(granges(atac_subset))

atac_subset
saveRDS(atac_subset, 'results/singlecell/atac_preQC.rds')

# see the genomic ranges associated with each feature in the object.
granges(atac_subset)
# GRanges object with 230904 ranges and 0 metadata columns:
#   seqnames            ranges strand
# <Rle>         <IRanges>  <Rle>
#   ...      ...               ...    ...
#   -------
#   seqinfo: 24 sequences from an unspecified genome; no seqlengths

atac_subset<-readRDS('results/singlecell/atac_preQC.rds')

######################################################
###               Drop the Doublets                ###
######################################################

table(atac_subset@meta.data$scDblFinder_score<=0.9)
table(atac_subset@meta.data$scDblFinder_score<=0.95)

atac_subset <- subset(atac_subset, subset = scDblFinder_score <= 0.95)
atac_subset #22123-4092=18031
######################################################
###               (2) computing QC                 ###
######################################################
DefaultAssay(atac_subset) <- "ATAC"

# ------------------- refer to Signac - Analyzing PBMC scATAC-seq https://stuartlab.org/signac/articles/pbmc_vignette -----
# compute TSS enrichment score per cell
# Transcriptional start site (TSS) enrichment score: The ENCODE project has defined
# an ATAC-seq targeting score based on the ratio of fragments centered at the TSS to fragments in TSS-flanking regions
atac_subset <- TSSEnrichment(atac_subset)

# compute nucleosome signal score per cell
# quantify the approximate ratio of mononucleosomal to nucleosome-free fragments
# Nucleosome banding pattern: The histogram of DNA fragment sizes (determined from the
# paired-end sequencing reads) should exhibit a strong nucleosome banding pattern corresponding to the length of DNA wrapped around a single nucleosome
atac_subset <- NucleosomeSignal(atac_subset)

# add blacklist ratio -
# Find the fraction of counts per cell that overlap a given set of genomic ranges
# Ratio reads in genomic blacklist regions:
# The ENCODE project has provided a list of blacklist regions, representing reads
# which are often associated with artefactual signal.
atac_subset$blacklist_fraction <- FractionCountsInRegion(
  object = atac_subset,
  assay = 'ATAC',
  regions = blacklist_hg19)

jpeg('results/singlecell/01_DensityScatter_ATAC.jpg',width=9, height = 7, units = 'in', res = 300)
DensityScatter(atac_subset, x = 'nCount_ATAC', y = 'TSS.enrichment', log_x = TRUE, quantiles = TRUE)
dev.off()

jpeg('results/singlecell/02_histogram_nucleosome_signal_NS1.jpg',width=12, height = 7, units = 'in', res = 300)
atac_subset$nucleosome_group <- ifelse(atac_subset$nucleosome_signal > 1, 'NS > 1', 'NS < 1')
FragmentHistogram(object = atac_subset, group.by = 'nucleosome_group')
dev.off()

jpeg('results/singlecell/03_histogram_nucleosome_signal_NS1.2.jpg',width=12, height = 7, units = 'in', res = 300)
atac_subset$nucleosome_group <- ifelse(atac_subset$nucleosome_signal > 1.2, 'NS > 1.2', 'NS < 1.2')
FragmentHistogram(object = atac_subset, group.by = 'nucleosome_group')
dev.off()

jpeg('results/singlecell/04_histogram_nucleosome_signal_NS1.1.jpg',width=12, height = 7, units = 'in', res = 300)
atac_subset$nucleosome_group <- ifelse(atac_subset$nucleosome_signal > 1.1, 'NS > 1.1', 'NS < 1.1')
FragmentHistogram(object = atac_subset, group.by = 'nucleosome_group')
dev.off()

jpeg('results/singlecell/05_VlnPlot_before_filter.jpg',width=14, height = 6, units = 'in', res = 300)
VlnPlot(
  object = atac_subset,
  features = c('nCount_RNA','nCount_ATAC', 'TSS.enrichment','blacklist_fraction',
               'nucleosome_signal'),
  pt.size = 0,
  ncol = 5)
dev.off()

# filter out low quality cells
atac_subset_qc <- subset(
  x = atac_subset,
  subset = nCount_ATAC < 100000 &
    nCount_RNA < 25000 &
    nCount_ATAC > 1000 &
    nCount_RNA > 1000 &
    nucleosome_signal < 1.2 &
    TSS.enrichment > 1 &
    blacklist_fraction < 0.02)
atac_subset_qc

saveRDS(atac_subset_qc, 'results/singlecell/atac_postQC.rds')

# postQC
atac_subset_qc<-readRDS('results/singlecell/atac_ATAC_RNAlabeled.rds')
Final_cell_barcode<-fread('results/singlecell/final_cell_barcodes.csv',header=TRUE)
atac_subset_qc <- subset(atac_subset_qc , cells = Final_cell_barcode$x)
DefaultAssay(atac_subset_qc)<-'ATAC'

DefaultAssay(atac_subset)<-'ATAC'
jpeg('results/singlecell/05_VlnPlot_after_filter.jpg',width=16, height = 10, units = 'in', res = 300)
p1<-VlnPlot(atac_subset,features = 'nCount_RNA', ncol=1,pt.size=0) + geom_hline(yintercept = 25000, color='red4') + theme(legend.position = 'none')
p2<-VlnPlot(atac_subset,features = 'nCount_ATAC', ncol=1,pt.size=0) + geom_hline(yintercept = 100000, color='red4')+ theme(legend.position = 'none')
p3<-VlnPlot(atac_subset,features = 'TSS.enrichment', ncol=1,pt.size=0) + geom_hline(yintercept = 1, color='red4')+ theme(legend.position = 'none')
p4<-VlnPlot(atac_subset,features = 'blacklist_fraction', ncol=1,pt.size=0) + geom_hline(yintercept = 0.02, color='red4')+ theme(legend.position = 'none')
p5<-VlnPlot(atac_subset,features = 'nucleosome_signal', ncol=1,pt.size=0) + geom_hline(yintercept = 1.2, color='red4')+ theme(legend.position = 'none')
p6<-VlnPlot(
  object = atac_subset_qc,
  features = c('nCount_RNA','nCount_ATAC', 'TSS.enrichment',
               'blacklist_fraction','nucleosome_signal'),
  pt.size = 0,
  ncol = 5)
design<-"ABCDE
         FFFFF"
wrap_plots(A=p1,B=p2,C=p3,D=p4,E=p5,F=p6,ncol = 5,nrow = 2, design = design)
dev.off()

DefaultAssay(atac_subset_qc)<-'RNA'
Idents(atac_subset_qc)<-'RNA'
jpeg('results/singlecell/05_VlnPlot_after_filter.jpg',width=16, height = 6, units = 'in', res = 300)
VlnPlot(
  object = atac_subset_qc,
  features = c('nCount_RNA','nCount_ATAC', 'TSS.enrichment',
               'blacklist_fraction','nucleosome_signal'),
  pt.size = 0,
  ncol = 5)
dev.off()
