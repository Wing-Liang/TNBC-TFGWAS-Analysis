## 07_coaccessibility_cicero.R
## Peak-peak co-accessibility (Cicero) on the discovery ATAC data; the connections and
## cis-co-accessible networks are stored as Signac Links on the ATAC object.
##
## Cicero is run on the ATAC UMAP coordinates, in four chromosome groups to bound memory
## (hg19 chr1-5, chr6-10, chr11-16, chr17-22/X/Y); the links of all groups are combined.
##
## Inputs   results/singlecell/atac_ATAC_RNAlabeled.rds        (05_atac_umap_figures.R)
## Outputs  results/singlecell/atac_links_general_all.rds      ATAC object with Links()
## Run from the repository root:  Rscript scripts/01_singlecell_annotation/07_coaccessibility_cicero.R

TNBC_ROOT <- Sys.getenv("TNBC_ROOT", ".")

library(Signac)
library(Seurat)
library(SeuratWrappers)
library(cicero)        # remotes::install_github("cole-trapnell-lab/cicero-release", ref = "monocle3")
library(monocle3)
library(BSgenome.Hsapiens.UCSC.hg19)

setwd(TNBC_ROOT)

atac_subset_qc <- readRDS("results/singlecell/atac_ATAC_RNAlabeled.rds")
DefaultAssay(atac_subset_qc) <- "ATAC"

## Cicero object on the ATAC UMAP
atac.cds    <- as.cell_data_set(x = atac_subset_qc)
atac.cicero <- make_cicero_cds(atac.cds, reduced_coordinates = reducedDims(atac.cds)$UMAP.ATAC)

## chromosome lengths (hg19), standard chromosomes only
seqinfo(atac_subset_qc) <- seqinfo(BSgenome.Hsapiens.UCSC.hg19)
genome <- seqlengths(atac_subset_qc)
genome <- genome[names(genome) %in% standardChromosomes(granges(atac_subset_qc))]

CHROM_GROUPS <- list(1:5, 6:10, 11:16, 17:24)
links <- lapply(CHROM_GROUPS, function(idx) {
  g <- genome[idx]
  genome.df <- data.frame(chr = names(g), length = g)
  conns <- run_cicero(atac.cicero, genomic_coords = genome.df, sample_num = 100)
  ccans <- generate_ccans(conns)
  ConnectionsToLinks(conns = conns, ccans = ccans)
})
links <- do.call(c, links)

Links(atac_subset_qc) <- links
saveRDS(atac_subset_qc, "results/singlecell/atac_links_general_all.rds")
