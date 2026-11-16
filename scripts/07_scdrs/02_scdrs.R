# ==============================================================================
# 02_scdrs.R
# MAGMA gene-level GWAS statistics and scDRS cell-level disease scores per subtype.
# Produces: Fig S10
# Inputs  : results/singlecell/atac_RNA_annotated.rds
#           results/singlecell/final_cell_barcodes.csv
#           results/gwas/bcac_gwas_filtered_TNBC.csv
#           results/gwas/bcac_gwas_filtered_LumA.csv
#           results/gwas/bcac_gwas_filtered_LumB.csv
#           ...
# Outputs : results/scdrs/rna_breastcancer.h5ad
#           results/scdrs/covariates.cov
#           results/scdrs/scDRS_results_all_BCAC.csv
# Run     : Rscript scripts/07_scdrs/02_scdrs.R   (from the repository root)
# ==============================================================================
## ---- paths (set the environment variables; see README) ----
TNBC_ROOT <- Sys.getenv("TNBC_ROOT", ".")
MAGMA_DIR <- Sys.getenv("MAGMA_DIR", "data/MAGMA")
Sys.setenv(MAGMA_DIR = MAGMA_DIR)   # used by the scdrs shell commands below
library(reticulate)
py_config()
# reinstall in bash
# pip install numpy==1.26.4 --force-reinstall
# pip install --force-reinstall scdrs anndata scanpy pandas scipy==1.14.1
# pip install --upgrade scdrs
use_python(Sys.which("python"), required = TRUE)
py_run_string("import sys; print(sys.version)")
py_run_string("import scdrs; print(scdrs.__version__)")
py_run_string("
import numpy as np
import scipy
import scdrs
import anndata
print('NumPy version:', np.__version__)
print('SciPy version:', scipy.__version__)
print('scDRS version:', scdrs.__version__)
print('AnnData version:', anndata.__version__)
")
#
# NOTE: Python version was forced by use_python() function
set.seed(202510)

##########################################
#             Prepare input files        #
##########################################
library(Seurat)
library(SeuratDisk)
library(data.table)

# Load your Seurat object
setwd(TNBC_ROOT)
rna <- readRDS('results/singlecell/atac_RNA_annotated.rds')
Final_cell_barcode<-fread('results/singlecell/final_cell_barcodes.csv',header=TRUE)
rna <- subset(rna , cells = Final_cell_barcode$x)
DefaultAssay(rna)<-'RNA'
rna[['RNA']]<-as(rna[['RNA']], Class="Assay5")

rna@assays$RNA

# Assay (v5) data with 32738 features for 16451 cells
# First 10 features:
#   MIR1302-10, FAM138A, OR4F5, RP11-34P13.7, RP11-34P13.8, AL627309.1,
# RP11-34P13.14, RP11-34P13.9, AP006222.2, RP4-669L17.10
# Layers:
#   counts, data, scale.data

Layers(rna, assay = "RNA")

head(GetAssayData(rna, assay = "RNA", layer = "data")[, 1:5])

colnames(rna@meta.data)

names(rna@graphs)

# Ensure Python environment with anndata and scipy is available
# Assuming you're using the scdrs_env conda environment

# (optional) step 1: Slim down a Seurat object. So you get raw counts, lognorm counts
seu = DietSeurat(
  rna,
  counts = TRUE, # so, raw counts save to adata.layers['counts']
  data = TRUE, # so, log1p counts save to adata.X when scale.data = False, else adata.layers['data']
  scale.data = FALSE, # set to false
  features = rownames(rna), # export all genes, not just top highly variable genes
  assays = "RNA",
  dimreducs = c("pca","umap.rna"),
  graphs = c("SCT_nn", "SCT_snn"), # to RNA_nn -> distances, RNA_snn -> connectivities
  misc = TRUE
)

# Step 2. Write Seurat to h5ad
# for seurat v5, need to JoinLayer first
DefaultAssay(seu) = "RNA"
seu <- JoinLayers(seu)
# single modality
MuDataSeurat::WriteH5AD(seu, "results/scdrs/rna_breastcancer.h5ad", assay="RNA")

# read .h5ad in python
# import anndata
# # -------------Read the .h5ad file--------------------
# # -----------Inspect the AnnData object---------------------
# # View observation metadata (e.g., cell annotations)
# # -------View variable metadata (e.g., gene names)----------------------
# #----------- Access the data matrix (sparse or dense)
# # -------------Check available layers----------------
# # ----------------Access embeddings (e.g., UMAP coordinates)-----------------------
# if "X_umap" in adata.obsm:

# GWAS data
library(dplyr)
library(data.table)
# TNBC
TNBC_gwas <- fread("results/gwas/bcac_gwas_filtered_TNBC.csv", header = TRUE)
TNBC_gwas <- TNBC_gwas %>% select(chrom, pos, rsid, p)
names(TNBC_gwas)<-c('CHR','BP','SNP','P')
TNBC_gwas$N <-106377
TNBC_gwas$SNP <- ifelse(TNBC_gwas$SNP == "", paste0("Chr", TNBC_gwas$CHR, "_", TNBC_gwas$BP), TNBC_gwas$SNP)
fwrite(TNBC_gwas, file.path(MAGMA_DIR, "input_gwas/TNBC.pval"), sep = "\t", row.names = FALSE)
# LumA
LumA_gwas <- fread("results/gwas/bcac_gwas_filtered_LumA.csv", header = TRUE)
LumA_gwas <- LumA_gwas %>% select(chrom, pos, rsid, p)
names(LumA_gwas)<-c('CHR','BP','SNP','P')
LumA_gwas$N <-148877
LumA_gwas$SNP <- ifelse(LumA_gwas$SNP == "", paste0("Chr", LumA_gwas$CHR, "_", LumA_gwas$BP), LumA_gwas$SNP)
fwrite(LumA_gwas, file.path(MAGMA_DIR, "input_gwas/LumA.pval"), sep = "\t", row.names = FALSE)
# LumB
LumB_gwas <- fread("results/gwas/bcac_gwas_filtered_LumB.csv", header = TRUE)
LumB_gwas <- LumB_gwas %>% select(chrom, pos, rsid, p)
names(LumB_gwas)<-c('CHR','BP','SNP','P')
LumB_gwas$N <-105277
LumB_gwas$SNP <- ifelse(LumB_gwas$SNP == "", paste0("Chr", LumB_gwas$CHR, "_", LumB_gwas$BP), LumB_gwas$SNP)
fwrite(LumB_gwas, file.path(MAGMA_DIR, "input_gwas/LumB.pval"), sep = "\t", row.names = FALSE)
# LumB_HER2Neg
LumB_HER2Neg_gwas <- fread("results/gwas/bcac_gwas_filtered_LumB_HER2Neg.csv", header = TRUE)
LumB_HER2Neg_gwas <- LumB_HER2Neg_gwas %>% select(chrom, pos, rsid, p)
names(LumB_HER2Neg_gwas)<-c('CHR','BP','SNP','P')
LumB_HER2Neg_gwas$N <-105277
LumB_HER2Neg_gwas$SNP <- ifelse(LumB_HER2Neg_gwas$SNP == "",
                                paste0("Chr", LumB_HER2Neg_gwas$CHR, "_", LumB_HER2Neg_gwas$BP), LumB_HER2Neg_gwas$SNP)
fwrite(LumB_HER2Neg_gwas, file.path(MAGMA_DIR, "input_gwas/LumB_HER2Neg.pval"), sep = "\t", row.names = FALSE)
# HER2
HER2_gwas <- fread("results/gwas/bcac_gwas_filtered_HER2.csv", header = TRUE)
HER2_gwas <- HER2_gwas %>% select(chrom, pos, rsid, p)
names(HER2_gwas)<-c('CHR','BP','SNP','P')
HER2_gwas$N <-105277
HER2_gwas$SNP <- ifelse(HER2_gwas$SNP == "",
                                paste0("Chr", HER2_gwas$CHR, "_", HER2_gwas$BP), HER2_gwas$SNP)
fwrite(HER2_gwas, file.path(MAGMA_DIR, "input_gwas/HER2.pval"), sep = "\t", row.names = FALSE)
# MAGMA instruction
# https://github.com/martinjzhang/scDRS/issues/2

# # Step 2: make gene annotation file for MAGMA using the following command, this only needs to be done
# # once for different GWAS summary statistics, and the file will be saved to out/step1.genes.annot
# cd $MAGMA_DIR
# mkdir -p out/step1
# $magma_v1.10/magma \
# --annotate window=10,10 \
# --snp-loc g1000_eur/g1000_eur.bim \
# --gene-loc NCBI37.1/NCBI37.3.gene.loc \
# --out out/step1

# magma_v1.10/magma \
# --annotate window=10,10 \
# --snp-loc g1000_eur/g1000_eur.bim \
# --gene-loc NCBI38/NCBI38.gene.loc \
# --out out/step1_hg38

# in Bash
#./magma_v1.10/magma --version
# magma_v1.10/magma --bfile g1000_eur/g1000_eur --pval input_gwas/TNBC.pval use='SNP,P' ncol='N' --gene-annot out/step1.genes.annot --out out/TNBC
# magma_v1.10/magma --bfile g1000_eur/g1000_eur --pval input_gwas/LumA.pval use='SNP,P' ncol='N' --gene-annot out/step1.genes.annot --out out/LumA
# magma_v1.10/magma --bfile g1000_eur/g1000_eur --pval input_gwas/LumB.pval use='SNP,P' ncol='N' --gene-annot out/step1.genes.annot --out out/LumB
# magma_v1.10/magma --bfile g1000_eur/g1000_eur --pval input_gwas/LumB_HER2Neg.pval use='SNP,P' ncol='N' --gene-annot out/step1.genes.annot --out out/LumB_HER2Neg
# magma_v1.10/magma --bfile g1000_eur/g1000_eur --pval input_gwas/HER2.pval use='SNP,P' ncol='N' --gene-annot out/step1.genes.annot --out out/HER2
# check output
# head out/LumA.genes.out

# Step 3: Convert Ensembl IDs to HGNC Symbols
library(biomaRt)
library(data.table)
library(dplyr)
# ------------TNBC---------------
magma_genes <- fread(file.path(MAGMA_DIR, "out/TNBC.genes.out"), header = TRUE)
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
gene_map <- getBM(
  attributes = c("entrezgene_id", "hgnc_symbol"),
  filters = "entrezgene_id",
  values = magma_genes$GENE,
  mart = ensembl
)
# Remove empty or NA HGNC symbols
gene_map <- gene_map[!is.na(gene_map$hgnc_symbol) & gene_map$hgnc_symbol != "", ]
# Merge with MAGMA output
magma_genes <- merge(magma_genes, gene_map, by.x = "GENE", by.y = "entrezgene_id")
# Select relevant columns
magma_genes <- magma_genes[, c("hgnc_symbol", "ZSTAT")]
#18189 genes
colnames(magma_genes) <- c("GENE", "ZSTAT")
magma_genes_unique <- magma_genes %>%
  group_by(GENE) %>%
  summarise(ZSTAT = ZSTAT[which.max(abs(ZSTAT))]) %>%
  ungroup()
# 18113 genes
magma_genes_unique <- magma_genes_unique[order(-magma_genes_unique$ZSTAT),]
fwrite(magma_genes_unique, file.path(MAGMA_DIR, "out/TNBC_hgnc.genes.out"), sep = "\t", row.names = FALSE)

# --------------LumA --------------------
magma_genes <- fread(file.path(MAGMA_DIR, "out/LumA.genes.out"), header = TRUE)
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
gene_map <- getBM(
  attributes = c("entrezgene_id", "hgnc_symbol"),
  filters = "entrezgene_id",
  values = magma_genes$GENE,
  mart = ensembl
)
# Remove empty or NA HGNC symbols
gene_map <- gene_map[!is.na(gene_map$hgnc_symbol) & gene_map$hgnc_symbol != "", ]
# Merge with MAGMA output
magma_genes <- merge(magma_genes, gene_map, by.x = "GENE", by.y = "entrezgene_id")
# Select relevant columns
magma_genes <- magma_genes[, c("hgnc_symbol", "ZSTAT")]
#18189 genes
colnames(magma_genes) <- c("GENE", "ZSTAT")
magma_genes_unique <- magma_genes %>%
  group_by(GENE) %>%
  summarise(ZSTAT = ZSTAT[which.max(abs(ZSTAT))]) %>%
  ungroup()
# 18113 genes
magma_genes_unique <- magma_genes_unique[order(-magma_genes_unique$ZSTAT),]
fwrite(magma_genes_unique, file.path(MAGMA_DIR, "out/LumA_hgnc.genes.out"), sep = "\t", row.names = FALSE)
# --------------LumB --------------------
magma_genes <- fread(file.path(MAGMA_DIR, "out/LumB.genes.out"), header = TRUE)
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
gene_map <- getBM(
  attributes = c("entrezgene_id", "hgnc_symbol"),
  filters = "entrezgene_id",
  values = magma_genes$GENE,
  mart = ensembl
)
# Remove empty or NA HGNC symbols
gene_map <- gene_map[!is.na(gene_map$hgnc_symbol) & gene_map$hgnc_symbol != "", ]
# Merge with MAGMA output
magma_genes <- merge(magma_genes, gene_map, by.x = "GENE", by.y = "entrezgene_id")
# Select relevant columns
magma_genes <- magma_genes[, c("hgnc_symbol", "ZSTAT")]
#18189 genes
colnames(magma_genes) <- c("GENE", "ZSTAT")
magma_genes_unique <- magma_genes %>%
  group_by(GENE) %>%
  summarise(ZSTAT = ZSTAT[which.max(abs(ZSTAT))]) %>%
  ungroup()
# 18113 genes
magma_genes_unique <- magma_genes_unique[order(-magma_genes_unique$ZSTAT),]
fwrite(magma_genes_unique, file.path(MAGMA_DIR, "out/LumB_hgnc.genes.out"), sep = "\t", row.names = FALSE)
# --------------LumB_HER2Neg --------------------
magma_genes <- fread(file.path(MAGMA_DIR, "out/LumB_HER2Neg.genes.out"), header = TRUE)
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
gene_map <- getBM(
  attributes = c("entrezgene_id", "hgnc_symbol"),
  filters = "entrezgene_id",
  values = magma_genes$GENE,
  mart = ensembl
)
# Remove empty or NA HGNC symbols
gene_map <- gene_map[!is.na(gene_map$hgnc_symbol) & gene_map$hgnc_symbol != "", ]
# Merge with MAGMA output
magma_genes <- merge(magma_genes, gene_map, by.x = "GENE", by.y = "entrezgene_id")
# Select relevant columns
magma_genes <- magma_genes[, c("hgnc_symbol", "ZSTAT")]
#18189 genes
colnames(magma_genes) <- c("GENE", "ZSTAT")
magma_genes_unique <- magma_genes %>%
  group_by(GENE) %>%
  summarise(ZSTAT = ZSTAT[which.max(abs(ZSTAT))]) %>%
  ungroup()
# 18113 genes
magma_genes_unique <- magma_genes_unique[order(-magma_genes_unique$ZSTAT),]
fwrite(magma_genes_unique, file.path(MAGMA_DIR, "out/LumB_HER2Neg_hgnc.genes.out"), sep = "\t", row.names = FALSE)
# --------------HER2 --------------------
magma_genes <- fread(file.path(MAGMA_DIR, "out/HER2.genes.out"), header = TRUE)
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
gene_map <- getBM(
  attributes = c("entrezgene_id", "hgnc_symbol"),
  filters = "entrezgene_id",
  values = magma_genes$GENE,
  mart = ensembl
)
# Remove empty or NA HGNC symbols
gene_map <- gene_map[!is.na(gene_map$hgnc_symbol) & gene_map$hgnc_symbol != "", ]
# Merge with MAGMA output
magma_genes <- merge(magma_genes, gene_map, by.x = "GENE", by.y = "entrezgene_id")
# Select relevant columns
magma_genes <- magma_genes[, c("hgnc_symbol", "ZSTAT")]
#18189 genes
colnames(magma_genes) <- c("GENE", "ZSTAT")
magma_genes_unique <- magma_genes %>%
  group_by(GENE) %>%
  summarise(ZSTAT = ZSTAT[which.max(abs(ZSTAT))]) %>%
  ungroup()
# 18113 genes
magma_genes_unique <- magma_genes_unique[order(-magma_genes_unique$ZSTAT),]
fwrite(magma_genes_unique, file.path(MAGMA_DIR, "out/HER2_hgnc.genes.out"), sep = "\t", row.names = FALSE)
# Step 4: Create .gs File with scdrs munge-gs
library(reticulate)
py_config()  # Verify Python environment

#TNBC
munge_command <- sprintf("scdrs munge-gs --zscore-file ${MAGMA_DIR}/out/TNBC_hgnc.genes.out --weight zscore --out-file results/scdrs/TNBC.gs")
py_run_string(sprintf("
import subprocess
subprocess.run('%s', shell=True, check=True)
", munge_command))
# Print info for the first 3 traits and first 10 genes
# Traits               ['ZSTAT']
# MDM4                 [7.9462]
# DDA1                 [7.6079]
# ABHD8                [7.3853]
# NR2F6                [7.186]
# USHBP1               [7.1567]
# TERT                 [7.0939]
# ANO8                 [7.0776]
# MRTFA                [6.712]
# POGLUT3              [6.4748]
# TCF7L2               [6.4471]
# --zscore-file have values above 1 or below 0. Seems fine.
gs <- fread("results/scdrs/TNBC.gs", sep = "\t", header = TRUE)
head(gs)

#LumA
munge_command <- sprintf("scdrs munge-gs --zscore-file ${MAGMA_DIR}/out/LumA_hgnc.genes.out --weight zscore --out-file results/scdrs/LumA.gs")
py_run_string(sprintf("
import subprocess
subprocess.run('%s', shell=True, check=True)
", munge_command))
# --zscore-file loaded: n_gene=18113, n_trait=1 (sys_time=0.1s)
# Print info for the first 3 traits and first 10 genes
# Traits               ['ZSTAT']
# RAD51B               [11.307]
# ADAM29               [9.1602]
# FGFR2                [8.863]
# ITPR1                [8.5996]
# EBF1                 [8.3456]
# TOX3                 [8.1856]
# SLC4A7               [8.0254]
# ESR1                 [7.9621]
# SKIDA1               [7.911]
# ZNF365               [7.9002]
# --zscore-file have values above 1 or below 0. Seems fine.
gs <- fread("results/scdrs/LumA.gs", sep = "\t", header = TRUE)
head(gs)
#LumB
munge_command <- sprintf("scdrs munge-gs --zscore-file ${MAGMA_DIR}/out/LumB_hgnc.genes.out --weight zscore --out-file results/scdrs/LumB.gs")
py_run_string(sprintf("
import subprocess
subprocess.run('%s', shell=True, check=True)
", munge_command))
# --zscore-file loaded: n_gene=18113, n_trait=1 (sys_time=0.0s)
# Print info for the first 3 traits and first 10 genes
# Traits               ['ZSTAT']
# TOX3                 [7.069]
# FGFR2                [6.4526]
# MAP3K1               [6.2435]
# TNNT3                [6.2319]
# LSP1                 [6.2238]
# PEX14                [5.9878]
# SETD9                [5.7454]
# MIER3                [5.5424]
# DFFA                 [5.487]
# ZNF365               [5.2939]
# --zscore-file have values above 1 or below 0. Seems fine.
gs <- fread("results/scdrs/LumB.gs", sep = "\t", header = TRUE)
head(gs)
#LumB_HER2Neg
munge_command <- sprintf("scdrs munge-gs --zscore-file ${MAGMA_DIR}/out/LumB_HER2Neg_hgnc.genes.out --weight zscore --out-file results/scdrs/LumB_HER2Neg.gs")
py_run_string(sprintf("
import subprocess
subprocess.run('%s', shell=True, check=True)
", munge_command))
# Print info for the first 3 traits and first 10 genes
# Traits               ['ZSTAT']
# FGFR2                [7.6897]
# NEK10                [7.1685]
# SETD9                [6.5138]
# MAP3K1               [6.2869]
# MIER3                [6.2541]
# RCCD1                [6.1071]
# PRC1                 [5.9555]
# TOX3                 [5.8855]
# ANXA10               [5.7754]
# ZBTB38               [5.101]
# --zscore-file have values above 1 or below 0. Seems fine.
gs <- fread("results/scdrs/LumB_HER2Neg.gs", sep = "\t", header = TRUE)
#HER2
munge_command <- sprintf("scdrs munge-gs --zscore-file ${MAGMA_DIR}/out/HER2_hgnc.genes.out --weight zscore --out-file results/scdrs/HER2.gs")
py_run_string(sprintf("
import subprocess
subprocess.run('%s', shell=True, check=True)
", munge_command))
# Print info for the first 3 traits and first 10 genes
# Traits               ['ZSTAT']
# PEX14                [6.3813]
# TOX3                 [6.0492]
# DFFA                 [5.5405]
# LGR6                 [5.0345]
# IKBIP                [4.4799]
# APAF1                [4.2251]
# CENPS-CORT           [3.9205]
# CENPS                [3.9193]
# SLC25A3              [3.8686]
# SYNE1                [3.7143]
# --zscore-file have values above 1 or below 0. Seems fine.
gs <- fread("results/scdrs/HER2.gs", sep = "\t", header = TRUE)
#-------- Validate Genes ----------
library(reticulate)
library(data.table)
py_run_string("import pip; pip.main(['install', 'anndata'])")
py_run_string("import anndata; print(anndata.__version__)")
py_run_string("
import anndata
adata = anndata.read_h5ad('results/scdrs/rna_breastcancer.h5ad')
print('Gene names:', list(adata.var_names[:10]))
print('Layers:', list(adata.layers.keys()))
print('Obs columns:', list(adata.obs.columns))
print('Var columns:', list(adata.var.columns))
")
# Gene names: ['MIR1302-10', 'FAM138A', 'OR4F5', 'RP11-34P13.7', 'RP11-34P13.8', 'AL627309.1', 'RP11-34P13.14', 'RP11-34P13.9', 'AP006222.2', 'RP4-669L17.10']
# Layers: ['counts']
# Obs columns: ['orig.ident', 'nCount_RNA', 'nFeature_RNA', 'scDblFinder_class', 'scDblFinder_score', 'nCount_ATAC', 'nFeature_ATAC', 'TSS.enrichment', 'TSS.percentile', 'nucleosome_signal', 'nucleosome_percentile', 'blacklist_fraction', 'nucleosome_group', 'nCount_SCT', 'nFeature_SCT', 'SCT_snn_res.0.06', 'seurat_clusters', 'celltype_RNA']
# Var columns: []

py_run_string("
import anndata
adata = anndata.read_h5ad('results/scdrs/rna_breastcancer.h5ad')
print('Counts sample:', adata.layers['counts'][:100, :100].toarray().sum())
")

gs <- fread("results/scdrs/TNBC.gs", sep = "\t", header = TRUE)
gs_string <- paste(strsplit(gs$GENESET, ",")[[1]], collapse=",")
py_run_string(sprintf("
import anndata
import numpy as np
adata = anndata.read_h5ad('results/scdrs/rna_breastcancer.h5ad')
gs = '%s'
gs_genes = [x.split(':')[0] for x in gs.split(',')]
h5ad_genes = adata.var_names
overlap = len(np.intersect1d(gs_genes, h5ad_genes)) / len(gs_genes)
print('Gene overlap:', overlap)
print('Matched genes sample:', np.intersect1d(gs_genes, h5ad_genes)[:10])
", gs_string))
# Gene overlap: 0.923
# Matched genes sample: ['AADAC' 'ABAT' 'ABCB4' 'ABCB8' 'ABHD15' 'ABHD8' 'ACADVL' 'ACAT1' 'ACBD4'
#                        'ACSBG1']

# Run scdrs compute-score
# Create a Covariate File
py_run_string("
import anndata
import pandas as pd
adata = anndata.read_h5ad('results/scdrs/rna_breastcancer.h5ad')
cov_df = adata.obs[['nCount_RNA', 'celltype_RNA']].reset_index()
cov_df.columns = ['CELL', 'NCOUNT_RNA', 'CLUSTER']
# Convert CLUSTER to dummy variables
cov_df = pd.get_dummies(cov_df, columns=['CLUSTER'], prefix='CLUSTER')
cov_df.to_csv('results/scdrs/covariates.cov', sep='\t', index=False)
")
cov_df<-read.table("results/scdrs/covariates.cov", sep = "\t", header = TRUE)

# compute-score -
# ----------------TNBC------------------
h5ad_file <- "results/scdrs/rna_breastcancer.h5ad"
gs_file <- "results/scdrs/TNBC.gs"
out_dir <- "results/scdrs/output_TNBC/"
dir.create(out_dir, recursive = TRUE)

cli_command <- sprintf(
  "scdrs compute-score --h5ad-file %s --h5ad-species human --gs-file %s --gs-species human --flag-filter-data True --flag-raw-count True --n-ctrl 1000 --flag-return-ctrl-raw-score False --flag-return-ctrl-norm-score True --out-folder %s",
  h5ad_file, gs_file,  out_dir
)
system(cli_command)
# Loading data:
#   --h5ad-file loaded: n_cell=16451, n_gene=18916 (sys_time=3.5s)
# First 3 cells: ['ACAGCGGGTAACGTGA-1', 'ACAGCGGGTCGGTTTG-1', 'ACAGCGGGTTCCTCAC-1']
# First 5 genes: ['RP11-34P13.7', 'AL627309.1', 'AP006222.2', 'RP4-669L17.10', 'RP5-857K21.4']
# --gs-file loaded: n_trait=1 (sys_time=3.5s)
# Print info for first 3 traits:
#   First 3 elements for 'ZSTAT': ['MDM4', 'DDA1', 'ABHD8'], [7.9462, 7.6079, 7.3853]
#
# Preprocessing:
#
#   Computing scDRS score:
#   Computing control scores: 100%|██████████| 1000/1000 [01:11<00:00, 13.98it/s]

list.files("results/scdrs/output_TNBC/")

# ----------------LumA------------------
h5ad_file <- "results/scdrs/rna_breastcancer.h5ad"
gs_file <- "results/scdrs/LumA.gs"
out_dir <- "results/scdrs/output_LumA/"
dir.create(out_dir, recursive = TRUE)

cli_command <- sprintf(
  "scdrs compute-score --h5ad-file %s --h5ad-species human --gs-file %s --gs-species human --flag-filter-data True --flag-raw-count True --n-ctrl 1000 --flag-return-ctrl-raw-score False --flag-return-ctrl-norm-score True --out-folder %s",
  h5ad_file, gs_file,  out_dir
)
system(cli_command)

list.files("results/scdrs/output_LumA/")

# ----------------LumB------------------
h5ad_file <- "results/scdrs/rna_breastcancer.h5ad"
gs_file <- "results/scdrs/LumB.gs"
out_dir <- "results/scdrs/output_LumB/"
dir.create(out_dir, recursive = TRUE)

cli_command <- sprintf(
  "scdrs compute-score --h5ad-file %s --h5ad-species human --gs-file %s --gs-species human --flag-filter-data True --flag-raw-count True --n-ctrl 1000 --flag-return-ctrl-raw-score False --flag-return-ctrl-norm-score True --out-folder %s",
  h5ad_file, gs_file,  out_dir
)
system(cli_command)

list.files("results/scdrs/output_LumB/")

# ----------------LumB_HER2Neg------------------
h5ad_file <- "results/scdrs/rna_breastcancer.h5ad"
gs_file <- "results/scdrs/LumB_HER2Neg.gs"
out_dir <- "results/scdrs/output_LumB_HER2Neg/"
dir.create(out_dir, recursive = TRUE)

cli_command <- sprintf(
  "scdrs compute-score --h5ad-file %s --h5ad-species human --gs-file %s --gs-species human --flag-filter-data True --flag-raw-count True --n-ctrl 1000 --flag-return-ctrl-raw-score False --flag-return-ctrl-norm-score True --out-folder %s",
  h5ad_file, gs_file,  out_dir
)
system(cli_command)

list.files("results/scdrs/output_LumB_HER2Neg/")

# ----------------HER2------------------
h5ad_file <- "results/scdrs/rna_breastcancer.h5ad"
gs_file <- "results/scdrs/HER2.gs"
out_dir <- "results/scdrs/output_HER2/"
dir.create(out_dir, recursive = TRUE)

cli_command <- sprintf(
  "scdrs compute-score --h5ad-file %s --h5ad-species human --gs-file %s --gs-species human --flag-filter-data True --flag-raw-count True --n-ctrl 1000 --flag-return-ctrl-raw-score False --flag-return-ctrl-norm-score True --out-folder %s",
  h5ad_file, gs_file,  out_dir
)
system(cli_command)

list.files("results/scdrs/output_HER2/")

# ---------- Check Metadata for --corr-analysis--------------
library(reticulate)
py_run_string("
import anndata
adata = anndata.read_h5ad('results/scdrs/rna_breastcancer.h5ad')
print('Obs columns:', list(adata.obs.columns))
")
# Obs columns: ['orig.ident', 'nCount_RNA', 'nFeature_RNA', 'scDblFinder_class', 'scDblFinder_score',
# 'nCount_ATAC', 'nFeature_ATAC', 'TSS.enrichment', 'TSS.percentile', 'nucleosome_signal', 'nucleosome_percentile',
# 'blacklist_fraction', 'nucleosome_group', 'nCount_SCT', 'nFeature_SCT', 'SCT_snn_res.0.06', 'seurat_clusters', 'celltype_RNA']

library(reticulate)
py_run_string("
import numpy as np
print('NumPy version:', np.__version__)
")

####################################
# Run scdrs perform-downstream    #
####################################
# pip install numpy==1.26.4 --force-reinstall
# pip install --force-reinstall scdrs anndata scanpy pandas scipy==1.14.1
# --------------TNBC----------------
h5ad_file <- "results/scdrs/rna_breastcancer.h5ad"
score_file <- "results/scdrs/output_TNBC/ZSTAT.full_score.gz"
out_dir <- "results/scdrs/output_TNBC/"

cli_downstream <- sprintf(
  "scdrs perform-downstream  --h5ad-file %s  --score-file %s  --out-folder %s  --group-analysis celltype_RNA  --gene-analysis --flag-filter-data True --flag-raw-count True",
  h5ad_file, score_file, out_dir
)
system(cli_downstream)
# ******************************************************************************
#   Call: scdrs perform-downstream \
# --h5ad-file results/scdrs/rna_breastcancer.h5ad \
# --score-file results/scdrs/output_TNBC/ZSTAT.full_score.gz \
# --out-folder results/scdrs/output_TNBC/ \
# --group-analysis celltype_RNA\
# --gene-analysis True \
# --flag-filter-data True \
# --flag-raw-count True \
# --knn-n-neighbors 15 \
# --knn-n-pcs 20
#
# --h5ad-file loaded: n_cell=16451, n_gene=18916 (sys_time=3.5s)
# First 3 cells: ['ACAGCGGGTAACGTGA-1', 'ACAGCGGGTCGGTTTG-1', 'ACAGCGGGTTCCTCAC-1']
# First 5 genes: ['RP11-34P13.7', 'AL627309.1', 'AP006222.2', 'RP4-669L17.10', 'RP5-857K21.4']
# Score file folder: results/scdrs/output_TNBC
# Find 1 score files: ZSTAT.full_score.gz
# --score-file loaded: n_trait=1 (sys_time=5.2s)
# Print info for the first trait, first 3 cells, and first 8 columns:
#   Trait                ZSTAT
# Cells                ['ACAGCGGGTAACGTGA-1', 'ACAGCGGGTCGGTTTG-1', 'ACAGCGGGTTCCTCAC-1']
# raw_score            [0.0634044  0.12215771 0.17082897]
# norm_score           [-1.2693849   0.39579168  0.38265002]
# mc_pval              [0.9080919  0.34565434 0.35164836]
# pval                 [0.9013484  0.33748809 0.34219012]
# nlog10_pval          [0.04510729 0.47174153 0.46573254]
# zscore               [-1.2892731   0.41932833  0.40649322]
# ctrl_norm_score_0    [-1.0855682  -0.86132044 -1.7804041 ]
# ctrl_norm_score_1    [-0.7236874  -1.4097202  -0.38075945]
#
# Performing scDRS group-analysis
# `connectivities` not found in `adata.obsp`; run `sc.pp.neighbors` first
# Finish group-analysis for ZSTAT (sys_time=128.7s)
#
# Performing scDRS gene-analysis
# Finish gene-analysis for ZSTAT (sys_time=129.1s)

# Results
library(data.table)
group_results_TNBC <- fread("results/scdrs/output_TNBC/ZSTAT.scdrs_group.celltype_RNA", sep = "\t", header = TRUE)
group_results_TNBC

# --------------LumA----------------
h5ad_file <- "results/scdrs/rna_breastcancer.h5ad"
score_file <- "results/scdrs/output_LumA/ZSTAT.full_score.gz"
out_dir <- "results/scdrs/output_LumA/"

cli_downstream <- sprintf(
  "scdrs perform-downstream  --h5ad-file %s  --score-file %s  --out-folder %s  --group-analysis celltype_RNA  --gene-analysis --flag-filter-data True --flag-raw-count True",
  h5ad_file, score_file, out_dir
)
system(cli_downstream)
# Results
library(data.table)
group_results_LumA <- fread("results/scdrs/output_LumA/ZSTAT.scdrs_group.celltype_RNA", sep = "\t", header = TRUE)
group_results_LumA

# --------------LumB----------------
h5ad_file <- "results/scdrs/rna_breastcancer.h5ad"
score_file <- "results/scdrs/output_LumB/ZSTAT.full_score.gz"
out_dir <- "results/scdrs/output_LumB/"

cli_downstream <- sprintf(
  "scdrs perform-downstream  --h5ad-file %s  --score-file %s  --out-folder %s  --group-analysis celltype_RNA  --gene-analysis --flag-filter-data True --flag-raw-count True",
  h5ad_file, score_file, out_dir
)
system(cli_downstream)
# Results
library(data.table)
group_results_LumB <- fread("results/scdrs/output_LumB/ZSTAT.scdrs_group.celltype_RNA", sep = "\t", header = TRUE)
group_results_LumB

# --------------LumB_HER2Neg----------------
h5ad_file <- "results/scdrs/rna_breastcancer.h5ad"
score_file <- "results/scdrs/output_LumB_HER2Neg/ZSTAT.full_score.gz"
out_dir <- "results/scdrs/output_LumB_HER2Neg/"

cli_downstream <- sprintf(
  "scdrs perform-downstream  --h5ad-file %s  --score-file %s  --out-folder %s  --group-analysis celltype_RNA  --gene-analysis --flag-filter-data True --flag-raw-count True",
  h5ad_file, score_file, out_dir
)
system(cli_downstream)
# Results
library(data.table)
group_results_LumB_HER2Neg <- fread("results/scdrs/output_LumB_HER2Neg/ZSTAT.scdrs_group.celltype_RNA", sep = "\t", header = TRUE)
group_results_LumB_HER2Neg

# --------------HER2----------------
h5ad_file <- "results/scdrs/rna_breastcancer.h5ad"
score_file <- "results/scdrs/output_HER2/ZSTAT.full_score.gz"
out_dir <- "results/scdrs/output_HER2/"

cli_downstream <- sprintf(
  "scdrs perform-downstream  --h5ad-file %s  --score-file %s  --out-folder %s  --group-analysis celltype_RNA  --gene-analysis --flag-filter-data True --flag-raw-count True",
  h5ad_file, score_file, out_dir
)
system(cli_downstream)
# Results
library(data.table)
group_results_HER2 <- fread("results/scdrs/output_HER2/ZSTAT.scdrs_group.celltype_RNA", sep = "\t", header = TRUE)
group_results_HER2

#combine the results
group_results_TNBC<-group_results_TNBC[,c(1:2,4)]
names(group_results_TNBC)[3]<-'TNBC'
group_results_LumA<-group_results_LumA[,c(1:2,4)]
names(group_results_LumA)[3]<-'LumA'
group_results_LumB<-group_results_LumB[,c(1:2,4)]
names(group_results_LumB)[3]<-'LumB'
group_results_LumB_HER2Neg<-group_results_LumB_HER2Neg[,c(1:2,4)]
names(group_results_LumB_HER2Neg)[3]<-'LumB_HER2Neg'
group_results_HER2<-group_results_HER2[,c(1:2,4)]
names(group_results_HER2)[3]<-'HER2'
all<-left_join(group_results_TNBC,group_results_LumA, by=c('group','n_cell'))
all<-left_join(all,group_results_LumB, by=c('group','n_cell'))
all<-left_join(all,group_results_LumB_HER2Neg, by=c('group','n_cell'))
all<-left_join(all,group_results_HER2, by=c('group','n_cell'))
all_scDRs_result<-all
fwrite(all_scDRs_result,'results/scdrs/scDRS_results_all_BCAC.csv')

# n_ctrl (num):Number of control gene sets (1,000, as specified by --n-ctrl 1000 in compute-score).
# Used to compute Monte Carlo p-values (mcp) and z-scores (mcz).

# assoc_mcp (num): Monte Carlo p-value for association, testing whether the scDRS scores for the cell type are significantly associated with the TNBC trait (ZSTAT).
# Range: 0 to 1 (lower = stronger association).
# Adjusted for multiple comparisons across control gene sets.

# assoc_mcz (num):Monte Carlo z-score for association, indicating the strength and direction of the association.
# Positive: Higher scDRS scores (stronger TNBC association).
# Negative: Lower scDRS scores (weaker association).
# Magnitude reflects effect size.

# hetero_mcp (num): Monte Carlo p-value for heterogeneity, testing whether scDRS scores vary significantly within the cell type (spatial or clustering effects).
# Lower values indicate greater heterogeneity (e.g., scores cluster in specific cells).

# hetero_mcz (num): Monte Carlo z-score for heterogeneity.
# Positive: Greater heterogeneity.
# Negative: Less heterogeneity.

# n_fdr_0.05 (num): Number of cells in the group with scDRS scores significant at FDR < 0.05.
# Indicates how many cells are strongly associated with TNBC.
# n_fdr_0.1 (num):Number of cells with FDR < 0.1.
# n_fdr_0.2 (num):Number of cells with FDR < 0.2.
