# ==============================================================================
# 01_differential_motif_activity.R
# Differential chromVAR motif activity between cell types.
# Produces: input to Fig 3F
# Inputs  : results/singlecell/final_cell_barcodes.csv
#           results/singlecell/atac_motif.rds
#           results/singlecell/jaspar2022_motifs_727.csv
# Outputs : results/singlecell/Motif_DiffAct_LHr_LS_Bs.csv
#           results/singlecell/Motif_DiffAct_LHr_LS_Bs_unique.csv
#           results/singlecell/Motif_DiffAct_other_celltype.csv
#           results/singlecell/Motif_DiffAct_other_celltype_unique.csv
# Run     : Rscript scripts/04_celltype_enrichment/01_differential_motif_activity.R   (from the repository root)
# ==============================================================================
## ---- paths (set the environment variables; see README) ----
TNBC_ROOT <- Sys.getenv("TNBC_ROOT", ".")
# Load required packages
library(lme4)      # For generalized mixed models
library(dplyr)     # For data manipulation
library(ggplot2)   # For visualization (optional)

# Load required packages
library(Seurat)        # For scATAC-seq object handling
library(Signac)        # For scATAC-seq utilities
library(dplyr)         # For data manipulation
library(GenomicRanges) # For genomic overlaps
library(motifmatchr)   # For motif matching
library(TFBSTools)     # For motif handling
library(data.table)

setwd(TNBC_ROOT)
#############################################
#         heatmap of motif enrichment       #
##############################################
Final_cell_barcode<-fread('results/singlecell/final_cell_barcodes.csv',header=TRUE)

#ATAC
atac_subset_qc<-readRDS('results/singlecell/atac_motif.rds')
atac_subset_qc <- subset(atac_subset_qc , cells = Final_cell_barcode$x)
DefaultAssay(atac_subset_qc)<-'chromvar'

# Find differential markers within LHr, LS, Bs
atac_LHr_LS_Bs <- subset(atac_subset_qc, idents=c('LHr','LS-1','LS-2','Bs'))
DAMotif_LHr_LS_Bs <- FindAllMarkers(
  object = atac_LHr_LS_Bs,
  only.pos = TRUE,
  mean.fxn = rowMeans,
  test.use = "LR",
  fc.name = "avg_diff",
  latent.vars = "nCount_ATAC"
)

motif_id_list_final<-fread('results/singlecell/jaspar2022_motifs_727.csv',header=TRUE)

subset.DAMotif_LHr_LS_Bs <- DAMotif_LHr_LS_Bs[(DAMotif_LHr_LS_Bs$p_val_adj < 0.05 &
                                                         DAMotif_LHr_LS_Bs$avg_diff > 0), ]

subset.DAMotif_LHr_LS_Bs <- subset.DAMotif_LHr_LS_Bs %>% rename('gene'='motif')
subset.DAMotif_LHr_LS_Bs <- subset.DAMotif_LHr_LS_Bs[subset.DAMotif_LHr_LS_Bs$motif %in% motif_id_list_final$motif,]
subset.DAMotif_LHr_LS_Bs  <- left_join(subset.DAMotif_LHr_LS_Bs ,motif_id_list_final,
                                            by=c('motif'))
head(subset.DAMotif_LHr_LS_Bs)

write.csv(subset.DAMotif_LHr_LS_Bs, 'results/singlecell/Motif_DiffAct_LHr_LS_Bs.csv',row.names = FALSE)

library(dplyr)

# Count how many clusters each motif appears in
motif_cluster_count <- subset.DAMotif_LHr_LS_Bs %>%
  group_by(motif) %>%
  summarize(cluster_count = n_distinct(cluster), .groups = 'drop')

# Filter motifs that appear in only one cluster
unique_motifs <- motif_cluster_count %>%
  filter(cluster_count == 1)

# Join back to original data to get full rows for those motifs
subset.DAMotif_unique <- subset.DAMotif_LHr_LS_Bs %>%
  filter(motif %in% unique_motifs$motif)

subset.DAMotif_unique
write.csv(subset.DAMotif_unique, 'results/singlecell/Motif_DiffAct_LHr_LS_Bs_unique.csv',row.names = FALSE)

##### Other cell type

# Find differential markers within LHr, LS, Bs
atac_other_celltype <- subset(atac_subset_qc, idents=c('Fb','Endo','T','B','Myl','Pc'))
DAMotif_other_celltype <- FindAllMarkers(
  object = atac_other_celltype,
  only.pos = TRUE,
  mean.fxn = rowMeans,
  test.use = "LR",
  fc.name = "avg_diff",
  latent.vars = "nCount_ATAC"
)

motif_id_list_final<-fread('results/singlecell/jaspar2022_motifs_727.csv',header=TRUE)

subset.DAMotif_other_celltype <- DAMotif_other_celltype[(DAMotif_other_celltype$p_val_adj < 0.05 &
                                                 DAMotif_other_celltype$avg_diff > 0), ]

subset.DAMotif_other_celltype <- subset.DAMotif_other_celltype %>% rename('gene'='motif')
subset.DAMotif_other_celltype <- subset.DAMotif_other_celltype[subset.DAMotif_other_celltype$motif %in% motif_id_list_final$motif,]
subset.DAMotif_other_celltype  <- left_join(subset.DAMotif_other_celltype ,motif_id_list_final,
                                       by=c('motif'))
head(subset.DAMotif_other_celltype)

write.csv(subset.DAMotif_other_celltype, 'results/singlecell/Motif_DiffAct_other_celltype.csv',row.names = FALSE)

# Count how many clusters each motif appears in
motif_cluster_count <- subset.DAMotif_other_celltype %>%
  group_by(motif) %>%
  summarize(cluster_count = n_distinct(cluster), .groups = 'drop')

# Filter motifs that appear in only one cluster
unique_motifs <- motif_cluster_count %>%
  filter(cluster_count == 1)

# Join back to original data to get full rows for those motifs
subset.DAMotif_unique <- subset.DAMotif_other_celltype %>%
  filter(motif %in% unique_motifs$motif)

head(subset.DAMotif_unique)
write.csv(subset.DAMotif_unique, 'results/singlecell/Motif_DiffAct_other_celltype_unique.csv',row.names = FALSE)
