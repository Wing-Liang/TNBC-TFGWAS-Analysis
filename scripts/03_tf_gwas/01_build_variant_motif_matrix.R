# ==============================================================================
# 01_build_variant_motif_matrix.R
# Join BCAC GWAS variants to LD blocks and build the variant x peak x motif indicator matrix.
# Produces: Table S3; input to 02_tfgwas_model.R
# Inputs  : results/singlecell/atac_motif.rds
#           results/singlecell/final_cell_barcodes.csv
#           results/ld_blocks/ld_blocks_hg38_ALL.csv
#           results/gwas/bcac_gwas_filtered.csv
#           results/gwas/bcac_gwas_ldblock.csv
#           results/tfgwas/variant_peak_motif_subtypes_ALL.csv
#           results/singlecell/motif_matrix_allPeaks.csv
#           ...
# Outputs : results/gwas/bcac_gwas_ldblock.csv
#           results/tfgwas/variant_peak_motif_subtypes_ALL.csv
#           results/tfgwas/LDblock_summary.csv
#           results/tfgwas/motif_summary_combined.csv
# Run     : Rscript scripts/03_tf_gwas/01_build_variant_motif_matrix.R   (from the repository root)
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

#########################################
#         matched peak and motif        #
#########################################
# --- Step 1: Load atac Data ---
setwd(TNBC_ROOT)
atac <- readRDS('results/singlecell/atac_motif.rds')
# atac
Final_cell_barcode<-fread('results/singlecell/final_cell_barcodes.csv',header=TRUE)
atac <- subset(atac , cells = Final_cell_barcode$x)
DefaultAssay(atac)<-'ATAC'

# # Load GWAS data
LDblock<-fread("results/ld_blocks/ld_blocks_hg38_ALL.csv", header = TRUE)
names(LDblock)<-c('Chr','Pos','chrpos37','pos38','block')

BCAC_filt <- fread('results/gwas/bcac_gwas_filtered.csv', header = TRUE)
BCAC_filt_select <- BCAC_filt %>% select(chr.Onco, Position.Onco,Triple_Neg_P_meta, Triple_Neg_X2_meta,
                                         Luminal_A_P_meta, Luminal_A_X2_meta, Luminal_B_P_meta, Luminal_B_X2_meta,
                                          Luminal_B_HER2Neg_P_meta,Luminal_B_HER2Neg_X2_meta,
                                         HER2_P_meta, HER2_X2_meta, chrpos37)
names(BCAC_filt_select)<- c('Chr','Pos','TNBC_pval','TNBC_chisquare',
                            "LuminalA_pval" , "LuminalA_chisquare","LuminalB_P_meta","LuminalB_chisquare",
                            "LuminalB_HER2Neg_pval","LuminalB_HER2Neg_chisquare",
                            "HER2_pval", "HER2_chisquare",'chrpos37')
BCAC_filt_select$Chr <- paste0('chr',BCAC_filt_select$Chr)
BCAC_filt_select$chrpos37 <- paste0('chr',BCAC_filt_select$chrpos37)
BCAC_filt_select <- BCAC_filt_select[Chr!='chr23',]
BCAC_filt_select<-left_join(BCAC_filt_select, LDblock, by=c('Chr','Pos','chrpos37'))
table(is.na(BCAC_filt_select$block))
fwrite(BCAC_filt_select, 'results/gwas/bcac_gwas_ldblock.csv', row.names = FALSE)

BCAC_filt<-fread('results/gwas/bcac_gwas_ldblock.csv', header = TRUE)

# --- Step 3: Extract Peaks and Cell Types from Seurat Object ---
# Check available assays
names(atac@assays)

# --- Step 4: Define TFij Using Motif Information ---
# Get motif positions (binary matrix of motif presence in peaks)
motif.obj <- SeuratObject::GetAssayData(
  object = atac[['ATAC']], slot = "motifs"
)
motif_positions<-Signac::GetMotifData(object = motif.obj)
motif_names <- colnames(motif_positions)

# Map selected motifs to indices

# Convert variants and peaks to GRanges
library(regioneR)
library(IRanges)
variants_gr <- toGRanges(data.frame(chr=BCAC_filt$Chr,start=BCAC_filt$Pos,end=BCAC_filt$Pos))
peaks_gr <- granges(atac)

# --- Step 5: Create Dataset with TFij per Cell Type ---
dataset <- BCAC_filt

# Find peaks overlapping variants
peak_overlaps <- findOverlaps(variants_gr, peaks_gr)

variant_peak_map <- rep(0, nrow(dataset))
variant_peak_map[queryHits(peak_overlaps)] <- subjectHits(peak_overlaps)  # Index of overlapping peak

# Define TFij for each selected motif
for (tf_idx in seq_along(motif_names)) {
  motif_id <- motif_names[tf_idx]
  tf_name <- motif_id
  motif_presence <- motif_positions[, tf_idx]  # Binary vector for this motif
  tfij <- rep(0, nrow(dataset))

  # Assign TFij: 1 if variant overlaps a peak with the motif, 0 otherwise
  for (i in 1:nrow(dataset)) {
    peak_idx <- variant_peak_map[i]
    if (peak_idx > 0 && motif_presence[peak_idx] == 1) {
      tfij[i] <- 1
    }
  }
  dataset[[tf_name]] <- tfij
}

# Add chromatin feature (binary: 1 if variant in any peak, 0 otherwise)
chromatin_name <- "Chromatin"
dataset[[chromatin_name]] <- ifelse(variant_peak_map > 0, 1, 0)

# Save dataset
fwrite(dataset, "results/tfgwas/variant_peak_motif_subtypes_ALL.csv", row.names = FALSE)

##########--- verify the dataset from motif_matrix_df ----------###############
library(dplyr)
library(tidyr)
library(stringr)
dataset<-fread("results/tfgwas/variant_peak_motif_subtypes_ALL.csv", header = TRUE)
a<-dataset[dataset$`MA0112.3`==1,] %>% select(Chr, Pos, `MA0112.3`, starts_with("Chromatin"))
motif_matrix_df <- fread("results/singlecell/motif_matrix_allPeaks.csv", header = TRUE)
b<-motif_matrix_df %>% filter(MA0112.3==TRUE, ) %>% select(MA0112.3,peak)
head(a)
head(b)
# Step 1: Separate peak into Chr, Start, End
b_clean <- b %>%
  filter(MA0112.3 == TRUE) %>%
  separate(peak, into = c("Chr", "Start", "End"), sep = "-", convert = TRUE)

# Step 2: Inner join by Chr and check if Pos falls within Start and End
a_within_peaks <- a %>%
  inner_join(b_clean, by = "Chr") %>%
  filter(Pos >= Start & Pos <= End)

# Step 3: Check if all positions in 'a' are covered
a_ids <- a %>% mutate(id = paste(Chr, Pos, sep = "_"))
a_matched_ids <- a_within_peaks %>% mutate(id = paste(Chr, Pos, sep = "_"))

unmatched <- setdiff(a_ids$id, a_matched_ids$id)

# Output results
if (length(unmatched) == 0) {
  cat("✅ All positions in dataset 'a' are covered by peaks in dataset 'b'\n")
} else {
  cat("❌ The following positions in 'a' are not found within any peak in 'b':\n")
  print(a_ids %>% filter(id %in% unmatched))
}
# ✅ All positions in dataset 'a' are covered by peaks in dataset 'b'

##########-----------------------------------------------------###############
# LD block ranges
LDblock<-fread("results/ld_blocks/ld_blocks_hg38_ALL.csv", header = TRUE)
names(LDblock)<-c('Chr','Pos','chrpos37','pos38','block')
head(LDblock)

# Load required package
library(data.table)

# Assuming LDblock is your data.table
# If it's not already a data.table, convert it: LDblock <- as.data.table(LDblock)

# Summarize LD blocks by Chr and block, calculating min and max for Pos (Build 37) and pos38 (Build 38)
summary_LD <- LDblock[, .(
  min_pos37 = min(Pos, na.rm = TRUE),
  max_pos37 = max(Pos, na.rm = TRUE),
  min_pos38 = min(pos38, na.rm = TRUE),
  max_pos38 = max(pos38, na.rm = TRUE),
  snp_count = .N
), by = .(Chr, block)]

# Print the summarized result
head(summary_LD)

fwrite(summary_LD, "results/tfgwas/LDblock_summary.csv", row.names = FALSE)

#########################
# summarize the number of motifs in gwas and ATAC PEAKS
# Load required libraries
library(data.table)
library(dplyr)

# Read the datasets
dataset <- fread("results/tfgwas/variant_peak_motif_subtypes_ALL.csv", header = TRUE)
motif_matrix_df <- fread("results/singlecell/motif_matrix_allPeaks.csv", header = TRUE)

# Identify motif columns (assuming they start with "MA" in both datasets)
motif_cols <- grep("^MA[0-9]", names(dataset), value = TRUE)
cat("Total number of motif columns in dataset:", length(motif_cols), "\n")

# Verify motif columns in motif_matrix_df
motif_cols_matrix <- grep("^MA[0-9]", names(motif_matrix_df), value = TRUE)
cat("Total number of motif columns in motif_matrix_df:", length(motif_cols_matrix), "\n")

# Ensure motif columns match between datasets
common_motifs <- intersect(motif_cols, motif_cols_matrix)
cat("Number of common motifs between datasets:", length(common_motifs), "\n")

# Initialize a data frame to store summary results
summary_df <- data.frame(
  Motif = character(),
  Dataset_Count = integer(),
  Dataset_Unique_Chroms = integer(),
  Dataset_Chroms = character(),
  Matrix_Count = integer(),
  Matrix_Unique_Chroms = integer(),
  Matrix_Chroms = character(),
  stringsAsFactors = FALSE
)

# Summarize each motif
for (motif in common_motifs) {
  # --- Summarize dataset (motif == 1) ---
  motif_data <- dataset[dataset[[motif]] == 1, ]
  total_count_dataset <- nrow(motif_data)
  unique_chroms_dataset <- length(unique(motif_data$Chr))
  chrom_dist_dataset <- paste(names(table(motif_data$Chr)), collapse = ",")

  # --- Summarize motif_matrix_df (motif == TRUE) ---
  motif_matrix_data <- motif_matrix_df[motif_matrix_df[[motif]] == TRUE, ]
  total_count_matrix <- nrow(motif_matrix_data)
  # Extract chromosome from peak column (assuming format chrX-start-end)
  motif_matrix_data$chrom <- sub("-.*", "", motif_matrix_data$peak)
  unique_chroms_matrix <- length(unique(motif_matrix_data$chrom))
  chrom_dist_matrix <- paste(names(table(motif_matrix_data$chrom)), collapse = ",")

  # Append to summary data frame
  summary_df <- rbind(summary_df, data.frame(
    Motif = motif,
    Dataset_Count = total_count_dataset,
    Dataset_Unique_Chroms = unique_chroms_dataset,
    Dataset_Chroms = chrom_dist_dataset,
    Matrix_Count = total_count_matrix,
    Matrix_Unique_Chroms = unique_chroms_matrix,
    Matrix_Chroms = chrom_dist_matrix
  ))

  # Detailed console output for each motif
  cat("\nMotif:", motif, "\n")
  cat("--- Dataset (variant_peak_motif_subtypes_ALL.csv) ---\n")
  cat("Total peaks with", motif, "= 1:", total_count_dataset, "\n")
  cat("Unique chromosomes:", unique_chroms_dataset, "\n")
  cat("Chromosome distribution:\n")
  print(table(motif_data$Chr))

  # Summarize p-values and chi-square values for dataset
  if (total_count_dataset > 0) {
    pval_cols <- grep("_pval$", names(dataset), value = TRUE)
    chisq_cols <- grep("_chisquare$", names(dataset), value = TRUE)
    cat("Summary of p-values for peaks with", motif, "= 1:\n")
    for (pval in pval_cols) {
      cat(pval, ":\n")
      print(summary(motif_data[[pval]]))
    }
    cat("Summary of chi-square values for peaks with", motif, "= 1:\n")
    for (chisq in chisq_cols) {
      cat(chisq, ":\n")
      print(summary(motif_data[[chisq]]))
    }
  } else {
    cat("No peaks with", motif, "= 1 in dataset, skipping p-value and chi-square summaries.\n")
  }

  cat("--- Motif Matrix (motif_matrix_allPeaks.csv) ---\n")
  cat("Total peaks with", motif, "= TRUE:", total_count_matrix, "\n")
  cat("Unique chromosomes:", unique_chroms_matrix, "\n")
  cat("Chromosome distribution:\n")
  print(table(motif_matrix_data$chrom))
}

motif_id_list_final<-fread('results/singlecell/jaspar2022_motifs_727.csv', header = TRUE)
names(motif_id_list_final)[1]<-'Motif'
summary_df<-left_join(motif_id_list_final, summary_df, by='Motif')

# Save summary to a CSV file
fwrite(summary_df, "results/tfgwas/motif_summary_combined.csv")
