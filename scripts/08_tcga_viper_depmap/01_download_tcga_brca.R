## 01_download_tcga_brca.R
## Download TCGA-BRCA RNA-seq (STAR - Counts) from the Genomic Data Commons and prepare the sample
## metadata with PAM50 subtypes.
##
## The analysis used GDC Data Release 44.0 (October 29, 2025); a later download may differ.
## Metadata are restricted to female patients, and adjacent normal tissue ("Solid Tissue Normal")
## is labelled NAT in the PAM50 column. Before restricting to white patients (done in the
## downstream scripts), primary tumours number Basal 197, Her2 82, LumA 571, LumB 209, Normal 40.
##
## Outputs  results/tcga/TCGA_BRCA_data_all.rds          SummarizedExperiment, 60,660 genes x 1,231 samples
##          results/tcga/TCGA_BRCA_PAM50_metadata.rds     colData of female samples, PAM50 + NAT
## Run from the repository root:  Rscript scripts/08_tcga_viper_depmap/01_download_tcga_brca.R

TNBC_ROOT <- Sys.getenv("TNBC_ROOT", ".")

library(TCGAbiolinks)
library(SummarizedExperiment)
library(dplyr)

setwd(TNBC_ROOT)
dir.create("results/tcga", recursive = TRUE, showWarnings = FALSE)

query <- GDCquery(
  project = "TCGA-BRCA",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts"
)
GDCdownload(query)
data <- GDCprepare(query)
metadata(data)$data_release
saveRDS(data, file = "results/tcga/TCGA_BRCA_data_all.rds")

metadata <- colData(data) %>% as.data.frame()
table(metadata$sample_type, metadata$paper_BRCA_Subtype_PAM50)
metadata <- metadata %>% filter(gender %in% "female")
metadata$paper_BRCA_Subtype_PAM50 <- ifelse(metadata$sample_type == "Solid Tissue Normal", "NAT",
                                            metadata$paper_BRCA_Subtype_PAM50)
saveRDS(metadata, file = "results/tcga/TCGA_BRCA_PAM50_metadata.rds")
