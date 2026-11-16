## 03_build_aracne_regulon.R
## ARACNe-BRCA regulon (aracne.networks::regulonbrca) with regulators and targets converted from
## Entrez IDs to the gene symbols of the TCGA-BRCA expression data.
##
## Each Entrez ID is mapped to its org.Hs.eg.db symbol when that symbol is present in the TCGA
## data; otherwise through its Ensembl ID to the TCGA gene name. Targets that map to the same
## symbol are de-duplicated (first kept). Unlike a regulator-only mapping, targets that are not
## themselves regulators are retained.
##
## Inputs   results/tcga/TCGA_BRCA_data_all.rds                    (01_download_tcga_brca.R)
## Outputs  results/tf_activity/viper/aracne_brca_regulon_V2_symbol.rds   (used by 04_viper_tf_activity.R)
## Run from the repository root:  Rscript scripts/08_tcga_viper_depmap/03_build_aracne_regulon.R

TNBC_ROOT <- Sys.getenv("TNBC_ROOT", ".")

suppressMessages({
  library(SummarizedExperiment); library(viper); library(aracne.networks); library(org.Hs.eg.db)
})

setwd(TNBC_ROOT)
OUT <- "results/tf_activity/viper"
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

se_full <- readRDS("results/tcga/TCGA_BRCA_data_all.rds")

data(regulonbrca, package = "aracne.networks")
all_entrez <- unique(c(names(regulonbrca), unlist(lapply(regulonbrca, function(r) names(r$tfmode)))))
e2a <- AnnotationDbi::select(org.Hs.eg.db, keys = all_entrez, keytype = "ENTREZID", columns = c("ENSEMBL", "SYMBOL"))
e2a <- e2a[!duplicated(e2a$ENTREZID), ]
map_ens    <- setNames(e2a$ENSEMBL, e2a$ENTREZID)
map_native <- setNames(toupper(e2a$SYMBOL), e2a$ENTREZID)

universe <- toupper(unique(rowData(se_full)$gene_name))
ens2sym  <- setNames(toupper(rowData(se_full)$gene_name), sub("[.].*$", "", rowData(se_full)$gene_id))
ens2sym  <- ens2sym[!duplicated(names(ens2sym))]
e2s <- function(ids) {
  nat <- map_native[ids]; ok <- !is.na(nat) & nat %in% universe
  ifelse(ok, nat, unname(ens2sym[map_ens[ids]]))
}

rs <- e2s(names(regulonbrca)); keep <- !is.na(rs)
reg <- regulonbrca[keep]; names(reg) <- rs[keep]
if (any(duplicated(names(reg))))
  cat("WARNING residual regulator-name collisions:", paste(unique(names(reg)[duplicated(names(reg))]), collapse = ", "), "\n")
reg <- lapply(reg, function(r) {
  tg <- e2s(names(r$tfmode)); k <- !is.na(tg) & !duplicated(tg)
  list(tfmode = setNames(r$tfmode[k], tg[k]), likelihood = r$likelihood[k])
})
class(reg) <- "regulon"

saveRDS(reg, file.path(OUT, "aracne_brca_regulon_V2_symbol.rds"))
cat("ARACNe-BRCA regulon:", length(reg), "regulators\n")
