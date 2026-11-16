# TCGA-BRCA expression, TF activity and DepMap

Differential expression (DESeq2) of the 61 TNBC-specific TF genes across PAM50 subtypes and adjacent normal tissue in TCGA-BRCA, inferred TF activity (msVIPER, ARACNe-BRCA regulons), and CRISPR gene effects in TNBC cell lines (DepMap Public 25Q3).

Run the scripts in the order listed, from the repository root (see the main README).

| Script | What it does | Produces |
|---|---|---|
| `01_download_tcga_brca.R` | Download TCGA-BRCA RNA-seq from the GDC and prepare PAM50 metadata. | Fig 4, Tables S9, S10, S14 |
| `02_tcga_deseq2_tf_genes.R` | TCGA-BRCA differential expression (DESeq2) of the 61 TF genes. | Fig 4A, Table S9 |
| `03_build_aracne_regulon.R` | ARACNe-BRCA regulon with gene-symbol identifiers. | Fig 4A-B |
| `04_run_viper.sh` | Run 04_viper_tf_activity.R with the reported settings. | Fig 4A-B, Table S10 |
| `04_viper_tf_activity.R` | msVIPER activity of the 61 TFs (ARACNe-BRCA regulons), TCGA-BRCA. | Fig 4A-B, Table S10 |
| `05_expression_activity_panel.R` | TF gene expression versus inferred TF activity. | Fig 4A-B |
| `06_table_tf_activity.py` | Format the TF-activity table. | Table S10 |
| `07_depmap_tf_genes.R` | CRISPR gene effects of the 61 TF genes in TNBC cell lines (DepMap 25Q3). | Fig 4C, Table S11 |
| `08_tcga_mdm4_rps18.R` | MDM4 and RPS18 expression across PAM50 subtypes, TCGA-BRCA. | Table S14 |
