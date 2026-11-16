# TF-GWAS model

Builds the variant x peak x motif matrix and fits, for each of 727 JASPAR 2022 motifs, `chi-square ~ motif-in-peak + (1 | LD block)` for overall breast cancer and five intrinsic subtypes.

Run the scripts in the order listed, from the repository root (see the main README).

| Script | What it does | Produces |
|---|---|---|
| `01_build_variant_motif_matrix.R` | Join BCAC GWAS variants to LD blocks and build the variant x peak x motif indicator matrix. | Table S3; input to 02_tfgwas_model.R |
| `02_tfgwas_model.R` | TF-GWAS linear mixed model per motif and subtype. | Table S4, Fig 3A |
| `03_build_overall_bc_gwas.R` | Join the BCAC overall breast-cancer GWAS to the variant table (GRCh38). | Table S4 (overall BC) |
| `04_tfgwas_overall_bc.R` | TF-GWAS model for overall breast cancer, discovery cohort. | Table S4 (overall BC) |
| `05_manhattan_plot.R` | TNBC Manhattan plot highlighting variants in peaks with a TNBC-specific motif. | Fig 3C |
| `manhattan_functions.R` | Plotting functions sourced by 05_manhattan_plot.R. | helper |
