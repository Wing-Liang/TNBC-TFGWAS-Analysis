# motifbreakR and GTEx

Predicts allele-specific TF-binding changes for genome-wide significant TNBC promoter variants (motifbreakR), tests colocalisation of TNBC GWAS and GTEx v8 breast cis-eQTL signals, and plots GTEx v11 genotype-expression associations.

Run the scripts in the order listed, from the repository root (see the main README).

| Script | What it does | Produces |
|---|---|---|
| `01_annotate_variant_regions.R` | Genomic-region annotation of the GWAS variants (ChIPseeker). | Fig 5A, Table S12 |
| `02_motifbreakR_setup.R` | Select genome-wide significant TNBC variants for motifbreakR. | input to 03_motifbreakR_tnbc.R |
| `03_motifbreakR_tnbc.R` | Predict allele-specific TF-binding changes (motifbreakR, JASPAR 2022) for TNBC variants. | Table S12, Fig 5A |
| `04_coloc_gtex_v8.R` | Colocalisation of TNBC GWAS with GTEx v8 breast cis-eQTL for promoter genes. | Table S13 |
| `05_coloc_harmonization_audit.R` | Audit allele harmonisation between the GWAS and eQTL data. | Table S13 (check) |
| `06_coloc_gtex_v8_harmonized.R` | Colocalisation with strand-harmonised alleles (reported analysis). | Table S13 |
| `07_locus_plots.R` | LocusCompare and regional association plots per locus. | Fig 5C |
| `08_gtex_v11_genotype_violins.R` | Gene expression by genotype, GTEx v11 breast tissue. | Fig 5D-H |
| `09_motifbreakR_motif_plots.R` | motifbreakR sequence-logo plots for promoter variants. | Fig 5D-H |
| `10_mdm4_rps18_celltype_dotplot.R` | MDM4 and RPS18 expression by cell type, discovery cohort. | Fig S11 |
