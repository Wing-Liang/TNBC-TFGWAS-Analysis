# Figures and tables → scripts

Folder names refer to `scripts/`. Within a folder, run scripts in numeric order (see each folder's README).
Steps that only prepare inputs (downloads, quality control, motif annotation, GWAS filtering) are
listed in the folder READMEs; the tables below name the scripts that draw each panel or write each table.

## Main figures

| Figure | Content | Folder | Script(s) |
|---|---|---|---|
| Fig 1 | Study workflow (schematic) | — | drawn by hand |
| Fig 2A–B, D | scRNA-seq UMAP, marker feature plots, cell counts | `01_singlecell_annotation` | `01_import_geo_data.R`, `08_marker_plots.R`, `09_cell_counts.R` |
| Fig 2C, E | scATAC-seq UMAP, peak annotation | `01_singlecell_annotation` | `02_call_peaks_macs2.R`, `03_qc_doublets.R`, `04_annotate_cell_types.R`, `06_wnn_integration.R`, `05_atac_umap_figures.R` |
| Fig 2F | Top-5 markers, RNA and gene activity | `01_singlecell_annotation` | `07_markers_and_da_peaks.R` |
| Fig 3A–B | UpSet of subtype motif sets; set × cell-type AddModuleScore heatmap | `04_celltype_enrichment` | `05_motif_sets_upset_heatmap.R` |
| Fig 3C | TNBC Manhattan plot, motif-overlapping SNPs | `03_tf_gwas` | `05_manhattan_plot.R` (sources `manhattan_functions.R`) |
| Fig 3D | TNBC-specific module-score UMAP | `04_celltype_enrichment` | `06_module_score_umap.R`, `07_module_score_umap_plot.R` |
| Fig 3E | Kruskal–Wallis / Dunn boxplots | `04_celltype_enrichment` | `03_subtype_specific_motifs.R` |
| Fig 3F | 57 TNBC-specific motifs in cell-type differentially accessible peaks | `04_celltype_enrichment` | `01_differential_motif_activity.R`, `03_subtype_specific_motifs.R` |
| Fig 4A–B | TF gene expression (DESeq2) and inferred activity (msVIPER), TCGA-BRCA | `08_tcga_viper_depmap` | `01_download_tcga_brca.R`, `02_tcga_deseq2_tf_genes.R`, `03_build_aracne_regulon.R`, `04_viper_tf_activity.R` (run by `04_run_viper.sh`), `05_expression_activity_panel.R` |
| Fig 4C | DepMap 25Q3 gene effects, 61 TF genes | `08_tcga_viper_depmap` | `07_depmap_tf_genes.R` |
| Fig 5A | motifbreakR summary, MDM4 / RPS18 promoter variants | `09_motifbreakR_gtex_coloc` | `01_annotate_variant_regions.R`, `02_motifbreakR_setup.R`, `03_motifbreakR_tnbc.R` |
| Fig 5B | Evidence schematic at MDM4 | — | drawn by hand |
| Fig 5C | LocusCompare and regional plot, MDM4 | `09_motifbreakR_gtex_coloc` | `07_locus_plots.R` |
| Fig 5D–H | GTEx v11 genotype–expression violins | `09_motifbreakR_gtex_coloc` | `08_gtex_v11_genotype_violins.R` |
| Fig 5D–H | motifbreakR motif-disruption plots | `09_motifbreakR_gtex_coloc` | `09_motifbreakR_motif_plots.R` |

## Supplementary figures

| Figure | Content | Folder | Script(s) |
|---|---|---|---|
| Fig S1 | Sample selection and filtering flowchart | — | drawn by hand |
| Fig S2 | Classic TF chromVAR scores across epithelial cells | `04_celltype_enrichment` | `04_classic_tf_chromvar.R` |
| Fig S3, S4, S7, S8 | Luminal / HER2 motifs; module-score boxplots; luminal-overlap motifs | `04_celltype_enrichment` | `02_motif_enrichment_by_subtype.R`, `03_subtype_specific_motifs.R` |
| Fig S5 | Unique motifs by TF family | `04_celltype_enrichment` | `05_motif_sets_upset_heatmap.R` |
| Fig S6 | Bubble plots by motif family | `04_celltype_enrichment` | `05_motif_sets_upset_heatmap.R` |
| Fig S9 | Dissociation-stress and AP-1 scores and AP-1 motif activity, both cohorts | `06_dissociation_ap1` | `01_geneset_scores.R`, `02_dissociation_ap1_scores.R`, `03_geneset_score_umaps.R`, `04_replication_atac_umap.R`, `05_ap1_motif_activity_umap.R`, `06_celltype_composition.R`, `07_celltype_composition_plot.py` |
| Fig S10 | scDRS | `07_scdrs` | `01_prepare_gwas_inputs.R`, `02_scdrs.R`, `03_scdrs_figure.R` |
| Fig S11 | MDM4 / RPS18 expression across cell types | `09_motifbreakR_gtex_coloc` | `10_mdm4_rps18_celltype_dotplot.R` |

## Supplementary tables

| Table | Content | Folder | Script(s) |
|---|---|---|---|
| S1 | Peak annotation by cell type | `01_singlecell_annotation` | `07_markers_and_da_peaks.R` |
| S2 | LD blocks (1,336, pyrho EUR, GRCh38) | `02_gwas_ld_blocks` | `01_filter_gwas.R`, `02_assign_ld_blocks.R` |
| S3 | SNP × open-chromatin counts per motif | `03_tf_gwas` | `01_build_variant_motif_matrix.R` |
| S4 | TF-GWAS results, discovery, subtypes and overall BC | `03_tf_gwas` | `02_tfgwas_model.R`, `03_build_overall_bc_gwas.R`, `04_tfgwas_overall_bc.R` |
| S5 | Dunn's test of module scores | `04_celltype_enrichment` | `03_subtype_specific_motifs.R` |
| S6 | Canonical marker detection, both cohorts | `05_replication` | `22_marker_detection_by_cohort.R` |
| S7 | TF-GWAS results, replication | `05_replication` | `06b_motif_matrix.R`, `07a_split_motif_batches.R`, `07b_tfgwas_model.R`, `08_tfgwas_combine.R`, `09a`–`09e` (overall BC), `14_motif_master_table.py` |
| S8 | Replication by TF family | `05_replication` | `13_tf_family_models.R`, `11_tf_family_annotation.R`, `15`–`22` |
| S9 | DESeq2, 61 TF genes, TCGA-BRCA | `08_tcga_viper_depmap` | `02_tcga_deseq2_tf_genes.R` |
| S10 | TF activity (msVIPER), 61 TF genes | `08_tcga_viper_depmap` | `04_viper_tf_activity.R`, `04_run_viper.sh`, `06_table_tf_activity.py` |
| S11 | DepMap gene effects, TF genes | `08_tcga_viper_depmap` | `07_depmap_tf_genes.R` |
| S12 | motifbreakR results | `09_motifbreakR_gtex_coloc` | `03_motifbreakR_tnbc.R` |
| S13 | Colocalisation, GTEx v8 breast (p12 = 5×10⁻⁵, ±300 kb) | `09_motifbreakR_gtex_coloc` | `04_coloc_gtex_v8.R`, `05_coloc_harmonization_audit.R`, `06_coloc_gtex_v8_harmonized.R` |
| S14 | MDM4 / RPS18 differential expression, TCGA-BRCA | `08_tcga_viper_depmap` | `08_tcga_mdm4_rps18.R` |
