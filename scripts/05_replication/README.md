# Replication cohort

Repeats the TF-GWAS in an independent 10x Multiome cohort of normal breast tissue from 15 donors (GEO GSE244585, GRCh38): data download, QC, peaks, motif matrix, model fitting, comparison with discovery, and grouping of motifs by TF family. These scripts run inside a `replication/` working folder (see the main README).

Run the scripts in the order listed, from the repository root (see the main README).

| Script | What it does | Produces |
|---|---|---|
| `00_check_h5_info.py` | List the contents of a GEO HDF5 file (replication cohort). | data check |
| `00_check_h5_ok.py` | Check that a GEO HDF5 file opens (replication cohort). | data check |
| `00_describe_data.R` | Describe the replication data: cells, fragments and pools. | data check |
| `01_download_geo_data.sh` | Download GSE244585 pools 9, 11 and 12 (European ancestry) from GEO. | replication input |
| `02_recover_pool11_fragments.sh` | Recover the truncated pool-11 fragment file deposited on GEO. | replication input |
| `03a_extract_cellxgene_data.py` | Extract cell barcodes, author cell types and raw RNA counts from the CELLxGENE object. | replication input |
| `03b_extract_cellxgene_embeddings.py` | Extract embeddings and cell metadata from the CELLxGENE object. | replication input |
| `03c_combine_fragments.sbatch` | Combine the three pools' ATAC fragments (CELLxGENE cells, autosomes). | replication input |
| `03d_call_peaks_macs2.sbatch` | MACS2 peak calling on the three pools. | replication peaks |
| `04_build_combined_object.R` | Build the replication Seurat object (RNA and ATAC peaks). | replication object |
| `05a_qc.R` | Replication quality control (median +/- 3 MAD) and doublet scoring. | replication cell set |
| `05b_annotate_rna.R` | RNA embedding, author cell types and marker genes. | replication object |
| `05c_annotate_atac.R` | ATAC gene activity and ATAC figures. | replication object |
| `05d_wnn.R` | Import the published weighted-nearest-neighbour embedding. | replication object |
| `06a_liftover_gwas_hg38.R` | Fill missing GRCh38 coordinates in the BCAC variant table. | input to 06b_motif_matrix.R |
| `06b_motif_matrix.R` | JASPAR 2022 motif annotation and the replication variant x peak x motif matrix. | input to 07b_tfgwas_model.R |
| `07a_split_motif_batches.R` | Split the replication variant x motif matrix into motif batches for the model array. | Table S7 |
| `07b_tfgwas_model.R` | Replication TF-GWAS model, one batch of motifs per array task. | Table S7 |
| `07b_tfgwas_model.sbatch` | SLURM array wrapper for 07b_tfgwas_model.R. | Table S7 |
| `08_tfgwas_combine.R` | Combine replication batches and compare with the discovery cohort. | Table S7 |
| `08_tfgwas_combine.sbatch` | SLURM wrapper for 08_tfgwas_combine.R. | Table S7 |
| `09a_overall_bc_annotation.R` | Attach the overall-BC chi-square to the replication motif batches. | Table S7 |
| `09b_split_discovery_overall_bc.R` | Split the discovery matrix into motif batches with the overall-BC chi-square. | Table S4 |
| `09c_tfgwas_overall_bc.R` | TF-GWAS model for overall breast cancer (replication; discovery with OVERALLBC_COHORT=discovery). | Tables S4 and S7 (overall BC) |
| `09c_tfgwas_overall_bc.sbatch` | SLURM array wrapper for 09c_tfgwas_overall_bc.R, replication cohort. | Table S7 (overall BC) |
| `09d_tfgwas_overall_bc_discovery.sbatch` | SLURM array wrapper for 09c_tfgwas_overall_bc.R, discovery cohort. | Table S4 (overall BC) |
| `09e_combine_overall_bc.R` | Combine the overall-BC batch results per cohort and add BH-FDR. | Tables S4, S7 (overall BC) |
| `10_tf_family_annotation.R` | Annotate motifs with TFClass transcription-factor families. | reference/group_TFfamily.csv |
| `11_peaks_per_motif.R` | Count peaks carrying each motif in both cohorts. | Table S8 (input) |
| `12_tf_family_models.R` | TF-GWAS at motif, TF-family and JASPAR-cluster level. | Table S8 |
| `12_tf_family_models.sbatch` | SLURM wrapper for 12_tf_family_models.R. | Table S8 |
| `13_motif_master_table.py` | Per-motif table: peak and variant counts, significance and subtype-uniqueness flags, both cohorts. | Tables S7-S8 |
| `14_cohort_matrix_check.R` | Structural check of the replication matrix against discovery. | data check |
| `15_summary_tables.py` | Summary counts of significant and subtype-unique units per level. | Table S8 |
| `16_cross_cohort_overlap.py` | Overlap of significant units between cohorts. | Table S8 |
| `17_list_shared_groups.py` | List TF families significant in both cohorts. | Table S8 |
| `18_list_unique_motifs.py` | List subtype-unique motifs in both cohorts. | Table S8 |
| `19_marker_detection_by_cohort.R` | Canonical marker detection by cell type in both cohorts. | Table S6 |
