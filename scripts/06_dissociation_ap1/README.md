# Dissociation stress and AP-1

Scores cells for a dissociation-stress signature and an AP-1 gene set, and compares AP-1 motif activity and cell-type composition between the two cohorts.

Run the scripts in the order listed, from the repository root (see the main README).

| Script | What it does | Produces |
|---|---|---|
| `01_geneset_scores.R` | Cell-cycle and dissociation-stress gene-set scores, both cohorts. | Fig S9 |
| `02_dissociation_ap1_scores.R` | Score cells for the AP-1 subset of the dissociation-stress signature. | Fig S9 |
| `03_ap1_score_plots.R` | Violin and scatter plots of the AP-1 score. | Fig S9 |
| `04_geneset_score_umaps.R` | UMAPs of dissociation-stress and AP-1 scores, both cohorts. | Fig S9 |
| `05_geneset_score_violins.R` | Violins of dissociation-stress scores, both cohorts. | Fig S9 |
| `06_replication_atac_umap.R` | ATAC UMAP of the replication cohort. | Fig S9 |
| `07_ap1_motif_activity_umap.R` | AP-1 motif activity on the UMAP, both cohorts. | Fig S9 |
| `08_celltype_composition.R` | Cell-type composition of the two cohorts. | Fig S9 |
| `09_celltype_composition_plot.py` | Plot of cell-type composition, both cohorts. | Fig S9 |
