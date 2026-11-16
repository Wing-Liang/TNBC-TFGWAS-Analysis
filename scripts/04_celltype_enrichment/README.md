# Subtype-specific motifs and cell-type enrichment

Identifies subtype-specific motifs and scores their accessibility by cell type (AddModuleScore on chromVAR deviations), with Kruskal-Wallis / Dunn tests and enrichment in cell-type differentially accessible peaks.

Run the scripts in the order listed, from the repository root (see the main README).

| Script | What it does | Produces |
|---|---|---|
| `01_differential_motif_activity.R` | Differential chromVAR motif activity between cell types. | input to Fig 3F |
| `02_motif_enrichment_by_subtype.R` | Motif enrichment in cell-type differentially accessible peaks and module scores, per subtype. | Figs S3, S4, S7, S8 |
| `03_subtype_specific_motifs.R` | Subtype-specific motifs: module scores, Kruskal-Wallis / Dunn tests and enrichment heatmaps. | Fig 3E-F, Table S5 |
| `04_classic_tf_chromvar.R` | chromVAR scores of classic mammary transcription factors across epithelial cell types. | Fig S2 |
| `05_motif_sets_upset_heatmap.R` | Motif sets per subtype, cell-type enrichment and module scores. | Fig 3A-B, Figs S5-S6 |
| `06_module_score_umap.R` | UMAP of the TNBC-specific motif module score. | Fig 3D |
| `07_module_score_umap_plot.R` | Draw the TNBC-specific module-score UMAP. | Fig 3D |
