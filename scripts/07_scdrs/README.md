# scDRS

Single-cell disease relevance scores (scDRS 1.0.3) for breast cancer subtypes, using MAGMA gene-level GWAS statistics and the discovery scRNA-seq data.

Run the scripts in the order listed, from the repository root (see the main README).

| Script | What it does | Produces |
|---|---|---|
| `01_prepare_gwas_inputs.R` | Per-subtype BCAC GWAS tables for MAGMA. | Fig S10 |
| `02_scdrs.R` | MAGMA gene-level GWAS statistics and scDRS cell-level disease scores per subtype. | Fig S10 |
| `03_scdrs_figure.R` | scDRS cell-type enrichment dot plots. | Fig S10 |
