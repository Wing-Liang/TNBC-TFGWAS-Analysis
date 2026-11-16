# GWAS filtering and LD blocks

Filters the BCAC intrinsic-subtype GWAS summary statistics (Zhang et al. 2020) and computes per-subtype chi-squares; variants (GRCh37) are then lifted to GRCh38 and assigned to the 1,336 European pyrho LD blocks, with variants that fail liftover or fall between blocks taking the block of the nearest assigned variant.

Run the scripts in the order listed, from the repository root (see the main README).

| Script | What it does | Produces |
|---|---|---|
| `01_filter_gwas.R` | Filter the BCAC subtype GWAS (control allele frequency 0.01-0.99) and compute per-subtype chi-squares. | input to Tables S2-S4 |
| `02_assign_ld_blocks.R` | Lift GWAS variants to GRCh38 and assign each to a European LD block. | Table S2 |
