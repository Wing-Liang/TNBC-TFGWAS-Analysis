## 03_scdrs_figure.R
## Cell-type enrichment of breast-cancer subtype GWAS signal by scDRS (Fig S10).
##
## Inputs   results/scdrs/scDRS_results_all_BCAC.csv     (02_scdrs.R; one row per cell type,
##                                                       one column of P values per subtype)
## Outputs  results/scdrs/dot_plot_scDRS_<subtype>.png
## Run from the repository root:  Rscript scripts/07_scdrs/03_scdrs_figure.R

TNBC_ROOT <- Sys.getenv("TNBC_ROOT", ".")

library(data.table)
library(ggplot2)

setwd(TNBC_ROOT)

all_scDRs_result <- fread("results/scdrs/scDRS_results_all_BCAC.csv", header = TRUE)
names(all_scDRs_result)[1] <- "Cell_type"
all_scDRs_result$Cell_type <- c("B", "Myl", "Pc", "T", "Bs", "Endo", "Fb", "LHr", "LS-1", "LS-2")

celltype_order <- c("LS-1", "LS-2", "LHr", "Bs", "Fb", "Endo", "T", "B", "Myl", "Pc")
all_scDRs_result <- all_scDRs_result[match(celltype_order, all_scDRs_result$Cell_type), ]

create_dot_plot <- function(subtype) {
  plot_data <- data.frame(
    Cell_type = celltype_order,
    value     = -log10(all_scDRs_result[[subtype]])
  )
  p <- ggplot(plot_data, aes(x = Cell_type, y = value)) +
    geom_point(size = 3, color = "#d47264") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
    theme_classic() +
    theme(
      axis.title = element_text(size = 10),
      legend.position = "bottom",
      plot.title = element_text(hjust = 0.5, size = 10)
    ) +
    labs(x = "Cell Type", y = "-log10(P value)") +
    scale_x_discrete(limits = celltype_order) +
    scale_y_continuous(limits = c(0, 8))
  ggsave(file.path("results/scdrs", paste0("dot_plot_scDRS_", subtype, ".png")), p,
         width = 5, height = 4, dpi = 150)
  p
}

for (subtype in c("TNBC", "LumA", "LumB", "LumB_HER2Neg", "HER2")) create_dot_plot(subtype)
