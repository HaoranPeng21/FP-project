#!/usr/bin/env Rscript
# Explicit setup command; analysis scripts never install packages or download data.
cran <- c("ape", "circlize", "cowplot", "data.table", "dplyr", "ggnewscale", "ggplot2",
  "ggrepel", "ggsci", "jsonlite", "patchwork", "readr", "scales", "tibble", "tidyr", "vegan")
missing <- cran[!vapply(cran, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) install.packages(missing, repos = "https://cloud.r-project.org")
bioc <- c("ComplexHeatmap", "ggtree")
missing <- bioc[!vapply(bioc, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) {
  if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager", repos = "https://cloud.r-project.org")
  BiocManager::install(missing, ask = FALSE, update = FALSE)
}
