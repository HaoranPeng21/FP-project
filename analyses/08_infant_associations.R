# Figure 3D. Selected coefficients from the prespecified adjusted longitudinal models.
# Keep the upstream full-family q-values; do not correct only these displayed cells.

# 1. Set these two private directories (or edit the paths directly).
data_dir <- Sys.getenv("FP_DATA_DIR", "/path/to/private/data/08_infant_associations")
out_dir <- Sys.getenv("FP_OUTPUT_DIR", "/path/to/private/results/08_infant_associations")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
font_family <- "sans"

# 2. Packages.
library(ggplot2)
library(scales)

# 3. Read the inputs, analyze, and draw the figures.
d <- read.delim(file.path(data_dir, "infant_models.tsv"), check.names = FALSE, quote = "")

stopifnot(all(c("heatmap_row", "heatmap_column", "display_effect", "q_value", "endpoint") %in% names(d)))

d$log_effect <- log(d$display_effect)

d$symbol <- ifelse(!is.na(d$q_value) & d$q_value < 0.05, "*", ifelse(tolower(as.character(d$raw_p_displayed)) %in%
    c("true", "t", "1"), "†", ""))

d$endpoint <- factor(d$endpoint, levels = c("detection", "positive_abundance"), labels = c("Detection", "Abundance growth"))

p <- ggplot(d, aes(heatmap_column, heatmap_row, fill = log_effect)) + geom_tile(colour = "white") + geom_text(aes(label = symbol),
    size = 3) + facet_grid(endpoint ~ ., scales = "free", space = "free") + scale_fill_gradient2(low = "#6485A6",
    mid = "white", high = "#CD5B48", midpoint = 0, name = "log(effect ratio)") + labs(x = NULL, y = NULL, fill = NULL, colour = NULL) + (ggplot2::theme_bw(base_size = 9,
    base_family = "sans") + ggplot2::theme(panel.grid.minor = ggplot2::element_blank(), legend.key = ggplot2::element_blank())) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), panel.grid = element_blank())

ggplot2::ggsave(filename = file.path(out_dir, "Figure_3D.pdf"), plot = p, width = 5, height = 5, device = grDevices::cairo_pdf,
    limitsize = FALSE)
