# Figure 3A. All 14 species passing overall stool prevalence >5% remain visible.
# Colour uses abundance among detected samples; bubble area uses prevalence.

# 1. Set these two private directories (or edit the paths directly).
data_dir <- Sys.getenv("FP_DATA_DIR", "/path/to/private/data/06_human_distribution")
out_dir <- Sys.getenv("FP_OUTPUT_DIR", "/path/to/private/results/06_human_distribution")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
font_family <- "sans"

# 2. Packages.
library(ggplot2)
library(dplyr)
library(scales)

# 3. Read the inputs, analyze, and draw the figures.
d <- read.delim(file.path(data_dir, "human_strata.tsv"), check.names = FALSE, quote = "")

order <- read.delim(file.path(data_dir, "human_species_order.tsv"), check.names = FALSE, quote = "")

stopifnot(all(c("n_positive", "n_samples", "mean_positive_abundance_pct", "group", "panel") %in% names(d)))

stopifnot(length(unique(order$species)) == 14, max(abs(100 * d$n_positive/d$n_samples - d$prevalence_pct)) < 1e-08)

d$species_label <- factor(sub("^Faecalibacterium ", "F. ", sub("^Faecalibacterium prausnitzii_F$", "Faecalibacterium langellae",
    sub("^Faecalibacterium sp900539945$", "Faecalibacterium harmsenii", sub("^s__", "", d$species)))), levels = rev(sub("^Faecalibacterium ",
    "F. ", sub("^Faecalibacterium prausnitzii_F$", "Faecalibacterium langellae", sub("^Faecalibacterium sp900539945$",
        "Faecalibacterium harmsenii", sub("^s__", "", order$species[order(order$display_rank)]))))))
d$panel <- factor(d$panel, levels = c("Healthy continent", "Adult disease", "Healthy age"))
d <- d[order(d$panel, d$group_order), ]

d$group_label <- factor(paste0(d$group, " (n=", format(d$n_samples, big.mark = ",", trim = TRUE), ")"), levels = unique(paste0(d$group,
    " (n=", format(d$n_samples, big.mark = ",", trim = TRUE), ")")))
d$colour_value <- log2(d$mean_positive_abundance_pct)
d$colour_value[d$n_positive == 0] <- NA_real_

p <- ggplot(d, aes(group_label, species_label)) + geom_point(aes(size = prevalence_pct, fill = colour_value), shape = 21,
    colour = "#8F9BAD", stroke = 0.25) + scale_size_area(max_size = 5, limits = c(0, 100), name = "Prevalence (%)") +
    scale_fill_gradient2(low = "#6485A6", mid = "white", high = "#CD5B48", midpoint = 0, limits = c(-3, 3), oob = scales::squish,
        na.value = "white", name = "log2 abundance (vs 1%)") + facet_grid(~panel, scales = "free_x", space = "free_x") +
    labs(x = NULL, y = NULL, fill = NULL, colour = NULL) + (ggplot2::theme_bw(base_size = 9, base_family = "sans") + ggplot2::theme(panel.grid.minor = ggplot2::element_blank(),
    legend.key = ggplot2::element_blank())) + theme(axis.text.x = element_text(angle = 60, hjust = 0), axis.text.y = element_text(face = "italic"),
    legend.position = "bottom") + scale_x_discrete(position = "top")

ggplot2::ggsave(filename = file.path(out_dir, "Figure_3A.pdf"), plot = p, width = 12, height = 6, device = grDevices::cairo_pdf,
    limitsize = FALSE)
