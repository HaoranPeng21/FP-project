# Figure 3B-C. The supplied figure uses the ten-species two-part analysis.
# C is arithmetic mean abundance among detected samples, with infant-bootstrap CIs.
# Stars retain binomial-logit detection / log-abundance Gaussian GEE results.


data_dir <- Sys.getenv("FP_DATA_DIR", "/path/to/private/data/07_infant_trajectories")
out_dir <- Sys.getenv("FP_OUTPUT_DIR", "/path/to/private/results/07_infant_trajectories")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
font_family <- "sans"
library(ggplot2)
library(ggrepel)
library(patchwork)

d <- read.delim(file.path(data_dir, "infant_summary.tsv"), check.names = FALSE, quote = "")

tests <- read.delim(file.path(data_dir, "infant_tests.tsv"), check.names = FALSE, quote = "")

pal <- read.delim(file.path(data_dir, "infant_palette.tsv"), check.names = FALSE, quote = "")
stopifnot(nrow(d) == 70, length(unique(d$species)) == 10, nrow(tests) == 20, setequal(tests$species,
    d$species), max(abs(d$prevalence - d$detected_samples/d$n_samples)) < 1e-08)
times <- c("W2", "M1", "M2", "M3", "M6", "M9", "M12")
d$x <- match(d$timepoint, times)
cols <- setNames(pal$display_color, pal$species)
endpoints <- c("detection", "positive_abundance")
mean_columns <- c("prevalence", "mean_abundance_if_detected")
lower_columns <- c("prevalence_ci_low", "mean_abundance_if_detected_ci_low")
upper_columns <- c("prevalence_ci_high", "mean_abundance_if_detected_ci_high")
y_titles <- c("Prevalence (%)", "Mean abundance among detected samples (%)")
plots <- vector("list", 2)

for (i in 1:2) {
    endpoint <- endpoints[i]
    y <- mean_columns[i]
    lo <- lower_columns[i]
    hi <- upper_columns[i]
    title <- y_titles[i]
    z <- d
    scale <- if (endpoint == "detection")
        100
    else 1
    z$y <- z[[y]] * scale
    z$lo <- z[[lo]] * scale
    z$hi <- z[[hi]] * scale
    lab <- z[z$x == 7, ]
    tt <- tests[tests$endpoint == endpoint, ]
    symbols <- tt$significance[match(lab$species, tt$species)]
    symbols[is.na(symbols)] <- ""
    lab$label <- paste(sub("^Faecalibacterium ", "F. ", sub("^Faecalibacterium prausnitzii_F$", "Faecalibacterium langellae",
        sub("^Faecalibacterium sp900539945$", "Faecalibacterium harmsenii", sub("^s__", "", lab$species)))),
        symbols)
    plots[[i]] <- ggplot(z, aes(x, y, colour = species, group = species)) + geom_errorbar(aes(ymin = lo,
        ymax = hi), width = 0.05, alpha = 0.4, linewidth = 0.3) + geom_line(linewidth = 0.5) + geom_point(size = 1) +
        ggrepel::geom_text_repel(data = lab, aes(label = label), nudge_x = 0.3, direction = "y", hjust = 0,
            size = 2.8, fontface = "italic", segment.colour = "grey65", min.segment.length = 0, max.overlaps = Inf,
            seed = 1, box.padding = 0.1) + scale_colour_manual(values = cols) + scale_x_continuous(breaks = 1:7,
        labels = c("2w", "1m", "2m", "3m", "6m", "9m", "12m"), limits = c(1, 10.5)) + labs(y = title, x = NULL, fill = NULL, colour = NULL) +
        (ggplot2::theme_bw(base_size = 9, base_family = "sans") + ggplot2::theme(panel.grid.minor = ggplot2::element_blank(),
            legend.key = ggplot2::element_blank())) + theme(legend.position = "none")
}

ggsave(file.path(out_dir, "Figure_3B.pdf"), plots[[1]], width = 6, height = 4.5, device = cairo_pdf)

ggsave(file.path(out_dir, "Figure_3C.pdf"), plots[[2]], width = 6, height = 4.5, device = cairo_pdf)

ggsave(file.path(out_dir, "Figure_3BC.pdf"), plots[[1]] + plots[[2]], width = 12, height = 4.5, device = cairo_pdf)
