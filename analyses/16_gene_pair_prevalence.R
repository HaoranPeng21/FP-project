# Figure S6C-F. Shared-family prevalence for highlighted pairs 1, 3, 4 and 7.

# 1. Set these two private directories (or edit the paths directly).
data_dir <- Sys.getenv("FP_DATA_DIR", "/path/to/private/data/16_gene_pair_prevalence")
out_dir <- Sys.getenv("FP_OUTPUT_DIR", "/path/to/private/results/16_gene_pair_prevalence")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
font_family <- "sans"

# 2. Packages.
library(data.table)
library(ggplot2)
library(patchwork)
library(scales)

# 3. Read the inputs, analyze, and draw the figures.
highlights <- read.delim(file.path(data_dir, "gene_highlights.tsv"), check.names = FALSE, quote = "")

map <- read.delim(file.path(data_dir, "gene_label_map.tsv"), check.names = FALSE, quote = "")
if (!"species_name" %in% names(map)) {
  name_column <- grep("standardized_name$", names(map), value = TRUE)
  stopifnot(length(name_column) == 1)
  names(map)[names(map) == name_column] <- "species_name"
}
map$short <- sub("^s__Faecalibacterium ", "F. ", map$species_name)
target <- highlights[match(c(1, 3, 4, 7), highlights$key_id), ]

all <- as.data.frame(data.table::fread(file.path(data_dir, "gene_pair_prevalence.tsv"), select = c("species_i",
    "species_j", "gene", "prevalence_i", "prevalence_j")))
plots <- lapply(seq_len(nrow(target)), function(i) {
    a <- map$source_label[match(target$species_i_short[i], map$short)]
    b <- map$source_label[match(target$species_j_short[i], map$short)]
    z <- all[all$species_i == a & all$species_j == b, c("gene", "prevalence_i", "prevalence_j")]
    rev <- all[all$species_i == b & all$species_j == a, c("gene", "prevalence_j", "prevalence_i")]
    names(rev) <- names(z)
    z <- rbind(z, rev)
    z <- z[is.finite(z$prevalence_i) & is.finite(z$prevalence_j) & z$prevalence_i > 0 & z$prevalence_j > 0, ]
    stopifnot(nrow(z) > 0, !anyDuplicated(z$gene))
    z$prevalence_i <- pmin(z$prevalence_i, 1)
    z$prevalence_j <- pmin(z$prevalence_j, 1)
    xr <- range(z$prevalence_i)
    yr <- range(z$prevalence_j)
    if (diff(xr) == 0)
        xr <- c(0, 1)
    if (diff(yr) == 0)
        yr <- c(0, 1)
    bx <- seq(xr[1], xr[2], length.out = 101)
    by <- seq(yr[1], yr[2], length.out = 101)
    z$ix <- cut(z$prevalence_i, bx, include.lowest = TRUE, labels = FALSE)
    z$iy <- cut(z$prevalence_j, by, include.lowest = TRUE, labels = FALSE)
    counts <- as.data.frame(table(factor(z$ix, levels = 1:100), factor(z$iy, levels = 1:100)))
    counts$x <- (bx[-1] + bx[-length(bx)])[as.integer(counts$Var1)]/2
    counts$y <- (by[-1] + by[-length(by)])[as.integer(counts$Var2)]/2
    p <- ggplot(counts, aes(x, y, fill = log1p(Freq))) + geom_tile() + scale_fill_gradientn(colours = c("white",
        "#C9879E", "#834C7D", "#2D1E3E"), limits = log1p(c(0, 10000)), breaks = log1p(c(0, 1, 10, 100, 1000, 10000)),
        labels = c(0, 1, 10, 100, 1000, 10000), name = "count (log1p)") + labs(x = target$species_i_short[i], y = target$species_j_short[i], fill = NULL, colour = NULL) +
        (ggplot2::theme_bw(base_size = 8, base_family = "sans") + ggplot2::theme(panel.grid.minor = ggplot2::element_blank(),
            legend.key = ggplot2::element_blank())) + coord_fixed(ratio = diff(xr)/diff(yr), expand = FALSE) +
        theme(legend.position = "right")
    top <- ggplot(z, aes(prevalence_i)) + geom_histogram(aes(y = after_stat(log1p(count))), bins = 100, fill = "#834C7D") +
        coord_cartesian(xlim = xr) + theme_void()
    right <- ggplot(z, aes(prevalence_j)) + geom_histogram(aes(y = after_stat(log1p(count))), bins = 100, fill = "#834C7D") +
        coord_flip(xlim = yr) + theme_void()
    write.table(data.frame(key_id = target$key_id[i], shared_families = nrow(z)), file.path(out_dir, paste0("pair_",
        target$key_id[i], "_summary.tsv")), sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
    (top + patchwork::plot_spacer())/(p + right) + patchwork::plot_layout(widths = c(4, 1), heights = c(1, 4))
})

ggplot2::ggsave(filename = file.path(out_dir, "Figure_S6CDEF.pdf"), plot = patchwork::wrap_plots(plots, nrow = 1,
    guides = "collect"), width = 14, height = 3.8, device = grDevices::cairo_pdf, limitsize = FALSE)
