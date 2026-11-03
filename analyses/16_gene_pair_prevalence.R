# Figure S6C-F. Shared-family prevalence for highlighted pairs 1, 3, 4 and 7.
run <- function() {
  need("data.table", "ggplot2", "patchwork", "scales")
  highlights <- read_input("gene_highlights"); map <- read_input("gene_label_map")
  map$short <- sub("^s__Faecalibacterium ", "F. ", map$final84_standardized_name)
  target <- highlights[match(c(1, 3, 4, 7), highlights$key_id), ]
  all <- as.data.frame(data.table::fread(input("gene_pair_prevalence"), select = c("species_i", "species_j", "gene", "prevalence_i", "prevalence_j")))
  plots <- lapply(seq_len(nrow(target)), function(i) {
    a <- map$source_label[match(target$species_i_short[i], map$short)]
    b <- map$source_label[match(target$species_j_short[i], map$short)]
    z <- all[all$species_i == a & all$species_j == b, c("gene", "prevalence_i", "prevalence_j")]
    rev <- all[all$species_i == b & all$species_j == a, c("gene", "prevalence_j", "prevalence_i")]
    names(rev) <- names(z); z <- rbind(z, rev)
    z <- z[is.finite(z$prevalence_i) & is.finite(z$prevalence_j) & z$prevalence_i > 0 & z$prevalence_j > 0, ]
    stopifnot(nrow(z) > 0, !anyDuplicated(z$gene))
    z$prevalence_i <- pmin(z$prevalence_i, 1); z$prevalence_j <- pmin(z$prevalence_j, 1)
    xr <- range(z$prevalence_i); yr <- range(z$prevalence_j)
    if (diff(xr) == 0) xr <- c(0, 1)
    if (diff(yr) == 0) yr <- c(0, 1)
    bx <- seq(xr[1], xr[2], length.out = 101); by <- seq(yr[1], yr[2], length.out = 101)
    z$ix <- cut(z$prevalence_i, bx, include.lowest = TRUE, labels = FALSE)
    z$iy <- cut(z$prevalence_j, by, include.lowest = TRUE, labels = FALSE)
    counts <- as.data.frame(table(factor(z$ix, levels = 1:100), factor(z$iy, levels = 1:100)))
    counts$x <- (bx[-1] + bx[-length(bx)])[as.integer(counts$Var1)] / 2
    counts$y <- (by[-1] + by[-length(by)])[as.integer(counts$Var2)] / 2
    p <- ggplot(counts, aes(x, y, fill = log1p(Freq))) + geom_tile() +
      scale_fill_gradientn(colours = c("white", "#C9879E", "#834C7D", "#2D1E3E"),
        limits = log1p(c(0, 10000)), breaks = log1p(c(0, 1, 10, 100, 1000, 10000)),
        labels = c(0, 1, 10, 100, 1000, 10000), name = "count (log1p)") +
      labs(x = target$species_i_short[i], y = target$species_j_short[i]) + figure_theme(8) +
      coord_fixed(ratio = diff(xr) / diff(yr), expand = FALSE) + theme(legend.position = "right")
    top <- ggplot(z, aes(prevalence_i)) + geom_histogram(aes(y = after_stat(log1p(count))), bins = 100, fill = "#834C7D") +
      coord_cartesian(xlim = xr) + theme_void()
    right <- ggplot(z, aes(prevalence_j)) + geom_histogram(aes(y = after_stat(log1p(count))), bins = 100, fill = "#834C7D") +
      coord_flip(xlim = yr) + theme_void()
    write_result(data.frame(key_id = target$key_id[i], shared_families = nrow(z)), paste0("pair_", target$key_id[i], "_summary.tsv"))
    (top + patchwork::plot_spacer()) / (p + right) + patchwork::plot_layout(widths = c(4, 1), heights = c(1, 4))
  })
  save_plot(patchwork::wrap_plots(plots, nrow = 1, guides = "collect"), "Figure_S6CDEF", 14, 3.8)
}
