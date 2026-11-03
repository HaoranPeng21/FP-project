# Figure 3B-C. The supplied figure uses the ten-species two-part analysis.
# C is arithmetic mean abundance among detected samples, with infant-bootstrap CIs.
# Stars retain binomial-logit detection / log-abundance Gaussian GEE results.
run <- function() {
  need("ggplot2", "ggrepel", "patchwork")
  d <- read_input("infant_summary"); tests <- read_input("infant_tests"); pal <- read_input("infant_palette")
  stopifnot(nrow(d) == 70, length(unique(d$species)) == 10,
    nrow(tests) == 20, setequal(tests$species, d$species),
    max(abs(d$prevalence - d$detected_samples / d$n_samples)) < 1e-8)
  times <- c("W2", "M1", "M2", "M3", "M6", "M9", "M12")
  d$x <- match(d$timepoint, times)
  cols <- setNames(pal$display_color, pal$species)
  make <- function(endpoint, y, lo, hi, title) {
    z <- d; scale <- if (endpoint == "detection") 100 else 1
    z$y <- z[[y]] * scale; z$lo <- z[[lo]] * scale; z$hi <- z[[hi]] * scale
    lab <- z[z$x == 7, ]; tt <- tests[tests$endpoint == endpoint, ]
    symbols <- tt$significance[match(lab$species, tt$species)]; symbols[is.na(symbols)] <- ""
    lab$label <- paste(short_species(lab$species), symbols)
    ggplot(z, aes(x, y, colour = species, group = species)) +
      geom_errorbar(aes(ymin = lo, ymax = hi), width = .05, alpha = .4, linewidth = .3) +
      geom_line(linewidth = .5) + geom_point(size = 1) +
      ggrepel::geom_text_repel(data = lab, aes(label = label), nudge_x = .3, direction = "y", hjust = 0,
        size = 2.8, fontface = "italic", segment.colour = "grey65", min.segment.length = 0,
        max.overlaps = Inf, seed = 1, box.padding = .1) +
      scale_colour_manual(values = cols) + scale_x_continuous(breaks = 1:7, labels = c("2w", "1m", "2m", "3m", "6m", "9m", "12m"), limits = c(1, 10.5)) +
      labs(x = NULL, y = title) + figure_theme(9) + theme(legend.position = "none")
  }
  b <- make("detection", "prevalence", "prevalence_ci_low", "prevalence_ci_high", "Prevalence (%)")
  c <- make("positive_abundance", "mean_abundance_if_detected", "mean_abundance_if_detected_ci_low",
    "mean_abundance_if_detected_ci_high", "Mean abundance among detected samples (%)")
  save_plot(b, "Figure_3B", 6, 4.5); save_plot(c, "Figure_3C", 6, 4.5)
  save_plot(b + c, "Figure_3BC", 12, 4.5)
}
