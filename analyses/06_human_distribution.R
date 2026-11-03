# Figure 3A. All 14 species passing overall stool prevalence >5% remain visible.
# Colour uses abundance among detected samples; bubble area uses prevalence.
run <- function() {
  need("ggplot2", "dplyr", "scales")
  d <- read_input("human_strata"); order <- read_input("human_species_order")
  check_columns(d, c("n_positive", "n_samples", "mean_positive_abundance_pct", "group", "panel"))
  stopifnot(length(unique(order$species)) == 14,
    max(abs(100 * d$n_positive / d$n_samples - d$prevalence_pct)) < 1e-8)
  d$species_label <- factor(short_species(d$species), levels = rev(short_species(order$species[order(order$display_rank)])))
  d$panel <- factor(d$panel, levels = c("Healthy continent", "Adult disease", "Healthy age"))
  d <- d[order(d$panel, d$group_order), ]
  d$group_label <- factor(paste0(d$group, " (n=", format(d$n_samples, big.mark = ",", trim = TRUE), ")"),
    levels = unique(paste0(d$group, " (n=", format(d$n_samples, big.mark = ",", trim = TRUE), ")")))
  d$colour_value <- log2(d$mean_positive_abundance_pct)
  d$colour_value[d$n_positive == 0] <- NA_real_
  p <- ggplot(d, aes(group_label, species_label)) +
    geom_point(aes(size = prevalence_pct, fill = colour_value), shape = 21, colour = "#8F9BAD", stroke = .25) +
    scale_size_area(max_size = 5, limits = c(0, 100), name = "Prevalence (%)") +
    scale_fill_gradient2(low = "#6485A6", mid = "white", high = "#CD5B48", midpoint = 0,
      limits = c(-3, 3), oob = scales::squish, na.value = "white", name = "log2 abundance (vs 1%)") +
    facet_grid(~panel, scales = "free_x", space = "free_x") + labs(x = NULL, y = NULL) + figure_theme(9) +
    theme(axis.text.x = element_text(angle = 60, hjust = 0), axis.text.y = element_text(face = "italic"),
      legend.position = "bottom") + scale_x_discrete(position = "top")
  save_plot(p, "Figure_3A", 12, 6)
}
