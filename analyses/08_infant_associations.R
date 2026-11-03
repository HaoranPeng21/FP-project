# Figure 3D. Selected coefficients from the prespecified adjusted longitudinal models.
# Keep the upstream full-family q-values; do not correct only these displayed cells.
run <- function() {
  need("ggplot2", "scales")
  d <- read_input("infant_models")
  check_columns(d, c("heatmap_row", "heatmap_column", "display_effect", "q_value", "endpoint"))
  d$log_effect <- log(d$display_effect)
  d$symbol <- ifelse(!is.na(d$q_value) & d$q_value < .05, "*", ifelse(as_bool(d$raw_p_displayed), "\u2020", ""))
  d$endpoint <- factor(d$endpoint, levels = c("detection", "positive_abundance"),
    labels = c("Detection", "Abundance growth"))
  p <- ggplot(d, aes(heatmap_column, heatmap_row, fill = log_effect)) +
    geom_tile(colour = "white") + geom_text(aes(label = symbol), size = 3) +
    facet_grid(endpoint ~ ., scales = "free", space = "free") +
    scale_fill_gradient2(low = "#6485A6", mid = "white", high = "#CD5B48", midpoint = 0, name = "log(effect ratio)") +
    labs(x = NULL, y = NULL) + figure_theme(9) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), panel.grid = element_blank())
  save_plot(p, "Figure_3D", 5, 5)
}
