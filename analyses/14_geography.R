# Figure S5. Country-level descriptive prevalence; no city-level substitution.
run <- function() {
  need("ggplot2", "patchwork", "ggnewscale", "jsonlite", "scales")
  source(file.path(repo_root, "R", "world_map.R"), local = TRUE)
  counts <- read_input("world_counts")
  polygons <- read_admin0_polygons(input("world_boundaries"))
  a <- build_country_choropleth_panel(polygons, counts[c("iso3", "count")],
    title = "Metagenome stool samples by country", subtitle = "", font_family = font_family)
  d <- read_input("geography_country_plot")
  stopifnot(length(unique(d$country)) == 31, length(unique(d$species)) == 14,
    all(d$country_total_samples >= 50), all(d$species_overall_prevalence > .05),
    max(abs(d$prevalence - d$detected_samples / d$country_total_samples)) < 1e-8)
  species <- unique(d[c("species", "species_phylogeny_order")]); species <- species[order(species$species_phylogeny_order), ]
  d <- d[order(d$country_order), ]
  d$country_label <- factor(paste0(d$country, " (n=", d$country_total_samples, ")"),
    levels = unique(paste0(d$country, " (n=", d$country_total_samples, ")")))
  d$species_label <- factor(short_species(d$species), levels = rev(short_species(species$species)))
  d$continent <- factor(d$continent, levels = c("Africa", "Asia", "Europe", "North America", "South America", "Oceania"))
  b <- ggplot(d, aes(country_label, species_label)) +
    geom_point(aes(size = prevalence * 100, fill = prevalence * 100), shape = 21, colour = "#6F7378", stroke = .2) +
    scale_size_area(max_size = 3.7, limits = c(0, 100), name = "Detection rate (%)") +
    scale_fill_gradient(low = "white", high = "#6485A6", limits = c(0, 100), guide = "none") +
    facet_grid(~continent, scales = "free_x", space = "free_x") +
    scale_x_discrete(position = "top") + labs(x = NULL, y = NULL) + figure_theme(8) +
    theme(axis.text.x = element_text(angle = 65, hjust = 0), axis.text.y = element_text(face = "italic"), legend.position = "right")
  save_plot(a / b + patchwork::plot_layout(heights = c(1.15, 1)), "Figure_S5", 12, 10)
}
