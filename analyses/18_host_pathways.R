# Figure S4. Retained nominal species-specific Fisher screen (not an FDR discovery set).
# Global q-values come from the upstream pathway/EC/transporter family, not these rows.
run <- function() {
  need("ggplot2", "patchwork")
  d <- read_input("host_pathways")
  stopifnot(all(d$feature_kind == "pathway"), all(d$fisher_p < .05),
    !anyDuplicated(d[c("species_id", "feature_id")]))
  pairs <- list(c("skani_group_48", "skani_group_17"), c("skani_group_45", "skani_group_52"), c("skani_group_59", "skani_group_37"))
  panel <- function(id) {
    z <- d[d$species_id == id, ]; z <- z[order(z$prevalence_diff, z$fisher_p, z$pathway_name), ]
    z$label <- paste0(z$pathway_name, ifelse(z$fisher_q_global < .05, " *", ""))
    z$label <- factor(z$label, levels = unique(z$label))
    ggplot(z, aes(y = label)) + geom_segment(aes(x = 100 * prevalence_human, xend = 100 * prevalence_animal, yend = label), colour = "grey65") +
      geom_point(aes(x = 100 * prevalence_animal, colour = "Animal"), size = 1.5) +
      geom_point(aes(x = 100 * prevalence_human, colour = "Human"), shape = 1, size = 1.5) +
      scale_colour_manual(values = c(Animal = "#E64B35", Human = "#4DBBD5")) +
      scale_x_continuous(limits = c(0, 100), breaks = seq(0, 100, 25)) +
      labs(title = short_species(z$display_name[1]), subtitle = sprintf("A %d | H %d", z$n_animal[1], z$n_human[1]),
        x = "Pathway prevalence (%)", y = NULL, colour = NULL) + figure_theme(6) +
      theme(legend.position = "bottom", panel.grid.major.y = element_blank())
  }
  columns <- lapply(pairs, function(p) patchwork::wrap_plots(lapply(p, panel), ncol = 1,
    heights = vapply(p, function(s) sum(d$species_id == s), numeric(1)) + 3))
  save_plot(patchwork::wrap_plots(columns, nrow = 1, guides = "collect"), "Figure_S4", 15, 10)
}
