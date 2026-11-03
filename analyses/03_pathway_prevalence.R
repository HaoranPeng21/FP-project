# Figure 2A: 19 species, 16 pathways; one BH family of 5,256 Fisher tests.
run <- function() {
  need("ComplexHeatmap", "circlize", "grid")
  sp <- read_input("species_scope"); sp <- sp[order(sp$figure_order), ]
  prev <- read_input("pathway_prevalence")
  tests <- read_input("pathway_tests")
  selected <- read_input("pathway_selection")
  stopifnot(nrow(tests) == 5256, nrow(selected) == 14)
  # Recompute each test from genome counts, then correct the complete test family.
  p <- vapply(seq_len(nrow(tests)), function(i) {
    z <- tests[i, ]; fisher.test(matrix(c(z$k_species, z$n_species - z$k_species,
      z$k_reference, z$n_reference - z$k_reference), nrow = 2))$p.value
  }, numeric(1))
  q <- p.adjust(p, "BH")
  stopifnot(max(abs(p - tests$p_value), na.rm = TRUE) < 1e-8,
    max(abs(q - tests$q_value), na.rm = TRUE) < 1e-8)
  limits <- ifelse(tests$species_name %in% c("F. longum", "F. taiwanense", "F. duncaniae"), 3L, 1L)
  chosen <- tests[as_bool(tests$eligible) & !is.na(tests$rank_within_species) & tests$rank_within_species <= limits, ]
  stopifnot(nrow(chosen) == 24, setequal(unique(chosen$representative_feature_id), selected$feature_id))
  ids <- c("AC-BUT-COATRANS", "PWY0-1314", selected$feature_id)
  stopifnot(length(ids) == 16, !anyDuplicated(ids))
  m <- matrix(NA_real_, nrow(sp), length(ids), dimnames = list(sp$display_name, ids))
  sig <- matrix(FALSE, nrow(m), ncol(m))
  for (j in seq_along(ids)) {
    ix <- match(paste(sp$display_name, ids[j]), paste(prev$display_name, prev$feature_id))
    m[, j] <- 100 * prev$prevalence[ix]
    ti <- match(paste(sp$display_name, ids[j]), paste(tests$species, tests$feature_id))
    sig[, j] <- !is.na(ti) & q[ti] < .05 & abs(tests$delta[ti]) >= .25
  }
  stopifnot(!anyNA(m), nrow(m) == 19)
  labels <- prev$pathway_name[match(ids, prev$feature_id)]
  labels[ids == "AC-BUT-COATRANS"] <- "Acetyl-CoA fermentation to butanoate"
  labels[ids == "PWY0-1314"] <- "Fructose degradation"
  labels[ids == "PWY-7395"] <- "D-GalN and GalNAc degradation"
  labels[ids == "PWY-6386"] <- "UDP-MurNAc-pentapeptide biosynthesis II"
  neq <- c(1, 1, selected$n_equivalent_pathways)
  labels <- ifelse(neq > 1, paste0(labels, " (+", neq - 1, ")"), labels)
  ht <- ComplexHeatmap::Heatmap(m, name = "Pathway prevalence (%)",
    col = circlize::colorRamp2(c(0, 25, 50, 75, 100), c("white", "#E9EFF5", "#C7D6E4", "#91ACC4", "#6485A6")),
    cluster_rows = FALSE, cluster_columns = TRUE,
    column_split = factor(c(rep("Shared", 2), rep("Differential", 14)), levels = c("Shared", "Differential")),
    cluster_column_slices = FALSE, clustering_distance_columns = "pearson", clustering_method_columns = "average",
    row_labels = sp$short_name, row_names_side = "left", column_labels = labels, column_names_rot = 45,
    row_names_gp = grid::gpar(fontsize = 8, fontface = "italic"), column_names_gp = grid::gpar(fontsize = 7),
    rect_gp = grid::gpar(col = "white"), cell_fun = function(j, i, x, y, w, h, fill) {
      if (isTRUE(sig[i, j])) grid::grid.text("*", x, y, gp = grid::gpar(fontsize = 9))
    })
  grDevices::cairo_pdf(file.path(out_dir, "Figure_2A.pdf"), width = 8.8, height = 6, family = font_family)
  ComplexHeatmap::draw(ht, heatmap_legend_side = "bottom"); grDevices::dev.off()
  write_result(data.frame(tests, recomputed_p = p, recomputed_q = q), "pathway_tests.tsv")
}
