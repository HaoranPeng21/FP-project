# Figure 2B. Colours are within-host species composition; numbers are genome counts.
run <- function() {
  need("ComplexHeatmap", "circlize", "grid")
  source(file.path(repo_root, "R", "trees.R"), local = TRUE)
  sp <- read_input("species_scope"); sp <- sp[order(sp$figure_order), ]
  d <- read_input("host_counts")
  counts <- t(as.matrix(d[sp$display_name])); colnames(counts) <- d$host
  stopifnot(sum(counts) == 15811, !anyNA(counts), all(colSums(counts) > 0))
  percent <- sweep(counts, 2, colSums(counts), "/") * 100
  stopifnot(max(abs(colSums(percent) - 100)) < 1e-8)
  ht <- ComplexHeatmap::Heatmap(percent, name = "Within-host composition (%)",
    col = circlize::colorRamp2(c(0, 25, 50, 75, 100), c("white", "#E9EFF5", "#C7D6E4", "#91ACC4", "#6485A6")),
    cluster_rows = FALSE, cluster_columns = FALSE, row_labels = sp$short_name, row_names_side = "left",
    top_annotation = ComplexHeatmap::HeatmapAnnotation(Host = d$host,
      col = list(Host = palette_for(d$host, host_palette)), show_annotation_name = FALSE, show_legend = FALSE),
    column_labels = paste0(d$host, " (n=", colSums(counts), ")"), column_names_rot = 45,
    row_names_gp = grid::gpar(fontsize = 8, fontface = "italic"), column_names_gp = grid::gpar(fontsize = 7),
    rect_gp = grid::gpar(col = "white"), cell_fun = function(j, i, x, y, w, h, fill) {
      if (counts[i, j] > 0) grid::grid.text(counts[i, j], x, y, gp = grid::gpar(fontsize = 7))
    })
  grDevices::cairo_pdf(file.path(out_dir, "Figure_2B.pdf"), width = 7.5, height = 5.5, family = font_family)
  ComplexHeatmap::draw(ht, heatmap_legend_side = "bottom"); grDevices::dev.off()
}
