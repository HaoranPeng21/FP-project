species_palette <- c(
  "s__Faecalibacterium prausnitzii" = "black", "s__Faecalibacterium tardum" = "#1f78b4",
  "s__Faecalibacterium duncaniae" = "#e41a1c", "s__Faecalibacterium taiwanense" = "#ff7f00",
  "s__Faecalibacterium prausnitzii_E" = "#4daf4a", "s__Faecalibacterium prausnitzii_F" = "#984ea3",
  "s__Faecalibacterium longum" = "#a65628", "s__Faecalibacterium hattorii" = "#66c2a5",
  "s__Faecalibacterium prausnitzii_I" = "#fc8d62", "s__Faecalibacterium prausnitzii_L" = "#e5c494",
  "s__Faecalibacterium prausnitzii_M" = "#b3de69", "s__Faecalibacterium wellingii" = "#D7BDE2",
  "s__Faecalibacterium butyricigenerans" = "#FDDB76", "s__Faecalibacterium gallinarum" = "#88CEE6",
  "s__Faecalibacterium intestinale" = "#E89DA0", "s__Faecalibacterium avium" = "#fccde5",
  "s__Faecalibacterium faecipullorum" = "#8dd3c7", "s__Faecalibacterium intestinigallinarum" = "#ffffb3",
  "s__Faecalibacterium gallistercoris" = "#bebada", "s__Faecalibacterium faecigallinarum" = "#fb8072",
  "s__Faecalibacterium intestinipullorum" = "#fdb462", "s__Faecalibacterium intestinavium" = "#bc80bd",
  "s__Faecalibacterium sp" = "#ccebc5", "s__Faecalibacterium sp900539945" = "#3C5488",
  "Candidate clades (14)" = "#bdbdbd")
host_palette <- c("Human gut" = "#1b9e77", "Wastewater" = "#8dd3c7", "Dog" = "#1f78b4",
  "Cat" = "#87CEFA", "Mongoose" = "#ffff99", "Rodent" = "#e6ab02", "Cattle" = "#e31a1c",
  "Goat" = "#fb9a99", "Primate" = "#6a3d9a", "Whale" = "#5398CD", "Pig" = "#fdbf6f",
  "Tamandua" = "#b15928", "Chicken" = "#cab2d6", "Raccoon" = "#F8766D", "Sea Lion" = "#00BA38", "Seal" = "#619CFF")
continent_palette <- c(Africa = "#FFCE30", Antarctica = "#E64B35", Asia = "#d62728",
  Europe = "#1f77b4", "North America" = "#ff7f0e", "South America" = "#2ca02c", Oceania = "#9467bd")
host1_palette <- c(Animal = "#d95f02", Human = "#1b9e77", Wastewater = "#ffb000")
palette_for <- function(values, palette) {
  missing <- sort(setdiff(unique(na.omit(values)), names(palette)))
  if (length(missing)) palette <- c(palette, setNames(grDevices::hcl.colors(length(missing), "Dark 3"), missing))
  palette
}
# Draw annular sectors using the angles of the plotted tips. Branch lengths are preserved.
fan_tree <- function(tr, meta, rings, palettes, title = "", subtitle = "", legends = TRUE) {
  meta <- meta[match(tr$tip.label, meta$label), , drop = FALSE]
  stopifnot(identical(tr$tip.label, meta$label), !anyDuplicated(meta$label))
  tr <- ape::ladderize(tr)
  meta <- meta[match(tr$tip.label, meta$label), , drop = FALSE]
  radius <- max(ape::node.depth.edgelength(tr)[seq_along(tr$tip.label)])
  stopifnot(is.finite(radius), radius > 0)
  tr$edge.length <- tr$edge.length / radius
  old <- par(mar = c(1, 1, 4.5, if (legends) 10 else 1), family = font_family, xpd = NA)
  on.exit(par(old))
  ape::plot.phylo(tr, type = "fan", show.tip.label = FALSE, edge.width = .25,
    x.lim = c(-1.28, 1.28), y.lim = c(-1.28, 1.28), no.margin = FALSE)
  xy <- get("last_plot.phylo", envir = ape::.PlotPhyloEnv)
  n <- length(tr$tip.label); theta <- atan2(xy$yy[seq_len(n)], xy$xx[seq_len(n)])
  idx <- order(theta); theta <- theta[idx]; half <- pi / n
  legend_y <- 1.2
  for (j in seq_along(rings)) {
    values <- as.character(meta[[rings[j]]])[idx]
    if (rings[j] == "haoran_species_group") values[grepl("sp\\. Haoran_", values)] <- "Candidate clades (14)"
    pal <- palette_for(values, palettes[[j]])
    cols <- unname(pal[values]); cols[is.na(cols)] <- "#eeeeee"
    r0 <- 1.01 + (j - 1) * .075; r1 <- r0 + .06
    for (i in seq_len(n)) {
      a <- seq(theta[i] - half, theta[i] + half, length.out = 4)
      polygon(c(r0 * cos(a), r1 * cos(rev(a))), c(r0 * sin(a), r1 * sin(rev(a))), col = cols[i], border = NA)
    }
    if (legends) {
      used <- names(pal)[names(pal) %in% values]
      heading <- if (rings[j] == "haoran_species_group") "Species" else if (rings[j] == "Host1") "Host" else rings[j]
      geometry <- legend(1.32, legend_y, legend = short_species(used), fill = pal[used],
        bty = "n", cex = .6, title = heading, title.cex = .7, xjust = 0)
      legend_y <- legend_y - geometry$rect$h - .08
    }
  }
  if ("Type_isolate" %in% names(meta)) {
    iso <- which(!is.na(meta$Type_isolate) & meta$Type_isolate == "Isolate")
    points(xy$xx[iso], xy$yy[iso], pch = 16, col = "#d62728", cex = .22)
  }
  title(main = title, line = 2.2, cex.main = .85)
  if (nzchar(subtitle)) mtext(subtitle, side = 3, line = .3, cex = .6)
  invisible(length(tr$tip.label))
}
