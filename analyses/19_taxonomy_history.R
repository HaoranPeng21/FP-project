# Figure S1. Public taxonomic-history schematic plus catalogue-derived release labels.
# Timeline text is transcribed from the supplied figure; it is not participant data.
run <- function() {
  need("ggplot2", "patchwork")
  timeline <- data.frame(
    x = c(0, 1.4, 2.5, 2.5, 3.8, 5, 6.2, 7.4),
    y = c(1, -1, 1, -1, -1, 1, -1, 1),
    year = c("2002", "2018", "2021", "2021", "2022", "2024", "2025", "2026"),
    label = c("Duncan et al.\nF. prausnitzii\nATCC 27768\nGenus established",
      "Fitzgerald et al.\nGenome-scale diversity\n'F. moorei' proposed; not a separate\naccepted species in this schematic",
      "Zou et al.\nF. longum (CM04-06)\nF. butyricigenerans (AF52-21)",
      "Gilroy et al.\nSeven chicken-gut candidatus taxa\navium; faecigallinarum; faecipullorum\ngallistercoris; intestinavium\nintestinigallinarum; intestinipullorum",
      "Sakamoto et al.\nF. duncaniae\nF. hattorii\nF. gallinarum",
      "Liou; Plomp & Harmsen\nF. taiwanense (HLW78)\nF. wellingii (HTF-F)",
      "Ang; Hitch et al.\nF. langellae (CNCM I-4541)\nF. tardum (CLA-AA-H175)\nF. intestinale (CLA-AA-H281)",
      "Pan et al.\nF. harmsenii* (HTF-05B)\nGTDB R232: sp900539945"))
  a <- ggplot(timeline, aes(x, 0)) + geom_hline(yintercept = 0, colour = "#6485A6") +
    geom_segment(aes(xend = x, yend = y * .65), colour = "#91ACC4") +
    geom_point(shape = 21, fill = "white", colour = "#6485A6", size = 2) +
    geom_text(aes(y = y * .8, label = year), fontface = "bold", size = 3) +
    geom_text(aes(y = y * 1.35, label = label), size = 2.15, lineheight = 1) +
    coord_cartesian(ylim = c(-2.3, 2.3), clip = "off") + theme_void() +
    labs(title = "Cultivated species descriptions (elapsed time compressed)") + theme(plot.title = element_text(size = 10, face = "bold"))
  catalog <- read.csv(input("catalogue"), check.names = FALSE)
  catalog <- catalog[catalog$ID %in% read_input("analysis_set")$ID, ]
  focal <- paste("s__Faecalibacterium", c("prausnitzii", "longum", "butyricigenerans", "duncaniae", "hattorii", "gallinarum", "taiwanense", "wellingii", "prausnitzii_F", "tardum", "intestinale", "sp900539945"))
  catalog$taxon <- sub("\\s*\\(.*\\)$", "", catalog$Haoran_species)
  release_labels <- do.call(rbind, lapply(focal, function(sp) {
    z <- catalog[catalog$taxon == sp & !is.na(catalog$taxon), c("Species", "Species_v220", "Species_v226")]
    z <- z[complete.cases(z), ]; key <- apply(z, 1, paste, collapse = "|")
    modal <- names(sort(table(key), decreasing = TRUE))[1]
    values <- strsplit(modal, "|", fixed = TRUE)[[1]]
    final <- if (sp == "s__Faecalibacterium sp900539945") "sp900539945 / F. harmsenii*" else short_species(sp)
    data.frame(taxon = sp, release = 1:4, label = c(sub("^s__Faecalibacterium ", "", values), final))
  }))
  nodes <- unique(release_labels[c("release", "label")])
  nodes$y <- NA_real_
  for (j in 1:4) {
    ix <- which(nodes$release == j)
    # Final-identity order follows the supplied schematic; earlier releases use lexical order.
    ii <- if (j == 4) seq_along(ix) else order(nodes$label[ix])
    nodes$y[ix[ii]] <- rev(seq_along(ix))
  }
  release_labels$y <- nodes$y[match(paste(release_labels$release, release_labels$label), paste(nodes$release, nodes$label))]
  edge <- do.call(rbind, lapply(focal, function(sp) {
    z <- release_labels[release_labels$taxon == sp, ]
    data.frame(taxon = sp, x = 1:3, xend = 2:4, y = z$y[1:3], yend = z$y[2:4])
  }))
  b <- ggplot() + geom_curve(data = edge, aes(x = x + .12, y, xend = xend - .12, yend = yend, colour = taxon),
    curvature = .12, linewidth = 1, alpha = .5) +
    geom_label(data = nodes, aes(release, y, label = label), size = 2.7, fill = "white", linewidth = .15) +
    scale_x_continuous(breaks = 1:4, labels = c("R207 | 2022", "R220 | 2024", "R226 | 2025", "R232 | 2026\nCultivated identity"), position = "top", expand = expansion(mult = .12)) +
    theme_void() + theme(axis.text.x = element_text(size = 9), legend.position = "none")
  save_plot(a / b + patchwork::plot_layout(heights = c(1, 1.35)), "Figure_S1", 13, 11)
  write_result(release_labels, "taxonomy_release_labels.tsv")
}
