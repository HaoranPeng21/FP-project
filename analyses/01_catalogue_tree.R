# Figure 1A. Full 17,103-tip bac120 phylogeny with species, host and continent rings.
run <- function() {
  need("ape")
  source(file.path(repo_root, "R", "trees.R"), local = TRUE)
  tr <- read_phylogeny(input("bac120_tree"))
  meta <- read.csv(input("catalogue_tips"), check.names = FALSE)
  stopifnot(nrow(meta) == 17103, !anyDuplicated(meta$label), all(meta$label %in% tr$tip.label))
  tr <- ape::keep.tip(tr, meta$label)
  meta$Continent <- gsub("_", " ", meta$Continent)
  meta$Continent[meta$Continent %in% c("East Africa", "South Africa", "West Africa")] <- "Africa"
  meta$Host1[tolower(meta$Host1) == "wastewater" & !is.na(meta$Host1)] <- "Wastewater"
  grDevices::cairo_pdf(file.path(out_dir, "Figure_1A.pdf"), width = 11.7, height = 10, family = font_family)
  fan_tree(tr, meta, c("haoran_species_group", "Host1", "Continent"),
    list(species_palette, host1_palette, continent_palette))
  grDevices::dev.off()
}
