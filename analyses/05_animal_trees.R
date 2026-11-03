# Figure 2C: four mixed-host subtrees. Figure S3: the complete selected animal overview.
run <- function() {
  need("ape")
  source(file.path(repo_root, "R", "trees.R"), local = TRUE)
  tr <- read_phylogeny(input("bac120_tree"))
  overview <- read.csv(input("animal_tips"), check.names = FALSE)
  overview <- overview[!is.na(overview$Host) & overview$Host != "Bioreactor Sludge", ]
  stopifnot(nrow(overview) == 998, all(overview$label %in% tr$tip.label))
  grDevices::cairo_pdf(file.path(out_dir, "Figure_S3.pdf"), width = 11.7, height = 9, family = font_family)
  fan_tree(ape::keep.tip(tr, overview$label), overview, c("haoran_species_group", "Host"), list(species_palette, host_palette))
  grDevices::dev.off()
  tips <- read.csv(input("animal_species_tips"), check.names = FALSE)
  species <- paste("s__Faecalibacterium", c("longum", "duncaniae", "taiwanense", "prausnitzii"))
  grDevices::cairo_pdf(file.path(out_dir, "Figure_2C.pdf"), width = 12, height = 3.5, family = font_family)
  par(mfrow = c(1, 4))
  for (sp in species) {
    z <- tips[tips$haoran_species == sp & !is.na(tips$haoran_species), ]
    stopifnot(nrow(z) > 1, all(z$label %in% tr$tip.label))
    fan_tree(ape::keep.tip(tr, z$label), z, c("haoran_species_group", "Host"), list(species_palette, host_palette),
      short_species(sp), paste0("n = ", nrow(z), "; ", paste(sort(unique(z$Host)), collapse = ", ")), legends = FALSE)
  }
  grDevices::dev.off()
}
