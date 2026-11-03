# Figure S7. Maast trees, patristic-distance PERMANOVA, PERMDISP and Fitch permutation.
# Analysis includes known continents with >=5 genomes; BH is across four species per test.
run <- function() {
  need("ape", "vegan")
  source(file.path(repo_root, "R", "trees.R"), local = TRUE)
  manifest <- read.csv(input("maast_tips"), check.names = FALSE)
  manifest$label <- clean_tip(manifest$label)
  manifest$Continent <- gsub("_", " ", trimws(manifest$Continent))
  manifest$Continent[manifest$Continent %in% c("", "NA", "Unknown", "not provided")] <- NA_character_
  manifest$Continent[manifest$Continent %in% c("South Africa", "West Africa")] <- "Africa"
  species <- paste("s__Faecalibacterium", c("longum", "taiwanense", "duncaniae", "prausnitzii"))
  nperm <- 1999L
  results <- vector("list", 4); plots <- vector("list", 4)
  for (i in seq_along(species)) {
    message("  PERMANOVA/Fitch: ", short_species(species[i]))
    z <- manifest[manifest$Host1 == "Human" & manifest$Species == species[i], ]
    stopifnot(length(unique(z$tree_file)) == 1, !anyDuplicated(z$label))
    tr <- read_phylogeny(input("maast_trees", unique(z$tree_file)))
    stopifnot(all(z$label %in% tr$tip.label), all(is.finite(tr$edge.length)), all(tr$edge.length >= 0))
    tr <- ape::keep.tip(tr, z$label)
    known <- z[!is.na(z$Continent), ]
    eligible <- names(which(table(known$Continent) >= 5))
    meta <- known[known$Continent %in% eligible, ]
    tree <- ape::keep.tip(tr, meta$label)
    meta <- meta[match(tree$tip.label, meta$label), ]
    meta$Continent <- factor(meta$Continent, levels = sort(unique(meta$Continent)))
    stopifnot(identical(tree$tip.label, meta$label), nlevels(meta$Continent) >= 2)
    distance <- as.dist(ape::cophenetic.phylo(tree))
    set.seed(7300 + i)
    per <- vegan::adonis2(distance ~ Continent, data = meta, permutations = nperm, by = NULL)
    dispersion <- vegan::betadisper(distance, meta$Continent, type = "median")
    set.seed(8300 + i)
    dis <- vegan::permutest(dispersion, permutations = nperm, pairwise = FALSE)
    # Precompute postorder children and represent state sets as integer bit masks.
    # This is the same Fitch union/intersection score as the production implementation.
    po <- ape::reorder.phylo(tree, "postorder")
    parents <- unique(po$edge[, 1]); children <- split(po$edge[, 2], po$edge[, 1])
    state <- as.integer(2^(as.integer(meta$Continent) - 1))
    score <- function(tip_state) {
      masks <- integer(length(tree$tip.label) + tree$Nnode)
      masks[seq_along(tip_state)] <- tip_state
      changes <- 0L
      for (parent in parents) {
        kids <- children[[as.character(parent)]]
        current <- masks[kids[1]]
        if (length(kids) > 1) for (kid in kids[-1]) {
          intersection <- bitwAnd(current, masks[kid])
          if (intersection > 0) current <- intersection else {
            current <- bitwOr(current, masks[kid]); changes <- changes + 1L
          }
        }
        masks[parent] <- current
      }
      changes
    }
    observed <- score(state)
    set.seed(9300 + i)
    null <- vapply(seq_len(nperm), function(j) score(sample(state, replace = FALSE)), integer(1))
    results[[i]] <- data.frame(species = species[i], n_total = nrow(z), n_known = nrow(known), n_analyzed = nrow(meta),
      n_continents = nlevels(meta$Continent), permanova_F = per$F[1], permanova_R2 = per$R2[1], permanova_p = per$`Pr(>F)`[1],
      betadisper_F = unname(dis$tab[1, "F"]), betadisper_p = unname(dis$tab[1, "Pr(>F)"]),
      fitch_observed = observed, fitch_null_mean = mean(null), fitch_clustering_strength = 1 - observed / mean(null),
      fitch_p = (1 + sum(null <= observed)) / (nperm + 1), permutations = nperm)
    plots[[i]] <- list(tree = tr, metadata = z)
  }
  stats <- do.call(rbind, results)
  for (test in c("permanova", "betadisper", "fitch")) stats[[paste0(test, "_q")]] <- p.adjust(stats[[paste0(test, "_p")]], "BH")
  stopifnot(identical(as.integer(stats$n_analyzed), c(1251L, 631L, 503L, 203L)))
  write_result(stats, "continent_phylogenetic_tests.tsv")
  grDevices::cairo_pdf(file.path(out_dir, "Figure_S7.pdf"), width = 8.27, height = 11.69, family = font_family)
  par(mfrow = c(2, 2), oma = c(4, 0, 0, 0))
  for (i in seq_len(4)) {
    r <- stats[i, ]
    subtitle <- sprintf("n=%d; PERMANOVA R2=%.3f, q=%.4f\nFitch C=%.3f, q=%.4f; PERMDISP F=%.3f, q=%.4f",
      r$n_analyzed, r$permanova_R2, r$permanova_q, r$fitch_clustering_strength, r$fitch_q, r$betadisper_F, r$betadisper_q)
    fan_tree(plots[[i]]$tree, plots[[i]]$metadata, c("Continent", "Host1"), list(continent_palette, host1_palette),
      short_species(species[i]), subtitle, legends = FALSE)
  }
  par(fig = c(0, 1, 0, .07), new = TRUE, mar = c(0, 0, 0, 0))
  plot.new()
  legend("center", legend = names(continent_palette)[names(continent_palette) != "Antarctica"],
    fill = continent_palette[names(continent_palette) != "Antarctica"], ncol = 3, bty = "n", cex = .7)
  grDevices::dev.off()
}
