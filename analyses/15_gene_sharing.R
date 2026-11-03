# Figure S6A-B. 25-species pangenome subset; 300 dependent species pairs.
# ANY-presence repertoire Jaccard is compared with representative bac120 path distance.
run <- function() {
  need("ape", "ggtree", "ggplot2", "ggrepel", "patchwork", "scales")
  map <- read_input("gene_label_map"); reps <- read_input("gene_representatives")
  stopifnot(nrow(map) == 25, !anyDuplicated(map$species_id))
  map$tip <- clean_tip(reps$representative_ID[match(map$species_id, reps$species_id)])
  tr <- read_phylogeny(input("bac120_tree"))
  stopifnot(!anyNA(map$tip), all(map$tip %in% tr$tip.label))
  tr <- ape::keep.tip(tr, map$tip)
  tr$tip.label <- map$species_id[match(tr$tip.label, map$tip)]
  raw <- as.matrix(read.delim(input("gene_ratio"), row.names = 1, check.names = FALSE))
  stopifnot(setequal(rownames(raw), map$source_label))
  rownames(raw) <- map$species_id[match(rownames(raw), map$source_label)]
  colnames(raw) <- map$species_id[match(colnames(raw), map$source_label)]
  jaccard <- raw[sort(map$species_id), sort(map$species_id)]
  stopifnot(max(abs(jaccard - t(jaccard))) < 1e-12, all(diag(jaccard) == 1))
  distance <- ape::cophenetic.phylo(tr)[rownames(jaccard), colnames(jaccard)]
  upper <- upper.tri(distance)
  index <- which(upper, arr.ind = TRUE)
  pairs <- data.frame(species_i_id = rownames(distance)[index[, 1]], species_j_id = colnames(distance)[index[, 2]],
    phylo_distance = distance[upper], shared_jaccard = jaccard[upper])
  pairs$species_i_short <- short_species(map$final84_standardized_name[match(pairs$species_i_id, map$species_id)])
  pairs$species_j_short <- short_species(map$final84_standardized_name[match(pairs$species_j_id, map$species_id)])
  rho <- cor(pairs$phylo_distance, pairs$shared_jaccard, method = "spearman")
  stopifnot(nrow(pairs) == 300, abs(rho - (-.726553365)) < 1e-8)
  # QAP permutes species labels together on the rows and columns of the distance matrix.
  # Standardized off-diagonal ranks can be reused because permutations preserve their multiset.
  standard_rank <- function(x) {r <- rank(x) - mean(rank(x)); r / sqrt(sum(r^2))}
  ranked_distance <- matrix(0, 25, 25)
  ranked_distance[upper] <- standard_rank(distance[upper])
  ranked_distance <- ranked_distance + t(ranked_distance)
  sharing_rank <- standard_rank(jaccard[upper])
  nperm <- 99999L; set.seed(20260823)
  null <- vapply(seq_len(nperm), function(i) {
    permutation <- sample.int(25)
    sum(ranked_distance[permutation, permutation][upper] * sharing_rank)
  }, numeric(1))
  qap_p <- (1 + sum(abs(null) >= abs(rho))) / (nperm + 1)
  frozen <- read_input("gene_qap")
  stopifnot(nrow(frozen) == 1, frozen$n_pairs == 300, frozen$qap_permutations == nperm,
    abs(frozen$spearman_rho - rho) < 1e-12)
  write_result(data.frame(n_species = 25, n_pairs = 300, rho, qap_permutations = nperm,
    qap_seed = 20260823, qap_p_R = qap_p, qap_p_figure = frozen$qap_two_sided_p,
    random_generator = "R Mersenne-Twister"), "qap_test.tsv")
  # R and the earlier NumPy pipeline use different RNGs; the test and seed are explicit.
  h <- pairs[(pairs$phylo_distance >= .1 & pairs$shared_jaccard >= .1) | pairs$shared_jaccard >= .15, ]
  h$highlight_group <- ifelse(h$phylo_distance >= .1 & h$shared_jaccard >= .1, "Divergent yet similar", "Top overall similarity")
  h <- h[order(-h$shared_jaccard), ]; h$key_id <- seq_len(nrow(h))
  stopifnot(nrow(h) == 13)
  write_result(pairs, "species_pairs.tsv"); write_result(h, "highlight_pairs.tsv")
  tree_plot <- ggtree::ggtree(tr, linewidth = .3)
  td <- tree_plot$data[tree_plot$data$isTip, ]; labels <- setNames(short_species(map$final84_standardized_name), map$species_id)
  td$name <- labels[td$label]
  tree_plot <- tree_plot + geom_text(data = td, aes(x = x, y = y, label = name), hjust = 0, nudge_x = .005,
    size = 2.2, fontface = "italic") + xlim(0, max(td$x) + .18) + theme_void()
  ypos <- setNames(td$y, td$label)
  diamonds <- do.call(rbind, lapply(seq_len(nrow(pairs)), function(i) {
    y1 <- ypos[pairs$species_i_id[i]]; y2 <- ypos[pairs$species_j_id[i]]
    x <- abs(y2 - y1); y <- (y1 + y2) / 2
    data.frame(pair = i, x = c(x - 1, x, x + 1, x), y = c(y, y - .5, y, y + .5), value = pairs$shared_jaccard[i])
  }))
  colours <- c("white", "#D8B3C6", "#9A6293", "#2D1E3E")
  heat <- ggplot(diamonds, aes(x, y, group = pair, fill = value)) + geom_polygon(colour = "white", linewidth = .1) +
    scale_fill_gradientn(colours = colours, trans = scales::pseudo_log_trans(sigma = .001), guide = "none") +
    scale_y_continuous(limits = c(.5, 25.5), expand = c(0, 0)) + theme_void()
  tree_plot <- tree_plot + scale_y_continuous(limits = c(.5, 25.5), expand = c(0, 0))
  scatter <- ggplot(pairs, aes(phylo_distance, shared_jaccard)) + geom_point(aes(fill = shared_jaccard), shape = 21, size = 1.7, colour = "#7B6879") +
    geom_point(data = h, aes(colour = highlight_group), size = 2.3) +
    ggrepel::geom_text_repel(data = h, aes(label = key_id, colour = highlight_group), size = 2.5, seed = 1, max.overlaps = Inf) +
    scale_colour_manual(values = c("Divergent yet similar" = "#00A087", "Top overall similarity" = "#E64B35")) +
    scale_fill_gradientn(colours = colours, trans = scales::pseudo_log_trans(sigma = .001), name = "Jaccard") +
    labs(x = "Phylogenetic distance (bac120 path length)", y = "Gene-repertoire Jaccard",
      subtitle = sprintf("300 pairs; Spearman rho = %.3f; QAP P = %.1e", rho, frozen$qap_two_sided_p), colour = NULL) + figure_theme(8) + theme(legend.position = "bottom")
  key <- h
  key$label <- paste(key$key_id, key$species_i_short, "vs", key$species_j_short)
  key$y <- rev(seq_len(nrow(key)))
  key_plot <- ggplot(key, aes(0, y, label = label, colour = highlight_group)) + geom_text(hjust = 0, size = 2.2) +
    xlim(0, 1) + theme_void() + theme(legend.position = "none") +
    scale_colour_manual(values = c("Divergent yet similar" = "#00A087", "Top overall similarity" = "#E64B35"))
  save_plot(tree_plot + heat + scatter + key_plot + patchwork::plot_layout(widths = c(3, 2.4, 4, 3)), "Figure_S6AB", 15, 6)
}
