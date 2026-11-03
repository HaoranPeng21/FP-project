# Figure S7. Maast trees, patristic-distance PERMANOVA, PERMDISP and Fitch permutation.
# Analysis includes known continents with >=5 genomes; BH is across four species per test.

# 1. Set these two private directories (or edit the paths directly).
data_dir <- Sys.getenv("FP_DATA_DIR", "/path/to/private/data/17_maast_tree_permanova")
out_dir <- Sys.getenv("FP_OUTPUT_DIR", "/path/to/private/results/17_maast_tree_permanova")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
font_family <- "sans"

# 2. Packages.
library(ape)
library(vegan)

# 3. Read the inputs, analyze, and draw the figures.
species_palette <- c(`s__Faecalibacterium prausnitzii` = "black", `s__Faecalibacterium tardum` = "#1f78b4", `s__Faecalibacterium duncaniae` = "#e41a1c",
    `s__Faecalibacterium taiwanense` = "#ff7f00", `s__Faecalibacterium prausnitzii_E` = "#4daf4a", `s__Faecalibacterium prausnitzii_F` = "#984ea3",
    `s__Faecalibacterium longum` = "#a65628", `s__Faecalibacterium hattorii` = "#66c2a5", `s__Faecalibacterium prausnitzii_I` = "#fc8d62",
    `s__Faecalibacterium prausnitzii_L` = "#e5c494", `s__Faecalibacterium prausnitzii_M` = "#b3de69", `s__Faecalibacterium wellingii` = "#D7BDE2",
    `s__Faecalibacterium butyricigenerans` = "#FDDB76", `s__Faecalibacterium gallinarum` = "#88CEE6", `s__Faecalibacterium intestinale` = "#E89DA0",
    `s__Faecalibacterium avium` = "#fccde5", `s__Faecalibacterium faecipullorum` = "#8dd3c7", `s__Faecalibacterium intestinigallinarum` = "#ffffb3",
    `s__Faecalibacterium gallistercoris` = "#bebada", `s__Faecalibacterium faecigallinarum` = "#fb8072", `s__Faecalibacterium intestinipullorum` = "#fdb462",
    `s__Faecalibacterium intestinavium` = "#bc80bd", `s__Faecalibacterium sp` = "#ccebc5", `s__Faecalibacterium sp900539945` = "#3C5488",
    `Candidate clades (14)` = "#bdbdbd")

host_palette <- c(`Human gut` = "#1b9e77", Wastewater = "#8dd3c7", Dog = "#1f78b4", Cat = "#87CEFA", Mongoose = "#ffff99",
    Rodent = "#e6ab02", Cattle = "#e31a1c", Goat = "#fb9a99", Primate = "#6a3d9a", Whale = "#5398CD", Pig = "#fdbf6f",
    Tamandua = "#b15928", Chicken = "#cab2d6", Raccoon = "#F8766D", `Sea Lion` = "#00BA38", Seal = "#619CFF")

continent_palette <- c(Africa = "#FFCE30", Antarctica = "#E64B35", Asia = "#d62728", Europe = "#1f77b4", `North America` = "#ff7f0e",
    `South America` = "#2ca02c", Oceania = "#9467bd")
host1_palette <- c(Animal = "#d95f02", Human = "#1b9e77", Wastewater = "#ffb000")

manifest <- read.csv(file.path(data_dir, "maast_tips.csv"), check.names = FALSE)
manifest$label <- sub("\\.fa(\\.tsv)?$", "", gsub("^['\"]|['\"]$", "", manifest$label))
manifest$Continent <- gsub("_", " ", trimws(manifest$Continent))
manifest$Continent[manifest$Continent %in% c("", "NA", "Unknown", "not provided")] <- NA_character_
manifest$Continent[manifest$Continent %in% c("South Africa", "West Africa")] <- "Africa"
species <- paste("s__Faecalibacterium", c("longum", "taiwanense", "duncaniae", "prausnitzii"))
nperm <- 1999L
results <- vector("list", 4)
plots <- vector("list", 4)

for (i in seq_along(species)) {
    message("  PERMANOVA/Fitch: ", sub("^Faecalibacterium ", "F. ", sub("^Faecalibacterium prausnitzii_F$", "Faecalibacterium langellae",
        sub("^Faecalibacterium sp900539945$", "Faecalibacterium harmsenii", sub("^s__", "", species[i])))))
    z <- manifest[manifest$Host1 == "Human" & manifest$Species == species[i], ]
    stopifnot(length(unique(z$tree_file)) == 1, !anyDuplicated(z$label))
    {
        tr <- ape::read.tree(file.path(data_dir, "maast_trees", unique(z$tree_file)))
        tr$tip.label <- sub("\\.fa(\\.tsv)?$", "", gsub("^['\"]|['\"]$", "", tr$tip.label))
        stopifnot(!anyDuplicated(tr$tip.label))
    }
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
    per <- vegan::adonis2(distance ~ Continent, data = meta, permutations = nperm)
    dispersion <- vegan::betadisper(distance, meta$Continent, type = "median")
    set.seed(8300 + i)
    dis <- vegan::permutest(dispersion, permutations = nperm, pairwise = FALSE)
    po <- ape::reorder.phylo(tree, "postorder")
    parents <- unique(po$edge[, 1])
    children <- split(po$edge[, 2], po$edge[, 1])
    state <- as.integer(2^(as.integer(meta$Continent) - 1))
    score <- function(tip_state) {
        masks <- integer(length(tree$tip.label) + tree$Nnode)
        masks[seq_along(tip_state)] <- tip_state
        changes <- 0L
        for (parent in parents) {
            kids <- children[[as.character(parent)]]
            current <- masks[kids[1]]
            if (length(kids) > 1)
                for (kid in kids[-1]) {
                  intersection <- bitwAnd(current, masks[kid])
                  if (intersection > 0)
                    current <- intersection
                  else {
                    current <- bitwOr(current, masks[kid])
                    changes <- changes + 1L
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
        betadisper_F = unname(dis$tab[1, "F"]), betadisper_p = unname(dis$tab[1, "Pr(>F)"]), fitch_observed = observed,
        fitch_null_mean = mean(null), fitch_clustering_strength = 1 - observed/mean(null), fitch_p = (1 + sum(null <=
            observed))/(nperm + 1), permutations = nperm)
    plots[[i]] <- list(tree = tr, metadata = z)
}

stats <- do.call(rbind, results)

for (test in c("permanova", "betadisper", "fitch")) stats[[paste0(test, "_q")]] <- p.adjust(stats[[paste0(test,
    "_p")]], "BH")
stopifnot(identical(as.integer(stats$n_analyzed), c(1251L, 631L, 503L, 203L)))

write.table(stats, file.path(out_dir, "continent_phylogenetic_tests.tsv"), sep = "\t", quote = FALSE, row.names = FALSE,
    na = "NA")

grDevices::cairo_pdf(file.path(out_dir, "Figure_S7.pdf"), width = 8.27, height = 11.69, family = font_family)

par(mfrow = c(2, 2), oma = c(4, 0, 0, 0))

for (i in seq_len(4)) {
    r <- stats[i, ]
    subtitle <- sprintf("n=%d; PERMANOVA R2=%.3f, q=%.4f\nFitch C=%.3f, q=%.4f; PERMDISP F=%.3f, q=%.4f", r$n_analyzed,
        r$permanova_R2, r$permanova_q, r$fitch_clustering_strength, r$fitch_q, r$betadisper_F, r$betadisper_q)
    {
        plot_tr <- plots[[i]]$tree
        plot_meta <- plots[[i]]$metadata
        plot_rings <- c("Continent", "Host1")
        plot_palettes <- list(continent_palette, host1_palette)
        plot_title <- sub("^Faecalibacterium ", "F. ", sub("^Faecalibacterium prausnitzii_F$", "Faecalibacterium langellae",
            sub("^Faecalibacterium sp900539945$", "Faecalibacterium harmsenii", sub("^s__", "", species[i]))))
        plot_subtitle <- subtitle
        plot_legends <- FALSE
        plot_meta <- plot_meta[match(plot_tr$tip.label, plot_meta$label), , drop = FALSE]
        stopifnot(identical(plot_tr$tip.label, plot_meta$label), !anyDuplicated(plot_meta$label))
        plot_tr <- ape::ladderize(plot_tr)
        plot_meta <- plot_meta[match(plot_tr$tip.label, plot_meta$label), , drop = FALSE]
        radius <- max(ape::node.depth.edgelength(plot_tr)[seq_along(plot_tr$tip.label)])
        stopifnot(is.finite(radius), radius > 0)
        plot_tr$edge.length <- plot_tr$edge.length/radius
        old <- par(mar = c(1, 1, 4.5, if (plot_legends) 10 else 1), family = font_family, xpd = NA)
        ape::plot.phylo(plot_tr, type = "fan", show.tip.label = FALSE, edge.width = 0.25, x.lim = c(-1.28, 1.28),
            y.lim = c(-1.28, 1.28), no.margin = FALSE)
        xy <- get("last_plot.phylo", envir = ape::.PlotPhyloEnv)
        n <- length(plot_tr$tip.label)
        theta <- atan2(xy$yy[seq_len(n)], xy$xx[seq_len(n)])
        idx <- order(theta)
        theta <- theta[idx]
        half <- pi/n
        legend_y <- 1.2
        for (j in seq_along(plot_rings)) {
            values <- as.character(plot_meta[[plot_rings[j]]])[idx]
            if (plot_rings[j] == "haoran_species_group")
                values[grepl("sp\\. Haoran_", values)] <- "Candidate clades (14)"
            pal <- {
                palette <- plot_palettes[[j]]
                missing <- sort(setdiff(unique(na.omit(values)), names(palette)))
                if (length(missing))
                  palette <- c(palette, setNames(grDevices::hcl.colors(length(missing), "Dark 3"), missing))
                palette
            }
            cols <- unname(pal[values])
            cols[is.na(cols)] <- "#eeeeee"
            r0 <- 1.01 + (j - 1) * 0.075
            r1 <- r0 + 0.06
            for (i in seq_len(n)) {
                a <- seq(theta[i] - half, theta[i] + half, length.out = 4)
                polygon(c(r0 * cos(a), r1 * cos(rev(a))), c(r0 * sin(a), r1 * sin(rev(a))), col = cols[i], border = NA)
            }
            if (plot_legends) {
                used <- names(pal)[names(pal) %in% values]
                heading <- if (plot_rings[j] == "haoran_species_group")
                  "Species"
                else if (plot_rings[j] == "Host1")
                  "Host"
                else plot_rings[j]
                geometry <- legend(1.32, legend_y, legend = sub("^Faecalibacterium ", "F. ", sub("^Faecalibacterium prausnitzii_F$",
                  "Faecalibacterium langellae", sub("^Faecalibacterium sp900539945$", "Faecalibacterium harmsenii",
                    sub("^s__", "", used)))), fill = pal[used], bty = "n", cex = 0.6, title = heading, title.cex = 0.7,
                  xjust = 0)
                legend_y <- legend_y - geometry$rect$h - 0.08
            }
        }
        if ("Type_isolate" %in% names(plot_meta)) {
            iso <- which(!is.na(plot_meta$Type_isolate) & plot_meta$Type_isolate == "Isolate")
            points(xy$xx[iso], xy$yy[iso], pch = 16, col = "#d62728", cex = 0.22)
        }
        graphics::title(main = plot_title, line = 2.2, cex.main = 0.85)
        if (nzchar(plot_subtitle))
            mtext(plot_subtitle, side = 3, line = 0.3, cex = 0.6)
        par(old)
    }
}

par(fig = c(0, 1, 0, 0.07), new = TRUE, mar = c(0, 0, 0, 0))
plot.new()

legend("center", legend = names(continent_palette)[names(continent_palette) != "Antarctica"], fill = continent_palette[names(continent_palette) !=
    "Antarctica"], ncol = 3, bty = "n", cex = 0.7)
grDevices::dev.off()
