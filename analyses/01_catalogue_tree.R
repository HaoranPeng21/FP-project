# Figure 1A. Full 17,103-tip bac120 phylogeny with species, host and continent rings.

# 1. Set these two private directories (or edit the paths directly).
data_dir <- Sys.getenv("FP_DATA_DIR", "/path/to/private/data/01_catalogue_tree")
out_dir <- Sys.getenv("FP_OUTPUT_DIR", "/path/to/private/results/01_catalogue_tree")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
font_family <- "sans"

# 2. Packages.
library(ape)

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
tr <- ape::read.tree(file.path(data_dir, "bac120_tree.nwk"))
tr$tip.label <- sub("\\.fa(\\.tsv)?$", "", gsub("^['\"]|['\"]$", "", tr$tip.label))
stopifnot(!anyDuplicated(tr$tip.label))

meta <- read.csv(file.path(data_dir, "catalogue_tips.csv"), check.names = FALSE)

stopifnot(nrow(meta) == 17103, length(unique(meta$haoran_species)) == 83, !anyDuplicated(meta$label), all(meta$label %in% tr$tip.label))

tr <- ape::keep.tip(tr, meta$label)
meta$Continent <- gsub("_", " ", meta$Continent)
meta$Continent[meta$Continent %in% c("East Africa", "South Africa", "West Africa")] <- "Africa"
meta$Host1[tolower(meta$Host1) == "wastewater" & !is.na(meta$Host1)] <- "Wastewater"

grDevices::cairo_pdf(file.path(out_dir, "Figure_1A.pdf"), width = 11.7, height = 10, family = font_family)

plot_tr <- tr
plot_meta <- meta
plot_rings <- c("haoran_species_group", "Host1", "Continent")
plot_palettes <- list(species_palette, host1_palette, continent_palette)
plot_title <- ""
plot_subtitle <- ""
plot_legends <- TRUE
plot_meta <- plot_meta[match(plot_tr$tip.label, plot_meta$label), , drop = FALSE]
stopifnot(identical(plot_tr$tip.label, plot_meta$label), !anyDuplicated(plot_meta$label))
plot_tr <- ape::ladderize(plot_tr)
plot_meta <- plot_meta[match(plot_tr$tip.label, plot_meta$label), , drop = FALSE]
radius <- max(ape::node.depth.edgelength(plot_tr)[seq_along(plot_tr$tip.label)])
stopifnot(is.finite(radius), radius > 0)
plot_tr$edge.length <- plot_tr$edge.length/radius
old <- par(mar = c(1, 1, 4.5, if (plot_legends) 10 else 1), family = font_family, xpd = NA)

ape::plot.phylo(plot_tr, type = "fan", show.tip.label = FALSE, edge.width = 0.25, x.lim = c(-1.28, 1.28), y.lim = c(-1.28,
    1.28), no.margin = FALSE)
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
            "Faecalibacterium langellae", sub("^Faecalibacterium sp900539945$", "Faecalibacterium harmsenii", sub("^s__",
                "", used)))), fill = pal[used], bty = "n", cex = 0.6, title = heading, title.cex = 0.7, xjust = 0)
        legend_y <- legend_y - geometry$rect$h - 0.08
    }
}

if ("Type_isolate" %in% names(plot_meta)) {
    iso <- which(!is.na(plot_meta$Type_isolate) & plot_meta$Type_isolate == "Isolate")
    points(xy$xx[iso], xy$yy[iso], pch = 16, col = "#d62728", cex = 0.22)
}

graphics::title(main = plot_title, line = 2.2, cex.main = 0.85)
if (nzchar(plot_subtitle)) mtext(plot_subtitle, side = 3, line = 0.3, cex = 0.6)
par(old)
grDevices::dev.off()
