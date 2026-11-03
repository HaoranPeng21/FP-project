# Figure 2B. Colours are within-host species composition; numbers are genome counts.

# 1. Set these two private directories (or edit the paths directly).
data_dir <- Sys.getenv("FP_DATA_DIR", "/path/to/private/data/04_host_distribution")
out_dir <- Sys.getenv("FP_OUTPUT_DIR", "/path/to/private/results/04_host_distribution")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
font_family <- "sans"

# 2. Packages.
library(ComplexHeatmap)
library(circlize)
library(grid)

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

sp <- read.delim(file.path(data_dir, "species_scope.tsv"), check.names = FALSE, quote = "")
sp <- sp[order(sp$figure_order), ]

d <- read.delim(file.path(data_dir, "host_counts.tsv"), check.names = FALSE, quote = "")
counts <- t(as.matrix(d[sp$display_name]))
colnames(counts) <- d$host
stopifnot(sum(counts) == 15811, !anyNA(counts), all(colSums(counts) > 0))
percent <- sweep(counts, 2, colSums(counts), "/") * 100
stopifnot(max(abs(colSums(percent) - 100)) < 1e-08)

ht <- ComplexHeatmap::Heatmap(percent, name = "Within-host composition (%)", col = circlize::colorRamp2(c(0, 25,
    50, 75, 100), c("white", "#E9EFF5", "#C7D6E4", "#91ACC4", "#6485A6")), cluster_rows = FALSE, cluster_columns = FALSE,
    row_labels = sp$short_name, row_names_side = "left", top_annotation = ComplexHeatmap::HeatmapAnnotation(Host = d$host,
        col = list(Host = {
            values <- d$host
            palette <- host_palette
            missing <- sort(setdiff(unique(na.omit(values)), names(palette)))
            if (length(missing)) palette <- c(palette, setNames(grDevices::hcl.colors(length(missing), "Dark 3"),
                missing))
            palette
        }), show_annotation_name = FALSE, show_legend = FALSE), column_labels = paste0(d$host, " (n=", colSums(counts),
        ")"), column_names_rot = 45, row_names_gp = grid::gpar(fontsize = 8, fontface = "italic"), column_names_gp = grid::gpar(fontsize = 7),
    rect_gp = grid::gpar(col = "white"), cell_fun = function(j, i, x, y, w, h, fill) {
        if (counts[i, j] > 0)
            grid::grid.text(counts[i, j], x, y, gp = grid::gpar(fontsize = 7))
    })

grDevices::cairo_pdf(file.path(out_dir, "Figure_2B.pdf"), width = 7.5, height = 5.5, family = font_family)

ComplexHeatmap::draw(ht, heatmap_legend_side = "bottom")
grDevices::dev.off()
