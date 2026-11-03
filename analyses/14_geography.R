# Figure S5. Country-level descriptive prevalence; no city-level substitution.

data_dir <- Sys.getenv("FP_DATA_DIR", "/path/to/private/data/14_geography")
out_dir <- Sys.getenv("FP_OUTPUT_DIR", "/path/to/private/results/14_geography")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
font_family <- "sans"
library(ggplot2)
library(patchwork)
library(ggnewscale)
library(jsonlite)
library(scales)

counts <- read.delim(file.path(data_dir, "world_counts.tsv"), check.names = FALSE, quote = "")
geojson_path <- file.path(data_dir, "world_boundaries.geojson")
raw <- jsonlite::fromJSON(geojson_path, simplifyVector = FALSE)
features <- raw$features
rows <- vector("list", length(features) * 4L)
n_rows <- 0L
part_counter <- 0L

for (feature in features) {
    {
        iso3 <- ""
        for (field in c("ADM0_A3", "ISO_A3", "SOV_A3")) {
            value <- feature$properties[[field]]
            if (!is.null(value) && nchar(trimws(as.character(value))) == 3 && value != "-99") {
                iso3 <- as.character(value)
                break
            }
        }
    }
    geom <- feature$geometry
    if (is.null(geom))
        next
    polys <- if (geom$type == "Polygon")
        list(geom$coordinates)
    else geom$coordinates
    for (poly in polys) {
        part_counter <- part_counter + 1L
        for (ring_idx in seq_along(poly)) {
            ring <- poly[[ring_idx]]
            if (length(ring) < 3)
                next
            mat <- do.call(rbind, lapply(ring, function(pt) c(pt[[1]], pt[[2]])))
            n_rows <- n_rows + 1L
            rows[[n_rows]] <- data.frame(iso3 = iso3, part_id = part_counter, ring_id = ring_idx, lon = mat[,
                1], lat = mat[, 2])
        }
    }
}

polygons <- do.call(rbind, rows[seq_len(n_rows)])
.robinson_lat <- seq(0, 90, 5)
.robinson_x <- c(1, 0.9986, 0.9954, 0.99, 0.9822, 0.973, 0.96, 0.9427, 0.9216, 0.8962, 0.8679, 0.835,
    0.7986, 0.7597, 0.7186, 0.6732, 0.6213, 0.5722, 0.5322)

.robinson_y <- c(0, 0.062, 0.124, 0.186, 0.248, 0.31, 0.372, 0.434, 0.4958, 0.5571, 0.6176, 0.6769, 0.7346,
    0.7903, 0.8435, 0.8936, 0.9394, 0.9761, 1)
polys_raw <- polygons
counts <- counts[c("iso3", "count")]
title <- "Metagenome stool samples by country"
subtitle <- ""
tag <- NULL
ink <- "#202124"
muted <- "#5B6670"
palette_colours <- c("#FFFFFF", "#E8EFF4", "#C4D4E0", "#91ADC2", "#6485A6")
no_data_col <- "#D9DCE1"
boundary_col <- "#FFFFFF"
vmin <- 1
vmax <- 20000

inset_region <- list(xlim = c(-13, 32), ylim = c(35, 72), reference_lat = 53, label = "Europe", box = c(left = 0.035,
    bottom = 0.06, right = 0.3, top = 0.46))
count_map <- setNames(counts$count, counts$iso3)
polys_raw$count <- unname(count_map[polys_raw$iso3])
polys_raw$count_fill <- ifelse(!is.na(polys_raw$count) & polys_raw$count > 0, polys_raw$count, NA_real_)
main_xy <- {
    lon <- polys_raw$lon
    lat <- polys_raw$lat
    lat <- pmin(pmax(lat, -90), 90)
    abs_lat <- abs(lat)
    xcoef <- approx(.robinson_lat, .robinson_x, xout = abs_lat, rule = 2)$y
    ycoef <- approx(.robinson_lat, .robinson_y, xout = abs_lat, rule = 2)$y
    data.frame(x = 0.8487 * (lon * pi/180) * xcoef, y = 1.3523 * sign(lat) * ycoef)
}

main_data <- cbind(polys_raw[, c("iso3", "part_id", "ring_id", "count_fill")], main_xy)
xr <- c(min(main_data$x), max(main_data$x))
p_map <- ggplot() + {
    data <- main_data
    palette_colours <- palette_colours
    no_data_col <- no_data_col
    boundary_col <- boundary_col
    boundary_linewidth <- 0.12
    vmin <- vmin
    vmax <- vmax
    show_guide <- TRUE
    guide_fill <- if (show_guide) {
        guide_colorbar(title = "Samples per country (log scale)", title.position = "bottom", title.hjust = 0.5,
            barwidth = unit(150, "pt"), barheight = unit(6, "pt"), order = 1)
    }
    else {
        "none"
    }
    guide_no_data <- if (show_guide)
        guide_legend(order = 2)
    else "none"
    list(geom_polygon(data = data, aes(x = x, y = y, group = part_id, subgroup = ring_id, fill = count_fill),
        colour = NA), scale_fill_gradientn(colours = palette_colours, trans = "log10", limits = c(vmin,
        vmax), breaks = c(1, 10, 100, 1000, 10000), labels = label_comma(), na.value = no_data_col, oob = squish,
        guide = guide_fill), new_scale_fill(), geom_polygon(data = data[which(is.na(data$count_fill))[1],
        , drop = FALSE], aes(x = x, y = y, group = part_id, fill = "No country-resolved samples"), colour = NA),
        scale_fill_manual(values = c(`No country-resolved samples` = no_data_col), name = NULL, guide = guide_no_data),
        geom_polygon(data = data, aes(x = x, y = y, group = part_id, subgroup = ring_id), fill = NA,
            colour = boundary_col, linewidth = boundary_linewidth))
} + coord_fixed(ratio = 1, xlim = xr, ylim = c(-1.37, 1.37), expand = FALSE, clip = "off") + labs(title = title,
    subtitle = subtitle, x = NULL, y = NULL, fill = NULL, colour = NULL) + theme_void(base_family = font_family) + theme(plot.title = element_text(size = 10,
    face = "bold", hjust = 0, colour = ink, family = font_family, margin = margin(b = 2)), plot.subtitle = element_text(size = 7.5,
    hjust = 0, colour = muted, family = font_family, margin = margin(b = 4)), plot.title.position = "plot",
    legend.position = "bottom", legend.title = element_text(size = 7.5, colour = ink, family = font_family),
    legend.text = element_text(size = 7, colour = ink, family = font_family), legend.margin = margin(t = -6),
    plot.margin = margin(4, 4, 2, 4))
if (!is.null(tag)) {
    p_map <- p_map + labs(tag = tag, x = NULL, y = NULL, fill = NULL, colour = NULL)
}

inset_xy <- {
    lon <- polys_raw$lon
    lat <- polys_raw$lat
    reference_lat <- inset_region$reference_lat
    data.frame(x = lon * cos(reference_lat * pi/180), y = lat)
}

inset_data <- cbind(polys_raw[, c("iso3", "part_id", "ring_id", "count_fill")], inset_xy)
inset_xlim <- inset_region$xlim * cos(inset_region$reference_lat * pi/180)
p_inset <- ggplot() + {
    data <- inset_data
    palette_colours <- palette_colours
    no_data_col <- no_data_col
    boundary_col <- boundary_col
    boundary_linewidth <- 0.15
    vmin <- vmin
    vmax <- vmax
    show_guide <- FALSE
    guide_fill <- if (show_guide) {
        guide_colorbar(title = "Samples per country (log scale)", title.position = "bottom", title.hjust = 0.5,
            barwidth = unit(150, "pt"), barheight = unit(6, "pt"), order = 1)
    }
    else {
        "none"
    }
    guide_no_data <- if (show_guide)
        guide_legend(order = 2)
    else "none"
    list(geom_polygon(data = data, aes(x = x, y = y, group = part_id, subgroup = ring_id, fill = count_fill),
        colour = NA), scale_fill_gradientn(colours = palette_colours, trans = "log10", limits = c(vmin,
        vmax), breaks = c(1, 10, 100, 1000, 10000), labels = label_comma(), na.value = no_data_col, oob = squish,
        guide = guide_fill), new_scale_fill(), geom_polygon(data = data[which(is.na(data$count_fill))[1],
        , drop = FALSE], aes(x = x, y = y, group = part_id, fill = "No country-resolved samples"), colour = NA),
        scale_fill_manual(values = c(`No country-resolved samples` = no_data_col), name = NULL, guide = guide_no_data),
        geom_polygon(data = data, aes(x = x, y = y, group = part_id, subgroup = ring_id), fill = NA,
            colour = boundary_col, linewidth = boundary_linewidth))
} + annotate("text", x = inset_xlim[1] + 0.03 * diff(inset_xlim), y = inset_region$ylim[1] + 0.05 * diff(inset_region$ylim),
    label = inset_region$label, hjust = 0, vjust = 0, size = 6.5/.pt, colour = muted, family = font_family) +
    coord_fixed(ratio = 1, xlim = inset_xlim, ylim = inset_region$ylim, expand = FALSE, clip = "on") +
    theme_void(base_family = font_family) + theme(legend.position = "none", panel.background = element_rect(fill = "white",
    colour = NA), panel.border = element_rect(fill = NA, colour = "#56616D", linewidth = 0.5), plot.margin = margin(0,
    0, 0, 0))
box <- inset_region$box

a <- p_map + inset_element(p_inset, left = box[["left"]], bottom = box[["bottom"]], right = box[["right"]],
    top = box[["top"]], align_to = "panel", ignore_tag = TRUE)

d <- read.delim(file.path(data_dir, "geography_country_plot.tsv"), check.names = FALSE, quote = "")

stopifnot(length(unique(d$country)) == 31, length(unique(d$species)) == 14, all(d$country_total_samples >=
    50), all(d$species_overall_prevalence > 0.05), max(abs(d$prevalence - d$detected_samples/d$country_total_samples)) <
    1e-08)
species <- unique(d[c("species", "species_phylogeny_order")])
species <- species[order(species$species_phylogeny_order), ]
d <- d[order(d$country_order), ]

d$country_label <- factor(paste0(d$country, " (n=", d$country_total_samples, ")"), levels = unique(paste0(d$country,
    " (n=", d$country_total_samples, ")")))

d$species_label <- factor(sub("^Faecalibacterium ", "F. ", sub("^Faecalibacterium prausnitzii_F$", "Faecalibacterium langellae",
    sub("^Faecalibacterium sp900539945$", "Faecalibacterium harmsenii", sub("^s__", "", d$species)))),
    levels = rev(sub("^Faecalibacterium ", "F. ", sub("^Faecalibacterium prausnitzii_F$", "Faecalibacterium langellae",
        sub("^Faecalibacterium sp900539945$", "Faecalibacterium harmsenii", sub("^s__", "", species$species))))))

d$continent <- factor(d$continent, levels = c("Africa", "Asia", "Europe", "North America", "South America",
    "Oceania"))

b <- ggplot(d, aes(country_label, species_label)) + geom_point(aes(size = prevalence * 100, fill = prevalence *
    100), shape = 21, colour = "#6F7378", stroke = 0.2) + scale_size_area(max_size = 3.7, limits = c(0,
    100), name = "Detection rate (%)") + scale_fill_gradient(low = "white", high = "#6485A6", limits = c(0,
    100), guide = "none") + facet_grid(~continent, scales = "free_x", space = "free_x") + scale_x_discrete(position = "top") +
    labs(x = NULL, y = NULL, fill = NULL, colour = NULL) + (ggplot2::theme_bw(base_size = 8, base_family = "sans") + ggplot2::theme(panel.grid.minor = ggplot2::element_blank(),
    legend.key = ggplot2::element_blank())) + theme(axis.text.x = element_text(angle = 65, hjust = 0),
    axis.text.y = element_text(face = "italic"), legend.position = "right")

ggplot2::ggsave(filename = file.path(out_dir, "Figure_S5.pdf"), plot = a/b + patchwork::plot_layout(heights = c(1.15,
    1)), width = 12, height = 10, device = grDevices::cairo_pdf, limitsize = FALSE)
