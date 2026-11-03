# Robinson projection and Natural Earth polygon rendering. Boundaries are an external input.
suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(ggnewscale)
  library(jsonlite)
  library(scales)
})

.robinson_lat <- seq(0, 90, 5)
.robinson_x <- c(1.0000, 0.9986, 0.9954, 0.9900, 0.9822, 0.9730, 0.9600, 0.9427, 0.9216,
                  0.8962, 0.8679, 0.8350, 0.7986, 0.7597, 0.7186, 0.6732, 0.6213, 0.5722, 0.5322)
.robinson_y <- c(0.0000, 0.0620, 0.1240, 0.1860, 0.2480, 0.3100, 0.3720, 0.4340, 0.4958,
                  0.5571, 0.6176, 0.6769, 0.7346, 0.7903, 0.8435, 0.8936, 0.9394, 0.9761, 1.0000)

project_robinson <- function(lon, lat) {
  lat <- pmin(pmax(lat, -90), 90)
  abs_lat <- abs(lat)
  xcoef <- approx(.robinson_lat, .robinson_x, xout = abs_lat, rule = 2)$y
  ycoef <- approx(.robinson_lat, .robinson_y, xout = abs_lat, rule = 2)$y
  data.frame(x = 0.8487 * (lon * pi / 180) * xcoef, y = 1.3523 * sign(lat) * ycoef)
}

project_regional <- function(lon, lat, reference_lat) {
  data.frame(x = lon * cos(reference_lat * pi / 180), y = lat)
}

read_admin0_polygons <- function(geojson_path) {
  raw <- jsonlite::fromJSON(geojson_path, simplifyVector = FALSE)
  features <- raw$features
  rows <- vector("list", length(features) * 4L)
  n_rows <- 0L
  part_counter <- 0L
  feature_iso3 <- function(props) {
    for (field in c("ADM0_A3", "ISO_A3", "SOV_A3")) {
      value <- props[[field]]
      if (!is.null(value) && nchar(trimws(as.character(value))) == 3 && value != "-99") {
        return(as.character(value))
      }
    }
    ""
  }
  for (feature in features) {
    iso3 <- feature_iso3(feature$properties)
    geom <- feature$geometry
    if (is.null(geom)) next
    polys <- if (geom$type == "Polygon") list(geom$coordinates) else geom$coordinates
    for (poly in polys) {
      part_counter <- part_counter + 1L
      for (ring_idx in seq_along(poly)) {
        ring <- poly[[ring_idx]]
        if (length(ring) < 3) next
        mat <- do.call(rbind, lapply(ring, function(pt) c(pt[[1]], pt[[2]])))
        n_rows <- n_rows + 1L
        rows[[n_rows]] <- data.frame(
          iso3 = iso3, part_id = part_counter, ring_id = ring_idx,
          lon = mat[, 1], lat = mat[, 2]
        )
      }
    }
  }
  do.call(rbind, rows[seq_len(n_rows)])
}

.map_fill_layers <- function(data, palette_colours, no_data_col, boundary_col,
                              boundary_linewidth, vmin, vmax, show_guide) {
  guide_fill <- if (show_guide) {
    guide_colorbar(
      title = "Samples per country (log scale)", title.position = "bottom",
      title.hjust = 0.5, barwidth = unit(150, "pt"), barheight = unit(6, "pt"),
      order = 1
    )
  } else {
    "none"
  }
  guide_no_data <- if (show_guide) guide_legend(order = 2) else "none"
  list(
    geom_polygon(
      data = data,
      aes(x = x, y = y, group = part_id, subgroup = ring_id, fill = count_fill),
      colour = NA
    ),
    scale_fill_gradientn(
      colours = palette_colours, trans = "log10", limits = c(vmin, vmax),
      breaks = c(1, 10, 100, 1000, 10000), labels = label_comma(),
      na.value = no_data_col, oob = squish, guide = guide_fill
    ),
    new_scale_fill(),
    geom_polygon(
      data = data[which(is.na(data$count_fill))[1], , drop = FALSE],
      aes(x = x, y = y, group = part_id, fill = "No country-resolved samples"),
      colour = NA
    ),
    scale_fill_manual(
      values = c("No country-resolved samples" = no_data_col), name = NULL,
      guide = guide_no_data
    ),
    geom_polygon(
      data = data, aes(x = x, y = y, group = part_id, subgroup = ring_id),
      fill = NA, colour = boundary_col, linewidth = boundary_linewidth
    )
  )
}

build_country_choropleth_panel <- function(
  polys_raw, counts, title, subtitle, tag = NULL,
  font_family = "Arial", ink = "#202124", muted = "#5B6670",
  palette_colours = c("#FFFFFF", "#E8EFF4", "#C4D4E0", "#91ADC2", "#6485A6"),
  no_data_col = "#D9DCE1", boundary_col = "#FFFFFF",
  vmin = 1, vmax = 20000,
  inset_region = list(
    xlim = c(-13, 32), ylim = c(35, 72), reference_lat = 53, label = "Europe",
    box = c(left = 0.035, bottom = 0.06, right = 0.30, top = 0.46)
  )
) {
  count_map <- setNames(counts$count, counts$iso3)
  polys_raw$count <- unname(count_map[polys_raw$iso3])
  polys_raw$count_fill <- ifelse(!is.na(polys_raw$count) & polys_raw$count > 0, polys_raw$count, NA_real_)

  main_xy <- project_robinson(polys_raw$lon, polys_raw$lat)
  main_data <- cbind(polys_raw[, c("iso3", "part_id", "ring_id", "count_fill")], main_xy)
  xr <- c(min(main_data$x), max(main_data$x))

  p_map <- ggplot() +
    .map_fill_layers(main_data, palette_colours, no_data_col, boundary_col, 0.12, vmin, vmax, show_guide = TRUE) +
    coord_fixed(ratio = 1, xlim = xr, ylim = c(-1.37, 1.37), expand = FALSE, clip = "off") +
    labs(title = title, subtitle = subtitle) +
    theme_void(base_family = font_family) +
    theme(
      plot.title = element_text(size = 10, face = "bold", hjust = 0, colour = ink,
                                 family = font_family, margin = margin(b = 2)),
      plot.subtitle = element_text(size = 7.5, hjust = 0, colour = muted,
                                    family = font_family, margin = margin(b = 4)),
      plot.title.position = "plot",
      legend.position = "bottom",
      legend.title = element_text(size = 7.5, colour = ink, family = font_family),
      legend.text = element_text(size = 7, colour = ink, family = font_family),
      legend.margin = margin(t = -6),
      plot.margin = margin(4, 4, 2, 4)
    )
  if (!is.null(tag)) {
    p_map <- p_map + labs(tag = tag)
  }

  if (is.null(inset_region)) {
    return(p_map)
  }

  inset_xy <- project_regional(polys_raw$lon, polys_raw$lat, inset_region$reference_lat)
  inset_data <- cbind(polys_raw[, c("iso3", "part_id", "ring_id", "count_fill")], inset_xy)
  inset_xlim <- inset_region$xlim * cos(inset_region$reference_lat * pi / 180)

  p_inset <- ggplot() +
    .map_fill_layers(inset_data, palette_colours, no_data_col, boundary_col, 0.15, vmin, vmax, show_guide = FALSE) +
    annotate("text", x = inset_xlim[1] + 0.03 * diff(inset_xlim),
             y = inset_region$ylim[1] + 0.05 * diff(inset_region$ylim),
             label = inset_region$label, hjust = 0, vjust = 0, size = 6.5 / .pt,
             colour = muted, family = font_family) +
    coord_fixed(ratio = 1, xlim = inset_xlim, ylim = inset_region$ylim, expand = FALSE, clip = "on") +
    theme_void(base_family = font_family) +
    theme(
      legend.position = "none",
      panel.background = element_rect(fill = "white", colour = NA),
      panel.border = element_rect(fill = NA, colour = "#56616D", linewidth = 0.5),
      plot.margin = margin(0, 0, 0, 0)
    )

  box <- inset_region$box
  p_map + inset_element(
    p_inset, left = box[["left"]], bottom = box[["bottom"]],
    right = box[["right"]], top = box[["top"]], align_to = "panel",
    ignore_tag = TRUE
  )
}
