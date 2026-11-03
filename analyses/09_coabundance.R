# Figure 4. Final tree/ring layout from cohort direction counts and REML results.
# Colour is the positive fraction among nominally significant cohorts.
# Thin frames: meta P<0.05. Thick frames: global BH q<0.05.
run <- function() {
  need("ape", "cowplot", "dplyr", "ggnewscale", "ggplot2", "ggsci", "ggtree", "grid", "readr", "scales", "tidyr", "tibble")
  read_tsv <- readr::read_tsv
  tree_file <- input("coabundance_tree")
  summary_tsv <- input("coabundance_ratios")
  meta_file <- input("coabundance_taxonomy")
  focal_order_file <- input("coabundance_order")
  meta_stats_file <- input("coabundance_meta")
  out_pdf <- file.path(out_dir, "Figure_4.pdf")
  out_png <- file.path(out_dir, "Figure_4.png")
  min_sig <- 1L
  min_partner_faec <- 2L
  highlight_min_faec <- 3L
  CANVAS_W <- 247.2
  CANVAS_H <- 247.2
  TREE_CX <- 123.60
  TREE_CY <- 123.60
  R_UNIT <- 66.693
  TREE_SIDE <- R_UNIT / 0.4
  CORE_R <- 33.035
  RING_R0 <- 34.905
  RING_R1 <- 62.330
  FAM_R0 <- 63.577
  FAM_R1 <- 66.693
  LAB_R <- 68.251

  DENDRO_GAP <- 1.247
  DENDRO_H <- 6.233
  COLLAB_GAP <- 1.247

  ROTATE_ANGLE <- 89.25

  FONT_PT <- 4.99
  TIP_FONT_PT <- 4.99
  FONT_TITLE_PT <- 9
  PT <- ggplot2::.pt

  POS_HEX <- "#D47463"
  NEG_HEX <- "#7B97B3"
  NEUTRAL_HEX <- "#A7A7A7"
  GRID_HEX <- "#D8D8D8"
  TILE_BORDER <- "grey45"
  OTHER_FAMILY_HEX <- "#D9D9D9"

  SIZE_FLOOR <- 1 / 3
  size_frac <- function(n, max_n) {
    n <- pmax(as.numeric(n), 1)
    max_n <- max(as.numeric(max_n), 1)
    if (max_n <= 1) return(rep(1, length(n)))
    SIZE_FLOOR + (1 - SIZE_FLOOR) * sqrt((n - 1) / (max_n - 1))
  }

  is_faec <- function(x) grepl("^s__Faecalibacterium ", as.character(x))

  short_partner <- function(x) {
    x <- sub("^s__", "", as.character(x))
    bits <- strsplit(x, " ", fixed = TRUE)
    vapply(bits, function(z) {
      if (length(z) < 2) return(z[[1]])
      if (grepl("^sp[0-9]", z[[2]])) return(z[[1]])
      z[[length(z)]]
    }, character(1))
  }

  abbreviate_genus <- function(x) {
    x <- sub("^s__", "", as.character(x))
    bits <- strsplit(x, " ", fixed = TRUE)
    labels <- vapply(bits, function(z) {
      if (length(z) < 2) return(z[[1]])
      paste0(substr(z[[1]], 1, 1), ".", paste(z[-1], collapse = " "))
    }, character(1))
    collisions <- duplicated(labels) | duplicated(labels, fromLast = TRUE)
    labels[collisions] <- vapply(bits[collisions], function(z) {
      paste0(substr(z[[1]], 1, 2), ".", paste(z[-1], collapse = " "))
    }, character(1))
    labels[x == "Bacteroides thetaiotaomicron"] <- "B.theta"
    labels
  }

  focal_order <- trimws(readLines(focal_order_file, warn = FALSE))
  focal_order <- focal_order[nzchar(focal_order) & !startsWith(focal_order, "#")]
  n_focal <- length(focal_order)
  if (n_focal != 13L || anyDuplicated(focal_order) || !all(is_faec(focal_order))) {
    stop("Expected 13 unique Faecalibacterium species in focal order file")
  }
  focal_short <- sub("^s__Faecalibacterium ", "F. ", focal_order)

  meta <- read_tsv(meta_file, col_names = c("leaf", "taxonomy"), show_col_types = FALSE) |>
    mutate(
      species = sub("^.*;s__", "s__", taxonomy),
      family = ifelse(
        grepl("(^|;)f__[^;]+", taxonomy),
        sub("^.*;f__([^;]+).*$", "\\1", taxonomy),
        "Unclassified"
      )
    )

  tree_all <- read.tree(tree_file)
  tree_all$tip.label <- sub("^'(.*)'$", "\\1", tree_all$tip.label)
  leaf_map <- setNames(meta$species, meta$leaf)
  tree_all$tip.label <- ifelse(
    tree_all$tip.label %in% names(leaf_map),
    leaf_map[tree_all$tip.label],
    tree_all$tip.label
  )
  if ("root" %in% tree_all$tip.label) tree_all <- drop.tip(tree_all, "root")
  all_bacterial_tips <- tree_all$tip.label

  dt <- read.delim(summary_tsv, sep = "\t", header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
  dt$otu_a <- as.character(dt$otu_a)
  dt$otu_b <- as.character(dt$otu_b)
  dt$ratio_pos_among_sig <- suppressWarnings(as.numeric(dt$ratio_pos_among_sig))
  dt$n_sig <- suppressWarnings(as.integer(dt$n_sig))
  dt <- dt[dt$n_sig >= min_sig & is.finite(dt$ratio_pos_among_sig), , drop = FALSE]
  dt <- dt[
    (dt$otu_a %in% focal_order | dt$otu_b %in% focal_order) &
      !(is_faec(dt$otu_a) & !(dt$otu_a %in% focal_order)) &
      !(is_faec(dt$otu_b) & !(dt$otu_b %in% focal_order)),
    , drop = FALSE
  ]

  MAX_N_SIG <- max(dt$n_sig, na.rm = TRUE)

  dt$edge_key <- paste(pmin(dt$otu_a, dt$otu_b), pmax(dt$otu_a, dt$otu_b), sep = "||")
  if (nzchar(meta_stats_file)) {
    if (!file.exists(meta_stats_file)) stop("Meta-analysis table not found: ", meta_stats_file)
    meta_stats <- read_tsv(meta_stats_file, show_col_types = FALSE) |>
      transmute(
        edge_key,
        meta_r = as.numeric(meta_r_reml),
        meta_p = as.numeric(meta_p_reml),
        meta_q = as.numeric(meta_q_bh_global),
        meta_symbol = as.character(meta_symbol),
        meta_star = as.character(meta_star),
        meta_display = as.character(meta_display),
        meta_n_cohorts = as.integer(n_cohorts_used),
        meta_i2 = as.numeric(i2_percent)
      ) |>
      distinct(edge_key, .keep_all = TRUE)
    dt <- dt |> left_join(meta_stats, by = "edge_key")
  } else {
    dt$meta_r <- NA_real_
    dt$meta_p <- NA_real_
    dt$meta_q <- NA_real_
    dt$meta_symbol <- NA_character_
    dt$meta_star <- NA_character_
    dt$meta_display <- NA_character_
    dt$meta_n_cohorts <- NA_integer_
    dt$meta_i2 <- NA_real_
  }
  dt <- dt |>
    mutate(
      meta_nominal_symbol = case_when(
        is.finite(meta_p) & meta_p < 0.05 & meta_r > 0 ~ "+",
        is.finite(meta_p) & meta_p < 0.05 & meta_r < 0 ~ "−",
        TRUE ~ ""
      )
    )

  partner_edges <- bind_rows(
    dt |>
      filter(otu_a %in% focal_order, !is_faec(otu_b)) |>
      transmute(focal = otu_a, partner = otu_b, ratio = ratio_pos_among_sig, n_sig,
                meta_r, meta_p, meta_q, meta_symbol, meta_nominal_symbol, meta_star, meta_display,
                meta_n_cohorts, meta_i2),
    dt |>
      filter(otu_b %in% focal_order, !is_faec(otu_a)) |>
      transmute(focal = otu_b, partner = otu_a, ratio = ratio_pos_among_sig, n_sig,
                meta_r, meta_p, meta_q, meta_symbol, meta_nominal_symbol, meta_star, meta_display,
                meta_n_cohorts, meta_i2)
  ) |>
    distinct(focal, partner, .keep_all = TRUE)

  partner_summary <- partner_edges |>
    group_by(partner) |>
    summarise(
      n_faec = n_distinct(focal),
      all_positive = all(abs(ratio - 1) < 1e-12),
      all_negative = all(abs(ratio) < 1e-12),
      min_ratio = min(ratio),
      max_ratio = max(ratio),
      .groups = "drop"
    ) |>
    filter(n_faec >= min_partner_faec) |>
    mutate(
      in_bac120_tree = partner %in% all_bacterial_tips,
      highlight = case_when(
        n_faec >= highlight_min_faec & all_positive ~ "Consistently positive",
        n_faec >= highlight_min_faec & all_negative ~ "Consistently negative",
        TRUE ~ "Other / mixed"
      ),
      display_label = short_partner(partner)
    ) |>
    arrange(desc(n_faec), partner)

  missing_from_tree <- partner_summary |> filter(!in_bac120_tree) |> pull(partner)
  if (length(missing_from_tree)) {
    message("Partner species absent from bac120 tree: ", paste(missing_from_tree, collapse = ", "))
  }

  keep_partners <- partner_summary |> filter(in_bac120_tree) |> pull(partner)
  tree <- keep.tip(tree_all, keep_partners)
  tip_order <- tree$tip.label
  partner_edges <- partner_edges |> filter(partner %in% tip_order)

  family_map <- meta |>
    distinct(species, family) |>
    group_by(species) |>
    summarise(family = first(family), .groups = "drop")

  family_df0 <- tibble(species = tip_order) |>
    left_join(family_map, by = "species") |>
    mutate(family = coalesce(family, "Unclassified"))
  family_counts_all <- family_df0 |>
    count(family, name = "n_partner") |>
    arrange(desc(n_partner), family)
  n_family <- nrow(family_counts_all)
  major_family_counts <- family_counts_all |> filter(n_partner > 3)
  major_family_levels <- major_family_counts$family
  n_major_family <- length(major_family_levels)
  family_levels <- c(major_family_levels, "Other")

  family_df <- family_df0 |>
    mutate(
      family_group = if_else(family %in% major_family_levels, family, "Other"),
      family_group = factor(family_group, levels = family_levels)
    )
  family_group_counts <- family_df |>
    count(family_group, name = "n_partner", .drop = FALSE) |>
    arrange(match(family_group, family_levels))

  colour_pool <- unique(toupper(substr(c(
    ggsci::pal_npg("nrc")(10),
    ggsci::pal_nejm("default")(8),
    ggsci::pal_lancet("lanonc")(9),
    ggsci::pal_jama("default")(7),
    ggsci::pal_aaas("default")(10),
    ggsci::pal_d3("category20")(20),
    ggsci::pal_futurama("planetexpress")(12),
    ggsci::pal_simpsons("springfield")(16),
    ggsci::pal_igv("default")(51)
  ), 1, 7)))
  colour_pool <- setdiff(colour_pool, c("#FFFFFF", "#FFFFFE"))
  if (length(colour_pool) < n_major_family) {
    stop("Colour pool too small for ", n_major_family, " displayed families")
  }
  family_colors <- c(
    setNames(colour_pool[seq_len(n_major_family)], major_family_levels),
    Other = OTHER_FAMILY_HEX
  )

  p0 <- ggtree(tree, layout = "fan", open.angle = 90, size = 0.18)
  tree_data <- p0$data
  D <- max(tree_data$x[tree_data$isTip], na.rm = TRUE)
  mm2x <- function(mm) mm * D / CORE_R
  x2mm <- function(x) x * CORE_R / D
  X_MAX <- mm2x(R_UNIT)

  tip_y <- tree_data |> filter(isTip) |> select(label, y)

  ring_edges <- seq(RING_R0, RING_R1, length.out = n_focal + 1L)
  ring_mid <- (head(ring_edges, -1) + tail(ring_edges, -1)) / 2
  ring_w <- diff(ring_edges)[1]

  ring_df <- tidyr::expand_grid(species = tip_order, faec_idx = seq_len(n_focal)) |>
    mutate(faec_row = focal_order[faec_idx]) |>
    left_join(partner_edges, by = c("species" = "partner", "faec_row" = "focal")) |>
    left_join(tip_y, by = c("species" = "label")) |>
    mutate(
      xpos = mm2x(ring_mid[faec_idx]), xw = mm2x(ring_w),
      size_f = size_frac(n_sig, MAX_N_SIG)
    )

  family_plot_df <- family_df |>
    left_join(tip_y, by = c("species" = "label")) |>
    mutate(xpos = mm2x((FAM_R0 + FAM_R1) / 2), xw = mm2x(FAM_R1 - FAM_R0))

  tip_df <- tibble(label = tip_order) |>
    left_join(partner_summary, by = c("label" = "partner")) |>
    mutate(
      plot_label = abbreviate_genus(label),
      highlight = factor(
        highlight,
        levels = c("Consistently positive", "Consistently negative", "Other / mixed")
      )
    )
  message("longest radial label: ", max(nchar(tip_df$plot_label)), " chars")

  p_tree <- p0 +
    scale_x_continuous(limits = c(0, X_MAX), expand = c(0, 0), oob = scales::oob_keep) +
    theme_tree() +
    theme(
      legend.position = "none",
      plot.margin = margin(0, 0, 0, 0),
      plot.background = element_blank(),
      panel.background = element_blank()
    )
  p_tree <- rotate_tree(p_tree, ROTATE_ANGLE)   # opening centred on vertical

  p_tree <- p_tree +
    geom_tile(
      data = ring_df,
      aes(x = xpos, y = y, width = xw, height = 1),
      inherit.aes = FALSE,
      fill = "white", color = GRID_HEX, linewidth = 0.08
    ) +
    geom_tile(
      data = ring_df |> filter(is.finite(ratio)),
      aes(x = xpos, y = y, width = xw * size_f, height = 1 * size_f, fill = ratio),
      inherit.aes = FALSE,
      color = TILE_BORDER, linewidth = 0.08
    ) +
    geom_tile(
      data = ring_df |> filter(is.finite(meta_p), meta_p < .05),
      aes(x = xpos, y = y, width = xw * size_f, height = size_f),
      inherit.aes = FALSE, fill = NA, colour = "black", linewidth = .12
    ) +
    geom_tile(
      data = ring_df |> filter(is.finite(meta_q), meta_q < .05),
      aes(x = xpos, y = y, width = xw * size_f, height = size_f),
      inherit.aes = FALSE, fill = NA, colour = "black", linewidth = .28
    ) +
    scale_fill_gradient2(
      low = NEG_HEX, mid = "white", high = POS_HEX,
      midpoint = 0.5, limits = c(0, 1), oob = squish, name = "Positive ratio"
    ) +
    ggnewscale::new_scale_fill() +
    geom_tile(
      data = family_plot_df,
      aes(x = xpos, y = y, width = xw, height = 1, fill = family_group),
      inherit.aes = FALSE,
      color = "white", linewidth = 0.05
    ) +
    scale_fill_manual(values = family_colors, drop = FALSE, guide = "none")

  p_tree <- p_tree %<+% tip_df +
    geom_tippoint(aes(color = highlight), size = 0.5, show.legend = FALSE) +
    geom_tiplab2(
      aes(label = plot_label, color = highlight),
      offset = mm2x(LAB_R) - D,
      align = TRUE, linetype = 0, hjust = 0,
      size = TIP_FONT_PT / PT, show.legend = FALSE
    ) +
    scale_color_manual(
      values = c(
        "Consistently positive" = POS_HEX,
        "Consistently negative" = NEG_HEX,
        "Other / mixed" = NEUTRAL_HEX
      ),
      na.value = NEUTRAL_HEX
    )

  tree_grob <- ggplotGrob(p_tree)
  tree_grob$layout$clip[tree_grob$layout$name == "panel"] <- "off"

  ff <- dt |>
    filter(otu_a %in% focal_order, otu_b %in% focal_order) |>
    mutate(pair_key = paste(pmin(otu_a, otu_b), pmax(otu_a, otu_b), sep = "||")) |>
    distinct(pair_key, .keep_all = TRUE)
  ff_ratio <- setNames(ff$ratio_pos_among_sig, ff$pair_key)
  ff_n_sig <- setNames(ff$n_sig, ff$pair_key)
  ff_meta_r <- setNames(ff$meta_r, ff$pair_key)
  ff_meta_p <- setNames(ff$meta_p, ff$pair_key)
  ff_meta_q <- setNames(ff$meta_q, ff$pair_key)
  ff_meta_symbol <- setNames(ff$meta_symbol, ff$pair_key)
  ff_meta_nominal_symbol <- setNames(ff$meta_nominal_symbol, ff$pair_key)

  HM_W <- n_focal * ring_w
  HM_X0 <- TREE_CX + RING_R0
  col_centre <- HM_X0 + (seq_len(n_focal) - 0.5) * ring_w
  row_lo <- TREE_CY + head(ring_edges, -1)
  row_hi <- TREE_CY + tail(ring_edges, -1)

  heat_df <- tidyr::expand_grid(row_i = seq_len(n_focal), col_i = seq_len(n_focal)) |>
    mutate(
      row_species = focal_order[row_i],
      col_species = focal_order[col_i],
      pair_key = paste(pmin(row_species, col_species), pmax(row_species, col_species), sep = "||"),
      ratio = ifelse(row_i == col_i, NA_real_, as.numeric(ff_ratio[pair_key])),
      n_sig = ifelse(row_i == col_i, NA_integer_, as.integer(ff_n_sig[pair_key])),
      meta_r = ifelse(row_i == col_i, NA_real_, as.numeric(ff_meta_r[pair_key])),
      meta_p = ifelse(row_i == col_i, NA_real_, as.numeric(ff_meta_p[pair_key])),
      meta_q = ifelse(row_i == col_i, NA_real_, as.numeric(ff_meta_q[pair_key])),
      meta_symbol = ifelse(row_i == col_i, NA_character_, as.character(ff_meta_symbol[pair_key])),
      meta_nominal_symbol = ifelse(
        row_i == col_i, NA_character_, as.character(ff_meta_nominal_symbol[pair_key])
      ),
      xmin = HM_X0 + (col_i - 1) * ring_w,
      xmax = HM_X0 + col_i * ring_w,
      ymin = row_lo[row_i],
      ymax = row_hi[row_i],
      size_f = size_frac(n_sig, MAX_N_SIG),
      cx = (xmin + xmax) / 2, cy = (ymin + ymax) / 2,
      sxmin = cx - (xmax - xmin) / 2 * size_f, sxmax = cx + (xmax - xmin) / 2 * size_f,
      symin = cy - (ymax - ymin) / 2 * size_f, symax = cy + (ymax - ymin) / 2 * size_f
    )

  tree_focal <- keep.tip(tree_all, focal_order)
  order_is_phylo <- local({
    pos <- setNames(seq_along(focal_order), focal_order)
    pp <- ape::prop.part(tree_focal)
    labs <- attr(pp, "labels")
    all(vapply(pp, function(ix) all(diff(sort(pos[labs[ix]])) == 1), logical(1)))
  })
  if (!order_is_phylo) {
    warning("Ring order is not a rotation of the bac120 topology; dendrogram omitted")
  }

  dd <- ggtree(tree_focal, branch.length = "none")$data
  dd$y <- NA_real_
  dd$y[dd$isTip] <- match(dd$label[dd$isTip], focal_order)
  repeat {
    todo <- which(is.na(dd$y))
    if (!length(todo)) break
    progressed <- FALSE
    for (i in todo) {
      kids <- dd$y[dd$parent == dd$node[i] & dd$node != dd$node[i]]
      if (length(kids) && !anyNA(kids)) {
        dd$y[i] <- mean(range(kids))
        progressed <- TRUE
      }
    }
    if (!progressed) stop("Could not resolve dendrogram node coordinates")
  }

  DEN_Y0 <- TREE_CY + RING_R1 + DENDRO_GAP
  dmax <- max(dd$x)
  den_seg <- dd |>
    filter(node != parent) |>
    left_join(dd |> select(node, xp = x, yp = y), by = c("parent" = "node")) |>
    mutate(
      cx_node = HM_X0 + (y - 0.5) * ring_w,
      cx_par = HM_X0 + (yp - 0.5) * ring_w,
      cy_node = DEN_Y0 + (1 - x / dmax) * DENDRO_H,
      cy_par = DEN_Y0 + (1 - xp / dmax) * DENDRO_H
    )
  if (!order_is_phylo) den_seg <- den_seg[0, ]

  COLLAB_Y <- DEN_Y0 + DENDRO_H + COLLAB_GAP

  lab_df <- tibble(
    x = HM_X0 - 1.2,
    y = (row_lo + row_hi) / 2,
    label = focal_short
  )

  col_df <- tibble(x = col_centre, y = COLLAB_Y, label = focal_short)

  overlay <- ggplot() +
    geom_rect(
      data = heat_df, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
      fill = "white", color = GRID_HEX, linewidth = 0.16
    ) +
    geom_rect(
      data = heat_df |> filter(is.finite(ratio)),
      aes(xmin = sxmin, xmax = sxmax, ymin = symin, ymax = symax, fill = ratio),
      color = TILE_BORDER, linewidth = 0.16
    ) +
    geom_rect(
      data = heat_df |> filter(is.finite(meta_p), meta_p < .05),
      aes(xmin = sxmin, xmax = sxmax, ymin = symin, ymax = symax),
      fill = NA, colour = "black", linewidth = .2
    ) +
    geom_rect(
      data = heat_df |> filter(is.finite(meta_q), meta_q < .05),
      aes(xmin = sxmin, xmax = sxmax, ymin = symin, ymax = symax),
      fill = NA, colour = "black", linewidth = .4
    ) +
    geom_text(
      data = heat_df |> filter(!is.na(meta_nominal_symbol), meta_nominal_symbol != ""),
      aes(cx, cy, label = meta_nominal_symbol), size = FONT_PT / PT, colour = "black"
    ) +
    scale_fill_gradient2(
      low = NEG_HEX, mid = "white", high = POS_HEX,
      midpoint = 0.5, limits = c(0, 1), oob = squish, guide = "none"
    ) +
    geom_segment(
      data = lab_df,
      aes(x = x + 0.1, xend = HM_X0, y = y, yend = y),
      color = "grey55", linewidth = 0.2
    ) +
    geom_text(
      data = lab_df, aes(x = x - 0.5, y = y, label = label),
      hjust = 1, vjust = 0.5, size = FONT_PT / PT, color = "grey15"
    ) +
    geom_segment(
      data = den_seg,
      aes(x = cx_par, xend = cx_node, y = cy_node, yend = cy_node),
      color = "grey35", linewidth = 0.22
    ) +
    geom_segment(
      data = den_seg,
      aes(x = cx_par, xend = cx_par, y = cy_par, yend = cy_node),
      color = "grey35", linewidth = 0.22
    ) +
    geom_text(
      data = col_df, aes(x = x, y = y, label = label),
      angle = 90, hjust = 0, vjust = 0.5, size = FONT_PT / PT, color = "grey15"
    ) +
    scale_x_continuous(limits = c(0, CANVAS_W), expand = c(0, 0)) +
    scale_y_continuous(limits = c(0, CANVAS_H), expand = c(0, 0)) +
    theme_void() +
    theme(
      plot.margin = margin(0, 0, 0, 0),
      plot.background = element_blank(),
      panel.background = element_blank(),
      legend.position = "none"
    )

  legend_theme <- theme(
    legend.title = element_text(size = FONT_TITLE_PT, color = "grey15"),
    legend.text = element_text(size = FONT_PT, color = "grey15"),
    legend.background = element_blank(),
    legend.key = element_blank(),
    legend.margin = margin(0, 0, 0, 0),
    legend.justification = "left"
  )

  grab_legend <- function(p) {
    g <- cowplot::get_plot_component(p, "guide-box-right", return_all = TRUE)
    if (is.list(g) && !inherits(g, "grob")) g <- g[[1]]
    g
  }

  p_leg_ratio <- ggplot(tibble(x = 1, y = 1, v = 0.5), aes(x, y, fill = v)) +
    geom_tile() +
    scale_fill_gradient2(
      low = NEG_HEX, mid = "white", high = POS_HEX, midpoint = 0.5,
      limits = c(0, 1), breaks = c(0, 0.25, 0.5, 0.75, 1),
      name = "Positive ratio among\nnominally significant cohorts",
      guide = guide_colorbar(
        title.position = "top", direction = "horizontal",
        barwidth = unit(38, "mm"), barheight = unit(3.4, "mm"),
        frame.colour = "grey40", frame.linewidth = 0.2, ticks.colour = "grey40"
      )
    ) +
    theme_void() + legend_theme + theme(legend.position = "right")

  p_leg_dir <- ggplot(
    tibble(k = factor(c("Consistently positive", "Consistently negative", "Other / mixed"),
                      levels = c("Consistently positive", "Consistently negative", "Other / mixed"))),
    aes(x = 1, y = k, color = k)
  ) +
    geom_point(size = 2.2) +
    scale_color_manual(
      values = c("Consistently positive" = POS_HEX,
                 "Consistently negative" = NEG_HEX,
                 "Other / mixed" = NEUTRAL_HEX),
      labels = c("All positive", "All negative", "Mixed / other"),
      name = paste0("Partner direction (>=", highlight_min_faec, " Faecalibacterium)"),
      guide = guide_legend(title.position = "top", ncol = 1, override.aes = list(size = 2.2))
    ) +
    theme_void() + legend_theme + theme(legend.position = "right")

  fam_legend_labels <- paste0(family_levels, " (", family_group_counts$n_partner, ")")
  p_leg_fam <- ggplot(
    tibble(f = factor(family_levels, levels = family_levels)),
    aes(x = 1, y = f, fill = f)
  ) +
    geom_tile() +
    scale_fill_manual(
      values = family_colors, drop = FALSE, labels = fam_legend_labels,
      name = "Family of partner species (>3 species shown)",
      guide = guide_legend(
        title.position = "top", ncol = 2, byrow = FALSE,
        override.aes = list(color = NA), keywidth = unit(3.2, "mm"),
        keyheight = unit(3.2, "mm")
      )
    ) +
    theme_void() + legend_theme + theme(legend.position = "right")

  p_leg_size <- local({
    breaks_n <- sort(unique(c(1L, round(MAX_N_SIG / 3), round(2 * MAX_N_SIG / 3), MAX_N_SIG)))
    df <- tibble(
      x = seq_along(breaks_n),
      lab = as.character(breaks_n),
      pt_size = 4.4 * size_frac(breaks_n, MAX_N_SIG)
    )
    ggplot(df, aes(x = x, y = 1)) +
      geom_point(aes(size = pt_size), shape = 15, color = "grey50") +
      scale_size_identity() +
      geom_text(aes(label = lab), y = 0.35, size = FONT_PT / PT, color = "grey15") +
      annotate("text", x = 0.4, y = 1.9, label = "Cohorts significant (n_sig)",
               hjust = 0, size = FONT_TITLE_PT / PT, color = "grey15") +
      scale_x_continuous(limits = c(0, nrow(df) + 1)) +
      scale_y_continuous(limits = c(0, 2.2)) +
      theme_void()
  })

  leg_ratio <- grab_legend(p_leg_ratio)
  leg_dir <- grab_legend(p_leg_dir)
  leg_fam <- grab_legend(p_leg_fam)

  final_plot <- cowplot::ggdraw() +
    theme(plot.background = element_rect(fill = NA, color = NA)) +
    cowplot::draw_grob(
      tree_grob,
      x = (TREE_CX - TREE_SIDE / 2) / CANVAS_W,
      y = (TREE_CY - TREE_SIDE / 2) / CANVAS_H,
      width = TREE_SIDE / CANVAS_W,
      height = TREE_SIDE / CANVAS_H,
      clip = "off"
    ) +
    cowplot::draw_plot(overlay, x = 0, y = 0, width = 1, height = 1) +
    cowplot::draw_grob(leg_fam, x = .49, y = .86, width = .49, height = .13) +
    cowplot::draw_plot(p_leg_size, x = .02, y = .02, width = .34, height = .08) +
    cowplot::draw_grob(leg_dir, x = .02, y = .10, width = .43, height = .075) +
    cowplot::draw_grob(leg_ratio, x = .69, y = .035, width = .29, height = .08) +
    cowplot::draw_label("Thin frame: meta P<0.05; thick frame: global BH q<0.05", x = .75, y = .02, size = 5)

  dir.create(dirname(out_pdf), recursive = TRUE, showWarnings = FALSE)
  ggsave(out_pdf, final_plot, width = CANVAS_W / 25.4, height = CANVAS_H / 25.4,
         device = cairo_pdf, bg = "transparent")
  if (Sys.getenv("FP_PNG", "1") == "1") ggsave(out_png, final_plot,
    width = CANVAS_W / 25.4, height = CANVAS_H / 25.4, dpi = 300, bg = "white")
}
