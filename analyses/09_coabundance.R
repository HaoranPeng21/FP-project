# Figure 4. Final tree/ring layout from cohort direction counts and REML results.
# Colour is the positive fraction among nominally significant cohorts.
# Thin frames: meta P<0.05. Thick frames: global BH q<0.05.

# 1. Set these two private directories (or edit the paths directly).
data_dir <- Sys.getenv("FP_DATA_DIR", "/path/to/private/data/09_coabundance")
out_dir <- Sys.getenv("FP_OUTPUT_DIR", "/path/to/private/results/09_coabundance")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
font_family <- "sans"

# 2. Packages.
library(ape)
library(cowplot)
library(dplyr)
library(ggnewscale)
library(ggplot2)
library(ggsci)
library(ggtree)
library(grid)
library(readr)
library(scales)
library(tidyr)
library(tibble)

# 3. Read the inputs, analyze, and draw the figures.
read_tsv <- readr::read_tsv
tree_file <- file.path(data_dir, "coabundance_tree.nwk")
summary_tsv <- file.path(data_dir, "coabundance_ratios.tsv")
meta_file <- file.path(data_dir, "coabundance_taxonomy.txt")
focal_order_file <- file.path(data_dir, "coabundance_order.txt")
meta_stats_file <- file.path(data_dir, "coabundance_meta.tsv")
out_pdf <- file.path(out_dir, "Figure_4.pdf")
out_png <- file.path(out_dir, "Figure_4.png")
min_sig <- 1L
min_partner_faec <- 2L
highlight_min_faec <- 3L
CANVAS_W <- 247.2
CANVAS_H <- 247.2
TREE_CX <- 123.6
TREE_CY <- 123.6
R_UNIT <- 66.693
TREE_SIDE <- R_UNIT/0.4
CORE_R <- 33.035
RING_R0 <- 34.905
RING_R1 <- 62.33
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
SIZE_FLOOR <- 1/3
focal_order <- trimws(readLines(focal_order_file, warn = FALSE))
focal_order <- focal_order[nzchar(focal_order) & !startsWith(focal_order, "#")]
n_focal <- length(focal_order)

if (n_focal != 13L || anyDuplicated(focal_order) || !all(grepl("^s__Faecalibacterium ", as.character(focal_order)))) {
    stop("Expected 13 unique Faecalibacterium species in focal order file")
}

focal_short <- sub("^s__Faecalibacterium ", "F. ", focal_order)

meta <- mutate(read_tsv(meta_file, col_names = c("leaf", "taxonomy"), show_col_types = FALSE), species = sub("^.*;s__",
    "s__", taxonomy), family = ifelse(grepl("(^|;)f__[^;]+", taxonomy), sub("^.*;f__([^;]+).*$", "\\1", taxonomy),
    "Unclassified"))
tree_all <- read.tree(tree_file)
tree_all$tip.label <- sub("^'(.*)'$", "\\1", tree_all$tip.label)
leaf_map <- setNames(meta$species, meta$leaf)

tree_all$tip.label <- ifelse(tree_all$tip.label %in% names(leaf_map), leaf_map[tree_all$tip.label], tree_all$tip.label)

if ("root" %in% tree_all$tip.label) tree_all <- drop.tip(tree_all, "root")
all_bacterial_tips <- tree_all$tip.label

dt <- read.delim(summary_tsv, sep = "\t", header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
dt$otu_a <- as.character(dt$otu_a)
dt$otu_b <- as.character(dt$otu_b)
dt$ratio_pos_among_sig <- suppressWarnings(as.numeric(dt$ratio_pos_among_sig))
dt$n_sig <- suppressWarnings(as.integer(dt$n_sig))
dt <- dt[dt$n_sig >= min_sig & is.finite(dt$ratio_pos_among_sig), , drop = FALSE]

dt <- dt[(dt$otu_a %in% focal_order | dt$otu_b %in% focal_order) & !(grepl("^s__Faecalibacterium ", as.character(dt$otu_a)) &
    !(dt$otu_a %in% focal_order)) & !(grepl("^s__Faecalibacterium ", as.character(dt$otu_b)) & !(dt$otu_b %in%
    focal_order)), , drop = FALSE]
MAX_N_SIG <- max(dt$n_sig, na.rm = TRUE)
dt$edge_key <- paste(pmin(dt$otu_a, dt$otu_b), pmax(dt$otu_a, dt$otu_b), sep = "||")
if (nzchar(meta_stats_file)) {
    if (!file.exists(meta_stats_file))
        stop("Meta-analysis table not found: ", meta_stats_file)
    meta_stats <- distinct(transmute(read_tsv(meta_stats_file, show_col_types = FALSE), edge_key, meta_r = as.numeric(meta_r_reml),
        meta_p = as.numeric(meta_p_reml), meta_q = as.numeric(meta_q_bh_global), meta_symbol = as.character(meta_symbol),
        meta_star = as.character(meta_star), meta_display = as.character(meta_display), meta_n_cohorts = as.integer(n_cohorts_used),
        meta_i2 = as.numeric(i2_percent)), edge_key, .keep_all = TRUE)
    dt <- left_join(dt, meta_stats, by = "edge_key")
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

dt <- mutate(dt, meta_nominal_symbol = case_when(is.finite(meta_p) & meta_p < 0.05 & meta_r > 0 ~ "+", is.finite(meta_p) &
    meta_p < 0.05 & meta_r < 0 ~ "−", TRUE ~ ""))

partner_edges <- distinct(bind_rows(transmute(filter(dt, otu_a %in% focal_order, !grepl("^s__Faecalibacterium ",
    as.character(otu_b))), focal = otu_a, partner = otu_b, ratio = ratio_pos_among_sig, n_sig, meta_r, meta_p,
    meta_q, meta_symbol, meta_nominal_symbol, meta_star, meta_display, meta_n_cohorts, meta_i2), transmute(filter(dt,
    otu_b %in% focal_order, !grepl("^s__Faecalibacterium ", as.character(otu_a))), focal = otu_b, partner = otu_a,
    ratio = ratio_pos_among_sig, n_sig, meta_r, meta_p, meta_q, meta_symbol, meta_nominal_symbol, meta_star, meta_display,
    meta_n_cohorts, meta_i2)), focal, partner, .keep_all = TRUE)

partner_summary <- arrange(mutate(filter(summarise(group_by(partner_edges, partner), n_faec = n_distinct(focal),
    all_positive = all(abs(ratio - 1) < 1e-12), all_negative = all(abs(ratio) < 1e-12), min_ratio = min(ratio),
    max_ratio = max(ratio), .groups = "drop"), n_faec >= min_partner_faec), in_bac120_tree = partner %in% all_bacterial_tips,
    highlight = case_when(n_faec >= highlight_min_faec & all_positive ~ "Consistently positive", n_faec >= highlight_min_faec &
        all_negative ~ "Consistently negative", TRUE ~ "Other / mixed"), display_label = {
        x <- partner
        x <- sub("^s__", "", as.character(x))
        bits <- strsplit(x, " ", fixed = TRUE)
        vapply(bits, function(z) {
            if (length(z) < 2)
                return(z[[1]])
            if (grepl("^sp[0-9]", z[[2]]))
                return(z[[1]])
            z[[length(z)]]
        }, character(1))
    }), desc(n_faec), partner)
missing_from_tree <- pull(filter(partner_summary, !in_bac120_tree), partner)
if (length(missing_from_tree)) {
    message("Partner species absent from bac120 tree: ", paste(missing_from_tree, collapse = ", "))
}

keep_partners <- pull(filter(partner_summary, in_bac120_tree), partner)
tree <- keep.tip(tree_all, keep_partners)
tip_order <- tree$tip.label
partner_edges <- filter(partner_edges, partner %in% tip_order)

family_map <- summarise(group_by(distinct(meta, species, family), species), family = first(family), .groups = "drop")

family_df0 <- mutate(left_join(tibble(species = tip_order), family_map, by = "species"), family = coalesce(family,
    "Unclassified"))
family_counts_all <- arrange(count(family_df0, family, name = "n_partner"), desc(n_partner), family)
n_family <- nrow(family_counts_all)
major_family_counts <- filter(family_counts_all, n_partner > 3)
major_family_levels <- major_family_counts$family
n_major_family <- length(major_family_levels)
family_levels <- c(major_family_levels, "Other")

family_df <- mutate(family_df0, family_group = if_else(family %in% major_family_levels, family, "Other"), family_group = factor(family_group,
    levels = family_levels))

family_group_counts <- arrange(count(family_df, family_group, name = "n_partner", .drop = FALSE), match(family_group,
    family_levels))

colour_pool <- unique(toupper(substr(c((ggsci::pal_npg("nrc"))(10), (ggsci::pal_nejm("default"))(8), (ggsci::pal_lancet("lanonc"))(9),
    (ggsci::pal_jama("default"))(7), (ggsci::pal_aaas("default"))(10), (ggsci::pal_d3("category20"))(20), (ggsci::pal_futurama("planetexpress"))(12),
    (ggsci::pal_simpsons("springfield"))(16), (ggsci::pal_igv("default"))(51)), 1, 7)))
colour_pool <- setdiff(colour_pool, c("#FFFFFF", "#FFFFFE"))
if (length(colour_pool) < n_major_family) {
    stop("Colour pool too small for ", n_major_family, " displayed families")
}

family_colors <- c(setNames(colour_pool[seq_len(n_major_family)], major_family_levels), Other = OTHER_FAMILY_HEX)

p0 <- ggtree(tree, layout = "fan", open.angle = 90, size = 0.18)
tree_data <- p0$data
D <- max(tree_data$x[tree_data$isTip], na.rm = TRUE)
X_MAX <- R_UNIT * D/CORE_R
tip_y <- select(filter(tree_data, isTip), label, y)
ring_edges <- seq(RING_R0, RING_R1, length.out = n_focal + 1L)
ring_mid <- (head(ring_edges, -1) + tail(ring_edges, -1))/2
ring_w <- diff(ring_edges)[1]

ring_df <- mutate(left_join(left_join(mutate(tidyr::expand_grid(species = tip_order, faec_idx = seq_len(n_focal)),
    faec_row = focal_order[faec_idx]), partner_edges, by = c(species = "partner", faec_row = "focal")), tip_y,
    by = c(species = "label")), xpos = ring_mid[faec_idx] * D/CORE_R, xw = ring_w * D/CORE_R, size_f = if (max(MAX_N_SIG,
    1) <= 1) rep(1, length(n_sig)) else 1/3 + 2/3 * sqrt((pmax(as.numeric(n_sig), 1) - 1)/(max(as.numeric(MAX_N_SIG),
    1) - 1)))

family_plot_df <- mutate(left_join(family_df, tip_y, by = c(species = "label")), xpos = (FAM_R0 + FAM_R1)/2 * D/CORE_R,
    xw = (FAM_R1 - FAM_R0) * D/CORE_R)

tip_df <- mutate(left_join(tibble(label = tip_order), partner_summary, by = c(label = "partner")), plot_label = {
    x <- label
    x <- sub("^s__", "", as.character(x))
    bits <- strsplit(x, " ", fixed = TRUE)
    labels <- vapply(bits, function(z) {
        if (length(z) < 2)
            return(z[[1]])
        paste0(substr(z[[1]], 1, 1), ".", paste(z[-1], collapse = " "))
    }, character(1))
    collisions <- duplicated(labels) | duplicated(labels, fromLast = TRUE)
    labels[collisions] <- vapply(bits[collisions], function(z) {
        paste0(substr(z[[1]], 1, 2), ".", paste(z[-1], collapse = " "))
    }, character(1))
    labels[x == "Bacteroides thetaiotaomicron"] <- "B.theta"
    labels
}, highlight = factor(highlight, levels = c("Consistently positive", "Consistently negative", "Other / mixed")))

message("longest radial label: ", max(nchar(tip_df$plot_label)), " chars")

p_tree <- p0 + scale_x_continuous(limits = c(0, X_MAX), expand = c(0, 0), oob = scales::oob_keep) + theme_tree() +
    theme(legend.position = "none", plot.margin = margin(0, 0, 0, 0), plot.background = element_blank(), panel.background = element_blank())

p_tree <- rotate_tree(p_tree, ROTATE_ANGLE)

p_tree <- p_tree + geom_tile(data = ring_df, aes(x = xpos, y = y, width = xw, height = 1), inherit.aes = FALSE,
    fill = "white", color = GRID_HEX, linewidth = 0.08) + geom_tile(data = filter(ring_df, is.finite(ratio)), aes(x = xpos,
    y = y, width = xw * size_f, height = 1 * size_f, fill = ratio), inherit.aes = FALSE, color = TILE_BORDER, linewidth = 0.08) +
    geom_tile(data = filter(ring_df, is.finite(meta_p), meta_p < 0.05), aes(x = xpos, y = y, width = xw * size_f,
        height = size_f), inherit.aes = FALSE, fill = NA, colour = "black", linewidth = 0.12) + geom_tile(data = filter(ring_df,
    is.finite(meta_q), meta_q < 0.05), aes(x = xpos, y = y, width = xw * size_f, height = size_f), inherit.aes = FALSE,
    fill = NA, colour = "black", linewidth = 0.28) + scale_fill_gradient2(low = NEG_HEX, mid = "white", high = POS_HEX,
    midpoint = 0.5, limits = c(0, 1), oob = squish, name = "Positive ratio") + ggnewscale::new_scale_fill() + geom_tile(data = family_plot_df,
    aes(x = xpos, y = y, width = xw, height = 1, fill = family_group), inherit.aes = FALSE, color = "white", linewidth = 0.05) +
    scale_fill_manual(values = family_colors, drop = FALSE, guide = "none")

p_tree <- p_tree %<+% tip_df + geom_tippoint(aes(color = highlight), size = 0.5, show.legend = FALSE) + geom_tiplab2(aes(label = plot_label,
    color = highlight), offset = LAB_R * D/CORE_R - D, align = TRUE, linetype = 0, hjust = 0, size = TIP_FONT_PT/PT,
    show.legend = FALSE) + scale_color_manual(values = c(`Consistently positive` = POS_HEX, `Consistently negative` = NEG_HEX,
    `Other / mixed` = NEUTRAL_HEX), na.value = NEUTRAL_HEX)
tree_grob <- ggplotGrob(p_tree)
tree_grob$layout$clip[tree_grob$layout$name == "panel"] <- "off"

ff <- distinct(mutate(filter(dt, otu_a %in% focal_order, otu_b %in% focal_order), pair_key = paste(pmin(otu_a,
    otu_b), pmax(otu_a, otu_b), sep = "||")), pair_key, .keep_all = TRUE)
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

heat_df <- mutate(tidyr::expand_grid(row_i = seq_len(n_focal), col_i = seq_len(n_focal)), row_species = focal_order[row_i],
    col_species = focal_order[col_i], pair_key = paste(pmin(row_species, col_species), pmax(row_species, col_species),
        sep = "||"), ratio = ifelse(row_i == col_i, NA_real_, as.numeric(ff_ratio[pair_key])), n_sig = ifelse(row_i ==
        col_i, NA_integer_, as.integer(ff_n_sig[pair_key])), meta_r = ifelse(row_i == col_i, NA_real_, as.numeric(ff_meta_r[pair_key])),
    meta_p = ifelse(row_i == col_i, NA_real_, as.numeric(ff_meta_p[pair_key])), meta_q = ifelse(row_i == col_i,
        NA_real_, as.numeric(ff_meta_q[pair_key])), meta_symbol = ifelse(row_i == col_i, NA_character_, as.character(ff_meta_symbol[pair_key])),
    meta_nominal_symbol = ifelse(row_i == col_i, NA_character_, as.character(ff_meta_nominal_symbol[pair_key])),
    xmin = HM_X0 + (col_i - 1) * ring_w, xmax = HM_X0 + col_i * ring_w, ymin = row_lo[row_i], ymax = row_hi[row_i],
    size_f = if (max(MAX_N_SIG, 1) <= 1) rep(1, length(n_sig)) else 1/3 + 2/3 * sqrt((pmax(as.numeric(n_sig), 1) -
        1)/(max(as.numeric(MAX_N_SIG), 1) - 1)), cx = (xmin + xmax)/2, cy = (ymin + ymax)/2, sxmin = cx - (xmax -
        xmin)/2 * size_f, sxmax = cx + (xmax - xmin)/2 * size_f, symin = cy - (ymax - ymin)/2 * size_f, symax = cy +
        (ymax - ymin)/2 * size_f)
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
    if (!length(todo))
        break
    progressed <- FALSE
    for (i in todo) {
        kids <- dd$y[dd$parent == dd$node[i] & dd$node != dd$node[i]]
        if (length(kids) && !anyNA(kids)) {
            dd$y[i] <- mean(range(kids))
            progressed <- TRUE
        }
    }
    if (!progressed)
        stop("Could not resolve dendrogram node coordinates")
}

DEN_Y0 <- TREE_CY + RING_R1 + DENDRO_GAP
dmax <- max(dd$x)

den_seg <- mutate(left_join(filter(dd, node != parent), select(dd, node, xp = x, yp = y), by = c(parent = "node")),
    cx_node = HM_X0 + (y - 0.5) * ring_w, cx_par = HM_X0 + (yp - 0.5) * ring_w, cy_node = DEN_Y0 + (1 - x/dmax) *
        DENDRO_H, cy_par = DEN_Y0 + (1 - xp/dmax) * DENDRO_H)
if (!order_is_phylo) den_seg <- den_seg[0, ]
COLLAB_Y <- DEN_Y0 + DENDRO_H + COLLAB_GAP
lab_df <- tibble(x = HM_X0 - 1.2, y = (row_lo + row_hi)/2, label = focal_short)
col_df <- tibble(x = col_centre, y = COLLAB_Y, label = focal_short)

overlay <- ggplot() + geom_rect(data = heat_df, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax), fill = "white",
    color = GRID_HEX, linewidth = 0.16) + geom_rect(data = filter(heat_df, is.finite(ratio)), aes(xmin = sxmin,
    xmax = sxmax, ymin = symin, ymax = symax, fill = ratio), color = TILE_BORDER, linewidth = 0.16) + geom_rect(data = filter(heat_df,
    is.finite(meta_p), meta_p < 0.05), aes(xmin = sxmin, xmax = sxmax, ymin = symin, ymax = symax), fill = NA,
    colour = "black", linewidth = 0.2) + geom_rect(data = filter(heat_df, is.finite(meta_q), meta_q < 0.05), aes(xmin = sxmin,
    xmax = sxmax, ymin = symin, ymax = symax), fill = NA, colour = "black", linewidth = 0.4) + geom_text(data = filter(heat_df,
    !is.na(meta_nominal_symbol), meta_nominal_symbol != ""), aes(cx, cy, label = meta_nominal_symbol), size = FONT_PT/PT,
    colour = "black") + scale_fill_gradient2(low = NEG_HEX, mid = "white", high = POS_HEX, midpoint = 0.5, limits = c(0,
    1), oob = squish, guide = "none") + geom_segment(data = lab_df, aes(x = x + 0.1, xend = HM_X0, y = y, yend = y),
    color = "grey55", linewidth = 0.2) + geom_text(data = lab_df, aes(x = x - 0.5, y = y, label = label), hjust = 1,
    vjust = 0.5, size = FONT_PT/PT, color = "grey15") + geom_segment(data = den_seg, aes(x = cx_par, xend = cx_node,
    y = cy_node, yend = cy_node), color = "grey35", linewidth = 0.22) + geom_segment(data = den_seg, aes(x = cx_par,
    xend = cx_par, y = cy_par, yend = cy_node), color = "grey35", linewidth = 0.22) + geom_text(data = col_df,
    aes(x = x, y = y, label = label), angle = 90, hjust = 0, vjust = 0.5, size = FONT_PT/PT, color = "grey15") +
    scale_x_continuous(limits = c(0, CANVAS_W), expand = c(0, 0)) + scale_y_continuous(limits = c(0, CANVAS_H),
    expand = c(0, 0)) + theme_void() + theme(plot.margin = margin(0, 0, 0, 0), plot.background = element_blank(),
    panel.background = element_blank(), legend.position = "none")

legend_theme <- theme(legend.title = element_text(size = FONT_TITLE_PT, color = "grey15"), legend.text = element_text(size = FONT_PT,
    color = "grey15"), legend.background = element_blank(), legend.key = element_blank(), legend.margin = margin(0,
    0, 0, 0), legend.justification = "left")

p_leg_ratio <- ggplot(tibble(x = 1, y = 1, v = 0.5), aes(x, y, fill = v)) + geom_tile() + scale_fill_gradient2(low = NEG_HEX,
    mid = "white", high = POS_HEX, midpoint = 0.5, limits = c(0, 1), breaks = c(0, 0.25, 0.5, 0.75, 1), name = "Positive ratio among\nnominally significant cohorts",
    guide = guide_colorbar(title.position = "top", direction = "horizontal", barwidth = unit(38, "mm"), barheight = unit(3.4,
        "mm"), frame.colour = "grey40", frame.linewidth = 0.2, ticks.colour = "grey40")) + theme_void() + legend_theme +
    theme(legend.position = "right")

p_leg_dir <- ggplot(tibble(k = factor(c("Consistently positive", "Consistently negative", "Other / mixed"), levels = c("Consistently positive",
    "Consistently negative", "Other / mixed"))), aes(x = 1, y = k, color = k)) + geom_point(size = 2.2) + scale_color_manual(values = c(`Consistently positive` = POS_HEX,
    `Consistently negative` = NEG_HEX, `Other / mixed` = NEUTRAL_HEX), labels = c("All positive", "All negative",
    "Mixed / other"), name = paste0("Partner direction (>=", highlight_min_faec, " Faecalibacterium)"), guide = guide_legend(title.position = "top",
    ncol = 1, override.aes = list(size = 2.2))) + theme_void() + legend_theme + theme(legend.position = "right")

fam_legend_labels <- paste0(family_levels, " (", family_group_counts$n_partner, ")")

p_leg_fam <- ggplot(tibble(f = factor(family_levels, levels = family_levels)), aes(x = 1, y = f, fill = f)) + geom_tile() +
    scale_fill_manual(values = family_colors, drop = FALSE, labels = fam_legend_labels, name = "Family of partner species (>3 species shown)",
        guide = guide_legend(title.position = "top", ncol = 2, byrow = FALSE, override.aes = list(color = NA),
            keywidth = unit(3.2, "mm"), keyheight = unit(3.2, "mm"))) + theme_void() + legend_theme + theme(legend.position = "right")

p_leg_size <- local({
    breaks_n <- sort(unique(c(1L, round(MAX_N_SIG/3), round(2 * MAX_N_SIG/3), MAX_N_SIG)))
    df <- tibble(x = seq_along(breaks_n), lab = as.character(breaks_n), pt_size = 4.4 * if (max(MAX_N_SIG, 1) <=
        1)
        rep(1, length(breaks_n))
    else 1/3 + 2/3 * sqrt((pmax(as.numeric(breaks_n), 1) - 1)/(max(as.numeric(MAX_N_SIG), 1) - 1)))
    ggplot(df, aes(x = x, y = 1)) + geom_point(aes(size = pt_size), shape = 15, color = "grey50") + scale_size_identity() +
        geom_text(aes(label = lab), y = 0.35, size = FONT_PT/PT, color = "grey15") + annotate("text", x = 0.4,
        y = 1.9, label = "Cohorts significant (n_sig)", hjust = 0, size = FONT_TITLE_PT/PT, color = "grey15") +
        scale_x_continuous(limits = c(0, nrow(df) + 1)) + scale_y_continuous(limits = c(0, 2.2)) + theme_void()
})
leg_ratio <- cowplot::get_legend(p_leg_ratio)
leg_dir <- cowplot::get_legend(p_leg_dir)
leg_fam <- cowplot::get_legend(p_leg_fam)

final_plot <- cowplot::ggdraw() + theme(plot.background = element_rect(fill = NA, color = NA)) + cowplot::draw_grob(tree_grob,
    x = (TREE_CX - TREE_SIDE/2)/CANVAS_W, y = (TREE_CY - TREE_SIDE/2)/CANVAS_H, width = TREE_SIDE/CANVAS_W, height = TREE_SIDE/CANVAS_H,
    clip = "off") + cowplot::draw_plot(overlay, x = 0, y = 0, width = 1, height = 1) + cowplot::draw_grob(leg_fam,
    x = 0.49, y = 0.86, width = 0.49, height = 0.13) + cowplot::draw_plot(p_leg_size, x = 0.02, y = 0.02, width = 0.34,
    height = 0.08) + cowplot::draw_grob(leg_dir, x = 0.02, y = 0.1, width = 0.43, height = 0.075) + cowplot::draw_grob(leg_ratio,
    x = 0.69, y = 0.035, width = 0.29, height = 0.08) + cowplot::draw_label("Thin frame: meta P<0.05; thick frame: global BH q<0.05",
    x = 0.75, y = 0.02, size = 5)
dir.create(dirname(out_pdf), recursive = TRUE, showWarnings = FALSE)

ggsave(out_pdf, final_plot, width = CANVAS_W/25.4, height = CANVAS_H/25.4, device = cairo_pdf, bg = "transparent")

if (Sys.getenv("FP_PNG", "1") == "1") ggsave(out_png, final_plot, width = CANVAS_W/25.4, height = CANVAS_H/25.4,
    dpi = 300, bg = "white")
