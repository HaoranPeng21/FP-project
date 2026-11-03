# Figure 6. Frozen PC5 logistic/Firth gene-family results and genome association scans.
# Brackets show complete physical-neighbourhood carriage; stars use member-family q.
# Upstream model fitting and genome assembly are outside this data-to-figure entry point.

# 1. Set these two private directories (or edit the paths directly).
data_dir <- Sys.getenv("FP_DATA_DIR", "/path/to/private/data/12_within_species")
out_dir <- Sys.getenv("FP_OUTPUT_DIR", "/path/to/private/results/12_within_species")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
font_family <- "sans"

# 2. Packages.
library(ggplot2)
library(ggsci)
library(patchwork)
library(grid)

# 3. Read the inputs, analyze, and draw the figures.
main_data_dir <- file.path(data_dir, "within_main")
supp_data_dir <- file.path(data_dir, "within_supplementary")
variant_data_dir <- file.path(data_dir, "within_presence")
npg <- pal_npg("nrc")(10)
col_positive <- npg[[1]]
col_negative <- npg[[2]]
col_present <- "#6485A6"
col_absent <- "#ECEEF1"
col_rest <- "#9AA1A9"
col_background <- "#C9CDD2"
col_cluster <- c(npg[[3]], npg[[5]], npg[[4]])
col_single <- "#A8B0B8"
col_ink <- "#1A1A1A"
col_ink_soft <- "#4D4D4D"
col_rule <- "#7C848C"
PT_PANEL <- 11
PT_AXIS_TTL <- 9
PT_BODY <- 8
plot_theme <- {
    theme_classic(base_size = PT_BODY, base_family = font_family) + theme(text = element_text(family = font_family,
        colour = col_ink), axis.title = element_text(size = PT_AXIS_TTL, colour = col_ink), axis.text = element_text(size = PT_BODY,
        colour = col_ink_soft), plot.title = element_text(size = PT_PANEL, face = "bold", hjust = 0, margin = margin(b = 2.5)),
        plot.title.position = "plot", legend.title = element_text(size = PT_BODY), legend.text = element_text(size = PT_BODY),
        axis.line = element_line(linewidth = 0.32, colour = col_rule), axis.ticks = element_line(linewidth = 0.32,
            colour = col_rule), axis.ticks.length = unit(1.1, "pt"), plot.margin = margin(2, 3, 2, 3))
}

overview <- read.delim(file.path(variant_data_dir, "panel_a_presence_species_phenotype_class_counts.tsv"), check.names = FALSE,
    quote = "")

class_display <- c(HealthDisease = "Health/disease", MedicationUse = "Medication use", Stool = "Stool", MedicalMeasurement = "Medical measurement",
    CurrentExposure = "Current exposure", Anthropology = "Anthropology", Diet = "Diet", Metabolite = "Metabolite",
    EarlyLifeExposure = "Early-life exposure", ExperimentTechniqueSample = "Technical/sample")
d <- overview
class_levels <- unique(d$Class[order(d$class_order)])
species_levels <- unique(d$species[order(d$species_order)])
d$class_label <- factor(d$Class, levels = class_levels, labels = class_display[class_levels])

d$species_label <- factor(paste0("F. ", gsub("_", " ", d$species, fixed = TRUE)), levels = paste0("F. ", gsub("_",
    " ", species_levels, fixed = TRUE)))
class_palette <- setNames(d$class_color[match(class_levels, d$Class)], levels(d$class_label))
totals <- unique(d[c("species_label", "total_phenotypes", "species_order")])
totals <- totals[order(totals$species_order), ]

p_a <- ggplot(d, aes(species_label, n_phenotypes, fill = class_label)) + geom_col(width = 0.74, colour = "white",
    linewidth = 0.18) + geom_text(data = totals, aes(species_label, total_phenotypes + 0.42, label = total_phenotypes),
    inherit.aes = FALSE, family = font_family, size = PT_BODY/.pt, colour = col_ink_soft) + scale_fill_manual(values = class_palette,
    drop = TRUE) + scale_y_continuous(limits = c(0, max(totals$total_phenotypes) + 2.2), breaks = scales::breaks_pretty(4),
    expand = expansion(mult = c(0, 0))) + labs(title = "a   Associated phenotypes per species", y = "Phenotypes (q ≤ 0.05)",
    fill = "Phenotype class") + plot_theme + theme(panel.grid.major.y = element_line(colour = "#EDEFF1", linewidth = 0.26),
    axis.text.x = element_text(face = "italic", angle = 32, hjust = 1, vjust = 1, size = PT_BODY), legend.position = "inside",
    legend.position.inside = c(0.995, 1.02), legend.justification = c(1, 1), legend.title = element_text(size = PT_BODY,
        face = "bold"), legend.key.height = unit(3, "mm"), legend.key.width = unit(3.4, "mm"), legend.key.spacing.y = unit(0.4,
        "mm"), legend.margin = margin(0, 0, 0, 0), legend.background = element_blank(), plot.margin = margin(1,
        3, 1, 2)) + guides(fill = guide_legend(ncol = 2, byrow = FALSE))

case_specs <- list(list(letter = "b", id = "presence__BMI__taiwanense__negative", lead = "BMI", species = "F. taiwanense",
    direction = "negative", assoc = file.path(variant_data_dir, "presence_bmi_ranked_genome_associations.tsv.gz"),
    track = file.path(variant_data_dir, "presence_bmi_gene_track.tsv"), heat = file.path(variant_data_dir, "presence_bmi_gene_presence_absence.tsv.gz"),
    stats = file.path(variant_data_dir, "presence_bmi_gene_statistics.tsv"), modules = file.path(variant_data_dir,
        "presence_bmi_exact_module_prevalence.tsv"), layout = "main"), list(letter = "c", id = "presence__TricyclicAntidepressants__longum__positive",
    lead = "Tricyclic antidepressants", species = "F. longum", direction = "positive", assoc = file.path(main_data_dir,
        "three_case_ranked_genome_associations.tsv.gz"), track = file.path(main_data_dir, "three_case_gene_track.tsv"),
    heat = file.path(main_data_dir, "three_case_gene_presence_absence.tsv.gz"), stats = file.path(main_data_dir,
        "three_case_gene_statistics.tsv"), modules = file.path(main_data_dir, "three_case_exact_module_prevalence.tsv"),
    layout = "main"), list(letter = "d", id = "presence__PPI__duncaniae__positive", lead = "Proton-pump inhibitors",
    species = "F. duncaniae", direction = "positive", assoc = file.path(supp_data_dir, "supplementary_two_case_associations.tsv.gz"),
    track = file.path(supp_data_dir, "supplementary_two_case_gene_track.tsv"), heat = file.path(supp_data_dir,
        "supplementary_two_case_heatmap.tsv.gz"), stats = file.path(supp_data_dir, "supplementary_two_case_gene_stats.tsv"),
    modules = file.path(supp_data_dir, "supplementary_two_case_module_prevalence.tsv"), layout = "supplementary"))

cases <- vector("list", length(case_specs))

for (case_index in seq_along(case_specs)) {
    spec <- case_specs[[case_index]]
    a <- subset(read.delim(gzfile(spec$assoc), check.names = FALSE, quote = ""), candidate_id == spec$id)
    a$is_q_hit <- tolower(as.character(a$is_q_hit)) %in% c("true", "t", "1")
    a$q_phenotype <- as.numeric(a$q_phenotype)
    t <- subset(read.delim(spec$track, check.names = FALSE, quote = ""), candidate_id == spec$id)
    t$is_cluster <- tolower(as.character(t$is_cluster)) %in% c("true", "t", "1")
    h <- subset(read.delim(gzfile(spec$heat), check.names = FALSE, quote = ""), candidate_id == spec$id)
    h <- if (spec$layout == "supplementary") {
        h[tolower(as.character(h$is_q_hit)) %in% c("true", "t", "1"), , drop = FALSE]
    }
    else {
        h[h$genome_group == "q≤0.05 genome", , drop = FALSE]
    }
    s <- subset(read.delim(spec$stats, check.names = FALSE, quote = ""), candidate_id == spec$id)
    s$hit_frequency <- as.numeric(if (spec$layout == "supplementary") s$primary_frequency else s$hit_frequency)
    s$rest_frequency <- as.numeric(s$rest_frequency)
    m <- subset(read.delim(spec$modules, check.names = FALSE, quote = ""), candidate_id == spec$id)
    m$genome_group[m$genome_group == "comparison rest"] <- "species rest"
    cases[[case_index]] <- list(spec = spec, assoc = a, track = t, heat = h, stats = s, modules = m)
}

case_panels <- vector("list", length(cases))

for (case_index in seq_along(cases)) {
    cs <- cases[[case_index]]
    spec <- cs$spec
    a <- cs$assoc
    hit_col <- if (spec$direction == "positive")
        col_positive
    else col_negative
    sig <- a[a$q_phenotype <= 0.05, ]
    ns <- a[a$q_phenotype > 0.05, ]
    fdr_line <- mean(c(min(sig$mlog10p, na.rm = TRUE), max(ns$mlog10p, na.rm = TRUE)))
    sig$sign_label <- ifelse(sig$beta >= 0, "positive", "negative")
    n_direction_hits <- length(unique(a$feature_key[a$is_q_hit]))
    y_top <- max(a$mlog10p, na.rm = TRUE) * 1.14
    x_rng <- range(a$within_species_rank, na.rm = TRUE)
    p_manhattan <- ggplot(a, aes(within_species_rank, mlog10p)) + geom_point(colour = col_background, size = 0.48,
        alpha = 0.6, stroke = 0) + geom_hline(yintercept = fdr_line, colour = col_rule, linewidth = 0.36, linetype = "22") +
        geom_point(data = sig, aes(fill = sign_label), shape = 21, colour = "#2B2B2B", stroke = 0.24, size = 1.6) +
        scale_fill_manual(values = c(positive = col_positive, negative = col_negative), guide = "none") + scale_x_continuous(breaks = scales::breaks_pretty(4),
        limits = x_rng, expand = expansion(mult = c(0.03, 0.05)), labels = function(v) format(v, big.mark = ",",
            trim = TRUE)) + scale_y_continuous(limits = c(0, y_top), breaks = scales::breaks_pretty(4), expand = expansion(mult = c(0.01,
        0))) + labs(title = paste(spec$letter, spec$lead, spec$species), x = "Within-species ANI rank", y = expression(-log[10] *
        "(" * italic(P)[meta] * ")")) + plot_theme + theme(plot.margin = margin(1, 5, 2, 1))
    t <- merge(cs$track, cs$stats[c("gene_family", "hit_frequency", "rest_frequency")], by = "gene_family", all.x = TRUE,
        sort = FALSE)
    t <- t[order(t$plot_x), ]
    cluster_ids <- unique(t$module_id[t$is_cluster])
    cluster_palette <- setNames(rep(col_cluster, length.out = length(cluster_ids)), cluster_ids)
    t$fill_key <- ifelse(t$is_cluster, t$module_id, "standalone")
    arrow_palette <- c(cluster_palette, standalone = col_single)
    brackets <- do.call(rbind, lapply(cluster_ids, function(mid) {
        z <- t[t$module_id == mid, ]
        mp <- cs$modules[cs$modules$module_id == mid, ]
        q_row <- mp[mp$genome_group == "q≤0.05 genomes", ][1, ]
        r_row <- mp[mp$genome_group == "species rest", ][1, ]
        label <- if (nrow(mp) == 2) {
            paste0(sprintf("%.0f%%", 100 * (q_row$n_exact_physical_module/q_row$n_group_genomes)), " vs ", sprintf("%.0f%%",
                100 * (r_row$n_exact_physical_module/r_row$n_group_genomes)))
        }
        else {
            ""
        }
        data.frame(module_id = mid, xmin = min(z$plot_x), xmax = max(z$plot_x), xmid = mean(range(z$plot_x)), label = label)
    }))
    arrows <- {
        d <- t
        do.call(rbind, lapply(seq_len(nrow(d)), function(i) {
            x <- d$plot_x[[i]]
            xx <- if (d$strand[[i]] == "-") {
                c(x + 0.42, x - 0.14, x - 0.14, x - 0.44, x - 0.14, x - 0.14, x + 0.42)
            }
            else {
                c(x - 0.42, x + 0.14, x + 0.14, x + 0.44, x + 0.14, x + 0.14, x - 0.42)
            }
            data.frame(x = xx, y = c(-0.17, -0.17, -0.31, 0, 0.31, 0.17, 0.17), polygon_id = paste0(d$module_id[[i]],
                "_", d$gene_family[[i]]), fill_key = d$fill_key[[i]])
        }))
    }
    x_min <- min(t$plot_x) - 0.62
    x_max <- max(t$plot_x) + 0.62
    p_track <- ggplot() + geom_segment(data = brackets, aes(x = xmin - 0.42, xend = xmax + 0.42, y = 0, yend = 0,
        colour = module_id), linewidth = 0.85) + geom_polygon(data = arrows, aes(x, y, group = polygon_id, fill = fill_key),
        colour = "white", linewidth = 0.22) + geom_text(data = t, aes(plot_x - 0.12, 0.46, label = gene_label),
        size = PT_BODY/.pt, family = font_family, fontface = "italic", angle = 45, hjust = 0, colour = col_ink) +
        geom_segment(data = brackets, aes(x = xmin - 0.44, xend = xmax + 0.44, y = -0.6, yend = -0.6), colour = col_rule,
            linewidth = 0.3) + geom_segment(data = brackets, aes(x = xmin - 0.44, xend = xmin - 0.44, y = -0.6,
        yend = -0.46), colour = col_rule, linewidth = 0.3) + geom_segment(data = brackets, aes(x = xmax + 0.44,
        xend = xmax + 0.44, y = -0.6, yend = -0.46), colour = col_rule, linewidth = 0.3) + geom_text(data = brackets,
        aes(xmid, -0.94, label = label), size = PT_BODY/.pt, family = font_family, colour = col_ink_soft) + scale_fill_manual(values = arrow_palette,
        guide = "none") + scale_colour_manual(values = cluster_palette, guide = "none") + scale_x_continuous(limits = c(x_min,
        x_max), expand = c(0, 0)) + coord_cartesian(ylim = c(-1.15, 1.75), clip = "off") + labs() + theme_void(base_family = font_family,
        base_size = PT_BODY) + theme(plot.margin = margin(1, 4, 0, 2))
    prev <- t[c("plot_x", "gene_label", "hit_frequency", "rest_frequency")]
    prev$hit_pct <- 100 * prev$hit_frequency
    prev$rest_pct <- 100 * prev$rest_frequency
    p_prev <- ggplot(prev) + geom_segment(aes(x = plot_x, xend = plot_x, y = rest_pct, yend = hit_pct), colour = "#C6CBD1",
        linewidth = 0.55) + geom_point(aes(plot_x, rest_pct), colour = col_rest, size = 1.15, stroke = 0) + geom_point(aes(plot_x,
        hit_pct), colour = col_present, size = 1.5, stroke = 0) + scale_x_continuous(limits = c(x_min, x_max),
        expand = c(0, 0)) + scale_y_continuous(limits = c(-6, 106), breaks = c(0, 50, 100), expand = c(0, 0)) +
        labs(y = "Prevalence\n(%)") + plot_theme + theme(panel.grid.major.y = element_line(colour = "#EDEFF1",
        linewidth = 0.26), axis.line.x = element_blank(), axis.ticks.x = element_blank(), axis.text.x = element_blank(),
        axis.line.y = element_blank(), axis.title.y = element_text(size = PT_BODY, colour = col_ink_soft, lineheight = 0.95),
        plot.margin = margin(0, 4, 0, 2))
    h <- cs$heat
    rank_by_genome <- tapply(h$within_species_rank, h$genome, min)
    ordered_genomes <- names(sort(rank_by_genome))
    h$y <- -match(h$genome, ordered_genomes)
    h$present_label <- ifelse(h$present == 1, "Present", "Absent")
    n_rows <- length(ordered_genomes)
    row_height <- if (n_rows > 30)
        1
    else 0.94
    p_heat <- ggplot(h, aes(plot_x, y)) + geom_tile(aes(fill = present_label), width = 0.9, height = row_height,
        colour = NA) + scale_fill_manual(values = c(Absent = col_absent, Present = col_present), guide = "none") +
        scale_x_continuous(limits = c(x_min, x_max), expand = c(0, 0)) + scale_y_continuous(limits = c(-n_rows -
        0.5, -0.5), expand = c(0, 0)) + coord_cartesian(clip = "off") + labs(y = paste0(format(n_rows, big.mark = ","),
        " genomes")) + theme_void(base_family = font_family, base_size = PT_BODY) + theme(axis.title.y = element_text(size = PT_BODY,
        colour = col_ink_soft, angle = 90, margin = margin(r = 1.5)), plot.margin = margin(1, 4, 1, 2))
    right <- p_track/p_prev/p_heat + plot_layout(heights = c(1, 0.66, 1.06))
    case_panels[[case_index]] <- (p_manhattan | right) + plot_layout(widths = c(1.04, 1.2))
}

key_items <- rbind(data.frame(row = 1, glyph = "head", label = "Genome scan", fill = NA), data.frame(row = 1, glyph = "point",
    label = "Tested genome", fill = col_background), data.frame(row = 1, glyph = "point", label = "q ≤ 0.05, β > 0",
    fill = col_positive), data.frame(row = 1, glyph = "point", label = "q ≤ 0.05, β < 0", fill = col_negative),
    data.frame(row = 1, glyph = "dash", label = "q = 0.05 threshold", fill = col_rule), data.frame(row = 2, glyph = "head",
        label = "Gene panels", fill = NA), data.frame(row = 2, glyph = "point", label = "q ≤ 0.05 genomes", fill = col_present),
    data.frame(row = 2, glyph = "point", label = "Species rest", fill = col_rest), data.frame(row = 2, glyph = "tile",
        label = "Gene present", fill = col_present), data.frame(row = 2, glyph = "tile", label = "Gene absent",
        fill = col_absent))

key_items$width <- ifelse(key_items$glyph == "head", 1.62 * nchar(key_items$label) + 3.4, 1.52 * nchar(key_items$label) +
    6.6)
key_items$x <- unlist(lapply(split(key_items, key_items$row), function(d) {
    cumsum(c(0, head(d$width, -1)))
}), use.names = FALSE)
key_items$glyph_x <- key_items$x + 1.6
key_items$text_x <- ifelse(key_items$glyph == "head", key_items$x, key_items$x + 4.4)
key_items$y <- ifelse(key_items$row == 1, 1, 0)
key_width <- max(vapply(split(key_items, key_items$row), function(d) sum(d$width), numeric(1)))
key_points <- key_items[key_items$glyph == "point", ]
key_tiles <- key_items[key_items$glyph == "tile", ]
key_dashes <- key_items[key_items$glyph == "dash", ]
key_heads <- key_items[key_items$glyph == "head", ]
key_labels <- key_items[key_items$glyph != "head", ]

p_key <- ggplot() + geom_point(data = key_points, aes(glyph_x, y), shape = 21, fill = key_points$fill, colour = "#2B2B2B",
    stroke = 0.24, size = 1.6) + geom_tile(data = key_tiles, aes(glyph_x, y), fill = key_tiles$fill, colour = "#C4CAD1",
    linewidth = 0.2, width = 2.1, height = 0.44) + geom_segment(data = key_dashes, aes(x = glyph_x - 1.5, xend = glyph_x +
    1.5, y = y, yend = y), colour = col_rule, linewidth = 0.36, linetype = "22") + geom_text(data = key_heads,
    aes(text_x, y, label = label), hjust = 0, size = PT_BODY/.pt, family = font_family, fontface = "bold", colour = col_ink) +
    geom_text(data = key_labels, aes(text_x, y, label = label), hjust = 0, size = PT_BODY/.pt, family = font_family,
        colour = col_ink_soft) + scale_x_continuous(limits = c(0, key_width * 1.02), expand = c(0, 0)) + scale_y_continuous(limits = c(-0.65,
    1.65), expand = c(0, 0)) + theme_void(base_family = font_family, base_size = PT_BODY) + theme(plot.margin = margin(2,
    3, 0, 2))

figure <- p_a/wrap_elements(full = case_panels[[1]], clip = FALSE)/wrap_elements(full = case_panels[[2]], clip = FALSE)/wrap_elements(full = case_panels[[3]],
    clip = FALSE)/p_key + plot_layout(heights = c(0.96, 1.02, 1.02, 1.02, 0.15)) + plot_annotation(theme = theme(plot.margin = margin(4,
    5, 3, 4)))
page_width_mm <- 210
page_height_mm <- 297
pdf_path <- file.path(out_dir, "Figure_6.pdf")
png_path <- file.path(out_dir, "Figure_6.png")

ggsave(pdf_path, figure, width = page_width_mm, height = page_height_mm, units = "mm", device = cairo_pdf)

if (Sys.getenv("FP_PNG", "1") == "1") {
    ggsave(png_path, figure, width = page_width_mm, height = page_height_mm, units = "mm", dpi = 300, bg = "white")
}
