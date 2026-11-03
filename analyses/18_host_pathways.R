# Figure S4. Retained nominal species-specific Fisher screen (not an FDR discovery set).
# Global q-values come from the upstream pathway/EC/transporter family, not these rows.


data_dir <- Sys.getenv("FP_DATA_DIR", "/path/to/private/data/18_host_pathways")
out_dir <- Sys.getenv("FP_OUTPUT_DIR", "/path/to/private/results/18_host_pathways")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
font_family <- "sans"
library(ggplot2)
library(patchwork)

d <- read.delim(file.path(data_dir, "host_pathways.tsv"), check.names = FALSE, quote = "")
stopifnot(all(d$feature_kind == "pathway"), all(d$fisher_p < 0.05), !anyDuplicated(d[c("species_id",
    "feature_id")]))

pairs <- list(c("skani_group_48", "skani_group_17"), c("skani_group_45", "skani_group_52"), c("skani_group_59",
    "skani_group_37"))
panels <- list()

for (species_id in unlist(pairs)) {
    id <- species_id
    z <- d[d$species_id == id, ]
    z <- z[order(z$prevalence_diff, z$fisher_p, z$pathway_name), ]
    z$label <- paste0(z$pathway_name, ifelse(z$fisher_q_global < 0.05, " *", ""))
    z$label <- factor(z$label, levels = unique(z$label))
    panels[[species_id]] <- ggplot(z, aes(y = label)) + geom_segment(aes(x = 100 * prevalence_human,
        xend = 100 * prevalence_animal, yend = label), colour = "grey65") + geom_point(aes(x = 100 *
        prevalence_animal, colour = "Animal"), size = 1.5) + geom_point(aes(x = 100 * prevalence_human,
        colour = "Human"), shape = 1, size = 1.5) + scale_colour_manual(values = c(Animal = "#E64B35",
        Human = "#4DBBD5")) + scale_x_continuous(limits = c(0, 100), breaks = seq(0, 100, 25)) + labs(title = sub("^Faecalibacterium ",
        "F. ", sub("^Faecalibacterium prausnitzii_F$", "Faecalibacterium langellae", sub("^Faecalibacterium sp900539945$",
            "Faecalibacterium harmsenii", sub("^s__", "", z$display_name[1])))), subtitle = sprintf("A %d | H %d",
        z$n_animal[1], z$n_human[1]), x = "Pathway prevalence (%)", y = NULL, fill = NULL, colour = NULL) + (ggplot2::theme_bw(base_size = 6,
        base_family = "sans") + ggplot2::theme(panel.grid.minor = ggplot2::element_blank(), legend.key = ggplot2::element_blank())) +
        theme(legend.position = "bottom", panel.grid.major.y = element_blank())
}

columns <- vector("list", 3)

for (i in seq_along(pairs)) {
    heights <- table(factor(d$species_id, levels = pairs[[i]])) + 3
    columns[[i]] <- patchwork::wrap_plots(panels[pairs[[i]]], ncol = 1, heights = as.numeric(heights))
}

ggsave(file.path(out_dir, "Figure_S4.pdf"), patchwork::wrap_plots(columns, nrow = 1, guides = "collect"),
    width = 15, height = 10, device = cairo_pdf)
