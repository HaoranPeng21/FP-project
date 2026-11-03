# Figure S2A-C. Final-catalogue QC (17,103 genomes / 83 species), rarefaction and SGB mapping.

# 1. Set these two private directories (or edit the paths directly).
data_dir <- Sys.getenv("FP_DATA_DIR", "/path/to/private/data/13_catalogue_qc")
out_dir <- Sys.getenv("FP_OUTPUT_DIR", "/path/to/private/results/13_catalogue_qc")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
font_family <- "sans"

# 2. Packages.
library(ggplot2)
library(dplyr)
library(patchwork)

# 3. Read the inputs, analyze, and draw the figures.
catalog <- read.csv(file.path(data_dir, "catalogue.csv"), check.names = FALSE)

analysis_set <- read.delim(file.path(data_dir, "analysis_set.tsv"), check.names = FALSE)
stopifnot(nrow(analysis_set) == 17103, !anyDuplicated(analysis_set$ID),
          length(unique(analysis_set$Haoran_species)) == 83)
catalog <- catalog[match(analysis_set$ID, catalog$ID), ]
stopifnot(nrow(catalog) == 17103, !anyNA(catalog$ID),
          identical(catalog$ID, analysis_set$ID))
catalog$Host1[is.na(catalog$Host1) | !nzchar(catalog$Host1)] <- "Unknown"
counts <- as.data.frame(table(host = catalog$Host1))
stopifnot(sum(counts$Freq) == 17103)
a <- ggplot(catalog, aes(CheckM2.completeness, CheckM2.contamination, colour = Host1)) +
  geom_point(size = .45, alpha = .35) +
  labs(x = "CheckM2 completeness (%)", y = "CheckM2 contamination (%)", colour = "Host",
       title = "Final catalogue: 17,103 genomes / 83 species", fill = NULL) + theme_bw(base_size = 8)

rare <- read.delim(file.path(data_dir, "rarefaction.tsv"), check.names = FALSE, quote = "")
stopifnot(max(rare$genomes_sampled) == 17103, max(rare$q975_species) <= 83)

b <- ggplot(rare, aes(genomes_sampled, mean_species)) + geom_ribbon(aes(ymin = q025_species, ymax = q975_species),
    fill = "#D6EAF1") + geom_line(colour = "#4DBBD5") + geom_hline(yintercept = 83, linetype = 2) + labs(x = "Sampled genomes",
    y = "Observed species", title = "Species rarefaction", fill = NULL, colour = NULL) + (ggplot2::theme_bw(base_size = 8, base_family = "sans") +
    ggplot2::theme(panel.grid.minor = ggplot2::element_blank(), legend.key = ggplot2::element_blank()))

mapping <- read.delim(file.path(data_dir, "benchmark_mapping.tsv"), check.names = FALSE, quote = "")

mapping$category <- ifelse(grepl("ANI", mapping$mapping_basis, ignore.case = TRUE), "ANI >=95%", ifelse(!is.na(mapping$species_id) &
    nzchar(mapping$species_id), "GTDB-name link", "Unresolved"))
mapping <- mapping[order(mapping$category, mapping$sgb), ]

mapping$label <- paste(mapping$sgb, ifelse(mapping$category == "Unresolved", "unresolved", sub("^Faecalibacterium ",
    "F. ", sub("^Faecalibacterium prausnitzii_F$", "Faecalibacterium langellae", sub("^Faecalibacterium sp900539945$",
        "Faecalibacterium harmsenii", sub("^s__", "", mapping$project_species))))))
mapping$y <- rev(seq_len(nrow(mapping)))

c <- ggplot(mapping, aes(x = 0, y = y, label = label, colour = category)) + geom_text(hjust = 0, size = 2.5) +
    xlim(0, 1.3) + labs(title = "SGB mapping", x = NULL, y = NULL, fill = NULL, colour = NULL) + theme_void() + theme(legend.position = "bottom")

ggplot2::ggsave(filename = file.path(out_dir, "Figure_S2ABC.pdf"), plot = a + b + c + patchwork::plot_layout(widths = c(1,
    1, 1.2)), width = 12, height = 5.5, device = grDevices::cairo_pdf, limitsize = FALSE)

write.table(counts, file.path(out_dir, "quality_counts.tsv"), sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
