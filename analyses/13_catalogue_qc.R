# Figure S2A-C. Catalogue QC, frozen species-rarefaction estimates, and SGB mapping.
run <- function() {
  need("ggplot2", "dplyr", "patchwork")
  # The source panel compares two nested quality tiers over the full catalogue.
  catalog <- read.csv(input("catalogue"), check.names = FALSE)
  check_columns(catalog, c("Host1", "CheckM2.completeness", "CheckM2.contamination"))
  catalog$Host1[is.na(catalog$Host1) | !nzchar(catalog$Host1)] <- "Unknown"
  catalog$quality <- ifelse(catalog$CheckM2.completeness >= 80 & catalog$CheckM2.contamination <= 5,
    "Completeness >=80, contamination <=5", "Completeness >=50, contamination <=10")
  catalog <- catalog[catalog$CheckM2.completeness >= 50 & catalog$CheckM2.contamination <= 10 &
    !is.na(catalog$CheckM2.completeness) & !is.na(catalog$CheckM2.contamination), ]
  counts <- as.data.frame(table(host = catalog$Host1, quality = catalog$quality))
  a <- ggplot(counts, aes(host, Freq, fill = quality)) + geom_col() +
    scale_fill_manual(values = c("#b5dce2", "#819bb0")) + labs(x = NULL, y = "Genomes", title = "Genomes by source") +
    figure_theme(8) + theme(legend.position = "top", axis.text.x = element_text(angle = 45, hjust = 1))
  rare <- read_input("rarefaction")
  b <- ggplot(rare, aes(genomes_sampled, mean_species)) +
    geom_ribbon(aes(ymin = q025_species, ymax = q975_species), fill = "#D6EAF1") +
    geom_line(colour = "#4DBBD5") + geom_hline(yintercept = 83, linetype = 2) +
    labs(x = "Sampled genomes", y = "Observed species", title = "Species rarefaction") + figure_theme(8)
  mapping <- read_input("benchmark_mapping")
  mapping$category <- ifelse(grepl("ANI", mapping$mapping_basis, ignore.case = TRUE), "ANI >=95%",
    ifelse(!is.na(mapping$species_id) & nzchar(mapping$species_id), "GTDB-name link", "Unresolved"))
  mapping <- mapping[order(mapping$category, mapping$sgb), ]
  mapping$label <- paste(mapping$sgb, ifelse(mapping$category == "Unresolved", "unresolved", short_species(mapping$project_species)))
  mapping$y <- rev(seq_len(nrow(mapping)))
  c <- ggplot(mapping, aes(x = 0, y = y, label = label, colour = category)) +
    geom_text(hjust = 0, size = 2.5) + xlim(0, 1.3) + labs(title = "SGB mapping") + theme_void() +
    theme(legend.position = "bottom")
  save_plot(a + b + c + patchwork::plot_layout(widths = c(1, 1, 1.2)), "Figure_S2ABC", 12, 5.5)
  write_result(counts, "quality_counts.tsv")
}
