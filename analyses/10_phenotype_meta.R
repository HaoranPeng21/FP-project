# Figure 5A. Input: frozen per-phenotype random-effects meta-analysis results.
# Retain the original phenome-wide BH family before any display filtering.

# 1. Set these two private directories (or edit the paths directly).
data_dir <- Sys.getenv("FP_DATA_DIR", "/path/to/private/data/10_phenotype_meta")
out_dir <- Sys.getenv("FP_OUTPUT_DIR", "/path/to/private/results/10_phenotype_meta")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
font_family <- "sans"

# 2. Packages.
library(ComplexHeatmap)
library(circlize)
library(grid)

# 3. Read the inputs, analyze, and draw the figures.
d <- read.delim(file.path(data_dir, "phenotype_meta.tsv"), check.names = FALSE, quote = "")

scope <- read.delim(file.path(data_dir, "phenotype_scope.tsv"), check.names = FALSE, quote = "")

info <- read.delim(file.path(data_dir, "phenotype_info.tsv"), check.names = FALSE, quote = "")

phy <- read.delim(file.path(data_dir, "species_phylogeny_order.tsv"), check.names = FALSE, quote = "")
d$species <- sub("^.*\\|", "", d$feature)
keep <- scope$species_label[scope$full_prevalence > 0.05]
species <- phy$species[order(phy$phylogeny_order)]
species <- species[species %in% keep]
stopifnot(length(species) == 14, !anyDuplicated(d[c("phenotype", "species")]))
d <- d[!is.na(d$k_study) & d$k_study >= 2 & d$species %in% species, ]
ph <- aggregate(p_meta_fdr ~ phenotype, d, min)
ph <- ph[ph$p_meta_fdr < 0.05, ]
ph <- ph[order(ph$p_meta_fdr), ]
classes <- info$Class[match(ph$phenotype, info$UnifiedName)]
classes[is.na(classes) | classes == ""] <- "HealthDisease"
ii <- order(match(classes, unique(classes)))
ph <- ph[ii, ]
classes <- classes[ii]
mat <- matrix(NA_real_, nrow(ph), length(species), dimnames = list(ph$phenotype, species))
stars <- matrix("", nrow(mat), ncol(mat))
k <- tapply(d$k_study, d$phenotype, max, na.rm = TRUE)

for (i in seq_len(nrow(d))) {
    r <- match(d$phenotype[i], rownames(mat))
    c <- match(d$species[i], colnames(mat))
    if (!is.na(r) && !is.na(c)) {
        mat[r, c] <- d$beta[i]/d$se[i]
        stars[r, c] <- if (!is.na(d$p_meta_fdr[i]) && d$p_meta_fdr[i] < 0.05)
            "*"
        else ""
    }
}

class_colours <- c(Metabolite = "#B59CCE", Diet = "#B4C384", MedicalMeasurement = "#74B8A1", HealthDisease = "#E6797A",
    Anthropology = "#D9AA59", MedicationUse = "#EBE997", Stool = "#39A3B2", CurrentExposure = "#E07A5F", EarlyLifeExposure = "#C98BA8",
    ExperimentTechniqueSample = "#9AA7B6")
missing <- setdiff(unique(classes), names(class_colours))

if (length(missing)) class_colours <- c(class_colours, setNames(grDevices::hcl.colors(length(missing)), missing))

ht <- ComplexHeatmap::Heatmap(mat, name = "Meta-analysis z", col = circlize::colorRamp2(c(-5, 0, 5), c("#6485A6",
    "white", "#CD5B48")), na_col = "#eeeeee", cluster_rows = FALSE, cluster_columns = FALSE, row_labels = paste0(rownames(mat),
    " (k=", k[rownames(mat)], ")"), column_labels = sub("^Faecalibacterium ", "F. ", sub("^Faecalibacterium prausnitzii_F$",
    "Faecalibacterium langellae", sub("^Faecalibacterium sp900539945$", "Faecalibacterium harmsenii", sub("^s__",
        "", colnames(mat))))), column_names_rot = 45, column_names_gp = grid::gpar(fontsize = 8, fontface = "italic"),
    row_names_gp = grid::gpar(fontsize = 6), row_names_side = "right", left_annotation = ComplexHeatmap::rowAnnotation(Class = classes,
        col = list(Class = class_colours)), rect_gp = grid::gpar(col = "white", lwd = 0.2), cell_fun = function(j,
        i, x, y, w, h, fill) grid::grid.text(stars[i, j], x, y, gp = grid::gpar(fontsize = 7)))

grDevices::cairo_pdf(file.path(out_dir, "Figure_5A.pdf"), width = 7, height = 11, family = font_family)
ComplexHeatmap::draw(ht, heatmap_legend_side = "bottom")
grDevices::dev.off()

write.table(data.frame(phenotype = rownames(mat), Class = classes, k_max = k[rownames(mat)]), file.path(out_dir,
    "displayed_phenotypes.tsv"), sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
