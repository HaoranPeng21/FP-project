options(stringsAsFactors = FALSE)
extra_libs <- Sys.getenv("FP_R_LIBS", "")
if (nzchar(extra_libs)) .libPaths(c(strsplit(extra_libs, .Platform$path.sep, fixed = TRUE)[[1]], .libPaths()))
data_root <- Sys.getenv("FP_DATA_ROOT", "")
output_root <- Sys.getenv("FP_OUTPUT_DIR", "")
if (!nzchar(data_root) || !nzchar(output_root)) {
  stop("Set FP_DATA_ROOT and FP_OUTPUT_DIR; see README.md. Inputs are private and are not included.")
}
data_root <- normalizePath(data_root, mustWork = TRUE)
dir.create(output_root, recursive = TRUE, showWarnings = FALSE)
output_root <- normalizePath(output_root, mustWork = TRUE)
if (output_root == repo_root || startsWith(paste0(output_root, "/"), paste0(repo_root, "/"))) {
  stop("FP_OUTPUT_DIR must be outside this code-only repository.")
}
source(file.path(repo_root, "R", "inputs.R"))
input_overrides <- list()
config_file <- Sys.getenv("FP_CONFIG_FILE", "")
if (nzchar(config_file)) source(normalizePath(config_file, mustWork = TRUE))
used_inputs <- new.env(parent = emptyenv())
input <- function(key, ...) {
  if (!key %in% names(input_paths)) stop("Unknown input key: ", key)
  path <- if (key %in% names(input_overrides)) input_overrides[[key]] else file.path(data_root, input_paths[[key]])
  path <- file.path(path, ...)
  if (!file.exists(path)) stop("Missing private input '", key, "'. Set input_overrides in FP_CONFIG_FILE. Expected: ", path)
  path <- normalizePath(path, mustWork = TRUE)
  used_inputs[[paste(key, ..., sep = "/")]] <- path
  path
}
read_tsv <- function(path, ...) read.delim(path, check.names = FALSE, quote = "", ...)
read_input <- function(key, ...) read_tsv(input(key, ...))
need <- function(...) {
  pkgs <- c(...)
  missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) stop("Missing R packages: ", paste(missing, collapse = ", "), ". Run install_dependencies.R or set FP_R_LIBS.")
  invisible(lapply(pkgs, function(p) suppressPackageStartupMessages(library(p, character.only = TRUE))))
}
check_columns <- function(d, columns) {
  missing <- setdiff(columns, names(d))
  if (length(missing)) stop("Missing input columns: ", paste(missing, collapse = ", "))
  invisible(d)
}
as_bool <- function(x) tolower(as.character(x)) %in% c("true", "t", "1")
display_species <- function(x) {
  x <- sub("^s__", "", x)
  x[x == "Faecalibacterium sp900539945"] <- "Faecalibacterium harmsenii"
  x[x == "Faecalibacterium prausnitzii_F"] <- "Faecalibacterium langellae"
  x
}
short_species <- function(x) sub("^Faecalibacterium ", "F. ", display_species(x))
font_family <- Sys.getenv("FP_FONT", "sans")
figure_theme <- function(base_size = 9) {
  ggplot2::theme_bw(base_size = base_size, base_family = font_family) +
    ggplot2::theme(panel.grid.minor = ggplot2::element_blank(),
                   plot.title = ggplot2::element_text(face = "bold"),
                   legend.key = ggplot2::element_blank())
}
save_plot <- function(p, name, width = 8, height = 6) {
  ggplot2::ggsave(file.path(out_dir, paste0(name, ".pdf")), p, width = width, height = height,
                  device = grDevices::cairo_pdf, family = font_family, limitsize = FALSE)
  if (Sys.getenv("FP_PNG", "1") == "1") ggplot2::ggsave(file.path(out_dir, paste0(name, ".png")), p,
                  width = width, height = height, dpi = 300, bg = "white", limitsize = FALSE)
  invisible(p)
}
write_result <- function(d, name) write.table(d, file.path(out_dir, name), sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
clean_tip <- function(x) sub("\\.fa(\\.tsv)?$", "", gsub("^['\"]|['\"]$", "", x))
read_phylogeny <- function(path) {
  need("ape")
  tr <- ape::read.tree(path)
  tr$tip.label <- clean_tip(tr$tip.label)
  if (anyDuplicated(tr$tip.label)) stop("Duplicate tree-tip labels after normalization")
  tr
}
wilson <- function(k, n) {
  p <- k / n; z <- qnorm(.975); den <- 1 + z^2 / n
  centre <- (p + z^2 / (2 * n)) / den
  half <- z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2)) / den
  data.frame(lower = centre - half, upper = centre + half)
}
