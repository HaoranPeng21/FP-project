# Figure 5B-F. Four ABO/FUT2 groups, carriers-only within-genus CLR, HC3-adjusted models.
# CLR denominator is all 33 named Faecalibacterium rows, with pseudocount 1e-6.
run <- function() {
  need("data.table", "ggplot2", "patchwork", "ggsci")
  clean_sample <- function(x) {
    x <- sub("\\.(fastq|fq)(\\.gz)?$", "", basename(x), ignore.case = TRUE)
    for (pat in c("_kneaddata.*$", "_paired_clean[12]?$", "_paired$", "_clean[12]?$", "_clean$", "_trimmed$", "_R[12]$", "_L00[12]$", "_S[0-9]+$", "_[12]$")) x <- sub(pat, "", x, ignore.case = TRUE)
    x
  }
  normalise_species <- function(x) sub("\\([0-9]+\\)$", "", sub("^s__", "", x))
  # Stream the large abundance matrix in native R; retain only species-level rows.
  con <- file(input("abo_abundance"), "rt"); on.exit(close(con), add = TRUE)
  header <- readLines(con, n = 1); selected <- list()
  repeat {
    lines <- readLines(con, n = 100, warn = FALSE)
    if (!length(lines)) break
    feature <- sub("\t.*", "", lines)
    terminal <- sub("^.*\\|", "", feature)
    take <- grepl("^s__Faecalibacterium", terminal) & !grepl("t__", feature)
    if (any(take)) selected[[length(selected) + 1]] <- lines[take]
  }
  abundance <- data.table::fread(text = paste(c(header, unlist(selected)), collapse = "\n"), data.table = FALSE, check.names = FALSE, quote = "")
  stopifnot(nrow(abundance) == 33, names(abundance)[1] == "clade_name")
  species_all <- normalise_species(sub("^.*\\|", "", abundance[[1]]))
  mat <- as.matrix(abundance[-1]); storage.mode(mat) <- "double"
  samples <- clean_sample(names(abundance)[-1]); stopifnot(!anyDuplicated(samples), !anyNA(mat), all(mat >= 0))
  rownames(mat) <- species_all; colnames(mat) <- samples
  blood <- read.delim(input("abo_phenotype"), check.names = FALSE, na.strings = c("", "NA", "NaN"))
  check_columns(blood, c("sample", "Bloodtype", "FUT2_secretor_status")); stopifnot(!anyDuplicated(blood$sample))
  blood <- blood[blood$sample %in% samples, ]; mat <- mat[, match(blood$sample, samples), drop = FALSE]
  levels <- c("AAB_secretor", "AAB_nonsecretor", "BO_secretor", "BO_nonsecretor")
  blood$group <- NA_character_
  ok <- blood$Bloodtype %in% c("A", "AB", "B", "O") & blood$FUT2_secretor_status %in% c("secretor", "non-secretor")
  blood$group[ok] <- paste0(ifelse(blood$Bloodtype[ok] %in% c("A", "AB"), "AAB_", "BO_"),
    ifelse(blood$FUT2_secretor_status[ok] == "secretor", "secretor", "nonsecretor"))
  blood$group <- factor(blood$group, levels = levels)
  fields <- c("sample_id", "study", "age", "gender", "BMI", "number_reads")
  master <- data.table::fread(input("phenotype_master"), select = fields, data.table = FALSE, na.strings = c("", "NA", "NaN"))
  master <- master[master$sample_id %in% blood$sample, ]; master <- master[order(master$sample_id, master$study != "DMP"), ]
  master <- master[!duplicated(master$sample_id), ]; j <- match(blood$sample, master$sample_id)
  blood$age <- as.numeric(master$age[j]); blood$BMI <- as.numeric(master$BMI[j])
  sex <- tolower(trimws(master$gender[j])); blood$gender <- ifelse(sex %in% c("m", "male", "1", "1.0"), 1, ifelse(sex %in% c("f", "female", "0", "0.0"), 0, NA_real_))
  reads <- suppressWarnings(as.numeric(master$number_reads[j])); blood$log10_reads <- ifelse(is.finite(reads) & reads > 0, log10(reads), NA_real_)
  scope <- read_input("phenotype_scope")
  species <- unique(normalise_species(scope$species_label[scope$full_prevalence > .05]))
  stopifnot(length(species) == 14, all(species %in% rownames(mat)))
  clr <- sweep(log(mat + 1e-6), 2, colMeans(log(mat + 1e-6)), "-")
  long <- do.call(rbind, lapply(species, function(sp) {
    i <- which(rownames(mat) == sp); carrier <- mat[i, ] > 0 & !is.na(blood$group)
    data.frame(species = sp, CLR = clr[i, carrier], blood[carrier, c("group", "age", "gender", "BMI", "log10_reads")])
  }))
  robust <- function(d, pooled = FALSE) {
    d <- d[complete.cases(d), ]; d$group <- factor(d$group, levels = levels)
    if (length(unique(d$group)) < 4 || nrow(d) < 100) stop("Insufficient support for the four-group adjusted model")
    d$focal <- as.integer(d$group == levels[1])
    fit <- if (pooled) lm(CLR ~ focal + age + gender + BMI + log10_reads, d) else lm(CLR ~ group + age + gender + BMI + log10_reads, d)
    X <- model.matrix(fit); bread <- solve(crossprod(X)); omega <- residuals(fit)^2 / pmax(1 - hatvalues(fit), 1e-8)^2
    vc <- bread %*% crossprod(X, X * omega) %*% bread
    terms <- if (pooled) "focal" else paste0("group", levels[-1])
    beta <- coef(fit)[terms] * if (pooled) 1 else -1
    se <- sqrt(diag(vc)[terms])
    data.frame(term = terms, effect = beta, se, lower = beta - 1.96 * se, upper = beta + 1.96 * se,
      p = 2 * pnorm(abs(beta / se), lower.tail = FALSE), n = nrow(d))
  }
  contrasts <- pooled <- tests <- list(); comparisons <- combn(levels, 2, simplify = FALSE)
  for (sp in species) {
    z <- long[long$species == sp, ]; z$species <- NULL
    contrasts[[sp]] <- data.frame(species = sp, robust(z))
    pooled[[sp]] <- data.frame(species = sp, robust(z, TRUE))
    tests[[sp]] <- do.call(rbind, lapply(comparisons, function(g) data.frame(species = sp, group1 = g[1], group2 = g[2],
      p = wilcox.test(z$CLR[z$group == g[1]], z$CLR[z$group == g[2]], exact = FALSE, correct = TRUE)$p.value)))
  }
  contrasts <- do.call(rbind, contrasts); pooled <- do.call(rbind, pooled); tests <- do.call(rbind, tests)
  stopifnot(nrow(contrasts) == 42, nrow(pooled) == 14, nrow(tests) == 84)
  contrasts$q <- p.adjust(contrasts$p, "BH"); pooled$q <- p.adjust(pooled$p, "BH"); tests$q <- p.adjust(tests$p, "BH")
  write_result(contrasts, "adjusted_42_contrasts.tsv"); write_result(pooled, "adjusted_14_pooled_contrasts.tsv"); write_result(tests, "pairwise_84_tests.tsv")
  show <- normalise_species(read_input("abo_species")$species_label)
  gal <- read_input("galnac_prevalence"); gal <- gal[gal$species_gtdb_r226 %in% show, ]; gal <- gal[order(gal$pct_pathway_positive), ]
  ci <- wilson(gal$n_pathway_positive, gal$n_genomes_total)
  gal$label <- factor(short_species(gal$species_gtdb_r226), levels = short_species(gal$species_gtdb_r226))
  c <- ggplot(gal, aes(100 * n_pathway_positive / n_genomes_total, label)) + geom_col(fill = "#C6B58A", width = .7) +
    geom_errorbar(aes(xmin = 100 * ci$lower, xmax = 100 * ci$upper), orientation = "y", width = .2) +
    xlim(0, 100) + labs(title = "A-antigen utilization completeness", x = "Prevalence (%)", y = NULL) + figure_theme(8)
  forest <- pooled[pooled$species %in% show, ]; forest <- forest[order(forest$effect), ]
  forest$label <- factor(short_species(forest$species), levels = short_species(forest$species))
  d <- ggplot(forest, aes(effect, label)) + geom_vline(xintercept = 0, colour = "grey65", linetype = 2) +
    geom_errorbar(aes(xmin = lower, xmax = upper), orientation = "y", width = .2) +
    geom_point(aes(shape = q < .05), size = 2) + scale_shape_manual(values = c(`TRUE` = 16, `FALSE` = 1)) +
    labs(title = "A/AB secretors versus Others", x = "Adjusted CLR difference (HC3 95% CI)", y = NULL, shape = "BH q <0.05") + figure_theme(8)
  box <- function(sp) {
    z <- long[long$species == sp, ]; tt <- tests[tests$species == sp, ]
    span <- diff(range(z$CLR)); if (span == 0) span <- 1
    tt$y <- max(z$CLR) + seq_len(nrow(tt)) * span * .09
    tt$x1 <- match(tt$group1, levels); tt$x2 <- match(tt$group2, levels)
    tt$label <- ifelse(tt$q < .001, "***", ifelse(tt$q < .01, "**", ifelse(tt$q < .05, "*", "ns")))
    ggplot(z, aes(as.numeric(group), CLR, fill = group)) + geom_boxplot(outlier.shape = NA, width = .65, linewidth = .3) +
      geom_segment(data = tt, aes(x = x1, xend = x2, y = y, yend = y), inherit.aes = FALSE, linewidth = .2) +
      geom_text(data = tt, aes(x = (x1 + x2) / 2, y = y, label = label), inherit.aes = FALSE, vjust = -.1, size = 2.3) +
      scale_x_continuous(breaks = 1:4, labels = c("A/AB\nS", "A/AB\nNS", "B/O\nS", "B/O\nNS")) +
      scale_fill_manual(values = ggsci::pal_npg()(10)[1:4]) + labs(title = short_species(sp), x = NULL, y = "CLR abundance") +
      figure_theme(8) + theme(legend.position = "none")
  }
  e <- box("Faecalibacterium longum"); f <- box("Faecalibacterium butyricigenerans")
  save_plot((c + d) / (e + f), "Figure_5CDEF", 10, 10)
  # Figure 5B: fixed biochemical pathway schematic, containing no participant data.
  scheme <- data.frame(x = 0:6,
    substrate = c("A-antigen", "GalNAc", "GalNAc-6P", "GalN-6P", "Tagatose-6P", "Tagatose-1,6-BP", "DHAP + GAP"),
    enzyme = c("GH109", "agaF+V+C+D (PTS)", "nagA", "agaS", "lacC", "gatY/kbaY or gatZ/kbaZ", ""))
  b <- ggplot(scheme, aes(x, 0)) + geom_point(shape = 21, fill = "#E3B65F", size = 7) +
    geom_segment(data = scheme[1:6, ], aes(x = x + .12, xend = x + .88, yend = 0), arrow = grid::arrow(length = grid::unit(2, "mm")), colour = "#E64B35") +
    geom_text(aes(y = -.25, label = substrate), size = 3) +
    geom_text(data = scheme[1:6, ], aes(x = x + .5, y = .3, label = enzyme), size = 2.5) +
    geom_text(data = scheme[1:6, ], aes(x = x + .5, y = -.5, label = paste("Step", x)), size = 2.5) +
    coord_cartesian(ylim = c(-.65, .6), xlim = c(-.5, 6.5), clip = "off") + theme_void()
  save_plot(b, "Figure_5B", 13, 2.3)
}
