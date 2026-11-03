# Figure 1B-D and Figure S2D-E. Frozen strict83 run selection; no resampling.
# The supplied PDF's 'same21' caption is stale: its TP/TN counts identify shared24.


data_dir <- Sys.getenv("FP_DATA_DIR", "/path/to/private/data/02_benchmark")
out_dir <- Sys.getenv("FP_OUTPUT_DIR", "/path/to/private/results/02_benchmark")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
font_family <- "sans"
library(ggplot2)
library(patchwork)
library(tidyr)
library(dplyr)
library(scales)

all83 <- read.delim(file.path(data_dir, "benchmark_all83.tsv"), check.names = FALSE, quote = "")

shared <- read.delim(file.path(data_dir, "benchmark_shared24.tsv"), check.names = FALSE, quote = "")

runs <- read.delim(file.path(data_dir, "benchmark_runs.tsv"), check.names = FALSE, quote = "")

stopifnot(length(unique(runs$run_id)) == 901, setequal(all83$run_id, runs$run_id), setequal(shared$run_id,
    runs$run_id), nrow(all83) == 901 * 83 * 2, nrow(shared) == 901 * 24 * 2, !anyDuplicated(all83[c("run_id",
    "tool", "species_id")]), !anyDuplicated(shared[c("run_id", "tool", "species_id")]))
sylph <- all83[all83$tool == "Sylph", ]
s <- shared[shared$tool == "Sylph", ]
mp <- shared[shared$tool == "MetaPhlAn4", ]
sp <- mp[grepl("prausnitzii$", sub("^s__", "", mp$standardized_name)), ]
stopifnot(nrow(sp) == 901)
datasets <- list(all83 = sylph, shared24_Sylph = s, shared24_MetaPhlAn4 = mp, SGB15332 = sp)
metrics_tables <- vector("list", 4)

for (i in 1:4) {
    d <- datasets[[i]]
    scope <- if (i %in% c(2, 3))
        "shared24"
    else names(datasets)[i]
    truth <- tolower(as.character(d$truth_present)) %in% c("true", "t", "1")
    detected <- d$observed_percent > 0
    tp <- sum(truth & detected)
    fn <- sum(truth & !detected)
    fp <- sum(!truth & detected)
    tn <- sum(!truth & !detected)
    metrics_tables[[i]] <- data.frame(tool = d$tool[1], scope = scope, TP = tp, FN = fn, FP = fp, TN = tn,
        recall = tp/(tp + fn), precision = tp/(tp + fp), FPR = fp/(fp + tn), rho = cor(d$expected_percent[truth],
            d$observed_percent[truth], method = "spearman"), MAE = mean(abs(d$observed_percent[truth] -
            d$expected_percent[truth])))
}

full_metrics <- metrics_tables[[1]]
shared_metrics <- rbind(metrics_tables[[2]], metrics_tables[[3]])
focal_metrics <- metrics_tables[[4]]

stopifnot(full_metrics$TP == 4453, full_metrics$FN == 47, full_metrics$FP == 1, full_metrics$TN == 70282,
    shared_metrics$TP == c(1293, 1299), shared_metrics$FN == c(7, 1), shared_metrics$FP == c(0, 849),
    shared_metrics$TN == c(20324, 19475), focal_metrics$TP == 63, focal_metrics$FN == 0, focal_metrics$FP ==
        581, focal_metrics$TN == 257)

write.table(rbind(full_metrics, shared_metrics), file.path(out_dir, "benchmark_metrics.tsv"), sep = "\t",
    quote = FALSE, row.names = FALSE)
cols <- c(Sylph = "#4DBBD5", MetaPhlAn4 = "#00A087")
plot_order <- c(1, 4, 2, 3)

plot_titles <- c("Sylph · 83 species", "F. prausnitzii · SGB15332", "Sylph · shared24", "MetaPhlAn4 · shared24")

plots <- vector("list", 4)

for (i in 1:4) {
    d <- datasets[[plot_order[i]]]
    title <- plot_titles[i]
    m <- metrics_tables[[plot_order[i]]]
    proxy <- 2e-04
    proxy <- 2e-04
    x <- pmax(d$expected_percent, proxy)
    y <- pmax(d$observed_percent, proxy)
    plots[[i]] <- ggplot(data.frame(x, y), aes(x, y)) + geom_abline(slope = 1, intercept = 0, colour = "grey55",
        linetype = 2, linewidth = 0.3) + geom_point(colour = cols[[d$tool[1]]], size = 0.6, alpha = 0.35) +
        scale_x_log10(breaks = c(proxy, 0.001, 0.01, 0.1, 1, 10), labels = c("0", ".001", ".01", ".1",
            "1", "10")) + scale_y_log10(breaks = c(proxy, 0.001, 0.01, 0.1, 1, 10), labels = c("0", ".001",
        ".01", ".1", "1", "10")) + annotate("text", x = 0.03, y = 0.00035, label = sprintf("TP %s  FN %s\nFP %s  TN %s",
        m$TP, m$FN, m$FP, m$TN), size = 2.6) + labs(title = title, x = "Expected relative abundance (%)",
        y = "Reported relative abundance (%)") + coord_equal() + (ggplot2::theme_bw(base_size = 8, base_family = "sans") +
        ggplot2::theme(panel.grid.minor = ggplot2::element_blank(), legend.key = ggplot2::element_blank()))
}

b <- plots[[1]]
d <- plots[[2]]

mm <- tidyr::pivot_longer(shared_metrics, c("recall", "precision", "FPR", "rho", "MAE"), names_to = "metric",
    values_to = "value")
mm$metric <- factor(mm$metric, levels = c("recall", "precision", "FPR", "rho", "MAE"))
mm$tool <- factor(mm$tool, levels = c("Sylph", "MetaPhlAn4"))

c <- ggplot(mm, aes(metric, value, fill = tool)) + geom_col(position = "dodge", width = 0.75) + geom_text(aes(label = sprintf("%.3f",
    value)), position = position_dodge(0.75), angle = 90, hjust = -0.1, vjust = 0.5, size = 2.4) + scale_fill_manual(values = cols) +
    scale_y_continuous(limits = c(0, 1.14), expand = c(0, 0)) + scale_x_discrete(labels = c(recall = "Rec.",
    precision = "Prec.", FPR = "FPR", rho = "ρ", MAE = "MAE (pp)")) + labs(y = "Value") + (ggplot2::theme_bw(base_size = 8,
    base_family = "sans") + ggplot2::theme(panel.grid.minor = ggplot2::element_blank(), legend.key = ggplot2::element_blank())) +
    theme(legend.position = "top")

ggplot2::ggsave(filename = file.path(out_dir, "Figure_1B.pdf"), plot = b, width = 3.7, height = 3.4,
    device = grDevices::cairo_pdf, limitsize = FALSE)

ggplot2::ggsave(filename = file.path(out_dir, "Figure_1C.pdf"), plot = c, width = 4, height = 3.4, device = grDevices::cairo_pdf,
    limitsize = FALSE)

ggplot2::ggsave(filename = file.path(out_dir, "Figure_1D.pdf"), plot = d, width = 3.7, height = 3.4,
    device = grDevices::cairo_pdf, limitsize = FALSE)

ggplot2::ggsave(filename = file.path(out_dir, "Figure_1BCD.pdf"), plot = b + c + d, width = 11.7, height = 3.5,
    device = grDevices::cairo_pdf, limitsize = FALSE)

ggsave(file.path(out_dir, "Figure_S2DE.pdf"), plots[[3]] + plots[[4]], width = 8, height = 3.8, device = cairo_pdf)
