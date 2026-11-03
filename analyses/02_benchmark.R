# Figure 1B-D and Figure S2D-E. Frozen strict83 run selection; no resampling.
# The supplied PDF's 'same21' caption is stale: its TP/TN counts identify shared24.
run <- function() {
  need("ggplot2", "patchwork", "tidyr", "dplyr", "scales")
  all83 <- read_input("benchmark_all83")
  shared <- read_input("benchmark_shared24")
  runs <- read_input("benchmark_runs")
  stopifnot(length(unique(runs$run_id)) == 901,
    setequal(all83$run_id, runs$run_id), setequal(shared$run_id, runs$run_id),
    nrow(all83) == 901 * 83 * 2, nrow(shared) == 901 * 24 * 2,
    !anyDuplicated(all83[c("run_id", "tool", "species_id")]),
    !anyDuplicated(shared[c("run_id", "tool", "species_id")]))
  metrics <- function(d, scope) {
    truth <- as_bool(d$truth_present); detected <- d$observed_percent > 0
    tp <- sum(truth & detected); fn <- sum(truth & !detected)
    fp <- sum(!truth & detected); tn <- sum(!truth & !detected)
    data.frame(tool = d$tool[1], scope = scope, TP = tp, FN = fn, FP = fp, TN = tn,
      recall = tp / (tp + fn), precision = tp / (tp + fp), FPR = fp / (fp + tn),
      rho = cor(d$expected_percent[truth], d$observed_percent[truth], method = "spearman"),
      MAE = mean(abs(d$observed_percent[truth] - d$expected_percent[truth])))
  }
  shared_metrics <- do.call(rbind, lapply(split(shared, shared$tool), metrics, scope = "shared24"))
  sm <- shared_metrics[shared_metrics$tool == "Sylph", ]
  mm <- shared_metrics[shared_metrics$tool == "MetaPhlAn4", ]
  stopifnot(sm$TP == 1293, sm$FN == 7, sm$FP == 0, sm$TN == 20324,
    mm$TP == 1299, mm$FN == 1, mm$FP == 849, mm$TN == 19475)
  sylph <- all83[all83$tool == "Sylph", ]
  full_metrics <- metrics(sylph, "all83")
  stopifnot(full_metrics$TP == 4453, full_metrics$FN == 47, full_metrics$FP == 1, full_metrics$TN == 70282)
  write_result(rbind(full_metrics, shared_metrics), "benchmark_metrics.tsv")
  cols <- c(Sylph = "#4DBBD5", MetaPhlAn4 = "#00A087")
  scatter <- function(d, title) {
    m <- metrics(d, title); proxy <- 2e-4
    x <- pmax(d$expected_percent, proxy); y <- pmax(d$observed_percent, proxy)
    ggplot(data.frame(x, y), aes(x, y)) +
      geom_abline(slope = 1, intercept = 0, colour = "grey55", linetype = 2, linewidth = .3) +
      geom_point(colour = cols[[d$tool[1]]], size = .6, alpha = .35) +
      scale_x_log10(breaks = c(proxy, .001, .01, .1, 1, 10), labels = c("0", ".001", ".01", ".1", "1", "10")) +
      scale_y_log10(breaks = c(proxy, .001, .01, .1, 1, 10), labels = c("0", ".001", ".01", ".1", "1", "10")) +
      annotate("text", x = .03, y = .00035, label = sprintf("TP %s  FN %s\nFP %s  TN %s", m$TP, m$FN, m$FP, m$TN), size = 2.6) +
      labs(title = title, x = "Expected relative abundance (%)", y = "Reported relative abundance (%)") +
      coord_equal() + figure_theme(8)
  }
  b <- scatter(sylph, "Sylph · 83 groups")
  s <- shared[shared$tool == "Sylph", ]; mp <- shared[shared$tool == "MetaPhlAn4", ]
  sp <- mp[grepl("prausnitzii$", sub("^s__", "", mp$standardized_name)), ]
  stopifnot(nrow(sp) == 901)
  focal_metrics <- metrics(sp, "SGB15332")
  stopifnot(focal_metrics$TP == 63, focal_metrics$FN == 0, focal_metrics$FP == 581, focal_metrics$TN == 257)
  d <- scatter(sp, "F. prausnitzii · SGB15332")
  mm <- tidyr::pivot_longer(shared_metrics, c("recall", "precision", "FPR", "rho", "MAE"), names_to = "metric", values_to = "value")
  mm$metric <- factor(mm$metric, levels = c("recall", "precision", "FPR", "rho", "MAE"))
  mm$tool <- factor(mm$tool, levels = c("Sylph", "MetaPhlAn4"))
  c <- ggplot(mm, aes(metric, value, fill = tool)) + geom_col(position = "dodge", width = .75) +
    geom_text(aes(label = sprintf("%.3f", value)), position = position_dodge(.75), angle = 90, hjust = -.1, vjust = .5, size = 2.4) +
    scale_fill_manual(values = cols) + scale_y_continuous(limits = c(0, 1.14), expand = c(0, 0)) +
    scale_x_discrete(labels = c(recall = "Rec.", precision = "Prec.", FPR = "FPR", rho = "\u03c1", MAE = "MAE (pp)")) +
    labs(x = NULL, y = "Value", fill = NULL) + figure_theme(8) + theme(legend.position = "top")
  save_plot(b, "Figure_1B", 3.7, 3.4); save_plot(c, "Figure_1C", 4, 3.4); save_plot(d, "Figure_1D", 3.7, 3.4)
  save_plot(b + c + d, "Figure_1BCD", 11.7, 3.5)
  save_plot(scatter(s, "Sylph · shared24") + scatter(mp, "MetaPhlAn4 · shared24"), "Figure_S2DE", 8, 3.8)
}
