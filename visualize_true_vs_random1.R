#!/usr/bin/env Rscript

# 可视化分析：
# 1) true_vs_random_global_totals1.tsv
# 2) true_vs_random_permutation_summary1.tsv

args <- commandArgs(trailingOnly = TRUE)
out_dir <- ifelse(length(args) >= 1, args[1], "analysis_outputs")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

global_file <- "true_vs_random_global_totals1.tsv"
perm_file <- "true_vs_random_permutation_summary1.tsv"

message("Reading: ", global_file)
message("Reading: ", perm_file)

global_df <- read.delim(global_file, header = TRUE, sep = "\t", check.names = FALSE)
perm_df <- read.delim(perm_file, header = TRUE, sep = "\t", check.names = FALSE)

# ---------- 图1：全局 total 对比（true vs random mean + range） ----------
png(file.path(out_dir, "global_totals_comparison.png"), width = 1400, height = 900, res = 140)
par(mar = c(8, 5, 4, 2))
metrics <- as.character(global_df$metric)
x <- seq_along(metrics)
true_vals <- global_df$true_total
rand_mean <- global_df$random_total_mean
bar_centers <- barplot(
  rbind(true_vals, rand_mean),
  beside = TRUE,
  col = c("#1f77b4", "#ff7f0e"),
  names.arg = metrics,
  las = 2,
  ylab = "Total value",
  main = "Global totals: true vs random mean"
)
# 给 random mean 柱子加 min/max 误差线
rand_bars <- bar_centers[2, ]
arrows(
  x0 = rand_bars,
  y0 = global_df$random_total_min,
  x1 = rand_bars,
  y1 = global_df$random_total_max,
  angle = 90,
  code = 3,
  length = 0.05,
  lwd = 1.8
)
legend("topright", legend = c("true_total", "random_total_mean", "random_total_range"),
       fill = c("#1f77b4", "#ff7f0e", NA), border = c("black", "black", NA),
       lty = c(0, 0, 1), col = c(NA, NA, "black"), bty = "n")
dev.off()

# ---------- 从 permutation summary 自动识别 p-value 列 ----------
p_cols <- grep("_p_two_sided$", colnames(perm_df), value = TRUE)
if (length(p_cols) == 0) {
  stop("No *_p_two_sided columns found in permutation summary.")
}

# ---------- 图2：各指标双侧 p 值分布 ----------
n_plot <- length(p_cols)
nr <- ceiling(sqrt(n_plot))
nc <- ceiling(n_plot / nr)
png(file.path(out_dir, "permutation_pvalue_histograms.png"), width = 1500, height = 1100, res = 140)
par(mfrow = c(nr, nc), mar = c(4, 4, 3, 1))
for (pc in p_cols) {
  pv <- perm_df[[pc]]
  pv <- pv[is.finite(pv)]
  hist(
    pv,
    breaks = seq(0, 1, by = 0.05),
    col = "#6baed6",
    border = "white",
    main = gsub("_p_two_sided", "", pc),
    xlab = "two-sided p-value"
  )
  abline(v = 0.05, col = "red", lty = 2, lwd = 2)
}
dev.off()

# ---------- 图3：显著比例（p < 0.05） ----------
sig_rate <- sapply(p_cols, function(pc) {
  pv <- perm_df[[pc]]
  pv <- pv[is.finite(pv)]
  if (length(pv) == 0) return(NA_real_)
  mean(pv < 0.05)
})

png(file.path(out_dir, "permutation_significant_rate.png"), width = 1200, height = 800, res = 140)
par(mar = c(8, 5, 4, 2))
barplot(
  sig_rate,
  names.arg = gsub("_p_two_sided", "", p_cols),
  las = 2,
  col = "#31a354",
  ylim = c(0, max(sig_rate, na.rm = TRUE) * 1.25),
  ylab = "Proportion with p < 0.05",
  main = "Permutation summary: significant proportion by metric"
)
abline(h = 0.05, col = "red", lty = 2, lwd = 2)
text(x = seq_along(sig_rate), y = sig_rate, labels = sprintf("%.3f", sig_rate), pos = 3, cex = 0.9)
dev.off()

# ---------- 输出简要统计 ----------
summary_file <- file.path(out_dir, "analysis_summary.txt")
con <- file(summary_file, open = "wt")
writeLines("[Global totals fold-change: true_total / random_total_mean]", con)
for (i in seq_len(nrow(global_df))) {
  fc <- global_df$true_total[i] / global_df$random_total_mean[i]
  line <- sprintf("%s\tfold_change=%.4f\ttrue=%.6f\trandom_mean=%.6f\tp_two_sided=%s",
                  global_df$metric[i], fc, global_df$true_total[i], global_df$random_total_mean[i], global_df$p_two_sided[i])
  writeLines(line, con)
}
writeLines("", con)
writeLines("[Permutation two-sided p-value summary]", con)
for (pc in p_cols) {
  pv <- perm_df[[pc]]
  pv <- pv[is.finite(pv)]
  line <- sprintf("%s\tn=%d\tmedian=%.6f\tprop_p_lt_0.05=%.6f",
                  pc, length(pv), median(pv), mean(pv < 0.05))
  writeLines(line, con)
}
close(con)

message("Done. Outputs written to: ", out_dir)
