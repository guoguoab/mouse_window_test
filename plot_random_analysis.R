#!/usr/bin/env Rscript

required_pkgs <- c("ggplot2", "dplyr", "tidyr", "forcats", "scales")
missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs) > 0) {
  stop(
    sprintf(
      "缺少R包: %s\n请先安装后重试，例如 install.packages(c(%s))",
      paste(missing_pkgs, collapse = ", "),
      paste(sprintf('"%s"', missing_pkgs), collapse = ", ")
    )
  )
}

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(forcats)
  library(scales)
})

presence_file <- "true_loc_presence_by_sample.tsv"
summary_file <- "all_sample_mean_summary.txt"
out_dir <- "plots_random_summary"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ---------------------------
# 图1：真实merge region在随机样本中的聚合保留情况
# ---------------------------
if (file.exists(presence_file)) {
  presence_df <- read.delim(presence_file, sep = "\t", check.names = FALSE)

  id_cols <- c("true_file", "true_file_path", "chip", "loc")
  sample_cols <- setdiff(names(presence_df), id_cols)

  # 尝试自动识别二值列，避免混入其他注释列
  sample_cols <- sample_cols[vapply(presence_df[sample_cols], function(x) {
    all(is.na(x) | x %in% c(0, 1))
  }, logical(1))]

  if (length(sample_cols) == 0) {
    warning("未在 true_loc_presence_by_sample.tsv 中检测到0/1样本列，跳过图1。")
  } else {
    # 50%截断：统计每个真实loc在所有sample中的匹配比例
    loc_match_stats <- presence_df %>%
      mutate(match_rate = rowMeans(across(all_of(sample_cols)), na.rm = TRUE)) %>%
      mutate(match_group = ifelse(match_rate >= 0.5, ">=50%", "<50%"))

    match_count_df <- loc_match_stats %>%
      count(match_group, name = "true_region_count") %>%
      complete(match_group = c(">=50%", "<50%"), fill = list(true_region_count = 0L)) %>%
      mutate(percent = true_region_count / sum(true_region_count))

    write.table(
      match_count_df,
      file = file.path(out_dir, "plot1_threshold50_region_count.tsv"),
      sep = "\t", row.names = FALSE, quote = FALSE
    )

    sample_presence <- presence_df %>%
      summarise(across(all_of(sample_cols), ~ mean(.x, na.rm = TRUE))) %>%
      pivot_longer(cols = everything(), names_to = "sample", values_to = "presence_rate") %>%
      mutate(sample = fct_reorder(sample, presence_rate))

    p1 <- ggplot(sample_presence, aes(x = sample, y = presence_rate, fill = presence_rate)) +
      geom_col(width = 0.72, color = "grey25") +
      coord_flip() +
      scale_fill_gradientn(colors = c("#E8F1FA", "#77ADD6", "#1F78B4")) +
      scale_y_continuous(labels = percent_format(accuracy = 1), expand = expansion(mult = c(0, 0.03))) +
      labs(
        title = "真实 merge regions 在随机样本中的保留比例",
        subtitle = "每个sample中被命中的loc占比（命中=1）",
        x = NULL,
        y = "保留比例"
      ) +
      theme_minimal(base_size = 12) +
      theme(
        legend.position = "none",
        panel.grid.major.y = element_blank(),
        plot.title = element_text(face = "bold")
      )

    ggsave(file.path(out_dir, "plot1_true_loc_presence_rate_by_sample.png"), p1,
      width = 9, height = 7, dpi = 300
    )

    # 50%阈值统计图
    p1_threshold <- ggplot(match_count_df, aes(x = match_group, y = true_region_count, fill = match_group)) +
      geom_col(width = 0.62, color = "grey25") +
      geom_text(aes(label = paste0(true_region_count, " (", percent(percent, accuracy = 0.1), ")")),
                vjust = -0.25, size = 4.2) +
      scale_fill_manual(values = c(">=50%" = "#2C7FB8", "<50%" = "#F39C6B")) +
      labs(
        title = "真实区域按50%匹配率阈值分组",
        subtitle = ">=50%：在至少一半随机sample中可匹配；<50%：匹配不稳定",
        x = "匹配率分组",
        y = "真实区域数量"
      ) +
      theme_minimal(base_size = 12) +
      theme(
        legend.position = "none",
        plot.title = element_text(face = "bold")
      ) +
      expand_limits(y = max(match_count_df$true_region_count) * 1.14)

    ggsave(file.path(out_dir, "plot1c_threshold50_region_count.png"), p1_threshold,
      width = 8, height = 5.6, dpi = 300
    )

    # 附加热图：展示出现频率最高的loc在不同sample中的分布
    top_n <- min(80, nrow(presence_df))
    loc_order <- presence_df %>%
      mutate(hit_sum = rowSums(across(all_of(sample_cols)), na.rm = TRUE)) %>%
      arrange(desc(hit_sum)) %>%
      slice_head(n = top_n) %>%
      mutate(loc_label = paste0(chip, " | ", loc))

    heat_df <- loc_order %>%
      select(loc_label, all_of(sample_cols)) %>%
      pivot_longer(cols = all_of(sample_cols), names_to = "sample", values_to = "present") %>%
      mutate(
        present = factor(present, levels = c(0, 1)),
        loc_label = factor(loc_label, levels = rev(unique(loc_order$loc_label)))
      )

    p1_heat <- ggplot(heat_df, aes(x = sample, y = loc_label, fill = present)) +
      geom_tile(color = "white", size = 0.1) +
      scale_fill_manual(values = c("0" = "#F1F1F1", "1" = "#2C7FB8"), na.value = "#F8F8F8") +
      labs(
        title = sprintf("命中频率最高的前%d个loc跨sample分布", top_n),
        x = NULL,
        y = "chip | loc",
        fill = "是否命中"
      ) +
      theme_minimal(base_size = 11) +
      theme(
        axis.text.y = element_text(size = 6),
        panel.grid = element_blank(),
        plot.title = element_text(face = "bold")
      )

    ggsave(file.path(out_dir, "plot1b_true_loc_top_heatmap.png"), p1_heat,
      width = 12, height = 10, dpi = 300
    )
  }
} else {
  warning("未找到 true_loc_presence_by_sample.tsv，图1已跳过。\n可先运行: Rscript compare_true_locs_in_data1.R true_data data1 true_loc_presence_by_sample.tsv")
}

# ---------------------------
# 图2：随机区域参数汇总（all_sample_mean_summary）
# ---------------------------
if (!file.exists(summary_file)) {
  stop("未找到 all_sample_mean_summary.txt")
}

summary_df <- read.delim(summary_file, sep = "\t", check.names = FALSE)

metrics <- c(
  "mean_region_num", "mean_enrich_num", "mean_glut_num",
  "mean_gaba_num", "mean_glut_gaba_ratio"
)

missing_metrics <- setdiff(metrics, names(summary_df))
if (length(missing_metrics) > 0) {
  stop(sprintf("all_sample_mean_summary.txt 缺少字段: %s", paste(missing_metrics, collapse = ", ")))
}

summary_long <- summary_df %>%
  select(slide, subclass, layer, all_of(metrics)) %>%
  pivot_longer(cols = all_of(metrics), names_to = "metric", values_to = "value") %>%
  mutate(
    layer = fct_infreq(layer),
    metric = factor(metric, levels = metrics)
  )

p2 <- ggplot(summary_long, aes(x = layer, y = value, fill = layer)) +
  geom_boxplot(width = 0.7, outlier.alpha = 0.18, color = "grey30") +
  geom_jitter(width = 0.2, alpha = 0.22, size = 0.8, color = "grey20") +
  facet_wrap(~metric, scales = "free_y", ncol = 2,
             labeller = as_labeller(c(
               mean_region_num = "Mean region num",
               mean_enrich_num = "Mean enrich num",
               mean_glut_num = "Mean glut num",
               mean_gaba_num = "Mean gaba num",
               mean_glut_gaba_ratio = "Mean glut/gaba ratio"
             ))) +
  scale_fill_brewer(palette = "Set2") +
  labs(
    title = "随机情况下区域参数分布（按layer）",
    subtitle = "all_sample_mean_summary.txt",
    x = NULL,
    y = "数值"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 20, hjust = 1),
    plot.title = element_text(face = "bold")
  )

ggsave(file.path(out_dir, "plot2_random_region_parameter_distribution.png"), p2,
  width = 12, height = 10, dpi = 320
)

# 附加：均值热图（按 subclass × metric）
heat_df2 <- summary_long %>%
  group_by(subclass, metric) %>%
  summarise(value = mean(value, na.rm = TRUE), .groups = "drop") %>%
  group_by(subclass) %>%
  summarise(total = sum(value, na.rm = TRUE), .groups = "drop") %>%
  slice_max(total, n = 35) %>%
  inner_join(summary_long, by = "subclass") %>%
  group_by(subclass, metric) %>%
  summarise(value = mean(value, na.rm = TRUE), .groups = "drop") %>%
  mutate(subclass = fct_reorder(subclass, value, .fun = sum))

p2_heat <- ggplot(heat_df2, aes(x = metric, y = subclass, fill = value)) +
  geom_tile(color = "white", size = 0.2) +
  scale_fill_gradientn(colors = c("#F7FBFF", "#6BAED6", "#08306B"), trans = "sqrt") +
  labs(
    title = "Top subclass参数热图（均值）",
    x = NULL,
    y = "subclass",
    fill = "均值"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    axis.text.x = element_text(angle = 20, hjust = 1),
    axis.text.y = element_text(size = 7),
    panel.grid = element_blank(),
    plot.title = element_text(face = "bold")
  )

ggsave(file.path(out_dir, "plot2b_random_region_parameter_heatmap.png"), p2_heat,
  width = 11, height = 10, dpi = 320
)

message("绘图完成。输出目录: ", out_dir)
