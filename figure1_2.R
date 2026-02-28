library(readr)
library(dplyr)

# 读入文件（路径自己换）
df <- read.delim(
  "D:/GUOGUO_project/gsmap/R/articlefigure/mouse_class_cluster_total.txt",
  sep = "\t",
  stringsAsFactors = FALSE,
  check.names = FALSE
)
df <- df[df$cell_Neuron_type %in% c("Gaba","Glut"), ]
df <- df %>%
  mutate(label = row_number())

library(dplyr)
library(stringr)
library(purrr)

## ===============================
## 1. 判断一行 region 是否包含不同子 region
## ===============================
has_different_region <- function(region_str) {
  
  if (is.na(region_str) || region_str == "") return(FALSE)
  
  ## 按逗号拆分并去空格
  parts <- str_trim(unlist(strsplit(region_str, ",")))
  
  ## 每 3 个拼成一个完整 region
  regions <- map_chr(
    split(parts, ceiling(seq_along(parts) / 3)),
    ~ paste(.x, collapse = ", ")
  )
  
  ## 是否存在不同 region
  length(unique(regions)) > 1
}

## ===============================
## 2. 应用到整个数据框
## ===============================
df_diff_region <- df %>%
  mutate(has_diff_region = map_lgl(region, has_different_region)) %>%
  filter(has_diff_region) %>%
  select(-has_diff_region)

## ===============================
## 3. 查看结果
## ===============================
head(df_diff_region)
df<-df_diff_region
res <- df %>%
  group_by(layer, cell_Neuron_type) %>%
  summarise(
    cell_count = n(),
    .groups = "drop"
  )
library(dplyr)
library(tidyr)
library(ggplot2)

## 计算每个 layer 的总数
df_plot <- res %>%
  group_by(layer) %>%
  mutate(
    total_cell = sum(cell_count)
  ) %>%
  ungroup() %>%
  mutate(
    proportion = cell_count / total_cell
  )
df_use <- df_plot %>%
  filter(cell_Neuron_type %in% c("Gaba", "Glut"))
max_count <- max(df_use$cell_count)

df_use <- df_use %>%
  mutate(
    prop_scaled = proportion * max_count
  )
fill_colors <- c(
  Gaba = "#F4A582",   # 浅红
  Glut = "#92C5DE"    # 浅蓝
)
line_colors <- c(
  Gaba = "#B2182B",   # 深红
  Glut = "#2166AC"    # 深蓝
)
p <- ggplot(df_use, aes(x = layer)) +
  
  ## ========= 柱状图（数量，左右分） =========
geom_col(
  aes(
    y = cell_count,
    fill = cell_Neuron_type
  ),
  position = position_dodge(width = 0.8),
  width = 0.65,
  color = "black",
  linewidth = 0.25
) +
  
  ## ========= 折线（比例，居中） =========
geom_line(
  aes(
    y = prop_scaled,
    color = cell_Neuron_type,
    group = cell_Neuron_type
  ),
  linewidth = 1.1
) +
  
  ## ========= 点（比例，居中） =========
geom_point(
  aes(
    y = prop_scaled,
    color = cell_Neuron_type
  ),
  size = 2.6
) +
  
  ## ========= 双 Y 轴 =========
scale_y_continuous(
  name = "Region count",
  expand = expansion(mult = c(0, 0.05)),
  sec.axis = sec_axis(
    ~ . / max_count,
    name = "Proportion of total cells"
  )
) +
  
  ## ========= 配色 =========
scale_fill_manual(
  values = c(
    Gaba = "#F4A582",
    Glut = "#92C5DE"
  ),
  name = "Cell type"
) +
  scale_color_manual(
    values = c(
      Gaba = "#B2182B",
      Glut = "#2166AC"
    ),
    name = "Cell type"
  ) +
  
  ## ========= 论文风格 =========
theme_classic(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "right",
    axis.line = element_line(linewidth = 0.6),
    axis.ticks = element_line(linewidth = 0.6)
  ) +
  
  labs(
    x = "Cortical layer",
    title = "GABAergic and Glutamatergic cell abundance and proportion by layer"
  )

p


suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
})

## =========================================
## 1. 先给 df 标记是否为跨区域簇
## =========================================

df <- df %>%
  mutate(
    is_diff_region = label %in% df_diff_region$label
  )

## =========================================
## 2. 计算每个 slide × layer 的 OR
## =========================================

or_df <- df %>%
  group_by(slide, layer) %>%
  summarise(
    a = sum(is_diff_region),                # 跨区域簇（该 layer）
    A = n(),                                # 所有簇（该 layer）
    .groups = "drop"
  ) %>%
  left_join(
    df %>%
      group_by(slide) %>%
      summarise(
        total_inter = sum(is_diff_region),
        total_all   = n(),
        .groups = "drop"
      ),
    by = "slide"
  ) %>%
  mutate(
    b = A - a,
    c = total_inter - a,
    d = (total_all - total_inter) - b,
    
    ## Haldane correction
    odds_ratio = ((a + 0.5) * (d + 0.5)) /
      ((b + 0.5) * (c + 0.5))
  )

## =========================================
## 3. 计算 log(OR) 的残差
##    （用于直方图）
## =========================================

or_df <- or_df %>%
  mutate(
    log_or = log(odds_ratio),
    resid  = log_or - mean(log_or, na.rm = TRUE)
  )

## =========================================
## 4. 画 OR vs layer（每个 slide 一条）
## =========================================
or_layer_sum <- or_df %>%
  group_by(layer) %>%
  summarise(
    mean_log_or = mean(log_or, na.rm = TRUE),
    sd_log_or   = sd(log_or, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    mean_or = exp(mean_log_or),
    lower   = exp(mean_log_or - 2 * sd_log_or),
    upper   = exp(mean_log_or + 2 * sd_log_or)
  )
p_or <- ggplot(
  or_layer_sum,
  aes(x = layer, y = mean_or)
) +
  geom_errorbar(
    aes(ymin = lower, ymax = upper),
    width = 0.2,
    linewidth = 0.8,
    color = "black"
  ) +
  geom_point(
    size = 3,
    shape = 21,
    fill = "white",
    stroke = 1
  ) +
  geom_text(
    aes(label = sprintf("%.2f", mean_or)),
    vjust = -0.8,
    size = 3
  ) +
  geom_hline(
    yintercept = 1,
    linetype = "dashed",
    color = "red"
  ) +
  scale_y_log10() +
  labs(
    title="GABA",
    x = "Layer",
    y = "Odds ratio (log scale)"
  ) +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

print(p_or)


