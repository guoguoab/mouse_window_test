library(dplyr)
library(ggplot2)

## ===============================
## 1. 读入数据
## ===============================
df <- read.delim(
  "D:/GUOGUO_project/gsmap/R/articlefigure/layer_cell_count.txt",
  sep = "\t",
  stringsAsFactors = FALSE,
  check.names = FALSE
)

df$layer <- gsub("'", "", df$layer)

## ===============================
## 2. 按 layer 汇总
## ===============================
df_sum <- df %>%
  group_by(layer) %>%
  summarise(
    rate_mean = mean(rate),
    rate_sd   = sd(rate),
    b_sum     = sum(b),
    .groups = "drop"
  )

## ===============================
## 3. 绘图（美化版）
## ===============================
p <- ggplot(df_sum, aes(x = layer, y = rate_mean)) +
  
  ## 柱子
  geom_col(
    width = 0.6,
    fill  = "#4DBBD5",
    color = "black",
    linewidth = 0.6
  ) +
  
  ## 误差线
  geom_errorbar(
    aes(
      ymin = rate_mean - rate_sd,
      ymax = rate_mean + rate_sd
    ),
    width = 0.18,
    linewidth = 0.8
  ) +
  
  ## y = 1 参考线
  geom_hline(
    yintercept = 1,
    linetype = "dashed",
    linewidth = 0.8,
    color = "grey40"
  ) +
  
  ## b 总数：从柱子底部竖着标注
  geom_text(
    aes(label = b_sum),
    y = 0.02,
    angle = 90,
    vjust = 0,
    size = 4,
    color = "black"
  ) +
  
  ## 坐标与标签
  scale_y_continuous(
    limits = c(0, max(df_sum$rate_mean + df_sum$rate_sd) * 1.15),
    expand = expansion(mult = c(0, 0.05))
  ) +
  
  labs(
    x = "Cortical layer",
    y = "Mean rate",
    title = "Rate by Layer",
    subtitle = "Bars show mean ± SD; numbers indicate total cell count (b)"
  ) +
  
  ## 主题美化
  theme_classic(base_size = 14) +
  theme(
    plot.title      = element_text(face = "bold", size = 16),
    plot.subtitle   = element_text(size = 12, color = "grey30"),
    axis.title.x    = element_text(face = "bold"),
    axis.title.y    = element_text(face = "bold"),
    axis.text       = element_text(color = "black"),
    axis.line       = element_line(linewidth = 0.8),
    axis.ticks      = element_line(linewidth = 0.8)
  )

print(p)
