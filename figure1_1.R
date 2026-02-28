library(dplyr)
library(ggplot2)
library(stringr)
library(readr)
library(tidyr)
library(colorspace)

input_dir <- "D:/GUOGUO_project/gsmap/R/articlefigure/figure1"

files <- list.files(
  input_dir,
  pattern = "\\.txt$",
  full.names = TRUE
)

## -------------------------------------------------
## 1. 读入数据
## -------------------------------------------------
df_all <- lapply(files, function(f) {
  df <- read.delim(f)
  df$section_id <- tools::file_path_sans_ext(basename(f))
  df
}) %>%
  bind_rows()

## -------------------------------------------------
## 2. 转成长表
## -------------------------------------------------
df_long <- df_all %>%
  pivot_longer(
    cols = c(proportion_glut, proportion_gaba),
    names_to = "type",
    values_to = "proportion"
  ) %>%
  mutate(
    type = recode(
      type,
      proportion_glut = "Glut",
      proportion_gaba = "GABA"
    )
  )

## -------------------------------------------------
## 3. 固定 CCF 顺序
## -------------------------------------------------
region_order <- unique(df_long$ccf_region_name1)

df_long <- df_long %>%
  mutate(ccf_region_name1 = factor(ccf_region_name1, levels = region_order))

## -------------------------------------------------
## 4. 内部编号
## -------------------------------------------------
df_long <- df_long %>%
  group_by(ccf_region_name1, type) %>%
  mutate(section_index = as.numeric(factor(brain_section_label))) %>%
  ungroup()

## -------------------------------------------------
## 5. 区域宽度
## -------------------------------------------------
region_width <- df_long %>%
  group_by(ccf_region_name1, type) %>%
  summarise(n_section = max(section_index), .groups = "drop") %>%
  group_by(ccf_region_name1) %>%
  summarise(width = max(n_section), .groups = "drop") %>%
  mutate(offset = cumsum(lag(width, default = 0)))

## -------------------------------------------------
## 6. 全局 x
## -------------------------------------------------
df_long <- df_long %>%
  left_join(region_width, by = "ccf_region_name1") %>%
  mutate(x_global = offset + section_index)

## -------------------------------------------------
## 7. 背景数据
## -------------------------------------------------
bg_df <- region_width %>%
  mutate(
    xmin = offset + 0.5,
    xmax = offset + width + 0.5,
    xmid = offset + width / 2 + 0.5
  )

## -------------------------------------------------
## 8. 无限 CCF 调色（关键修改点）
## -------------------------------------------------
n_region <- length(region_order)

ccf_colors <- setNames(
  qualitative_hcl(
    n = n_region,
    palette = "Dynamic",  # 区分度最高
    l = 90,               # 亮度（越大越浅）
    c = 35                # 饱和度（控制柔和程度）
  ),
  region_order
)

## -------------------------------------------------
## 9. 绘图
## -------------------------------------------------
p <- ggplot(df_long, aes(x = x_global, y = proportion)) +
  
  geom_rect(
    data = bg_df,
    aes(
      xmin = xmin,
      xmax = xmax,
      ymin = -Inf,
      ymax = Inf,
      fill = ccf_region_name1
    ),
    inherit.aes = FALSE,
    alpha = 0.35
  ) +
  
  geom_vline(
    xintercept = bg_df$xmax,
    linetype = "dashed",
    linewidth = 0.3,
    color = "grey50"
  ) +
  
  geom_line(
    aes(
      color = type,
      group = interaction(layer, ccf_region_name1, type)
    ),
    linewidth = 1
  ) +
  
  facet_wrap(~ layer, ncol = 1, scales = "free_y") +
  
  scale_x_continuous(
    breaks = bg_df$xmid,
    labels = bg_df$ccf_region_name1,
    expand = c(0, 0)
  ) +
  
  scale_color_manual(
    values = c(Glut = "#2b6cb0", GABA = "#c53030"),
    name = NULL
  ) +
  
  scale_fill_manual(
    values = ccf_colors,
    guide = "none"
  ) +
  
  labs(
    x = "CCF region",
    y = "Proportion"
  ) +
  
  theme_classic(base_size = 11) +
  theme(
    legend.position = "top",
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold"),
    panel.spacing.y = unit(0.8, "lines")
  )

print(p)

