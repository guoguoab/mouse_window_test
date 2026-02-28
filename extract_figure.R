suppressPackageStartupMessages({
  library(dplyr)
  library(stringr)
  library(ggplot2)
})

# ===============================
# 1. 读取 region 汇总文件
# ===============================
region_df <- read.delim(
  "D:/GUOGUO_project/gsmap/R/articlefigure/mouse_class_cluster_total.txt",
  sep = "\t",
  stringsAsFactors = FALSE,
  check.names = FALSE
)
x1<-region_df[region_df$slide=='C57BL6J-1.079',]
# ===============================
# 2. 指定 slide / layer / region × class
# ===============================
slide_use <- "C57BL6J-2.031"
layer_use <- "layer 4"

merge_region_1 <- "regions_1"
class_1        <- "035 OB Eomes Ms4a15 Glut"

merge_region_2 <- "regions_1"
class_2        <- "039 OB Meis2 Thsd7b Gaba"

# ===============================
# 3. 精确匹配 region × class（允许多行）
# ===============================
df_sub <- region_df %>%
  filter(slide == slide_use, layer == layer_use) %>%
  filter((merge_regions == merge_region_1 & class == class_1) |
           (merge_regions == merge_region_2 & class == class_2))

extract_ids <- function(x) {
  x <- x[!is.na(x) & x != ""]
  if(length(x) == 0) return(character(0))
  unlist(lapply(x, function(xx) {
    formatC(as.numeric(str_split(xx, ",")[[1]]), digits = 10, format = "e")
  }))
}

get_region_ids <- function(df_rows) {
  c(
    extract_ids(df_rows$Glut_Neruon_cell_ids),
    extract_ids(df_rows$GABA_Neruon_cell_ids),
    extract_ids(df_rows$Non_Neruon_cell_ids)
  )
}

region_cell_id_1 <- get_region_ids(df_sub %>% filter(merge_regions == merge_region_1 & class == class_1))
region_cell_id_2 <- get_region_ids(df_sub %>% filter(merge_regions == merge_region_2 & class == class_2))

# ===============================
# 4. 读取 slide 对应的细胞坐标文件
# ===============================
cell_df <- read.delim(
  file.path("D:/GUOGUO_project/gsmap/R/articlefigure/data",
            paste0(slide_use, ".txt")),
  sep = "\t",
  stringsAsFactors = FALSE
) %>%
  mutate(cell_label_fmt = formatC(as.numeric(cell_label), digits = 10, format = "e"))

# ===============================
# 5. 图 1：灰色背景 + region1(红) / region2(蓝) + 手动图注
# ===============================
p1 <- ggplot(cell_df, aes(x = x, y = y)) +
  geom_point(color = "grey85", size = 0.5, alpha = 1) +  # 背景细胞
  
  geom_point(
    data = cell_df %>% filter(cell_label_fmt %in% region_cell_id_1),
    color = "red", size = 1
  ) +
  geom_point(
    data = cell_df %>% filter(cell_label_fmt %in% region_cell_id_2),
    color = "blue", size = 1
  ) +
  
  # 手动注释（论文风格）
  annotate(
    "text",
    x = Inf, y = Inf,
    label = paste0("Red: ", merge_region_1, "\n", class_1),
    hjust = 1.05, vjust = 2,
    color = "red",
    size = 3.5, fontface = "bold"
  ) +
  annotate(
    "text",
    x = Inf, y = Inf,
    label = paste0("Blue: ", merge_region_2, "\n", class_2),
    hjust = 1.05, vjust = 4,
    color = "blue",
    size = 3.5, fontface = "bold"
  ) +
  
  coord_equal() +
  theme_void() +
  ggtitle(paste(slide_use, layer_use, "Region Comparison")) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    plot.margin = margin(10, 120, 10, 10)  # 给右侧图注留空间
  )

print(p1)

# ===============================
# 6. 图 2：灰色背景 + subclass 着色（论文风格）
# ===============================
highlight_df <- cell_df %>%
  filter(cell_label_fmt %in% c(region_cell_id_1, region_cell_id_2))

p2 <- ggplot(cell_df, aes(x = x, y = y)) +
  geom_point(color = "grey85", size = 0.5, alpha = 1) +  # 背景细胞
  
  geom_point(
    data = highlight_df,
    aes(color = subclass),
    size = 1
  ) +
  
  coord_equal() +
  theme_void() +
  theme(
    legend.title = element_blank(),
    legend.position = "right",
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14)
  ) +
  ggtitle(paste(slide_use, layer_use, "Subclass Colored"))

print(p2)

