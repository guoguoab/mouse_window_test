df_allfiltered <- df_all %>%
  filter(
    layer%in%c("layer 2/3","layer 5","layer 6a","layer 6b"),
    a_enrich_class_cell_ids_num>=3&b_enrich_class_cell_ids_num>=3|a_enrich_class_cell_ids_num/(a_glutotal+a_gabatotal)>=0.8&b_enrich_class_cell_ids_num/(b_glutotal+b_gabatotal)>=0.8
  ) %>%
  group_by(slide, layer) %>%
  slice_max(
    order_by = overlap_a_in_b,
    n = 1,
    with_ties = FALSE
  ) %>%
  ungroup()
plot_param_df1<-plot_param_df[18,]

df <- read.csv("D:/GUOGUO_project/gsmap/R/analysefig/cell_metadata1.csv")

df<-df[df$brain_section_label=='Zhuang-ABCA-1.057',]
plot_param_dfx <- region_df1 %>%
  filter(
    a_cell_Neuron_type == "Gaba",
    b_cell_Neuron_type == "Glut")

plot_param_dfx1 <- region_df1 %>%
  filter(
    a_cell_Neuron_type == "Glut",
    b_cell_Neuron_type == "Gaba")
plot_param_dfxt <- region_dft %>%
  filter(
   cluster.1_cell_Neuron_type == "Gaba",
    cluster.2_cell_Neuron_type == "Glut")

plot_param_dfx1t <- region_dft %>%
  filter(
    cluster.1_cell_Neuron_type== "Glut",
    cluster.2_cell_Neuron_type == "Gaba")

plot_param_dfx11 <- plot_param_dfx1 %>%
  filter(
   overlap_b_in_a==1
  )

plot_param_dfxxx <- plot_param_dfx %>%
  filter(
    overlap_a_in_b==1
  )
df_all <- bind_rows(plot_param_dfx11, plot_param_dfxxx)
df_all <- df_all %>%
  filter(
    overlap_a_in_b==1
  )
region_df1 <- result_df %>%
  filter(
    overlap_a_in_b>0&overlap_b_in_a>0
  )
plot_param_dfx1 <- plot_param_dfx1 %>%
  filter(
    overlap_a_in_b>0
  )
plot_param_dfx <- plot_param_dfx %>%
  filter(
    overlap_a_in_b>0
  )
plot_param_dfx12 <- region_df %>%
  filter(
    cell_Neuron_type=="Glut")
plot_param_dfx123 <- region_df %>%
  filter(
    cell_Neuron_type=="Gaba")
extract_ids <- function(x) {
  x <- x[!is.na(x) & x != ""]
  if (length(x) == 0) return(character(0))
  unlist(lapply(x, function(xx) {
    formatC(
      as.numeric(str_split(xx, ",")[[1]]),
      digits = 10,
      format = "e"
    )
  }))
}

get_region_ids <- function(df_rows) {
  c(
    extract_ids(df_rows$Glut_Neruon_cell_ids),
    extract_ids(df_rows$GABA_Neruon_cell_ids),
    extract_ids(df_rows$Non_Neruon_cell_ids)
  )
}
out_dir <- "D:/GUOGUO_project/gsmap/R/articlefigure/overlap_plots4"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
plot_param_df<-df_allfiltered
savex<-list()
for (i in seq_len(nrow(plot_param_df))) {
  
  message("Plotting ", i, "/", nrow(plot_param_df))
  
  slide_use <- plot_param_df$slide[i]
  layer_use <- plot_param_df$layer[i]
  
  merge_region_1 <- plot_param_df$a_merge_region[i]
  class_1        <- plot_param_df$a_class[i]   # GABA
  
  merge_region_2 <- plot_param_df$b_merge_region[i]
  class_2        <- plot_param_df$b_class[i]   # GLUT
  
  # ===============================
  # 1. region 汇总表中取对应行
  # ===============================
  df_sub <- region_df %>%
    filter(
      slide == slide_use,
      layer == layer_use
    ) %>%
    filter(
      (merge_regions == merge_region_1 & class == class_1) |
        (merge_regions == merge_region_2 & class == class_2)
    )
  
  if (nrow(df_sub) == 0) next
  
  region_cell_id_1 <- get_region_ids(
    df_sub %>% filter(merge_regions == merge_region_1 & class == class_1)
  )
  region_cell_id_2 <- get_region_ids(
    df_sub %>% filter(merge_regions == merge_region_2 & class == class_2)
  )
  
  # ===============================
  # 2. 读取 cell 坐标
  # ===============================
  cell_df <- read.delim(
    file.path(
      "D:/GUOGUo_project/gsmap/R/articlefigure/data",
      paste0(slide_use, ".txt")
    ),
    sep = "\t",
    stringsAsFactors = FALSE
  ) %>%
    mutate(
      cell_label_fmt = formatC(
        as.numeric(cell_label),
        digits = 10,
        format = "e"
      )
    )
  cell_df <- cell_df %>%
    mutate(
      layer = str_trim(
        sapply(
          str_split(ccf_region_name, ","),
          function(x) tail(x, 1)
        )
      )
    )
  
  highlight_df <- cell_df %>%
    filter(cell_label_fmt %in% c(region_cell_id_1, region_cell_id_2))
  x<-cell_df %>% filter(cell_label_fmt %in% region_cell_id_2|cell_label_fmt %in% region_cell_id_1)
  if(length(unique(x$subclass))<=11){
    savex <- bind_rows(savex, plot_param_df[i,])
  # ===============================
  # 3. 图 1：GABA vs GLUT（论文风格）
  # ===============================
  if(plot_param_df[i,]$a_cell_Neuron_type=="Glut"){
  p1 <- ggplot(cell_df, aes(x = x, y = y)) +
    geom_point(color = "grey85", size = 0.2) +
    geom_point(
      data = cell_df %>% filter(cell_label_fmt %in% region_cell_id_1),
      color = "#1F77B4", size = 0.2
    ) +
    geom_point(
      data = cell_df %>% filter(cell_label_fmt %in% region_cell_id_2),
      color = "#D62728", size = 0.2
    ) +scale_y_reverse()+
    annotate(
      "text",
      x = Inf, y = Inf,
      label = paste0("Red (Gaba):\n", class_2),
      hjust = 1.05, vjust = -8.7,
      color = "#D62728",
      size = 2, fontface = "bold"
    ) +
    annotate(
      "text",
      x = Inf, y = Inf,
      label = paste0("Blue (GlUt):\n", class_1),
      hjust = 1.05, vjust =-6,
      color = "#1F77B4",
      size =2, fontface = "bold"
    ) +
    coord_equal(clip = "off") +
    theme_void() +
    ggtitle(paste(slide_use, layer_use)) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      plot.margin = margin(10, 120, 10, 10)
    )
  }
  else{
    
    p1 <- ggplot(cell_df, aes(x = x, y = y)) +
      geom_point(color = "grey85", size = 0.2) +
      geom_point(
        data = cell_df %>% filter(cell_label_fmt %in% region_cell_id_2),
        color = "#1F77B4", size = 0.2
      ) +
      geom_point(
        data = cell_df %>% filter(cell_label_fmt %in% region_cell_id_1),
        color = "#D62728", size = 0.2
      ) +scale_y_reverse()+
      annotate(
        "text",
        x = Inf, y = Inf,
        label = paste0("Red (Gaba):\n", class_1),
        hjust = 1.05, vjust = -8.7,
        color = "#D62728",
        size = 2, fontface = "bold"
      ) +
      annotate(
        "text",
        x = Inf, y = Inf,
        label = paste0("Blue (GlUt):\n", class_2),
        hjust = 1.05, vjust =-6,
        color = "#1F77B4",
        size =2, fontface = "bold"
      ) +
      coord_equal(clip = "off") +
      theme_void() +
      ggtitle(paste(slide_use, layer_use)) +
      theme(
        plot.title = element_text(hjust = 0.5, face = "bold"),
        plot.margin = margin(10, 120, 10, 10)
      )
    
  }
  # ===============================
  # 4. 图 2：subclass 着色
  # ===============================
  p2 <- ggplot(cell_df, aes(x = x, y = y)) +
    geom_point(color = "grey85", size = 0.4) +
    geom_point(
      data = highlight_df,
      aes(color = subclass),
      size = 0.2
    ) + scale_y_reverse()+
    coord_equal() +
    theme_void() +
    theme(
      legend.title = element_blank(),
      plot.title = element_text(hjust = 0.5, face = "bold")
    ) +
    ggtitle(paste(slide_use, layer_use, "Subclass"))
  # ===============================
  # 6. 图 3：按 layer 着色（整张切片）
  # ===============================
  p3 <- ggplot(cell_df, aes(x = x, y = y)) +
    geom_point(
      aes(color = layer),
      size = 0.2
    ) +
    scale_y_reverse() +
    coord_equal() +
    theme_void() +
    ggtitle(paste(slide_use, "Layer")) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      legend.title = element_blank(),
      legend.position = "right",
      legend.text = element_text(size = 8)
    )
  
  # ===============================
  # 5. 保存 PNG
  # ===============================
  ggsave(
    filename = file.path(
      out_dir,
      paste0(slide_use, "_", gsub("/", "_", layer_use), gsub("/", "_", class_2),"_",gsub("/", "_", class_1),"_Gaba_vs_Glut.png")
    ),
    plot = p1,
    width = 6, height = 6, dpi = 300
  )
  
  ggsave(
    filename = file.path(
      out_dir,
      paste0(slide_use, "_", gsub("/", "_", layer_use),gsub("/", "_", class_2),"_",gsub("/", "_", class_1), "_Subclass.png")
    ),
    plot = p2,
    width = 6, height = 6, dpi = 300
    
  )
  ggsave(
    filename = file.path(
      out_dir,
      paste0(
        slide_use, "_",
        gsub("/", "_", layer_use),
        "_Layer.png"
      )
    ),
    plot = p3,
    width = 6, height = 6, dpi = 300
  )
  }
}
install.packages("openxlsx")   # 如果还没装
library(openxlsx)
write.xlsx(
  plot_param_df,
  file = "D:/GUOGUO_project/gsmap/R/articlefigure/result_df.xlsx",
  rowNames = FALSE
)
write.xlsx(
  x,
  file = "D:/GUOGUO_project/gsmap/R/articlefigure/resultall_df.xlsx",
  rowNames = FALSE
)
write.table(
  plot_param_df,
  file = "D:/GUOGUO_project/gsmap/R/articlefigure/result_df.txt",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = TRUE
)
