suppressPackageStartupMessages({
  library(dplyr)
  library(stringr)
  library(tidyr)
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
region_dft <- read.delim(
  "D:/GUOGUO_project/gsmap/R/analyse/Table_2_total_cell.txt",
  sep = "\t",
  stringsAsFactors = FALSE,
  check.names = FALSE
)
glutlist<-unique(region_df$class)
glutlist<-glutlist[grepl("^(00[1-9]|01[0-9]|02[0-6])", glutlist)]
gabalist <- c(
  "046 Vip Gaba",
  "047 Sncg Gaba",
  "048 RHP-COA Ndnf Gaba",
  "049 Lamp5 Gaba",
  "050 Lamp5 Lhx6 Gaba",
  "051 Pvalb chandelier Gaba",
  "052 Pvalb Gaba",
  "053 Sst Gaba"
)
target_class <- c(glutlist, gabalist)

region_df <- region_df[
  region_df$class %in% target_class,
]


# ===============================
# 2. 函数：提取 cell_ids 并统一精度
# ===============================
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

get_all_ids <- function(row) {
  c(
    extract_ids(row$Glut_Neruon_cell_ids),
    extract_ids(row$GABA_Neruon_cell_ids),
    extract_ids(row$Non_Neruon_cell_ids)
  )
}

# ===============================
# 3. 初始化结果列表
# ===============================
result_list <- list()

# ===============================
# 4. 对每个 slide × layer 做两两比较
# ===============================
for (slide_i in unique(region_df$slide)) {
  
  df_slide <- region_df %>% filter(slide == slide_i)
  
  for (layer_i in unique(df_slide$layer)) {
    
    df_layer <- df_slide %>% filter(layer == layer_i)
    if (nrow(df_layer) < 2) next
    
    combn_indices <- combn(seq_len(nrow(df_layer)), 2)
    
    for (j in seq_len(ncol(combn_indices))) {
      
      idx_a <- combn_indices[1, j]
      idx_b <- combn_indices[2, j]
      
      row_a <- df_layer[idx_a, ]
      row_b <- df_layer[idx_b, ]
      
      ids_a <- get_all_ids(row_a)
      ids_b <- get_all_ids(row_b)
      
      n_a <- length(ids_a)
      n_b <- length(ids_b)
      
      # 限制条件：b 的 cell_id 数 ≥ a
      
      if (n_a == 0 || n_b == 0) next
      
      overlap_a_in_b <- sum(ids_a %in% ids_b) / n_a
      overlap_b_in_a <- sum(ids_b %in% ids_a) / n_b
      
      if (overlap_a_in_b >= 0|| overlap_b_in_a >= 0) {
        result_list[[length(result_list) + 1]] <- data.frame(
          ## slide / layer（整体）
          slide = slide_i,
          layer = layer_i,
          
          ## a / b 各自的 layer
          a_layer = row_a$layer,
          b_layer = row_b$layer,
          
          ## class
          a_class = row_a$class,
          b_class = row_b$class,
          a_glutotal=row_a$Glut_Neruon_cell_ids_num,
          a_gabatotal=row_a$GABA_Neruon_cell_ids_num,
          b_glutotal=row_b$Glut_Neruon_cell_ids_num,
          b_gabatotal=row_b$GABA_Neruon_cell_ids_num,
          ## region（原始 region）
          a_region = row_a$region,
          b_region = row_b$region,
          a_ei=row_a$Glut_Neruon_cell_ids_num/row_a$GABA_Neruon_cell_ids_num,
          b_ei=row_b$Glut_Neruon_cell_ids_num/row_b$GABA_Neruon_cell_ids_num,
          ## merge_regions（如果你还想保留）
          a_merge_region = row_a$merge_regions,
          b_merge_region = row_b$merge_regions,
          
          ## neuron type
          a_cell_Neuron_type = row_a$cell_Neuron_type,
          b_cell_Neuron_type = row_b$cell_Neuron_type,
          
          ## enrich 数量
          a_enrich_class_cell_ids_num = row_a$enrich_class_cell_ids_num,
          b_enrich_class_cell_ids_num = row_b$enrich_class_cell_ids_num,
          
          ## cell_id 数量
          a_cell_id_n = n_a,
          b_cell_id_n = n_b,
          
          ## overlap
          overlap_a_in_b = overlap_a_in_b,
          overlap_b_in_a = overlap_b_in_a,
          
          stringsAsFactors = FALSE
        )
        
      }
    }
  }
}

# ===============================
# 5. 汇总结果
# ===============================
result_df <- bind_rows(result_list)



# ===============================
# 6. 输出 CSV
# ===============================
write.csv(
  result_df,
  "D:/GUOGUO_project/gsmap/R/articlefigure/cellid_overlap_summary_filtered.csv",
  row.names = FALSE
)

message("Done! Found ", nrow(result_df), " overlapping pairs.")
