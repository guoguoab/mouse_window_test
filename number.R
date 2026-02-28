library(dplyr)
library(readr)
library(stringr)

## ===============================
## 1. 基础路径
## ===============================
base_dir <- "/picb/neurosys/chenrenrui/Mouse_ST/A molecularly defined and spatially resolved cell atlas of the whole mouse brain/subclass_10K/neocortex_sample_subclass/merge_region"

## ===============================
## 2. 遍历所有 slide (C57BL6J-*)
## ===============================
slide_dirs <- list.dirs(
  base_dir,
  recursive = FALSE,
  full.names = TRUE
)

slide_dirs <- slide_dirs[file.info(slide_dirs)$isdir]

## ===============================
## 3. 读取所有 cell_id.txt
## ===============================
df_list <- list()

for (slide_path in slide_dirs) {
  
  slide_name <- basename(slide_path)
  
  sample_dirs <- list.dirs(
    slide_path,
    recursive = FALSE,
    full.names = TRUE
  )
  
  sample_dirs <- sample_dirs[grepl("^sample_", basename(sample_dirs))]
  
  for (sample_path in sample_dirs) {
    
    sample_name <- basename(sample_path)
    
    txt_files <- list.files(
      sample_path,
      pattern = "cell_id\\.txt$",
      full.names = TRUE
    )
    
    if (length(txt_files) == 0) next
    
    for (f in txt_files) {
      
      tmp <- tryCatch(
        read.delim(f, stringsAsFactors = FALSE),
        error = function(e) NULL
      )
      
      if (is.null(tmp)) next
      
      ## ===============================
      ## ★ 核心：统一 enrich_subclass_cell_ids 为 character
      ## ===============================
      if ("enrich_subclass_cell_ids" %in% colnames(tmp)) {
        tmp$enrich_subclass_cell_ids <- as.character(tmp$enrich_subclass_cell_ids)
      }
      
      tmp$slide  <- slide_name
      tmp$sample <- sample_name
      
      df_list[[length(df_list) + 1]] <- tmp
      
    }
  }
}

## ===============================
## 4. 合并成总表
## ===============================
df_all <- bind_rows(df_list)

## ===============================
## 5. 类型转换（非常关键）
## ===============================
df_all <- df_all %>%
  mutate(
    enrich_subclass_cell_ids_num = as.numeric(enrich_subclass_cell_ids_num),
    Glut_Neruon_cell_ids_num     = as.numeric(Glut_Neruon_cell_ids_num),
    GABA_Neruon_cell_ids_num     = as.numeric(GABA_Neruon_cell_ids_num)
  )

## ===============================
## 6. sample 级别汇总
## ===============================
df_sample_level <- df_all %>%
  group_by(
    slide,
    sample,
    subclass,
    layer
  ) %>%
  summarise(
    region_num = n_distinct(region),
    enrich_num = mean(enrich_subclass_cell_ids_num, na.rm = TRUE),
    glut_num   = mean(Glut_Neruon_cell_ids_num, na.rm = TRUE),
    gaba_num   = mean(GABA_Neruon_cell_ids_num, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    glut_gaba_ratio = glut_num / gaba_num
  )

## ===============================
## 7. 500 sample 取平均（最终结果）
## ===============================
df_final <- df_sample_level %>%
  group_by(
    slide,
    subclass,
    layer
  ) %>%
  summarise(
    mean_region_num         = mean(region_num, na.rm = TRUE),
    mean_enrich_num         = mean(enrich_num, na.rm = TRUE),
    mean_glut_gaba_ratio    = mean(glut_gaba_ratio, na.rm = TRUE),
    mean_enrich_num_again   = mean(enrich_num, na.rm = TRUE), # 按你要求保留
    sample_n                = n(),
    .groups = "drop"
  )

## ===============================
## 8. 保存结果
## ===============================
write.table(
  df_final,
  file = "slide_subclass_layer_summary.txt",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
