library(dplyr)
library(readr)
library(stringr)

## ===============================
## 1. 基础路径（默认使用当前项目下 data1）
## ===============================
args <- commandArgs(trailingOnly = TRUE)
base_dir <- ifelse(length(args) >= 1, args[1], "data1")

if (!dir.exists(base_dir)) {
  stop("base_dir 不存在: ", base_dir)
}

## ===============================
## 2. 遍历所有 sample 目录（通用，不限制 sample 个数）
## ===============================
sample_dirs <- list.dirs(
  base_dir,
  recursive = FALSE,
  full.names = TRUE
)

sample_dirs <- sample_dirs[grepl("^sample_", basename(sample_dirs))]

if (length(sample_dirs) == 0) {
  stop("在 ", base_dir, " 下没有找到 sample_* 目录")
}

## ===============================
## 3. 读取所有 *cell_id.txt 文件
## ===============================
df_list <- list()

for (sample_path in sample_dirs) {
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

    if (is.null(tmp) || nrow(tmp) == 0) next

    needed_cols <- c(
      "slide", "sample", "subclass", "layer", "region",
      "enrich_subclass_cell_ids_num",
      "Glut_Neruon_cell_ids_num",
      "GABA_Neruon_cell_ids_num"
    )

    if (!all(needed_cols %in% colnames(tmp))) next

    ## 只保留后续统计必需列，并统一类型，避免 bind_rows 的类型冲突
    tmp <- tmp %>%
      transmute(
        slide = as.character(slide),
        sample = as.character(sample),
        subclass = as.character(subclass),
        layer = as.character(layer),
        region = as.character(region),
        enrich_subclass_cell_ids_num = as.numeric(enrich_subclass_cell_ids_num),
        Glut_Neruon_cell_ids_num = as.numeric(Glut_Neruon_cell_ids_num),
        GABA_Neruon_cell_ids_num = as.numeric(GABA_Neruon_cell_ids_num)
      )

    df_list[[length(df_list) + 1]] <- tmp
  }
}

if (length(df_list) == 0) {
  stop("没有读取到有效的 *cell_id.txt 数据")
}

## ===============================
## 4. 合并成总表
## ===============================
df_all <- bind_rows(df_list)

## ===============================
## 5. sample 级别汇总
## ===============================
df_sample_level <- df_all %>%
  group_by(slide, sample, subclass, layer) %>%
  summarise(
    region_num = n_distinct(region),
    enrich_num = mean(enrich_subclass_cell_ids_num, na.rm = TRUE),
    glut_num = mean(Glut_Neruon_cell_ids_num, na.rm = TRUE),
    gaba_num = mean(GABA_Neruon_cell_ids_num, na.rm = TRUE),
    glut_gaba_ratio = ifelse(gaba_num == 0, NA_real_, glut_num / gaba_num),
    .groups = "drop"
  )

## ===============================
## 6. 所有 sample 的平均值（最终结果）
## ===============================
df_final <- df_sample_level %>%
  group_by(slide, subclass, layer) %>%
  summarise(
    mean_region_num = mean(region_num, na.rm = TRUE),
    mean_enrich_num = mean(enrich_num, na.rm = TRUE),
    mean_glut_num = mean(glut_num, na.rm = TRUE),
    mean_gaba_num = mean(gaba_num, na.rm = TRUE),
    mean_glut_gaba_ratio = mean(glut_gaba_ratio, na.rm = TRUE),
    sample_n = n(),
    .groups = "drop"
  )

## ===============================
## 7. 保存结果
## ===============================
write.table(
  df_sample_level,
  file = "sample_level_summary.txt",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  df_final,
  file = "all_sample_mean_summary.txt",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

cat("完成：\n")
cat("- sample 级别结果: sample_level_summary.txt\n")
cat("- 所有 sample 平均结果: all_sample_mean_summary.txt\n")
