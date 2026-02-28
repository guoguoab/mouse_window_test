suppressPackageStartupMessages({
  library(dplyr)
  library(stringr)
})

# ===== 参数 =====
input_dir  <- "/picb/neurosys/fuwen1/mouse_cortex/cellwindowcortexall1"
output_dir <- "/picb/neurosys/fuwen1/mouse_cortex/figure1"

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# 找到所有 rdata 文件
rdata_files <- list.files(
  input_dir,
  pattern = "\\.rdata$",
  full.names = TRUE,
  ignore.case = TRUE
)

# ===== 主循环 =====
for (f in rdata_files) {
  
  # 清空环境，避免对象名冲突
  rm(list = ls())
  
  load(f)  # 假设只加载一个 data.frame
  
  # 自动找到 data.frame
  obj_name <- ls()[sapply(ls(), function(x) is.data.frame(get(x)))]
  if (length(obj_name) != 1) {
    warning(basename(f), " 中未检测到唯一 data.frame，跳过")
    next
  }
  
  df <- get(obj_name)
  
  # ===== 1. 提取脑区 =====
  df <- df %>%
    mutate(
      ccf_region_name1 = str_trim(str_split_fixed(ccf_region_name, ",", 2)[,1])
    )
  
  # ===== 2. 计算兴奋/抑制比例 =====
  res <- df %>%
    filter(cell_Neuron_type %in% c("Glut", "Gaba")) %>%
    group_by(
      brain_section_label,
      ccf_region_name1,
      layer
    ) %>%
    summarise(
      rate = sum(cell_Neuron_type == "Glut") /
        sum(cell_Neuron_type %in% c("Glut", "Gaba")),
      .groups = "drop"
    )
  
  # ===== 3. 输出 =====
  out_file <- file.path(
    output_dir,
    paste0(tools::file_path_sans_ext(basename(f)), ".txt")
  )
  
  write.table(
    res,
    file = out_file,
    sep = "\t",
    row.names = FALSE,
    quote = FALSE
  )
}
