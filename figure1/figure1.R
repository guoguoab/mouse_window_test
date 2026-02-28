suppressPackageStartupMessages({
  library(dplyr)
  library(stringr)
})

# ===== 参数 =====
input_dir  <- "/picb/neurosys/fuwen1/mouse_cortex/cellwindowcortexall1"
output_dir <- "/picb/neurosys/fuwen1/mouse_cortex/figure1_1"

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
  

  load(f)          # 假设加载出的对象是 file_tmp
  
  df <- file_tmp
  
  # ===== 1. 提取脑区 =====
  df <- df %>%
    mutate(
      ccf_region_name1 = str_trim(
        str_split_fixed(ccf_region_name, ",", 2)[, 1]
      )
    )
  
  # ===== 2. 计算 Glut / GABA 占比 =====
  res <- df %>%
    filter(cell_Neuron_type %in% c("Glut", "GABA")) %>%
    group_by(
      brain_section_label,
      ccf_region_name1,
      layer
    ) %>%
    summarise(
      n_glut = sum(cell_Neuron_type == "Glut"),
      n_gaba = sum(cell_Neuron_type == "GABA"),
      n_total = n_glut + n_gaba,
      proportion_glut = n_glut / n_total,
      proportion_gaba = n_gaba / n_total,
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

