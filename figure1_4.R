## ===============================
## 1. 设置 RData 文件目录
## ===============================
rdata_dir <- "/picb/neurosys/fuwen1/mouse_cortex/cellwindowcortexall1"   # 改成你的路径

files <- list.files(
  rdata_dir,
  pattern = "\\.RData$",
  full.names = TRUE
)

## ===============================
## 2. 初始化计数器
## ===============================
n_glut <- 0
n_gaba <- 0

## ===============================
## 3. 批量读取并统计
## ===============================
for (f in files) {
  
  ## 清空环境，避免变量污染

  load(f)
  
  if (!exists("file_tmp")) {
    warning(paste("file_tmp not found in", f))
    next
  }
  
  ## 统计
  n_glut <- n_glut + sum(file_tmp$cell_Neuron_type == "Glut", na.rm = TRUE)
  n_gaba <- n_gaba + sum(file_tmp$cell_Neuron_type == "GABA", na.rm = TRUE)
}

## ===============================
## 4. 输出结果
## ===============================
cat("Glut cell count :", n_glut, "\n")
cat("GABA cell count :", n_gaba, "\n")
cat("Total (Glut + GABA) :", n_glut + n_gaba, "\n")
