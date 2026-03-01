#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
})

args <- commandArgs(trailingOnly = TRUE)
true_dir <- if (length(args) >= 1) args[[1]] else "true_data"
data1_dir <- if (length(args) >= 2) args[[2]] else "data1"
out_file <- if (length(args) >= 3) args[[3]] else "true_loc_presence_by_sample.tsv"
tol_digits <- if (length(args) >= 4) as.integer(args[[4]]) else 8L

if (is.na(tol_digits) || tol_digits < 0) {
  stop("tol_digits 必须是非负整数")
}

normalize_num <- function(x, digits = 8L) {
  format(round(as.numeric(x), digits = digits), scientific = FALSE, trim = TRUE, nsmall = 0)
}

coord_loc_key <- function(xstart, ystart, digits = 8L) {
  paste0(normalize_num(xstart, digits), "_", normalize_num(ystart, digits))
}

extract_chip_from_true_filename <- function(path) {
  bn <- basename(path)
  sub("^((C57BL6J-[^_]+)).*$", "\\1", bn)
}

extract_chip_from_sample_filename <- function(path) {
  bn <- basename(path)
  # sample_1_C57BL6J-1.025_layer 1_xxx_merged_rgions_loc.txt
  sub("^sample_[^_]+_((C57BL6J-[^_]+)).*$", "\\1", bn)
}

find_loc_files <- function(root_dir) {
  # 兼容不同命名：*_loc.txt / *_rgions_loc.txt / *_regions_loc.txt，以及大小写 .TXT
  all_files <- list.files(root_dir, recursive = TRUE, full.names = TRUE)
  if (length(all_files) == 0) return(character(0))

  file_names <- basename(all_files)
  hit <- grepl(
    "(_loc\\.txt$)|(_rgions_loc\\.txt$)|(_regions_loc\\.txt$)",
    file_names,
    ignore.case = TRUE
  )
  unique(all_files[hit])
}

# 构建 data1 每个 sample 对每个 chip 的 loc 集合
build_sample_chip_loc_sets <- function(sample_dir, digits = 8L) {
  files <- find_loc_files(sample_dir)
  out <- list()

  if (length(files) == 0) {
    return(out)
  }

  for (f in files) {
    dt <- tryCatch(
      fread(f, sep = "\t", header = TRUE, fill = TRUE, showProgress = FALSE),
      error = function(e) NULL
    )

    if (is.null(dt) || nrow(dt) == 0) next

    chip <- NA_character_
    if ("chip" %in% names(dt)) {
      chip_vals <- unique(dt$chip[!is.na(dt$chip) & dt$chip != ""])
      if (length(chip_vals) > 0) chip <- chip_vals[[1]]
    }
    if (is.na(chip) || chip == "") {
      chip <- extract_chip_from_sample_filename(f)
    }

    if (is.na(chip) || chip == "" || identical(chip, basename(f))) next

    # 只按 loc 匹配。优先使用 loc / loc_1；若缺失则降级用 xstart_ystart 生成 loc key。
    loc_keys <- character(0)
    if ("loc" %in% names(dt)) {
      loc_keys <- as.character(dt$loc)
    } else if ("loc_1" %in% names(dt)) {
      loc_keys <- as.character(dt$loc_1)
      # sample 文件的 loc_1 可能是 layer 1_xstart_ystart，只取坐标部分
      loc_keys <- sub("^layer[^_]*_", "", loc_keys)
    } else if (all(c("xstart", "ystart") %in% names(dt))) {
      loc_keys <- coord_loc_key(dt$xstart, dt$ystart, digits = digits)
    }

    loc_keys <- loc_keys[!is.na(loc_keys) & loc_keys != ""]
    if (length(loc_keys) == 0) next

    if (is.null(out[[chip]])) {
      out[[chip]] <- unique(loc_keys)
    } else {
      out[[chip]] <- unique(c(out[[chip]], loc_keys))
    }
  }

  out
}

true_loc_files <- find_loc_files(true_dir)
if (length(true_loc_files) == 0) {
  stop(sprintf("在 %s 中没有找到 loc 文件（支持 *_loc.txt / *_rgions_loc.txt / *_regions_loc.txt）", true_dir))
}

sample_dirs <- list.dirs(data1_dir, recursive = FALSE, full.names = TRUE)
if (length(sample_dirs) == 0) {
  stop(sprintf("在 %s 中没有找到 sample 目录", data1_dir))
}

message(sprintf("发现 true_data loc 文件: %d", length(true_loc_files)))
message(sprintf("发现 sample 目录: %d", length(sample_dirs)))

sample_chip_locs <- vector("list", length(sample_dirs))
names(sample_chip_locs) <- basename(sample_dirs)
for (i in seq_along(sample_dirs)) {
  sname <- basename(sample_dirs[[i]])
  sample_chip_locs[[sname]] <- build_sample_chip_loc_sets(sample_dirs[[i]], digits = tol_digits)

  chip_count <- length(sample_chip_locs[[sname]])
  loc_count <- sum(vapply(sample_chip_locs[[sname]], length, integer(1)))
  message(sprintf("[%s] chip数: %d, loc总数: %d", sname, chip_count, loc_count))
}

results <- list()
idx <- 1L

for (f in true_loc_files) {
  dt <- tryCatch(
    fread(f, sep = "\t", header = TRUE, fill = TRUE, showProgress = FALSE),
    error = function(e) NULL
  )

  if (is.null(dt) || nrow(dt) == 0) next

  chip <- extract_chip_from_true_filename(f)

  # true 文件优先使用 loc；若缺失则用 xstart_ystart 生成 loc key
  if ("loc" %in% names(dt)) {
    dt[, loc_key := as.character(loc)]
  } else if (all(c("xstart", "ystart") %in% names(dt))) {
    dt[, loc_key := coord_loc_key(xstart, ystart, digits = tol_digits)]
  } else {
    warning(sprintf("跳过文件(缺少 loc 与坐标列): %s", f))
    next
  }

  dt <- dt[!is.na(loc_key) & loc_key != ""]
  if (nrow(dt) == 0) next

  # 用户要求：按每个 loc 检索；不考虑内部子 loc 组成
  uniq_locs <- unique(dt$loc_key)

  for (loc_name in uniq_locs) {
    row <- data.table(
      true_file = basename(f),
      true_file_path = f,
      chip = chip,
      loc = loc_name
    )

    for (sname in names(sample_chip_locs)) {
      chip_bucket <- sample_chip_locs[[sname]][[chip]]
      found <- !is.null(chip_bucket) && (loc_name %in% chip_bucket)
      row[[sname]] <- as.integer(found)
    }

    results[[idx]] <- row
    idx <- idx + 1L
  }
}

if (length(results) == 0) {
  stop("没有产生任何结果，请检查输入文件格式")
}

out_dt <- rbindlist(results, fill = TRUE)
setcolorder(
  out_dt,
  c("true_file", "true_file_path", "chip", "loc", sort(setdiff(names(out_dt), c("true_file", "true_file_path", "chip", "loc"))))
)

fwrite(out_dt, out_file, sep = "\t")
message(sprintf("完成: %s", out_file))
message(sprintf("记录数: %d", nrow(out_dt)))
