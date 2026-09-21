# enriched cell 数替换为 Z 值的统计逻辑

`compare_subclass_enriched_cell_zscore.py` 延续原脚本“真实值与 500 次随机
enrich 结果比较”的统计单位，但把纵轴统一成置换零分布的 Z 值。

对每个 subclass 分别执行：

1. 第 `b` 次随机置换内可能有多个 cluster，先计算该次置换的平均 enriched
   cell 数 `R_b`。因此 500 次置换产生 500 个彼此独立的零分布统计量，而不是
   把所有随机 cluster 混在一起造成伪重复。
2. 计算零分布均值 `mu_R = mean(R_b)` 和样本标准差
   `s_R = sd(R_b, ddof=1)`。
3. 真实数据的统计量是所有真实 cluster 的平均 enriched cell 数 `T`，转换为
   `Z_true = (T - mu_R) / s_R`。
4. 为了画出与真实 Z 可直接比较的随机箱线图，每次置换也用**同一组**参数转换：
   `Z_b = (R_b - mu_R) / s_R`。不能给每次置换单独计算均值和标准差，否则其 Z
   值会失去共同尺度。
5. P 值仍使用带 `+1` 校正的经验置换检验。例如富集方向为
   `P = (1 + count(Z_b >= Z_true)) / (B + 1)`；脚本同时输出减少方向和双侧 P。

## 与截图中 window Z / cluster mass 的区别

截图的方法需要**每个空间窗口的单侧 raw P 值和窗口邻接关系**：先计算
`z_i = Phi^-1(1 - p_i)`，再把相邻且通过 cluster-forming threshold 的窗口合并，
并求 `cluster_mass = sum(z_i)`。每次置换必须重新完成窗口检验、阈值化和聚类，
最后保存该次置换的最大 cluster mass。

当前 `*_merged_regions_table_cell_id.txt` 只有已经合并好的 cluster 及 enriched
cell 数，没有 window raw P、window 标识或邻接关系，所以不能从这些文件严谨地
反推出截图中的 cluster mass。本脚本实现的是在现有数据可支持范围内，与原统计
最接近的 **subclass enriched-cell Z 标准化置换比较**，不把它误称为 cluster
mass inference。

## 运行

```bash
python compare_subclass_enriched_cell_zscore.py \
  --true-dir true_data \
  --random-dir data1 \
  --workers 20
```

会输出 summary TSV、每次置换的 detail TSV，以及按 20 个 subclass 拆分的 SVG。
若某个 subclass 少于两个随机重复，或随机均值标准差为 0，Z 值无定义，脚本会
明确报告并跳过该 subclass。

## 更保守的“随机最大 Z”版本

原脚本不是只把真实 Z 与“随机 Z 的平均值”做一次二元比较：随机 Z 的均值按定义
接近 0，图中保留完整随机 Z 箱线图，并用全部随机 Z 计算经验 P 值。如果需要更
严格的直接门槛，可以运行：

```bash
python compare_subclass_enriched_cell_zscore_max.py \
  --true-dir true_data \
  --random-dir data1 \
  --workers 20
```

新脚本对每个 subclass 计算 `max(Z_b)`，直接判断
`Z_true > max(Z_b)`，并输出真实 Z、随机最大 Z、产生最大值的 sample、两者差值和
是否超过最大值。这个条件等价于没有任何一次已观测随机置换达到真实值；它比与
随机中心比较更保守，但并不产生比 `1 / (B + 1)` 更精细的 P 值。
