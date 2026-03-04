# true_vs_random 两个结果文件的快速解读

## 全局 totals（`true_vs_random_global_totals1.tsv`）
- `region_num`: true = 41016, random mean = 46682.264，fold-change = 0.8786（true 更低）。
- `glut_num`: true = 1488127.581366, random mean = 972570.889858，fold-change = 1.5301（true 更高）。
- `gaba_num`: true = 244806.208899, random mean = 156213.666545，fold-change = 1.5671（true 更高）。
- 三个指标的 `p_two_sided` 均为 0.00199601。

## permutation summary（`true_vs_random_permutation_summary1.tsv`）
样本行数：17815。

双侧 p 值显著比例（p < 0.05）：
- `region_num_p_two_sided`: 0.0593（n = 17815，median p = 0.6886）
- `glut_num_p_two_sided`: 0.2290（n = 17815，median p = 0.3233）
- `gaba_num_p_two_sided`: 0.2006（n = 17815，median p = 0.4004）
- `ei_ratio_p_two_sided`: 0.1713（n = 17166，median p = 0.2628）

## 可视化脚本
请运行：

```bash
Rscript visualize_true_vs_random1.R
```

输出文件（默认目录 `analysis_outputs/`）：
- `global_totals_comparison.png`
- `permutation_pvalue_histograms.png`
- `permutation_significant_rate.png`
- `analysis_summary.txt`
