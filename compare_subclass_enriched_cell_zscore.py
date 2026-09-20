#!/usr/bin/env python3
"""Compare true enriched-cell counts with a permutation-null Z distribution.

The independent null observations are the per-subclass means from each random
replicate (sample_1 ... sample_500), not all random clusters pooled together.
For subclass ``s`` the transformation is

    Z = (T - mean(R_1, ..., R_B)) / sample_sd(R_1, ..., R_B)

where T is the mean enriched-cell count in the true clusters and R_b is the
mean for random replicate b.  The same null mean and SD transform every R_b,
which makes the plotted random distribution directly comparable with true Z.
"""

from __future__ import annotations

import argparse
import math
import os
from pathlib import Path
from statistics import mean, median, stdev
from xml.sax.saxutils import escape

from compare_subclass_enriched_cell_size import (
    empirical_pvalues,
    load_random_parallel,
    load_true,
    quantile,
    significance,
    write_tsv,
)


def z_score(value: float, null_mean: float, null_sd: float) -> float:
    """Standardize *value* against a permutation-null mean and sample SD."""
    if not math.isfinite(null_sd) or null_sd <= 0:
        raise ValueError("随机重复的标准差必须大于 0，才能计算 Z 值")
    return (value - null_mean) / null_sd


def compute_results(
    true: dict[str, list[float]], random: dict[str, dict[str, list[float]]]
) -> tuple[list[dict[str, object]], list[dict[str, object]], list[str]]:
    """Return subclass summary, replicate detail, and skipped subclasses."""
    summary: list[dict[str, object]] = []
    detail: list[dict[str, object]] = []
    skipped: list[str] = []
    samples = sorted(random)

    for subclass in sorted(true):
        true_counts = true[subclass]
        true_mean = mean(true_counts)
        replicate_rows = [
            (sample, len(counts), mean(counts))
            for sample in samples
            if (counts := random[sample].get(subclass, []))
        ]
        # One value cannot estimate a null SD; a constant null has no defined Z.
        if len(replicate_rows) < 2:
            skipped.append(f"{subclass}: 少于 2 个随机重复")
            continue
        replicate_means = [row[2] for row in replicate_rows]
        null_mean = mean(replicate_means)
        null_sd = stdev(replicate_means)
        if null_sd == 0:
            skipped.append(f"{subclass}: 随机重复均值的标准差为 0")
            continue

        true_z = z_score(true_mean, null_mean, null_sd)
        random_z = [z_score(value, null_mean, null_sd) for value in replicate_means]
        p_gt, p_lt, p_two = empirical_pvalues(true_z, random_z)

        for (sample, cluster_n, replicate_mean), replicate_z in zip(
            replicate_rows, random_z
        ):
            detail.append({
                "subclass": subclass,
                "sample": sample,
                "random_cluster_n": cluster_n,
                "random_mean_enriched_cell_num": replicate_mean,
                "null_mean_enriched_cell_num": null_mean,
                "null_sd_enriched_cell_num": null_sd,
                "random_z": replicate_z,
                "true_mean_enriched_cell_num": true_mean,
                "true_z": true_z,
            })

        summary.append({
            "subclass": subclass,
            "true_cluster_n": len(true_counts),
            "true_mean_enriched_cell_num": true_mean,
            "true_median_enriched_cell_num": median(true_counts),
            "random_sample_n": len(random_z),
            "null_mean_enriched_cell_num": null_mean,
            "null_sd_enriched_cell_num": null_sd,
            "true_z": true_z,
            "random_z_mean": mean(random_z),
            "random_z_median": median(random_z),
            "random_z_q1": quantile(random_z, 0.25),
            "random_z_q3": quantile(random_z, 0.75),
            "random_z_min": min(random_z),
            "random_z_max": max(random_z),
            "p_one_sided_true_gt": p_gt,
            "p_one_sided_true_lt": p_lt,
            "p_two_sided": p_two,
        })
    return summary, detail, skipped


def draw_svg(
    rows: list[dict[str, object]], path: Path, subclasses_per_plot: int
) -> list[Path]:
    """Draw random Z boxplots and true-Z diamonds in manageable panels."""
    if not rows:
        raise RuntimeError("没有可计算 Z 值的 subclass，无法绘图")
    if subclasses_per_plot < 1:
        raise ValueError("每张图的 subclass 数必须至少为 1")
    rows = sorted(rows, key=lambda row: float(row["true_z"]), reverse=True)
    all_values = [0.0]
    for row in rows:
        all_values.extend((float(row["random_z_min"]), float(row["random_z_max"]), float(row["true_z"])))
    low, high = min(all_values), max(all_values)
    padding = max((high - low) * 0.1, 0.5)
    y_min, y_max = low - padding, high + padding
    panels = [rows[i:i + subclasses_per_plot] for i in range(0, len(rows), subclasses_per_plot)]
    written: list[Path] = []
    path.parent.mkdir(parents=True, exist_ok=True)

    for panel_number, panel in enumerate(panels, 1):
        width, height = max(1000, 105 + 70 * len(panel)), 700
        left, right, top, bottom = 85, 35, 70, 230
        chart_h, chart_w = height - top - bottom, width - left - right
        y = lambda value: top + chart_h * (y_max - float(value)) / (y_max - y_min)
        out = [
            f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}">',
            '<rect width="100%" height="100%" fill="white"/>',
            f'<text x="{width/2}" y="35" text-anchor="middle" font-size="24" font-family="sans-serif" font-weight="bold">True vs permutation-null enriched-cell Z by subclass (part {panel_number}/{len(panels)})</text>',
        ]
        for tick in range(6):
            value = y_min + (y_max - y_min) * tick / 5
            py = y(value)
            out += [
                f'<line x1="{left}" y1="{py:.1f}" x2="{width-right}" y2="{py:.1f}" stroke="#dddddd"/>',
                f'<text x="{left-8}" y="{py+5:.1f}" text-anchor="end" font-size="12" font-family="sans-serif">{value:.2f}</text>',
            ]
        if y_min <= 0 <= y_max:
            out.append(f'<line x1="{left}" y1="{y(0):.1f}" x2="{width-right}" y2="{y(0):.1f}" stroke="#777" stroke-dasharray="5,4"/>')
        step = chart_w / len(panel)
        for index, row in enumerate(panel):
            x = left + step * (index + 0.5)
            box_w = min(30, step * 0.55)
            low_z, q1, med, q3, high_z = (float(row[key]) for key in (
                "random_z_min", "random_z_q1", "random_z_median", "random_z_q3", "random_z_max"
            ))
            true_y = y(row["true_z"])
            out += [
                f'<line x1="{x:.1f}" y1="{y(low_z):.1f}" x2="{x:.1f}" y2="{y(high_z):.1f}" stroke="#d97820" stroke-width="2"/>',
                f'<rect x="{x-box_w/2:.1f}" y="{y(q3):.1f}" width="{box_w:.1f}" height="{max(1, y(q1)-y(q3)):.1f}" fill="#f5ad63" stroke="#b85c00"/>',
                f'<line x1="{x-box_w/2:.1f}" y1="{y(med):.1f}" x2="{x+box_w/2:.1f}" y2="{y(med):.1f}" stroke="#703500" stroke-width="2"/>',
                f'<polygon points="{x:.1f},{true_y-6:.1f} {x+6:.1f},{true_y:.1f} {x:.1f},{true_y+6:.1f} {x-6:.1f},{true_y:.1f}" fill="#2878b5"/>',
                f'<text x="{x:.1f}" y="{max(55, true_y-10):.1f}" text-anchor="middle" font-size="11" font-family="sans-serif">{significance(float(row["p_two_sided"]))}</text>',
                f'<text transform="translate({x+4:.1f},{top+chart_h+12}) rotate(55)" text-anchor="start" font-size="11" font-family="sans-serif">{escape(str(row["subclass"]))}</text>',
            ]
        out += [
            f'<text transform="translate(20,{top+chart_h/2}) rotate(-90)" text-anchor="middle" font-size="15" font-family="sans-serif">Z score of mean enriched-cell count</text>',
            f'<rect x="{width-325}" y="50" width="14" height="14" fill="#f5ad63" stroke="#b85c00"/><text x="{width-305}" y="62" font-size="12" font-family="sans-serif">random replicate Z</text>',
            f'<polygon points="{width-155},50 {width-149},56 {width-155},62 {width-161},56" fill="#2878b5"/><text x="{width-143}" y="61" font-size="12" font-family="sans-serif">true Z</text>',
            "</svg>",
        ]
        panel_path = path.with_name(f"{path.stem}_part_{panel_number:02d}{path.suffix}")
        panel_path.write_text("\n".join(out), encoding="utf-8")
        written.append(panel_path)
    return written


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="按 subclass 将 enriched cell num 转为置换零分布 Z 值并比较")
    parser.add_argument("--true-dir", type=Path, default=Path("true_data"))
    parser.add_argument("--random-dir", type=Path, default=Path("data1"), help="包含 sample_1 ... sample_500 的目录")
    parser.add_argument("--workers", type=int, default=min(20, os.cpu_count() or 1))
    parser.add_argument("--out-summary", type=Path, default=Path("subclass_enriched_cell_zscore_summary.tsv"))
    parser.add_argument("--out-detail", type=Path, default=Path("subclass_enriched_cell_zscore_random_detail.tsv"))
    parser.add_argument("--out-plot", type=Path, default=Path("subclass_enriched_cell_zscore_comparison.svg"))
    parser.add_argument("--subclasses-per-plot", type=int, default=20)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    true = load_true(args.true_dir)
    random, file_count = load_random_parallel(args.random_dir, args.workers)
    summary, detail, skipped = compute_results(true, random)
    write_tsv(args.out_summary, summary)
    write_tsv(args.out_detail, detail)
    panels = draw_svg(summary, args.out_plot, args.subclasses_per_plot)
    print(f"真实 subclass 数: {len(true)}")
    print(f"随机文件数: {file_count}; 随机 sample 数: {len(random)}")
    if len(random) != 500:
        print(f"警告: 检测到 {len(random)} 个随机 sample（预期 500）")
    print(f"成功比较 subclass 数: {len(summary)}; 跳过: {len(skipped)}")
    for reason in skipped:
        print(f"跳过: {reason}")
    print(f"输出: {args.out_summary}, {args.out_detail}")
    for panel in panels:
        print(f"输出: {panel}")


if __name__ == "__main__":
    main()
