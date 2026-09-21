#!/usr/bin/env python3
"""Compare the true enriched-cell Z score with the largest random Z score.

For each subclass, one random statistic is first calculated for every
permutation: the mean enriched-cell count among that permutation's clusters.
The true mean and all random means are standardized with the mean and sample
SD of those random statistics.  This conservative variant then compares
``true_z`` with ``max(random_z)`` rather than with the centre of the null
distribution.
"""

from __future__ import annotations

import argparse
import math
import os
from pathlib import Path
from statistics import mean, stdev
from xml.sax.saxutils import escape

from compare_subclass_enriched_cell_size import load_random_parallel, load_true, write_tsv
from compare_subclass_enriched_cell_zscore import z_score


def compute_results(
    true: dict[str, list[float]], random: dict[str, dict[str, list[float]]]
) -> tuple[list[dict[str, object]], list[str]]:
    """Return maximum-null comparisons and descriptions of skipped subclasses."""
    rows: list[dict[str, object]] = []
    skipped: list[str] = []
    samples = sorted(random)

    for subclass in sorted(true):
        replicate_rows = [
            (sample, mean(counts))
            for sample in samples
            if (counts := random[sample].get(subclass, []))
        ]
        if len(replicate_rows) < 2:
            skipped.append(f"{subclass}: 少于 2 个随机重复")
            continue

        replicate_means = [value for _, value in replicate_rows]
        null_mean = mean(replicate_means)
        null_sd = stdev(replicate_means)
        if not math.isfinite(null_sd) or null_sd <= 0:
            skipped.append(f"{subclass}: 随机重复均值的标准差为 0 或非有限值")
            continue

        true_mean = mean(true[subclass])
        true_z = z_score(true_mean, null_mean, null_sd)
        random_rows = [
            (sample, z_score(value, null_mean, null_sd))
            for sample, value in replicate_rows
        ]
        random_z_max = max(value for _, value in random_rows)
        max_samples = [sample for sample, value in random_rows if value == random_z_max]
        difference = true_z - random_z_max

        rows.append({
            "subclass": subclass,
            "true_cluster_n": len(true[subclass]),
            "random_sample_n": len(random_rows),
            "true_mean_enriched_cell_num": true_mean,
            "null_mean_enriched_cell_num": null_mean,
            "null_sd_enriched_cell_num": null_sd,
            "true_z": true_z,
            "random_z_max": random_z_max,
            "random_z_max_sample": ",".join(max_samples),
            "difference_true_z_minus_random_z_max": difference,
            "true_z_gt_random_z_max": true_z > random_z_max,
        })

    return rows, skipped


def draw_svg(
    rows: list[dict[str, object]], path: Path, subclasses_per_plot: int
) -> list[Path]:
    """Draw paired true-Z and maximum-random-Z markers in panelled SVGs."""
    if not rows:
        raise RuntimeError("没有可与随机最大 Z 值比较的 subclass，无法绘图")
    if subclasses_per_plot < 1:
        raise ValueError("每张图的 subclass 数必须至少为 1")

    rows = sorted(
        rows,
        key=lambda row: float(row["difference_true_z_minus_random_z_max"]),
        reverse=True,
    )
    values = [0.0]
    for row in rows:
        values.extend((float(row["true_z"]), float(row["random_z_max"])))
    low, high = min(values), max(values)
    padding = max((high - low) * 0.1, 0.5)
    y_min, y_max = low - padding, high + padding
    panels = [
        rows[index:index + subclasses_per_plot]
        for index in range(0, len(rows), subclasses_per_plot)
    ]
    path.parent.mkdir(parents=True, exist_ok=True)
    written: list[Path] = []

    for panel_number, panel in enumerate(panels, 1):
        width, height = max(1000, 105 + 70 * len(panel)), 700
        left, right, top, bottom = 85, 35, 70, 230
        chart_h, chart_w = height - top - bottom, width - left - right
        y = lambda value: top + chart_h * (y_max - float(value)) / (y_max - y_min)
        out = [
            f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}">',
            '<rect width="100%" height="100%" fill="white"/>',
            f'<text x="{width/2}" y="35" text-anchor="middle" font-size="24" font-family="sans-serif" font-weight="bold">True Z vs maximum permutation Z by subclass (part {panel_number}/{len(panels)})</text>',
        ]
        for tick in range(6):
            value = y_min + (y_max - y_min) * tick / 5
            py = y(value)
            out += [
                f'<line x1="{left}" y1="{py:.1f}" x2="{width-right}" y2="{py:.1f}" stroke="#dddddd"/>',
                f'<text x="{left-8}" y="{py+5:.1f}" text-anchor="end" font-size="12" font-family="sans-serif">{value:.2f}</text>',
            ]
        step = chart_w / len(panel)
        for index, row in enumerate(panel):
            x = left + step * (index + 0.5)
            true_y, max_y = y(row["true_z"]), y(row["random_z_max"])
            colour = "#198754" if bool(row["true_z_gt_random_z_max"]) else "#777777"
            out += [
                f'<line x1="{x:.1f}" y1="{true_y:.1f}" x2="{x:.1f}" y2="{max_y:.1f}" stroke="{colour}" stroke-width="2"/>',
                f'<circle cx="{x:.1f}" cy="{max_y:.1f}" r="5" fill="#d97820"/>',
                f'<polygon points="{x:.1f},{true_y-7:.1f} {x+7:.1f},{true_y:.1f} {x:.1f},{true_y+7:.1f} {x-7:.1f},{true_y:.1f}" fill="#2878b5"/>',
                f'<text x="{x:.1f}" y="{min(true_y, max_y)-11:.1f}" text-anchor="middle" font-size="12" font-family="sans-serif" fill="{colour}">{"PASS" if bool(row["true_z_gt_random_z_max"]) else ""}</text>',
                f'<text transform="translate({x+4:.1f},{top+chart_h+12}) rotate(55)" text-anchor="start" font-size="11" font-family="sans-serif">{escape(str(row["subclass"]))}</text>',
            ]
        out += [
            f'<text transform="translate(20,{top+chart_h/2}) rotate(-90)" text-anchor="middle" font-size="15" font-family="sans-serif">Z score of mean enriched-cell count</text>',
            f'<circle cx="{width-310}" cy="57" r="5" fill="#d97820"/><text x="{width-298}" y="62" font-size="12" font-family="sans-serif">maximum random Z</text>',
            f'<polygon points="{width-155},50 {width-148},57 {width-155},64 {width-162},57" fill="#2878b5"/><text x="{width-143}" y="62" font-size="12" font-family="sans-serif">true Z</text>',
            "</svg>",
        ]
        panel_path = path.with_name(
            f"{path.stem}_part_{panel_number:02d}{path.suffix}"
        )
        panel_path.write_text("\n".join(out), encoding="utf-8")
        written.append(panel_path)
    return written


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="按 subclass 将真实 Z 值与随机 Z 最大值比较")
    parser.add_argument("--true-dir", type=Path, default=Path("true_data"))
    parser.add_argument("--random-dir", type=Path, default=Path("data1"))
    parser.add_argument("--workers", type=int, default=min(20, os.cpu_count() or 1))
    parser.add_argument(
        "--out-summary",
        type=Path,
        default=Path("subclass_enriched_cell_zscore_max_summary.tsv"),
    )
    parser.add_argument(
        "--out-plot",
        type=Path,
        default=Path("subclass_enriched_cell_zscore_max_comparison.svg"),
    )
    parser.add_argument("--subclasses-per-plot", type=int, default=20)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    true = load_true(args.true_dir)
    random, file_count = load_random_parallel(args.random_dir, args.workers)
    rows, skipped = compute_results(true, random)
    write_tsv(args.out_summary, rows)
    panels = draw_svg(rows, args.out_plot, args.subclasses_per_plot)
    passed = sum(bool(row["true_z_gt_random_z_max"]) for row in rows)
    print(f"真实 subclass 数: {len(true)}")
    print(f"随机文件数: {file_count}; 随机 sample 数: {len(random)}")
    if len(random) != 500:
        print(f"警告: 检测到 {len(random)} 个随机 sample（预期 500）")
    print(f"成功比较 subclass 数: {len(rows)}; 真实 Z 超过随机最大 Z: {passed}; 跳过: {len(skipped)}")
    for reason in skipped:
        print(f"跳过: {reason}")
    print(f"输出: {args.out_summary}")
    for panel in panels:
        print(f"输出: {panel}")


if __name__ == "__main__":
    main()
