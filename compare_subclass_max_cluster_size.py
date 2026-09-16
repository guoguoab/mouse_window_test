#!/usr/bin/env python3
"""Compare true mean cluster sizes with random maximum cluster sizes by subclass.

Each row in a ``*_merged_regions_table_cell_id.txt`` file represents one
cluster.  ``total_cell_num`` is its size.  The mean true cluster size is compared
with the maximum cluster size in each random replicate (normally sample_1 ...
sample_500).  Random files are parsed in parallel because that is the expensive
input step.
"""

from __future__ import annotations

import argparse
import csv
import math
import os
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
from statistics import mean, median
from typing import Sequence
from xml.sax.saxutils import escape


def clean_subclass(value: str) -> str:
    return " ".join(value.strip().split())


def read_cluster_file(path: Path) -> tuple[str, dict[str, list[float]]]:
    """Return ``(sample, subclass -> cluster sizes)`` for one table."""
    values: dict[str, list[float]] = defaultdict(list)
    sample = ""
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:  # Empty result files contain no clusters.
            return path.parent.name, {}
        subclass_column = "subclass" if "subclass" in reader.fieldnames else "Class"
        if subclass_column not in reader.fieldnames:
            raise ValueError(f"文件缺少 subclass/Class 列: {path}")
        if "total_cell_num" not in reader.fieldnames:
            raise ValueError(f"文件缺少 total_cell_num 列: {path}")
        for row in reader:
            subclass = clean_subclass(row.get(subclass_column, "") or "")
            raw_size = (row.get("total_cell_num") or "").strip()
            if not subclass or not raw_size:
                continue
            try:
                size = float(raw_size)
            except ValueError:
                continue
            if not math.isfinite(size) or size < 0:
                continue
            values[subclass].append(size)
            if not sample:
                sample = (row.get("sample") or "").strip()
    if not sample:
        # Random file names conventionally begin with sample_N_.
        parts = path.name.split("_", 2)
        if len(parts) >= 2 and parts[0] == "sample" and parts[1].isdigit():
            sample = f"sample_{parts[1]}"
        else:
            sample = path.parent.name
    return sample, dict(values)


def find_cluster_files(directory: Path) -> list[Path]:
    files = sorted(directory.rglob("*_merged_regions_table_cell_id.txt"))
    if not files:
        raise FileNotFoundError(f"在 {directory} 中未找到 merged_regions_table_cell_id 文件")
    return files


def load_true(directory: Path) -> dict[str, list[float]]:
    combined: dict[str, list[float]] = defaultdict(list)
    for path in find_cluster_files(directory):
        _, by_subclass = read_cluster_file(path)
        for subclass, sizes in by_subclass.items():
            combined[subclass].extend(sizes)
    return dict(combined)


def load_random_parallel(
    directory: Path, workers: int
) -> tuple[dict[str, dict[str, list[float]]], int]:
    """Read random result files concurrently, then combine by sample/subclass."""
    files = find_cluster_files(directory)
    by_sample: dict[str, dict[str, list[float]]] = defaultdict(lambda: defaultdict(list))
    with ProcessPoolExecutor(max_workers=max(1, workers)) as executor:
        for sample, by_subclass in executor.map(read_cluster_file, files, chunksize=8):
            if not sample:
                continue
            for subclass, sizes in by_subclass.items():
                by_sample[sample][subclass].extend(sizes)
    return {sample: dict(values) for sample, values in by_sample.items()}, len(files)


def empirical_pvalues(
    true_value: float, random_values: Sequence[float]
) -> tuple[float, float, float]:
    n = len(random_values)
    if n == 0:
        return math.nan, math.nan, math.nan
    p_gt = (sum(value >= true_value for value in random_values) + 1) / (n + 1)
    p_lt = (sum(value <= true_value for value in random_values) + 1) / (n + 1)
    center = mean(random_values)
    distance = abs(true_value - center)
    p_two = (sum(abs(value - center) >= distance for value in random_values) + 1) / (n + 1)
    return p_gt, p_lt, p_two


def quantile(values: Sequence[float], fraction: float) -> float:
    ordered = sorted(values)
    if len(ordered) == 1:
        return ordered[0]
    position = (len(ordered) - 1) * fraction
    lower = math.floor(position)
    upper = math.ceil(position)
    if lower == upper:
        return ordered[lower]
    return ordered[lower] * (upper - position) + ordered[upper] * (position - lower)


def compute_results(
    true: dict[str, list[float]], random: dict[str, dict[str, list[float]]]
) -> tuple[list[dict[str, object]], list[dict[str, object]]]:
    detail: list[dict[str, object]] = []
    summary: list[dict[str, object]] = []
    samples = sorted(random)
    for subclass in sorted(true):
        true_sizes = true[subclass]
        true_mean = mean(true_sizes)
        replicate_maxima: list[float] = []
        for sample in samples:
            sizes = random[sample].get(subclass, [])
            if not sizes:
                continue
            sample_maximum = max(sizes)
            replicate_maxima.append(sample_maximum)
            detail.append({
                "subclass": subclass,
                "sample": sample,
                "random_cluster_n": len(sizes),
                "random_max_cluster_size": sample_maximum,
                "true_mean_cluster_size": true_mean,
                "difference_true_mean_minus_random_max": true_mean - sample_maximum,
            })
        if not replicate_maxima:
            continue
        p_gt, p_lt, p_two = empirical_pvalues(true_mean, replicate_maxima)
        random_mean_of_replicate_maxima = mean(replicate_maxima)
        summary.append({
            "subclass": subclass,
            "true_cluster_n": len(true_sizes),
            "true_mean_cluster_size": true_mean,
            "true_median_cluster_size": median(true_sizes),
            "random_sample_n": len(replicate_maxima),
            "random_mean_of_replicate_maxima": random_mean_of_replicate_maxima,
            "random_median_of_replicate_maxima": median(replicate_maxima),
            "random_q1": quantile(replicate_maxima, 0.25),
            "random_q3": quantile(replicate_maxima, 0.75),
            "random_min": min(replicate_maxima),
            "random_max": max(replicate_maxima),
            "difference_true_mean_minus_random_max_mean": true_mean - random_mean_of_replicate_maxima,
            "fold_change_true_mean_over_random_max_mean": true_mean / random_mean_of_replicate_maxima if random_mean_of_replicate_maxima else math.nan,
            "p_one_sided_true_gt": p_gt,
            "p_one_sided_true_lt": p_lt,
            "p_two_sided": p_two,
        })
    return summary, detail


def write_tsv(path: Path, rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        if not rows:
            return
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def significance(p: float) -> str:
    if p < 0.001:
        return "***"
    if p < 0.01:
        return "**"
    if p < 0.05:
        return "*"
    return "ns"


def draw_svg(rows: list[dict[str, object]], path: Path) -> None:
    """Draw random-replicate boxplots and true-value diamonds without dependencies."""
    if not rows:
        raise RuntimeError("真实数据与随机数据没有共同的 subclass，无法绘图")
    rows = sorted(rows, key=lambda row: float(row["true_mean_cluster_size"]), reverse=True)
    width = max(1000, 105 + 70 * len(rows))
    height = 700
    left, right, top, bottom = 85, 35, 70, 230
    chart_h = height - top - bottom
    chart_w = width - left - right
    maximum = max(
        max(float(row["random_max"]), float(row["true_mean_cluster_size"]))
        for row in rows
    ) * 1.12 or 1
    y = lambda value: top + chart_h * (1 - float(value) / maximum)
    out = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}">',
        '<rect width="100%" height="100%" fill="white"/>',
        f'<text x="{width/2}" y="35" text-anchor="middle" font-size="24" font-family="sans-serif" font-weight="bold">True mean vs random maximum cluster size by subclass</text>',
    ]
    for tick in range(6):
        value = maximum * tick / 5
        py = y(value)
        out += [
            f'<line x1="{left}" y1="{py:.1f}" x2="{width-right}" y2="{py:.1f}" stroke="#dddddd"/>',
            f'<text x="{left-8}" y="{py+5:.1f}" text-anchor="end" font-size="12" font-family="sans-serif">{value:.0f}</text>',
        ]
    step = chart_w / len(rows)
    for index, row in enumerate(rows):
        x = left + step * (index + 0.5)
        box_w = min(30, step * 0.55)
        low, q1, med, q3, high = (
            float(row[key])
            for key in (
                "random_min", "random_q1", "random_median_of_replicate_maxima",
                "random_q3", "random_max",
            )
        )
        out += [
            f'<line x1="{x:.1f}" y1="{y(low):.1f}" x2="{x:.1f}" y2="{y(high):.1f}" stroke="#d97820" stroke-width="2"/>',
            f'<rect x="{x-box_w/2:.1f}" y="{y(q3):.1f}" width="{box_w:.1f}" height="{max(1, y(q1)-y(q3)):.1f}" fill="#f5ad63" stroke="#b85c00"/>',
            f'<line x1="{x-box_w/2:.1f}" y1="{y(med):.1f}" x2="{x+box_w/2:.1f}" y2="{y(med):.1f}" stroke="#703500" stroke-width="2"/>',
        ]
        true_y = y(row["true_mean_cluster_size"])
        out += [
            f'<polygon points="{x:.1f},{true_y-6:.1f} {x+6:.1f},{true_y:.1f} {x:.1f},{true_y+6:.1f} {x-6:.1f},{true_y:.1f}" fill="#2878b5"/>',
            f'<text x="{x:.1f}" y="{max(55, true_y-10):.1f}" text-anchor="middle" font-size="11" font-family="sans-serif">{significance(float(row["p_two_sided"]))}</text>',
            f'<text transform="translate({x+4:.1f},{top+chart_h+12}) rotate(55)" text-anchor="start" font-size="11" font-family="sans-serif">{escape(str(row["subclass"]))}</text>',
        ]
    out += [
        f'<text transform="translate(20,{top+chart_h/2}) rotate(-90)" text-anchor="middle" font-size="15" font-family="sans-serif">Cluster size (total_cell_num)</text>',
        f'<rect x="{width-300}" y="50" width="14" height="14" fill="#f5ad63" stroke="#b85c00"/><text x="{width-280}" y="62" font-size="12" font-family="sans-serif">random replicate maxima</text>',
        f'<polygon points="{width-155},50 {width-149},56 {width-155},62 {width-161},56" fill="#2878b5"/><text x="{width-143}" y="61" font-size="12" font-family="sans-serif">true</text>',
        "</svg>",
    ]
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(out), encoding="utf-8")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="按 subclass 比较真实 cluster size 均值与随机 cluster size 最大值，并绘制 SVG"
    )
    parser.add_argument("--true-dir", type=Path, default=Path("true_data"))
    parser.add_argument(
        "--random-dir", type=Path, default=Path("data1"),
        help="包含 sample_1 ... sample_500 的目录",
    )
    parser.add_argument(
        "--workers", type=int, default=min(20, os.cpu_count() or 1),
        help="并行读取随机文件的进程数",
    )
    parser.add_argument(
        "--out-summary", type=Path,
        default=Path("subclass_max_cluster_size_summary.tsv"),
    )
    parser.add_argument(
        "--out-detail", type=Path,
        default=Path("subclass_max_cluster_size_random_detail.tsv"),
    )
    parser.add_argument(
        "--out-plot", type=Path,
        default=Path("subclass_max_cluster_size_comparison.svg"),
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    true = load_true(args.true_dir)
    random, file_count = load_random_parallel(args.random_dir, args.workers)
    summary, detail = compute_results(true, random)
    write_tsv(args.out_summary, summary)
    write_tsv(args.out_detail, detail)
    draw_svg(summary, args.out_plot)
    print(f"真实 subclass 数: {len(true)}")
    print(f"随机文件数: {file_count}; 随机 sample 数: {len(random)}")
    if len(random) != 500:
        print(f"警告: 检测到 {len(random)} 个随机 sample（预期 500）")
    print(f"成功比较 subclass 数: {len(summary)}")
    print(f"输出: {args.out_summary}, {args.out_detail}, {args.out_plot}")


if __name__ == "__main__":
    main()
