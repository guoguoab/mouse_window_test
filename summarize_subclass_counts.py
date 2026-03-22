#!/usr/bin/env python3
"""并行统计 data 目录下所有 txt 文件的 subclass 细胞数量并汇总。"""

from __future__ import annotations

import argparse
import csv
from collections import Counter
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
from typing import Iterable


DEFAULT_WORKERS = 20


def count_subclass_in_file(file_path: Path) -> Counter[str]:
    """读取单个 txt 文件并返回该文件的 subclass 计数。"""
    counts: Counter[str] = Counter()

    with file_path.open("r", encoding="utf-8", newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        if not reader.fieldnames or "subclass" not in reader.fieldnames:
            raise ValueError(f"文件缺少 subclass 列: {file_path}")

        for row in reader:
            subclass = (row.get("subclass") or "").strip()
            if subclass:
                counts[subclass] += 1

    print(f"[进度] 已处理文件: {file_path.name}，subclass 数: {len(counts)}")
    return counts


def merge_counters(counters: Iterable[Counter[str]]) -> Counter[str]:
    """合并多个 Counter。"""
    merged: Counter[str] = Counter()
    for counter in counters:
        merged.update(counter)
    return merged


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="并行读取 txt 文件，按 subclass 汇总细胞数量，并输出 'subclass ;number'。"
    )
    parser.add_argument("--data-dir", default="data", help="txt 文件所在目录（默认: data）")
    parser.add_argument(
        "--workers",
        type=int,
        default=DEFAULT_WORKERS,
        help=f"并行进程数（默认: {DEFAULT_WORKERS}）",
    )
    parser.add_argument(
        "--output",
        default="subclass_summary.txt",
        help="输出文件路径（默认: subclass_summary.txt）",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    data_dir = Path(args.data_dir)

    if not data_dir.exists() or not data_dir.is_dir():
        raise FileNotFoundError(f"数据目录不存在或不是目录: {data_dir}")

    txt_files = sorted(data_dir.glob("*.txt"))
    if not txt_files:
        raise FileNotFoundError(f"目录中未找到 txt 文件: {data_dir}")

    workers = max(1, args.workers)
    with ProcessPoolExecutor(max_workers=workers) as executor:
        counters = executor.map(count_subclass_in_file, txt_files)
        total_counts = merge_counters(counters)

    output_path = Path(args.output)
    with output_path.open("w", encoding="utf-8") as f:
        f.write("subclass ;number\n")
        for subclass, number in sorted(total_counts.items()):
            f.write(f"{subclass} ;{number}\n")

    print(f"已处理文件数: {len(txt_files)}")
    print(f"汇总 subclass 数: {len(total_counts)}")
    print(f"结果已输出到: {output_path}")


if __name__ == "__main__":
    main()
