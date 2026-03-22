#!/usr/bin/env python3
"""并行统计 data 目录下所有 txt 文件中每个 subclass 的细胞数量，并汇总输出。"""

from __future__ import annotations

import argparse
import csv
from collections import Counter
from concurrent.futures import ProcessPoolExecutor
from itertools import chain
from pathlib import Path
from typing import Iterable



def iter_aligned_rows(file_path: Path) -> Iterable[dict[str, str]]:
    """读取 tab 分隔文本；若表头少一列（常见于首列无列名）则自动向后对齐。"""
    with file_path.open("r", encoding="utf-8", newline="") as f:
        reader = csv.reader(f, delimiter="\t")

        try:
            header = next(reader)
        except StopIteration as exc:
            raise ValueError(f"空文件: {file_path}") from exc

        try:
            first_row = next(reader)
        except StopIteration:
            first_row = None

        if first_row is not None and len(first_row) == len(header) + 1:
            header = ["_unnamed_index"] + header

        rows_iter = reader if first_row is None else chain([first_row], reader)
        for row in rows_iter:
            if len(row) < len(header):
                row = row + [""] * (len(header) - len(row))
            elif len(row) > len(header):
                row = row[: len(header)]
            yield dict(zip(header, row))


def count_subclass_in_file(file_path: Path) -> Counter:
    """统计单个 txt 文件中的 subclass 计数。"""
    counter: Counter[str] = Counter()

    rows = iter_aligned_rows(file_path)
    first_row = next(rows, None)
    if first_row is None:
        return counter

    if "subclass" not in first_row:
        raise ValueError(f"文件缺少 subclass 列: {file_path}")

    for row in chain([first_row], rows):
        subclass = (row.get("subclass") or "").strip()
        if subclass:
            counter[subclass] += 1

    return counter


def merge_counters(counters: Iterable[Counter]) -> Counter:
    total: Counter[str] = Counter()
    for counter in counters:
        total.update(counter)
    return total


def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "并行读取指定目录内所有 .txt 文件，按 subclass 统计每个文件的细胞数，"
            "最终累加为总表并输出为 'subclass;number'。"
        )
    )
    parser.add_argument(
        "--input-dir",
        default="data",
        help="待统计的 txt 文件目录（默认: data）",
    )
    parser.add_argument(
        "--workers",
        type=int,
        default=20,
        help="并行进程数（默认: 20）",
    )
    parser.add_argument(
        "--output",
        default="subclass_total_summary.txt",
        help="输出汇总文件路径（默认: subclass_total_summary.txt）",
    )
    args = parser.parse_args()

    input_dir = Path(args.input_dir)
    if not input_dir.exists() or not input_dir.is_dir():
        raise FileNotFoundError(f"输入目录不存在或不是目录: {input_dir}")

    txt_files = sorted(input_dir.glob("*.txt"))
    if not txt_files:
        raise FileNotFoundError(f"目录中未找到 .txt 文件: {input_dir}")

    max_workers = max(1, args.workers)
    with ProcessPoolExecutor(max_workers=max_workers) as executor:
        file_counters = list(executor.map(count_subclass_in_file, txt_files))

    total_counter = merge_counters(file_counters)

    output_path = Path(args.output)
    with output_path.open("w", encoding="utf-8", newline="") as out_f:
        out_f.write("subclass;number\n")
        for subclass, number in sorted(total_counter.items(), key=lambda x: (-x[1], x[0])):
            out_f.write(f"{subclass};{number}\n")

    print(f"处理文件数: {len(txt_files)}")
    print(f"输出文件: {output_path}")


if __name__ == "__main__":
    main()
