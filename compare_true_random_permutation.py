#!/usr/bin/env python3
"""比较 sample_level_summary 与 true_data 的 region 数目，并计算置换检验风格 p 值。"""

from __future__ import annotations

import argparse
import csv
import glob
import os
from collections import Counter, defaultdict
from typing import Dict, Iterable, List, Tuple


Key = Tuple[str, str, str]  # (slide, subclass, layer)


def normalize_layer(layer: str) -> str:
    if layer is None:
        return ""
    value = layer.strip()
    value = value.replace("2&3", "2/3")
    value = value.replace("2.3", "2/3")
    return value


def clean_subclass(value: str) -> str:
    if value is None:
        return ""
    return " ".join(value.strip().split())


def read_tsv(path: str) -> Iterable[dict]:
    with open(path, "r", encoding="utf-8") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            yield row


def load_random_summary(path: str) -> Dict[Key, List[Tuple[str, float]]]:
    random_map: Dict[Key, List[Tuple[str, float]]] = defaultdict(list)

    for row in read_tsv(path):
        try:
            region_num = float(row["region_num"])
        except (KeyError, TypeError, ValueError):
            continue

        key = (
            row.get("slide", "").strip(),
            clean_subclass(row.get("subclass", "")),
            normalize_layer(row.get("layer", "")),
        )
        sample = row.get("sample", "").strip() or "NA"

        if not all(key):
            continue

        random_map[key].append((sample, region_num))

    return random_map


def choose_true_files(true_dir: str) -> List[str]:
    patt = os.path.join(true_dir, "*_merged_regions_table_cell_id.txt")
    files = sorted(glob.glob(patt))
    if files:
        return files

    # 兼容拼写 merged_rgions
    patt2 = os.path.join(true_dir, "*_merged_rgions_loc.txt")
    return sorted(glob.glob(patt2))


def load_true_counts(true_dir: str) -> Tuple[Dict[Key, List[float]], Counter]:
    true_files = choose_true_files(true_dir)
    if not true_files:
        raise FileNotFoundError(f"在 {true_dir} 中找不到可用真实文件")

    true_map: Dict[Key, List[float]] = defaultdict(list)
    skipped = Counter()

    for path in true_files:
        rows = list(read_tsv(path))
        if not rows:
            skipped["empty_file"] += 1
            continue

        first = rows[0]
        slide = first.get("slide", "").strip()
        layer = normalize_layer(first.get("layer", ""))
        subclass = clean_subclass(first.get("Class", first.get("subclass", "")))

        if not (slide and layer and subclass):
            skipped["missing_key_columns"] += 1
            continue

        regions = {
            r.get("region", "").strip()
            for r in rows
            if r.get("region", "").strip()
        }
        if not regions:
            skipped["missing_region"] += 1
            continue

        key = (slide, subclass, layer)
        true_map[key].append(float(len(regions)))

    return true_map, skipped


def compute_statistics(
    true_map: Dict[Key, List[float]],
    random_map: Dict[Key, List[Tuple[str, float]]],
) -> Tuple[List[dict], List[dict], Counter]:
    summary_rows: List[dict] = []
    detail_rows: List[dict] = []
    missing = Counter()

    for key in sorted(true_map.keys()):
        slide, subclass, layer = key
        true_values = true_map[key]
        true_region_num = sum(true_values) / len(true_values)

        random_values = random_map.get(key)
        if not random_values:
            missing["no_random_match"] += 1
            continue

        n = len(random_values)
        true_gt = 0
        random_ge = 0

        for sample, random_region_num in random_values:
            is_true_gt = int(true_region_num > random_region_num)
            if is_true_gt:
                true_gt += 1
            else:
                random_ge += 1

            detail_rows.append(
                {
                    "slide": slide,
                    "subclass": subclass,
                    "layer": layer,
                    "sample": sample,
                    "true_region_num": f"{true_region_num:.6f}",
                    "random_region_num": f"{random_region_num:.6f}",
                    "true_gt_random": is_true_gt,
                }
            )

        # 单侧：H1 true > random，越小越显著
        p_value = (random_ge + 1.0) / (n + 1.0)

        summary_rows.append(
            {
                "slide": slide,
                "subclass": subclass,
                "layer": layer,
                "true_region_num": f"{true_region_num:.6f}",
                "true_file_count": len(true_values),
                "random_sample_n": n,
                "true_gt_random_count": true_gt,
                "true_gt_random_ratio": f"{(true_gt / n):.6f}",
                "perm_p_value_one_sided": f"{p_value:.6g}",
            }
        )

    return summary_rows, detail_rows, missing


def write_tsv(path: str, rows: List[dict]) -> None:
    if not rows:
        with open(path, "w", encoding="utf-8") as f:
            f.write("")
        return

    with open(path, "w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="比较真实 region 数和随机 sample region_num，并输出置换检验风格统计"
    )
    parser.add_argument(
        "--sample-summary",
        default="sample_level_summary.txt",
        help="number.R 生成的 sample_level_summary.txt",
    )
    parser.add_argument(
        "--true-dir",
        default="true_data",
        help="真实数据目录（文件需包含 slide/layer/Class 或 subclass/region 列）",
    )
    parser.add_argument(
        "--out-summary",
        default="true_vs_random_permutation_summary.tsv",
        help="汇总输出文件",
    )
    parser.add_argument(
        "--out-detail",
        default="true_vs_random_permutation_detail.tsv",
        help="逐 sample 输出文件",
    )
    args = parser.parse_args()

    random_map = load_random_summary(args.sample_summary)
    if not random_map:
        raise RuntimeError(f"未从 {args.sample_summary} 读取到有效随机汇总")

    true_map, skipped = load_true_counts(args.true_dir)
    if not true_map:
        raise RuntimeError(f"未从 {args.true_dir} 读取到有效真实统计")

    summary_rows, detail_rows, missing = compute_statistics(true_map, random_map)

    write_tsv(args.out_summary, summary_rows)
    write_tsv(args.out_detail, detail_rows)

    print("完成")
    print(f"- 随机键数量: {len(random_map)}")
    print(f"- 真实键数量: {len(true_map)}")
    print(f"- 成功匹配键数量: {len(summary_rows)}")
    print(f"- 逐sample比较记录数: {len(detail_rows)}")
    if skipped:
        print(f"- 跳过真实文件统计: {dict(skipped)}")
    if missing:
        print(f"- 未匹配统计: {dict(missing)}")
    print(f"- 汇总输出: {args.out_summary}")
    print(f"- 详情输出: {args.out_detail}")


if __name__ == "__main__":
    main()
