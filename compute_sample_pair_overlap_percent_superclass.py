#!/usr/bin/env python3
"""按 superclass 统计 sample pair overlap 百分比并输出 CSV。

逻辑参考 plot_sample_overlap_permutation.py 中 sample_pair_overlap_percent.csv 计算部分：
1) 读取 overlap_results/*_cellid_overlap_summary_filtered.csv。
2) 将 a_class/b_class 映射到 superclass：
   - 优先使用 Merfish_brain_cell_type_subclass.txt 的 class -> superclass。
   - 若 superclass 为 none/空，则回退到 class -> cell_Neuron_type。
   - 若映射缺失，再回退到输入行中的 a_cell_Neuron_type/b_cell_Neuron_type。
3) 对所有出现过的 superclass 做 unique，生成 pair_type（A-B，含方向）。
4) 按原脚本相同规则统计 overlap 分箱百分比（0-20,20-40,40-60,60-80,80-100）。
5) 输出新的 CSV。
"""

from __future__ import annotations

import argparse
import csv
from multiprocessing import Pool
from pathlib import Path
from typing import Dict, Iterable, List, Set, Tuple

BINS = ["0–20%", "20–40%", "40–60%", "60–80%", "80–100%"]


def normalize_token(value: str | None) -> str:
    token = (value or "").strip()
    return token if token else "Unknown"


def is_none_like(value: str | None) -> bool:
    return (value or "").strip().lower() in {"", "none", "na", "nan", "null"}


def load_class_to_superclass(mapping_path: Path) -> Dict[str, str]:
    class2super: Dict[str, str] = {}
    with mapping_path.open("r", encoding="utf-8", newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            cls = normalize_token(row.get("class"))
            neuron_type = normalize_token(row.get("cell_Neuron_type"))
            superclass = normalize_token(row.get("superclass"))
            class2super[cls] = neuron_type if is_none_like(superclass) else superclass
    return class2super


def class_to_superclass(raw_class: str | None, raw_neuron: str | None, class2super: Dict[str, str]) -> str:
    cls = normalize_token(raw_class)
    if cls in class2super:
        return class2super[cls]

    neuron = normalize_token(raw_neuron)
    if neuron != "Unknown":
        return neuron
    return cls


def overlap_bin(x: float) -> str:
    if x < 0.2:
        return "0–20%"
    if x < 0.4:
        return "20–40%"
    if x < 0.6:
        return "40–60%"
    if x < 0.8:
        return "60–80%"
    return "80–100%"


def split_pair(pair_name: str) -> Tuple[str, str]:
    left, right = pair_name.split("-", 1)
    return left, right


def read_total_table(path: Path, class2super: Dict[str, str]) -> List[Dict[str, object]]:
    rows: List[Dict[str, object]] = []
    with path.open("r", encoding="utf-8", newline="") as fh:
        reader = csv.DictReader(fh)
        for row in reader:
            try:
                oa = float(row.get("overlap_a_in_b") or "")
                ob = float(row.get("overlap_b_in_a") or "")
            except ValueError:
                continue

            if oa <= 0 and ob <= 0:
                continue

            a_super = class_to_superclass(row.get("a_class"), row.get("a_cell_Neuron_type"), class2super)
            b_super = class_to_superclass(row.get("b_class"), row.get("b_cell_Neuron_type"), class2super)
            rows.append(
                {
                    "a_super": a_super,
                    "b_super": b_super,
                    "overlap_a_in_b": oa,
                    "overlap_b_in_a": ob,
                }
            )
    return rows


def pair_type_counts(total_rows: Iterable[Dict[str, object]], pair_order: List[str]) -> Dict[str, int]:
    counts = {pair_name: 0 for pair_name in pair_order}
    for row in total_rows:
        key = f"{row['a_super']}-{row['b_super']}"
        if key in counts:
            counts[key] += 1
    return counts


def collect_directional_values(total_rows: Iterable[Dict[str, object]], pair_name: str) -> List[float]:
    left_type, right_type = split_pair(pair_name)
    values: List[float] = []
    for row in total_rows:
        ab = f"{row['a_super']}-{row['b_super']}"
        a2b = float(row["overlap_a_in_b"])
        b2a = float(row["overlap_b_in_a"])

        if left_type == right_type:
            if ab == pair_name:
                values.extend([a2b, b2a])
            continue

        reverse_name = f"{right_type}-{left_type}"
        if ab == pair_name:
            values.append(a2b)
        elif ab == reverse_name:
            values.append(b2a)
    return values


def denominator_for_pair(pair_name: str, counts: Dict[str, int]) -> int:
    left_type, right_type = split_pair(pair_name)
    if left_type == right_type:
        return 2 * counts.get(pair_name, 0)
    reverse_name = f"{right_type}-{left_type}"
    return counts.get(pair_name, 0) + counts.get(reverse_name, 0)


def compute_distribution(values: List[float], denominator: int) -> Dict[str, float]:
    counts = {b: 0 for b in BINS}
    if denominator <= 0:
        return {b: 0.0 for b in BINS}

    for v in values:
        counts[overlap_bin(v)] += 1
    return {b: counts[b] * 100.0 / denominator for b in BINS}


def compute_sample_distribution(path: Path, class2super: Dict[str, str], pair_order: List[str]) -> Dict[str, Dict[str, float]]:
    total_rows = read_total_table(path, class2super)
    counts = pair_type_counts(total_rows, pair_order)

    dist_by_pair: Dict[str, Dict[str, float]] = {}
    for pair_name in pair_order:
        vals = collect_directional_values(total_rows, pair_name)
        den = denominator_for_pair(pair_name, counts)
        dist_by_pair[pair_name] = compute_distribution(vals, den)
    return dist_by_pair


def gather_unique_superclasses(sample_files: List[Path], class2super: Dict[str, str]) -> List[str]:
    tokens: Set[str] = set()
    for path in sample_files:
        with path.open("r", encoding="utf-8", newline="") as fh:
            reader = csv.DictReader(fh)
            for row in reader:
                a_super = class_to_superclass(row.get("a_class"), row.get("a_cell_Neuron_type"), class2super)
                b_super = class_to_superclass(row.get("b_class"), row.get("b_cell_Neuron_type"), class2super)
                tokens.add(a_super)
                tokens.add(b_super)
    return sorted(tokens)


def compute_task(task: Tuple[Path, Dict[str, str], List[str]]) -> Tuple[str, Dict[str, Dict[str, float]]]:
    path, class2super, pair_order = task
    sample_name = path.stem.replace("_cellid_overlap_summary_filtered", "")
    return sample_name, compute_sample_distribution(path, class2super, pair_order)


def write_sample_table(
    sample_dists: Dict[str, Dict[str, Dict[str, float]]],
    pair_order: List[str],
    out_csv: Path,
) -> None:
    out_csv.parent.mkdir(parents=True, exist_ok=True)
    with out_csv.open("w", encoding="utf-8", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=["sample", "pair_type", "overlap_bin", "percent"])
        writer.writeheader()
        for sample, dist_by_pair in sorted(sample_dists.items()):
            for pair_type in pair_order:
                for bin_name in BINS:
                    writer.writerow(
                        {
                            "sample": sample,
                            "pair_type": pair_type,
                            "overlap_bin": bin_name,
                            "percent": dist_by_pair[pair_type][bin_name],
                        }
                    )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Compute sample overlap percentage table grouped by superclass pairs.")
    parser.add_argument("--overlap-dir", type=Path, default=Path("overlap_results"))
    parser.add_argument("--sample-pattern", type=str, default="*_cellid_overlap_summary_filtered.csv")
    parser.add_argument("--mapping", type=Path, default=Path("Merfish_brain_cell_type_subclass.txt"))
    parser.add_argument("--out-csv", type=Path, default=Path("overlap_results/sample_pair_overlap_percent_superclass.csv"))
    parser.add_argument("--jobs", type=int, default=20)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    sample_files = sorted(args.overlap_dir.glob(args.sample_pattern))
    if not sample_files:
        raise SystemExit(f"No sample files matched: {args.overlap_dir / args.sample_pattern}")

    class2super = load_class_to_superclass(args.mapping)
    superclasses = gather_unique_superclasses(sample_files, class2super)
    pair_order = [f"{a}-{b}" for a in superclasses for b in superclasses]

    jobs = max(1, args.jobs)
    tasks = [(path, class2super, pair_order) for path in sample_files]
    sample_dists: Dict[str, Dict[str, Dict[str, float]]] = {}

    with Pool(processes=jobs) as pool:
        total = len(tasks)
        for idx, (sample_name, dist) in enumerate(pool.imap_unordered(compute_task, tasks), start=1):
            sample_dists[sample_name] = dist
            print(f"[{idx}/{total}] Processed sample: {sample_name}")

    write_sample_table(sample_dists, pair_order, args.out_csv)
    print(f"Unique superclasses ({len(superclasses)}): {', '.join(superclasses)}")
    print(f"Output CSV -> {args.out_csv}")


if __name__ == "__main__":
    main()
