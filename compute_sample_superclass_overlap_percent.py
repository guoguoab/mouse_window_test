#!/usr/bin/env python3
"""按 superclass 统计 sample overlap 百分比并输出 CSV。

逻辑参考 plot_sample_overlap_permutation.py 中 sample_pair_overlap_percent.csv 的计算：
1) 从 Merfish_brain_cell_type_subclass.txt 建立 class -> superclass 映射；
   - 若 superclass 为 none/空，则退回 cell_Neuron_type。
2) 对所有映射得到的 superclass 去重（unique），构建有向 pair_type（A-B）。
3) 读取 overlap_results/*_cellid_overlap_summary_filtered.csv，
   将 a_class/b_class 映射到 superclass 后，按 overlap 分箱统计百分比。
4) 输出新的 sample 百分比表。
"""

from __future__ import annotations

import argparse
import csv
from multiprocessing import Pool
from pathlib import Path
from typing import Dict, Iterable, List, Tuple

PairKey = Tuple[str, str]

BINS = ["0–20%", "20–40%", "40–60%", "60–80%", "80–100%"]


def normalize_label(value: str | None) -> str:
    token = (value or "").strip()
    return token if token else "Unknown"


def load_superclass_mapping(path: Path) -> Tuple[Dict[str, str], List[str]]:
    class_to_super: Dict[str, str] = {}
    unique_supers: set[str] = set()

    with path.open("r", encoding="utf-8", newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            cls = normalize_label(row.get("class"))
            if cls == "Unknown":
                continue

            superclass = normalize_label(row.get("superclass"))
            neuron_type = normalize_label(row.get("cell_Neuron_type"))

            if superclass.lower() == "none" or superclass == "Unknown":
                superclass = neuron_type

            superclass = normalize_label(superclass)
            class_to_super[cls] = superclass
            unique_supers.add(superclass)

    if not class_to_super:
        raise ValueError(f"No class mapping loaded from: {path}")

    return class_to_super, sorted(unique_supers)


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


def pair_type_counts(total_rows: Iterable[Dict[str, object]], pair_order: List[PairKey]) -> Dict[PairKey, int]:
    counts = {pair_name: 0 for pair_name in pair_order}
    for row in total_rows:
        key = (str(row['a_super']), str(row['b_super']))
        if key in counts:
            counts[key] += 1
    return counts


def collect_directional_values(total_rows: Iterable[Dict[str, object]], pair_name: PairKey) -> List[float]:
    left_type, right_type = pair_name
    values: List[float] = []
    for row in total_rows:
        ab = (str(row['a_super']), str(row['b_super']))
        a2b = float(row["overlap_a_in_b"])
        b2a = float(row["overlap_b_in_a"])

        if left_type == right_type:
            if ab == pair_name:
                values.extend([a2b, b2a])
            continue

        reverse_name = (right_type, left_type)
        if ab == pair_name:
            values.append(a2b)
        elif ab == reverse_name:
            values.append(b2a)
    return values


def denominator_for_pair(pair_name: PairKey, counts: Dict[PairKey, int]) -> int:
    left_type, right_type = pair_name
    if left_type == right_type:
        return 2 * counts.get(pair_name, 0)
    reverse_name = (right_type, left_type)
    return counts.get(pair_name, 0) + counts.get(reverse_name, 0)


def compute_distribution(values: List[float], denominator: int) -> Dict[str, float]:
    counts = {b: 0 for b in BINS}
    if denominator <= 0:
        return {b: 0.0 for b in BINS}

    for v in values:
        counts[overlap_bin(v)] += 1
    return {b: counts[b] * 100.0 / denominator for b in BINS}


def map_class_to_superclass(cls: str, class_to_super: Dict[str, str]) -> str:
    cls = normalize_label(cls)
    if cls in class_to_super:
        return class_to_super[cls]

    low = cls.lower()
    if "gaba" in low:
        return "Gaba"
    if "glut" in low:
        return "Glut"
    if "non" in low and "neuron" in low:
        return "NonNeuron"
    return cls


def read_total_table(path: Path, class_to_super: Dict[str, str]) -> List[Dict[str, object]]:
    rows: List[Dict[str, object]] = []
    with path.open("r", encoding="utf-8", newline="") as fh:
        reader = csv.DictReader(fh)
        for row in reader:
            a_class = normalize_label(row.get("a_class"))
            b_class = normalize_label(row.get("b_class"))
            a_super = map_class_to_superclass(a_class, class_to_super)
            b_super = map_class_to_superclass(b_class, class_to_super)

            try:
                oa = float(row.get("overlap_a_in_b") or "")
                ob = float(row.get("overlap_b_in_a") or "")
            except ValueError:
                continue

            # 与原逻辑保持一致：至少一个方向 > 0
            if oa <= 0 and ob <= 0:
                continue

            rows.append(
                {
                    "a_super": a_super,
                    "b_super": b_super,
                    "overlap_a_in_b": oa,
                    "overlap_b_in_a": ob,
                }
            )
    return rows


def compute_sample_distribution(
    path: Path,
    class_to_super: Dict[str, str],
    pair_order: List[PairKey],
) -> Dict[PairKey, Dict[str, float]]:
    total_rows = read_total_table(path, class_to_super)
    counts = pair_type_counts(total_rows, pair_order)
    dist_by_pair: Dict[PairKey, Dict[str, float]] = {}
    for pair_name in pair_order:
        vals = collect_directional_values(total_rows, pair_name)
        den = denominator_for_pair(pair_name, counts)
        dist_by_pair[pair_name] = compute_distribution(vals, den)
    return dist_by_pair


def compute_sample_distribution_task(
    task: Tuple[Path, Dict[str, str], List[PairKey]],
) -> Tuple[str, Dict[PairKey, Dict[str, float]]]:
    path, class_to_super, pair_order = task
    sample_name = path.stem.replace("_cellid_overlap_summary_filtered", "")
    return sample_name, compute_sample_distribution(path, class_to_super, pair_order)


def write_sample_table(
    sample_dists: Dict[str, Dict[PairKey, Dict[str, float]]],
    pair_order: List[PairKey],
    out_csv: Path,
) -> None:
    out_csv.parent.mkdir(parents=True, exist_ok=True)
    with out_csv.open("w", encoding="utf-8", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=["sample", "pair_type", "overlap_bin", "percent"])
        writer.writeheader()
        for sample, dist_by_pair in sorted(sample_dists.items()):
            for pair_name in pair_order:
                pair_label = f"{pair_name[0]}-{pair_name[1]}"
                for b in BINS:
                    writer.writerow(
                        {
                            "sample": sample,
                            "pair_type": pair_label,
                            "overlap_bin": b,
                            "percent": dist_by_pair[pair_name][b],
                        }
                    )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Compute sample overlap percentages grouped by superclass pair type.")
    parser.add_argument("--overlap-dir", type=Path, default=Path("overlap_results"))
    parser.add_argument("--mapping", type=Path, default=Path("Merfish_brain_cell_type_subclass.txt"))
    parser.add_argument("--out-csv", type=Path, default=Path("overlap_results/sample_pair_overlap_percent_superclass.csv"))
    parser.add_argument("--sample-pattern", type=str, default="*_cellid_overlap_summary_filtered.csv")
    parser.add_argument("--jobs", type=int, default=20, help="并行处理 sample 文件的进程数（默认 20）")
    return parser.parse_args()


def main() -> None:
    args = parse_args()

    class_to_super, unique_supers = load_superclass_mapping(args.mapping)
    pair_order: List[PairKey] = [(a, b) for a in unique_supers for b in unique_supers]

    sample_files = sorted(args.overlap_dir.glob(args.sample_pattern))
    if not sample_files:
        raise SystemExit(f"No sample files matched: {args.overlap_dir / args.sample_pattern}")

    tasks = [(path, class_to_super, pair_order) for path in sample_files]
    jobs = max(1, args.jobs)

    sample_dists: Dict[str, Dict[PairKey, Dict[str, float]]] = {}
    total = len(tasks)
    with Pool(processes=jobs) as pool:
        for idx, (sample_name, dist) in enumerate(pool.imap_unordered(compute_sample_distribution_task, tasks), start=1):
            sample_dists[sample_name] = dist
            print(f"[{idx}/{total}] Processed sample: {sample_name}")

    write_sample_table(sample_dists, pair_order, args.out_csv)

    print(f"Samples: {len(sample_dists)}")
    print(f"Unique superclasses: {len(unique_supers)}")
    print(f"Superclass pair types: {len(pair_order)}")
    print(f"Output -> {args.out_csv}")


if __name__ == "__main__":
    main()
