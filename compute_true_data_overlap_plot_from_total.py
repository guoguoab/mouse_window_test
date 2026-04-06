#!/usr/bin/env python3
"""基于真实值 total overlap 表按 superclass 计算 overlap 百分比并可选绘图。

实现逻辑对齐 compute_sample_superclass_overlap_percent.py：
1) 读取 Merfish_brain_cell_type_subclass.txt 建立 class -> superclass 映射；
   - 若 superclass 为 none/空，则回退到 cell_Neuron_type。
2) 从 total overlap 表读取 a_class/b_class，映射到 superclass。
3) 仅保留 overlap_a_in_b > 0 或 overlap_b_in_a > 0 的记录。
4) 对所有 superclass 的有向 pair_type（A-B）统计 overlap 分箱百分比。
5) 输出真实值 pair percent 表，并可选输出 SVG 图。
"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path
from typing import Dict, Iterable, List, Tuple


PairKey = Tuple[str, str]
BINS = ["0–20%", "20–40%", "40–60%", "60–80%", "80–100%"]


def normalize_label(value: str | None) -> str:
    token = (value or "").strip()
    return token if token else "Unknown"


def safe_float(value: str | None) -> float:
    try:
        return float(value or "")
    except ValueError:
        return 0.0


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


def load_total_table(path: Path, class_to_super: Dict[str, str]) -> List[Dict[str, object]]:
    rows: List[Dict[str, object]] = []
    with path.open("r", encoding="utf-8", newline="") as fh:
        reader = csv.DictReader(fh)
        for row in reader:
            a_class = normalize_label(row.get("a_class"))
            b_class = normalize_label(row.get("b_class"))
            a_super = map_class_to_superclass(a_class, class_to_super)
            b_super = map_class_to_superclass(b_class, class_to_super)
            rows.append(
                {
                    "a_super": a_super,
                    "b_super": b_super,
                    "overlap_a_in_b": safe_float(row.get("overlap_a_in_b")),
                    "overlap_b_in_a": safe_float(row.get("overlap_b_in_a")),
                }
            )
    return rows


def keep_pair(row: Dict[str, object]) -> bool:
    oa = row.get("overlap_a_in_b")
    ob = row.get("overlap_b_in_a")
    oa_num = oa if isinstance(oa, float) else 0.0
    ob_num = ob if isinstance(ob, float) else 0.0
    return (oa_num > 0) or (ob_num > 0)


def pair_type_counts(total_rows: Iterable[Dict[str, object]], pair_order: List[PairKey]) -> Dict[PairKey, int]:
    counts = {k: 0 for k in pair_order}
    for row in total_rows:
        key = (str(row.get("a_super", "")), str(row.get("b_super", "")))
        if key in counts:
            counts[key] += 1
    return counts


def collect_directional_values(total_rows: Iterable[Dict[str, object]], pair_name: PairKey) -> List[float]:
    values: List[float] = []
    left_type, right_type = pair_name

    for row in total_rows:
        ab = (str(row.get("a_super", "")), str(row.get("b_super", "")))
        oa = float(row.get("overlap_a_in_b", 0.0) or 0.0)
        ob = float(row.get("overlap_b_in_a", 0.0) or 0.0)

        if left_type == right_type:
            if ab == pair_name:
                values.append(oa)
                values.append(ob)
            continue

        reverse_name = (right_type, left_type)
        if ab == pair_name:
            values.append(oa)
        elif ab == reverse_name:
            values.append(ob)
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


def draw_plot(dist_by_pair: Dict[PairKey, Dict[str, float]], pair_order: List[PairKey], plot_path: Path) -> None:
    pair_labels = [f"{a}-{b}" for a, b in pair_order]
    colors = ["#d9d2ad", "#b7d4ca", "#57a9a5", "#3d80ad", "#294f72"]

    width = max(1600, 200 * len(pair_order) + 260)
    height = 620
    left, right, top, bottom = 90, 40, 100, 90
    chart_w = width - left - right
    chart_h = height - top - bottom
    y_max = 100.0

    def y_to_px(v: float) -> float:
        return top + chart_h - (v / y_max) * chart_h

    svg: List[str] = []
    svg.append(f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}">')
    svg.append('<rect width="100%" height="100%" fill="white"/>')
    svg.append('<text x="600" y="40" text-anchor="middle" font-size="28" font-weight="700">Sub Cluster Cell Overlap Percent</text>')

    svg.append(f'<line x1="{left}" y1="{top}" x2="{left}" y2="{top + chart_h}" stroke="#333" stroke-width="2"/>')
    svg.append(f'<line x1="{left}" y1="{top + chart_h}" x2="{left + chart_w}" y2="{top + chart_h}" stroke="#333" stroke-width="2"/>')
    for tick in [0, 20, 40, 60, 80, 100]:
        y = y_to_px(float(tick))
        svg.append(f'<line x1="{left - 6}" y1="{y:.2f}" x2="{left}" y2="{y:.2f}" stroke="#333" stroke-width="2"/>')
        svg.append(f'<text x="{left - 14}" y="{y + 5:.2f}" text-anchor="end" font-size="18" fill="#333">{tick}%</text>')

    legend_start_x, legend_y = 320, 64
    for i, b in enumerate(BINS):
        lx = legend_start_x + i * 150
        svg.append(f'<rect x="{lx}" y="{legend_y - 14}" width="28" height="16" fill="{colors[i]}" stroke="#333"/>')
        svg.append(f'<text x="{lx + 38}" y="{legend_y}" font-size="18" fill="#333">{b}</text>')

    group_w = chart_w / len(pair_order)
    bar_w = 30
    for g, pair in enumerate(pair_order):
        center = left + group_w * (g + 0.5)
        for i, b in enumerate(BINS):
            v = dist_by_pair[pair][b]
            x = center - (2.5 * bar_w) + i * (bar_w + 3)
            y = y_to_px(v)
            h = top + chart_h - y
            svg.append(
                f'<rect x="{x:.2f}" y="{y:.2f}" width="{bar_w}" height="{h:.2f}" '
                f'fill="{colors[i]}" stroke="#333" stroke-width="1.2"/>'
            )
            if v > 0:
                svg.append(
                    f'<text x="{x + bar_w/2:.2f}" y="{max(18, y - 6):.2f}" text-anchor="middle" '
                    f'font-size="16" font-weight="700" fill="#222">{v:.1f}%</text>'
                )
        svg.append(f'<text x="{center:.2f}" y="{top + chart_h + 40}" text-anchor="middle" font-size="24" fill="#333">{pair_labels[g]}</text>')

    svg.append('</svg>')

    plot_path.parent.mkdir(parents=True, exist_ok=True)
    out_svg = plot_path if plot_path.suffix.lower() == ".svg" else plot_path.with_suffix(".svg")
    out_svg.write_text("\n".join(svg), encoding="utf-8")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Build superclass overlap distribution from existing true-data total table.")
    parser.add_argument("--total-table", type=Path, default=Path("true_data_overlap_results/true_data_total_overlap_table.csv"))
    parser.add_argument("--mapping", type=Path, default=Path("Merfish_brain_cell_type_subclass.txt"))
    parser.add_argument("--output-dir", type=Path, default=Path("true_data_overlap_results"))
    parser.add_argument("--pair-percent-table", type=str, default="true_data_pair_overlap_percent_superclass.csv")
    parser.add_argument("--plot", type=str, default="true_data_superclass_overlap_percent.svg")
    parser.add_argument("--no-plot", action="store_true", help="仅输出 CSV，不绘制 SVG。")
    return parser.parse_args()


def main() -> None:
    args = parse_args()

    class_to_super, unique_supers = load_superclass_mapping(args.mapping)
    pair_order: List[PairKey] = [(a, b) for a in unique_supers for b in unique_supers]

    total_table = load_total_table(args.total_table, class_to_super)
    filtered_rows = [r for r in total_table if keep_pair(r)]

    type_counts = pair_type_counts(filtered_rows, pair_order)
    dist_by_pair: Dict[PairKey, Dict[str, float]] = {}
    out_pair_rows: List[Dict[str, object]] = []
    for pair_name in pair_order:
        values = collect_directional_values(filtered_rows, pair_name)
        denominator = denominator_for_pair(pair_name, type_counts)
        dist = compute_distribution(values, denominator)
        dist_by_pair[pair_name] = dist
        pair_label = f"{pair_name[0]}-{pair_name[1]}"
        reverse_name = (pair_name[1], pair_name[0])
        for b in BINS:
            out_pair_rows.append(
                {
                    "pair_type": pair_label,
                    "overlap_bin": b,
                    "percent": dist[b],
                    "directional_values_in_bin_denominator": denominator,
                    "directional_values_count": len(values),
                    "source_pair_count": type_counts.get(pair_name, 0),
                    "mixed_pair_total": type_counts.get(pair_name, 0) + type_counts.get(reverse_name, 0),
                }
            )

    out_pair = args.output_dir / args.pair_percent_table
    out_pair.parent.mkdir(parents=True, exist_ok=True)
    with out_pair.open("w", encoding="utf-8", newline="") as fh:
        writer = csv.DictWriter(
            fh,
            fieldnames=[
                "pair_type",
                "overlap_bin",
                "percent",
                "directional_values_in_bin_denominator",
                "directional_values_count",
                "source_pair_count",
                "mixed_pair_total",
            ],
        )
        writer.writeheader()
        writer.writerows(out_pair_rows)

    out_plot = args.output_dir / args.plot
    if not args.no_plot:
        draw_plot(dist_by_pair, pair_order, out_plot)

    print(f"Input total pairs: {len(total_table)} from {args.total_table}")
    print(f"Pairs kept (overlap_a_in_b > 0 OR overlap_b_in_a > 0): {len(filtered_rows)}")
    print(f"Unique superclasses: {len(unique_supers)}")
    print(f"Superclass pair types: {len(pair_order)}")
    print(f"Pair distribution table -> {out_pair}")
    if not args.no_plot:
        print(f"Plot -> {out_plot}")


if __name__ == "__main__":
    main()
