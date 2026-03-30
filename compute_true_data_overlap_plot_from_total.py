#!/usr/bin/env python3
"""Build pair-overlap statistics/plot from an existing total-overlap table.

This script skips total-table computation and only runs downstream analysis.
Extra AB-pair rule:
- Keep a pair only when overlap_a_in_b > 0 OR overlap_b_in_a > 0.
- If both are 0 (or invalid), the pair is excluded from statistics.
"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path
from typing import Dict, Iterable, List, Tuple


CELL_TYPES = ["Gaba", "Glut", "NonNeuron"]


def safe_float(value: str | None) -> float | None:
    if value is None or value == "":
        return None
    try:
        return float(value)
    except ValueError:
        return None


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


def load_total_table(path: Path) -> List[Dict[str, object]]:
    rows: List[Dict[str, object]] = []
    with path.open("r", encoding="utf-8", newline="") as fh:
        reader = csv.DictReader(fh)
        for row in reader:
            rows.append(
                {
                    "a_cell_Neuron_type": (row.get("a_cell_Neuron_type") or "").strip(),
                    "b_cell_Neuron_type": (row.get("b_cell_Neuron_type") or "").strip(),
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


def pair_type_counts(total_rows: Iterable[Dict[str, object]], pair_order: List[str]) -> Dict[str, int]:
    counts = {k: 0 for k in pair_order}
    for row in total_rows:
        a_type = str(row.get("a_cell_Neuron_type", ""))
        b_type = str(row.get("b_cell_Neuron_type", ""))
        key = f"{a_type}-{b_type}"
        if key in counts:
            counts[key] += 1
    return counts


def split_pair(pair_name: str) -> Tuple[str, str]:
    parts = pair_name.split("-", 1)
    if len(parts) != 2:
        return "", ""
    return parts[0], parts[1]


def collect_directional_values(total_rows: Iterable[Dict[str, object]], pair_name: str) -> List[float]:
    values: List[float] = []
    left_type, right_type = split_pair(pair_name)
    reverse_name = f"{right_type}-{left_type}"

    for row in total_rows:
        a_type = str(row.get("a_cell_Neuron_type", ""))
        b_type = str(row.get("b_cell_Neuron_type", ""))
        ab = f"{a_type}-{b_type}"
        oa = float(row.get("overlap_a_in_b", 0.0) or 0.0)
        ob = float(row.get("overlap_b_in_a", 0.0) or 0.0)

        if left_type == right_type:
            if ab == pair_name:
                values.append(oa)
                values.append(ob)
            continue

        if ab == pair_name:
            values.append(oa)
        elif ab == reverse_name:
            values.append(ob)
    return values


def denominator_for_pair(pair_name: str, counts: Dict[str, int]) -> int:
    left_type, right_type = split_pair(pair_name)
    reverse_name = f"{right_type}-{left_type}"
    if left_type == right_type:
        return 2 * counts.get(pair_name, 0)
    return counts.get(pair_name, 0) + counts.get(reverse_name, 0)


def compute_distribution(values: List[float], bins: List[str], denominator: int) -> Dict[str, float]:
    counts = {b: 0 for b in bins}
    if denominator <= 0:
        return {b: 0.0 for b in bins}

    for v in values:
        counts[overlap_bin(v)] += 1

    return {b: counts[b] * 100.0 / denominator for b in bins}


def draw_plot(dist_by_pair: Dict[str, Dict[str, float]], plot_path: Path) -> None:
    pair_order = [
        "Gaba-Gaba",
        "Gaba-Glut",
        "Gaba-NonNeuron",
        "Glut-Gaba",
        "Glut-Glut",
        "Glut-NonNeuron",
        "NonNeuron-Gaba",
        "NonNeuron-Glut",
        "NonNeuron-NonNeuron",
    ]
    pair_labels = [
        "GABA-GABA",
        "GABA-Glut",
        "GABA-NonNeuron",
        "Glut-GABA",
        "Glut-Glut",
        "Glut-NonNeuron",
        "NonNeuron-GABA",
        "NonNeuron-Glut",
        "NonNeuron-NonNeuron",
    ]
    bins = ["0–20%", "20–40%", "40–60%", "60–80%", "80–100%"]
    colors = ["#d9d2ad", "#b7d4ca", "#57a9a5", "#3d80ad", "#294f72"]

    width, height = 2100, 620
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
    for i, b in enumerate(bins):
        lx = legend_start_x + i * 150
        svg.append(f'<rect x="{lx}" y="{legend_y - 14}" width="28" height="16" fill="{colors[i]}" stroke="#333"/>')
        svg.append(f'<text x="{lx + 38}" y="{legend_y}" font-size="18" fill="#333">{b}</text>')

    group_w = chart_w / len(pair_order)
    bar_w = 30
    for g, pair in enumerate(pair_order):
        center = left + group_w * (g + 0.5)
        for i, b in enumerate(bins):
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
    parser = argparse.ArgumentParser(description="Build overlap plot from existing true-data total table.")
    parser.add_argument("--total-table", type=Path, default=Path("true_data_overlap_results/true_data_total_overlap_table.csv"))
    parser.add_argument("--output-dir", type=Path, default=Path("true_data_overlap_results"))
    parser.add_argument("--pair-percent-table", type=str, default="true_data_pair_overlap_percent_nonzero_pair.csv")
    parser.add_argument("--plot", type=str, default="true_data_subcluster_overlap_percent_nonzero_pair.svg")
    return parser.parse_args()


def main() -> None:
    args = parse_args()

    total_table = load_total_table(args.total_table)
    filtered_rows = [r for r in total_table if keep_pair(r)]

    bins = ["0–20%", "20–40%", "40–60%", "60–80%", "80–100%"]
    pair_order = [f"{a}-{b}" for a in CELL_TYPES for b in CELL_TYPES]

    type_counts = pair_type_counts(filtered_rows, pair_order)
    dist_by_pair: Dict[str, Dict[str, float]] = {}
    out_pair_rows: List[Dict[str, object]] = []
    for pair_name in pair_order:
        values = collect_directional_values(filtered_rows, pair_name)
        denominator = denominator_for_pair(pair_name, type_counts)
        dist = compute_distribution(values, bins, denominator)
        dist_by_pair[pair_name] = dist
        for b in bins:
            out_pair_rows.append(
                {
                    "pair_type": pair_name,
                    "overlap_bin": b,
                    "percent": dist[b],
                    "directional_values_in_bin_denominator": denominator,
                    "directional_values_count": len(values),
                    "source_pair_count": type_counts.get(pair_name, 0),
                    "mixed_pair_total": type_counts.get(pair_name, 0) + type_counts.get(f"{pair_name.split('-', 1)[1]}-{pair_name.split('-', 1)[0]}", 0),
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
    draw_plot(dist_by_pair, out_plot)

    print(f"Input total pairs: {len(total_table)} from {args.total_table}")
    print(f"Pairs kept (overlap_a_in_b > 0 OR overlap_b_in_a > 0): {len(filtered_rows)}")
    print(f"Pair distribution table -> {out_pair}")
    print(f"Plot -> {out_plot}")


if __name__ == "__main__":
    main()
