#!/usr/bin/env python3
"""Compute overlap summary table for true_data and draw pair-overlap percent plot.

Workflow:
1) Reproduce compute_data1_overlap.py logic on files in true_data/.
2) Save the combined pair table ("总表").
3) Build overlap-percentage distributions for ordered neuron-type pairs and draw a grouped bar plot.
"""

from __future__ import annotations

import argparse
import csv
import itertools
from pathlib import Path
from typing import Dict, Iterable, List, Tuple



def safe_int(value: str | None) -> int | None:
    if value is None or value == "":
        return None
    try:
        return int(float(value))
    except ValueError:
        return None


def safe_float(value: str | None) -> float | None:
    if value is None or value == "":
        return None
    try:
        return float(value)
    except ValueError:
        return None


def format_scientific_id(token: str) -> str | None:
    token = token.strip()
    if not token:
        return None
    try:
        return format(float(token), ".10e")
    except ValueError:
        return None


def extract_ids(text: str | None) -> List[str]:
    if not text:
        return []
    out: List[str] = []
    for token in text.split(","):
        normalized = format_scientific_id(token)
        if normalized is not None:
            out.append(normalized)
    return out


def get_all_ids(row: Dict[str, str]) -> List[str]:
    ids: List[str] = []
    ids.extend(extract_ids(row.get("Glut_Neruon_cell_ids")))
    ids.extend(extract_ids(row.get("GABA_Neruon_cell_ids")))
    ids.extend(extract_ids(row.get("Non_Neruon_cell_ids")))
    return ids


def class_value(row: Dict[str, str]) -> str:
    return (row.get("class") or row.get("Class") or row.get("subclass") or "").strip()


def read_rows(true_data_dir: Path) -> List[Dict[str, str]]:
    rows: List[Dict[str, str]] = []
    for path in sorted(true_data_dir.glob("*_merged_regions_table_cell_id.txt")):
        with path.open("r", encoding="utf-8", newline="") as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            for row in reader:
                row["__source_file"] = path.name
                rows.append(row)
    return rows


def calc_ei(glu: str | None, gaba: str | None) -> float | None:
    g = safe_float(glu)
    b = safe_float(gaba)
    if g is None or b is None or b == 0:
        return None
    return g / b


def load_type_mapping(path: Path) -> Dict[str, str]:
    mapping: Dict[str, str] = {}
    with path.open("r", encoding="utf-8", newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            cls = (row.get("class") or "").strip()
            neuron = (row.get("cell_Neuron_type") or "").strip()
            if cls and neuron:
                mapping[cls] = neuron
    return mapping


def normalize_type(value: str | None) -> str:
    token = (value or "").strip()
    if not token:
        return "Unknown"
    low = token.lower()
    if low == "gaba":
        return "Gaba"
    if low == "glut":
        return "Glut"
    if low == "nonneuron":
        return "NonNeuron"
    return token


def fill_type(raw_type: str | None, cls: str, mapping: Dict[str, str]) -> str:
    ntype = normalize_type(raw_type)
    if ntype != "Unknown":
        return ntype
    return normalize_type(mapping.get(cls, "Unknown"))


def compute_total_table(rows: List[Dict[str, str]], mapping: Dict[str, str]) -> List[Dict[str, object]]:
    groups: Dict[Tuple[str, str], List[Dict[str, str]]] = {}
    for row in rows:
        key = (row.get("slide", ""), row.get("layer", ""))
        groups.setdefault(key, []).append(row)

    out_rows: List[Dict[str, object]] = []
    for (slide, layer), group_rows in groups.items():
        if len(group_rows) < 2:
            continue

        for row_a, row_b in itertools.combinations(group_rows, 2):
            ids_a = get_all_ids(row_a)
            ids_b = get_all_ids(row_b)
            n_a = len(ids_a)
            n_b = len(ids_b)
            if n_a == 0 or n_b == 0:
                continue

            b_set = set(ids_b)
            a_set = set(ids_a)
            overlap_a_in_b = sum(i in b_set for i in ids_a) / n_a
            overlap_b_in_a = sum(i in a_set for i in ids_b) / n_b

            a_cls = class_value(row_a)
            b_cls = class_value(row_b)
            a_type = fill_type(row_a.get("cell_Neuron_type"), a_cls, mapping)
            b_type = fill_type(row_b.get("cell_Neuron_type"), b_cls, mapping)

            out_rows.append(
                {
                    "slide": slide,
                    "layer": layer,
                    "a_layer": row_a.get("layer", ""),
                    "b_layer": row_b.get("layer", ""),
                    "a_class": a_cls,
                    "b_class": b_cls,
                    "a_glutotal": safe_int(row_a.get("Glut_Neruon_cell_ids_num")),
                    "a_gabatotal": safe_int(row_a.get("GABA_Neruon_cell_ids_num")),
                    "b_glutotal": safe_int(row_b.get("Glut_Neruon_cell_ids_num")),
                    "b_gabatotal": safe_int(row_b.get("GABA_Neruon_cell_ids_num")),
                    "a_region": row_a.get("region", ""),
                    "b_region": row_b.get("region", ""),
                    "a_ei": calc_ei(row_a.get("Glut_Neruon_cell_ids_num"), row_a.get("GABA_Neruon_cell_ids_num")),
                    "b_ei": calc_ei(row_b.get("Glut_Neruon_cell_ids_num"), row_b.get("GABA_Neruon_cell_ids_num")),
                    "a_merge_region": row_a.get("merge_regions", ""),
                    "b_merge_region": row_b.get("merge_regions", ""),
                    "a_cell_Neuron_type": a_type,
                    "b_cell_Neuron_type": b_type,
                    "a_enrich_class_cell_ids_num": safe_int(
                        row_a.get("enrich_class_cell_ids_num") or row_a.get("enrich_subclass_cell_ids_num")
                    ),
                    "b_enrich_class_cell_ids_num": safe_int(
                        row_b.get("enrich_class_cell_ids_num") or row_b.get("enrich_subclass_cell_ids_num")
                    ),
                    "a_cell_id_n": n_a,
                    "b_cell_id_n": n_b,
                    "overlap_a_in_b": overlap_a_in_b,
                    "overlap_b_in_a": overlap_b_in_a,
                    "a_source_file": row_a.get("__source_file", ""),
                    "b_source_file": row_b.get("__source_file", ""),
                }
            )
    return out_rows


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


def pair_type_counts(total_rows: Iterable[Dict[str, object]]) -> Dict[str, int]:
    counts = {"Gaba-Gaba": 0, "Gaba-Glut": 0, "Glut-Gaba": 0, "Glut-Glut": 0}
    for row in total_rows:
        a_type = str(row.get("a_cell_Neuron_type", ""))
        b_type = str(row.get("b_cell_Neuron_type", ""))
        key = f"{a_type}-{b_type}"
        if key in counts:
            counts[key] += 1
    return counts


def collect_directional_values(total_rows: Iterable[Dict[str, object]], pair_name: str) -> List[float]:
    """Collect overlap values in the requested direction.

    Rules:
    - Gaba-Glut: use Gaba->Glut direction (a->b for Gaba-Glut rows, b->a for Glut-Gaba rows)
    - Glut-Gaba: use Glut->Gaba direction (a->b for Glut-Gaba rows, b->a for Gaba-Glut rows)
    - Gaba-Gaba / Glut-Glut: both directions are valid, so keep overlap_a_in_b and overlap_b_in_a.
    """
    values: List[float] = []
    for row in total_rows:
        a_type = str(row.get("a_cell_Neuron_type", ""))
        b_type = str(row.get("b_cell_Neuron_type", ""))
        ab = f"{a_type}-{b_type}"

        if pair_name in {"Gaba-Gaba", "Glut-Glut"} and ab == pair_name:
            values.append(float(row["overlap_a_in_b"]))
            values.append(float(row["overlap_b_in_a"]))
            continue

        if pair_name == "Gaba-Glut":
            if ab == "Gaba-Glut":
                values.append(float(row["overlap_a_in_b"]))
            elif ab == "Glut-Gaba":
                values.append(float(row["overlap_b_in_a"]))
        elif pair_name == "Glut-Gaba":
            if ab == "Glut-Gaba":
                values.append(float(row["overlap_a_in_b"]))
            elif ab == "Gaba-Glut":
                values.append(float(row["overlap_b_in_a"]))
    return values


def denominator_for_pair(pair_name: str, counts: Dict[str, int]) -> int:
    if pair_name == "Gaba-Gaba":
        return 2 * counts["Gaba-Gaba"]
    if pair_name == "Glut-Glut":
        return 2 * counts["Glut-Glut"]
    return counts["Gaba-Glut"] + counts["Glut-Gaba"]


def compute_distribution(values: List[float], bins: List[str], denominator: int) -> Dict[str, float]:
    counts = {b: 0 for b in bins}
    if denominator <= 0:
        return {b: 0.0 for b in bins}

    for v in values:
        counts[overlap_bin(v)] += 1

    return {b: counts[b] * 100.0 / denominator for b in bins}


def draw_plot(dist_by_pair: Dict[str, Dict[str, float]], plot_path: Path) -> None:
    """Draw grouped bar chart to SVG without third-party plotting dependencies."""
    pair_order = ["Gaba-Gaba", "Gaba-Glut", "Glut-Gaba", "Glut-Glut"]
    pair_labels = ["GABA-GABA", "GABA-Glut", "Glut-GABA", "Glut-Glut"]
    bins = ["0–20%", "20–40%", "40–60%", "60–80%", "80–100%"]
    colors = ["#d9d2ad", "#b7d4ca", "#57a9a5", "#3d80ad", "#294f72"]

    width, height = 1200, 620
    left, right, top, bottom = 90, 40, 100, 90
    chart_w = width - left - right
    chart_h = height - top - bottom
    y_max = 55.0

    def y_to_px(v: float) -> float:
        return top + chart_h - (v / y_max) * chart_h

    svg: List[str] = []
    svg.append(f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}">')
    svg.append('<rect width="100%" height="100%" fill="white"/>')
    svg.append('<text x="600" y="40" text-anchor="middle" font-size="28" font-weight="700">Sub Cluster Cell Overlap Percent</text>')

    # Axes and y ticks
    svg.append(f'<line x1="{left}" y1="{top}" x2="{left}" y2="{top + chart_h}" stroke="#333" stroke-width="2"/>')
    svg.append(f'<line x1="{left}" y1="{top + chart_h}" x2="{left + chart_w}" y2="{top + chart_h}" stroke="#333" stroke-width="2"/>')
    for tick in [0, 10, 20, 30, 40, 50]:
        y = y_to_px(float(tick))
        svg.append(f'<line x1="{left - 6}" y1="{y:.2f}" x2="{left}" y2="{y:.2f}" stroke="#333" stroke-width="2"/>')
        svg.append(f'<text x="{left - 14}" y="{y + 5:.2f}" text-anchor="end" font-size="18" fill="#333">{tick}%</text>')

    # Legend
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


def save_total_table(rows: List[Dict[str, object]], out_path: Path) -> None:
    fieldnames = [
        "slide",
        "layer",
        "a_layer",
        "b_layer",
        "a_class",
        "b_class",
        "a_glutotal",
        "a_gabatotal",
        "b_glutotal",
        "b_gabatotal",
        "a_region",
        "b_region",
        "a_ei",
        "b_ei",
        "a_merge_region",
        "b_merge_region",
        "a_cell_Neuron_type",
        "b_cell_Neuron_type",
        "a_enrich_class_cell_ids_num",
        "b_enrich_class_cell_ids_num",
        "a_cell_id_n",
        "b_cell_id_n",
        "overlap_a_in_b",
        "overlap_b_in_a",
        "a_source_file",
        "b_source_file",
    ]
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w", encoding="utf-8", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Compute true_data overlap total table and pair-overlap plot.")
    parser.add_argument("--data-dir", type=Path, default=Path("true_data"))
    parser.add_argument("--mapping", type=Path, default=Path("Merfish_brain_cell_type_subclass.txt"))
    parser.add_argument("--output-dir", type=Path, default=Path("true_data_overlap_results"))
    parser.add_argument("--total-table", type=str, default="true_data_total_overlap_table.csv")
    parser.add_argument("--pair-percent-table", type=str, default="true_data_pair_overlap_percent.csv")
    parser.add_argument("--plot", type=str, default="true_data_subcluster_overlap_percent.svg")
    return parser.parse_args()


def main() -> None:
    args = parse_args()

    rows = read_rows(args.data_dir)
    if not rows:
        raise SystemExit(f"No *_merged_regions_table_cell_id.txt files found in {args.data_dir}")

    mapping = load_type_mapping(args.mapping)
    total_table = compute_total_table(rows, mapping)

    out_total = args.output_dir / args.total_table
    save_total_table(total_table, out_total)

    bins = ["0–20%", "20–40%", "40–60%", "60–80%", "80–100%"]
    pair_order = ["Gaba-Gaba", "Gaba-Glut", "Glut-Gaba", "Glut-Glut"]

    type_counts = pair_type_counts(total_table)
    dist_by_pair: Dict[str, Dict[str, float]] = {}
    out_pair_rows: List[Dict[str, object]] = []
    for pair_name in pair_order:
        values = collect_directional_values(total_table, pair_name)
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
                    "mixed_pair_total": type_counts["Gaba-Glut"] + type_counts["Glut-Gaba"],
                }
            )

    out_pair = args.output_dir / args.pair_percent_table
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

    print(f"Total pairs: {len(total_table)} -> {out_total}")
    print(f"Pair distribution table -> {out_pair}")
    print(f"Plot -> {out_plot}")


if __name__ == "__main__":
    main()
