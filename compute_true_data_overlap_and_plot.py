#!/usr/bin/env python3
"""Compute overlap summary from true_data and plot overlap-percent distributions.

Workflow:
1. Read all `*_merged_regions_table_cell_id.txt` under `true_data/`.
2. Build a "total table" with pairwise overlaps (same overlap logic as
   `compute_data1_overlap.py`).
3. Fill `a_cell_Neuron_type` / `b_cell_Neuron_type` from source columns when
   available; otherwise map from `Merfish_brain_cell_type_subclass.txt` by class.
4. Compute percentage distributions by overlap bins for pair relations
   (GABA-GABA, GABA-Glut, Glut-GABA, Glut-Glut) and draw grouped bars.
"""

from __future__ import annotations

import argparse
import csv
import itertools
from collections import defaultdict
from pathlib import Path
from typing import Dict, Iterable, List, Tuple



PAIR_ORDER = ["GABA-GABA", "GABA-Glut", "Glut-GABA", "Glut-Glut"]
BIN_LABELS = ["0–20%", "20–40%", "40–60%", "60–80%", "80–100%"]
BIN_EDGES = [0.0, 0.2, 0.4, 0.6, 0.8, 1.000000001]  # include 1.0
BIN_COLORS = ["#d8d1ad", "#b8d5cc", "#5ca7a1", "#377cae", "#274f71"]


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
        norm = format_scientific_id(token)
        if norm is not None:
            out.append(norm)
    return out


def get_all_ids(row: Dict[str, str]) -> List[str]:
    ids: List[str] = []
    ids.extend(extract_ids(row.get("Glut_Neruon_cell_ids")))
    ids.extend(extract_ids(row.get("GABA_Neruon_cell_ids")))
    ids.extend(extract_ids(row.get("Non_Neruon_cell_ids")))
    return ids


def calc_ei(glu: str | None, gaba: str | None) -> float | None:
    g = safe_float(glu)
    b = safe_float(gaba)
    if g is None or b is None or b == 0:
        return None
    return g / b


def class_value(row: Dict[str, str]) -> str:
    return (row.get("class") or row.get("Class") or row.get("subclass") or "").strip()


def normalize_type(value: str | None) -> str:
    if not value:
        return ""
    v = value.strip().lower()
    if v == "gaba":
        return "GABA"
    if v == "glut":
        return "Glut"
    return value.strip()


def load_class_type_map(path: Path) -> Dict[str, str]:
    mapping: Dict[str, str] = {}
    with path.open("r", encoding="utf-8", newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            cls = (row.get("class") or "").strip()
            ctype = normalize_type(row.get("cell_Neuron_type"))
            if cls and ctype:
                mapping[cls] = ctype
    return mapping


def read_true_rows(true_data_dir: Path) -> List[Dict[str, str]]:
    rows: List[Dict[str, str]] = []
    for path in sorted(true_data_dir.glob("*_merged_regions_table_cell_id.txt")):
        with path.open("r", encoding="utf-8", newline="") as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            for row in reader:
                row["__source_file"] = path.name
                rows.append(row)
    return rows


def fill_cell_type(row: Dict[str, str], which: str, class_type_map: Dict[str, str]) -> str:
    direct_keys = [f"{which}_cell_Neuron_type", f"{which}_cell_neurontype", "cell_Neuron_type", "cell_neurontype"]
    for key in direct_keys:
        if key in row:
            val = normalize_type(row.get(key))
            if val:
                return val

    c = class_value(row)
    return normalize_type(class_type_map.get(c, ""))


def build_total_table(rows: List[Dict[str, str]], class_type_map: Dict[str, str]) -> List[Dict[str, object]]:
    groups: Dict[Tuple[str, str], List[Dict[str, str]]] = defaultdict(list)
    for row in rows:
        key = (row.get("slide", ""), row.get("layer", ""))
        groups[key].append(row)

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

            a_type = fill_cell_type(row_a, "a", class_type_map)
            b_type = fill_cell_type(row_b, "b", class_type_map)

            out_rows.append(
                {
                    "slide": slide,
                    "layer": layer,
                    "a_layer": row_a.get("layer", ""),
                    "b_layer": row_b.get("layer", ""),
                    "a_class": class_value(row_a),
                    "b_class": class_value(row_b),
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


def overlap_bin_idx(value: float) -> int:
    for i in range(len(BIN_EDGES) - 1):
        if BIN_EDGES[i] <= value < BIN_EDGES[i + 1]:
            return i
    return len(BIN_LABELS) - 1


def collect_directed_overlap_values(total_rows: Iterable[Dict[str, object]], x_type: str, y_type: str) -> List[float]:
    vals: List[float] = []
    for row in total_rows:
        a = normalize_type(str(row.get("a_cell_Neuron_type", "")))
        b = normalize_type(str(row.get("b_cell_Neuron_type", "")))

        if a == x_type and b == y_type:
            vals.append(float(row["overlap_a_in_b"]))
        if a == y_type and b == x_type:
            vals.append(float(row["overlap_b_in_a"]))
    return vals


def compute_pair_bin_percent(total_rows: List[Dict[str, object]]) -> Dict[str, List[float]]:
    pair_defs = {
        "GABA-GABA": ("GABA", "GABA"),
        "GABA-Glut": ("GABA", "Glut"),
        "Glut-GABA": ("Glut", "GABA"),
        "Glut-Glut": ("Glut", "Glut"),
    }

    result: Dict[str, List[float]] = {}
    for pair_name in PAIR_ORDER:
        x_type, y_type = pair_defs[pair_name]
        values = collect_directed_overlap_values(total_rows, x_type, y_type)
        counts = [0] * len(BIN_LABELS)
        for v in values:
            counts[overlap_bin_idx(v)] += 1
        total = len(values)
        if total == 0:
            result[pair_name] = [0.0] * len(BIN_LABELS)
        else:
            result[pair_name] = [c / total * 100 for c in counts]
    return result


def write_csv(path: Path, rows: List[Dict[str, object]]) -> None:
    if not rows:
        path.write_text("", encoding="utf-8")
        return
    fieldnames = list(rows[0].keys())
    with path.open("w", encoding="utf-8", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def write_pair_percent_csv(path: Path, pair_pct: Dict[str, List[float]]) -> None:
    with path.open("w", encoding="utf-8", newline="") as fh:
        writer = csv.writer(fh)
        writer.writerow(["pair_type", *BIN_LABELS])
        for pair in PAIR_ORDER:
            writer.writerow([pair, *[f"{v:.6f}" for v in pair_pct[pair]]])


def make_plot(pair_pct: Dict[str, List[float]], output_path: Path) -> None:
    """Write a lightweight grouped-bar SVG using only stdlib."""
    w, h = 1100, 620
    margin = {"left": 85, "right": 35, "top": 95, "bottom": 95}
    plot_w = w - margin["left"] - margin["right"]
    plot_h = h - margin["top"] - margin["bottom"]

    max_val = max(max(v) for v in pair_pct.values()) if pair_pct else 1.0
    y_max = max(60.0, (int(max_val / 10) + 1) * 10.0)

    def x_pos(pair_i: int, bin_i: int) -> tuple[float, float]:
        group_w = plot_w / len(PAIR_ORDER)
        bar_w = 28
        gap = 4
        total = len(BIN_LABELS) * bar_w + (len(BIN_LABELS) - 1) * gap
        x0 = margin["left"] + pair_i * group_w + (group_w - total) / 2 + bin_i * (bar_w + gap)
        return x0, bar_w

    def y_pos(value: float) -> float:
        return margin["top"] + plot_h * (1 - value / y_max)

    parts: List[str] = []
    parts.append(f'<svg xmlns="http://www.w3.org/2000/svg" width="{w}" height="{h}" viewBox="0 0 {w} {h}">')
    parts.append('<rect x="0" y="0" width="100%" height="100%" fill="#ffffff"/>')

    # Title
    parts.append('<text x="550" y="38" text-anchor="middle" font-size="34" font-weight="700">Sub Cluster Cell Overlap Percent</text>')

    # Legend
    lx = 300
    for i, label in enumerate(BIN_LABELS):
        x = lx + i * 130
        parts.append(f'<rect x="{x}" y="58" width="24" height="14" fill="{BIN_COLORS[i]}" stroke="#333" stroke-width="1"/>')
        parts.append(f'<text x="{x+34}" y="70" font-size="19" fill="#333">{label}</text>')

    # Axes
    x_left, x_right = margin["left"], margin["left"] + plot_w
    y_top, y_bottom = margin["top"], margin["top"] + plot_h
    parts.append(f'<line x1="{x_left}" y1="{y_bottom}" x2="{x_right}" y2="{y_bottom}" stroke="#222" stroke-width="2"/>')
    parts.append(f'<line x1="{x_left}" y1="{y_top}" x2="{x_left}" y2="{y_bottom}" stroke="#222" stroke-width="2"/>')

    # Y ticks
    for tick in range(0, int(y_max) + 1, 10):
        y = y_pos(float(tick))
        parts.append(f'<line x1="{x_left-6}" y1="{y}" x2="{x_left}" y2="{y}" stroke="#222" stroke-width="1"/>')
        parts.append(f'<text x="{x_left-12}" y="{y+6}" text-anchor="end" font-size="22" fill="#444">{tick}%</text>')

    # Bars and value labels
    for p_i, pair in enumerate(PAIR_ORDER):
        vals = pair_pct[pair]
        for b_i, val in enumerate(vals):
            x, bw = x_pos(p_i, b_i)
            y = y_pos(val)
            bh = y_bottom - y
            parts.append(
                f'<rect x="{x:.1f}" y="{y:.1f}" width="{bw}" height="{bh:.1f}" fill="{BIN_COLORS[b_i]}" stroke="#2f2f2f" stroke-width="1"/>'
            )
            parts.append(
                f'<text x="{x + bw/2:.1f}" y="{max(y - 6, y_top + 12):.1f}" text-anchor="middle" font-size="18" font-weight="700">{val:.1f}%</text>'
            )

        # group label
        group_w = plot_w / len(PAIR_ORDER)
        cx = margin["left"] + p_i * group_w + group_w / 2
        parts.append(f'<text x="{cx:.1f}" y="{y_bottom+40}" text-anchor="middle" font-size="29" fill="#333">{pair}</text>')

    parts.append('</svg>')
    output_path.write_text("\n".join(parts), encoding="utf-8")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Compute true_data overlap and plot pair overlap distributions.")
    parser.add_argument("--true-data-dir", type=Path, default=Path("true_data"))
    parser.add_argument("--mapping-file", type=Path, default=Path("Merfish_brain_cell_type_subclass.txt"))
    parser.add_argument("--out-total-table", type=Path, default=Path("true_data_total_overlap_table.csv"))
    parser.add_argument("--out-pair-percent", type=Path, default=Path("true_data_pair_overlap_percent.csv"))
    parser.add_argument("--out-plot", type=Path, default=Path("true_data_pair_overlap_percent.svg"))
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    class_type_map = load_class_type_map(args.mapping_file)
    rows = read_true_rows(args.true_data_dir)
    if not rows:
        raise SystemExit(f"No *_merged_regions_table_cell_id.txt found in {args.true_data_dir}")

    total_rows = build_total_table(rows, class_type_map)
    write_csv(args.out_total_table, total_rows)
    pair_pct = compute_pair_bin_percent(total_rows)
    write_pair_percent_csv(args.out_pair_percent, pair_pct)
    make_plot(pair_pct, args.out_plot)

    print(f"Rows read: {len(rows)}")
    print(f"Total pair rows: {len(total_rows)} -> {args.out_total_table}")
    print(f"Pair percent table -> {args.out_pair_percent}")
    print(f"Plot -> {args.out_plot}")


if __name__ == "__main__":
    main()
