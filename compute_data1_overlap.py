#!/usr/bin/env python3
"""Compute pairwise region overlap summaries for each sample in data1/.

This script ports the core logic in extract_figure1.R to Python and runs each
sample in parallel (default: 20 workers).
"""

from __future__ import annotations

import argparse
import csv
import itertools
import math
import re
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
from typing import Dict, Iterable, List, Sequence

GABA_LIST = {
    "046 Vip Gaba",
    "047 Sncg Gaba",
    "048 RHP-COA Ndnf Gaba",
    "049 Lamp5 Gaba",
    "050 Lamp5 Lhx6 Gaba",
    "051 Pvalb chandelier Gaba",
    "052 Pvalb Gaba",
    "053 Sst Gaba",
}
GLUT_PATTERN = re.compile(r"^(00[1-9]|01[0-9]|02[0-6])")


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
    result: List[str] = []
    for token in text.split(","):
        normalized = format_scientific_id(token)
        if normalized is not None:
            result.append(normalized)
    return result


def get_all_ids(row: Dict[str, str]) -> List[str]:
    ids: List[str] = []
    ids.extend(extract_ids(row.get("Glut_Neruon_cell_ids")))
    ids.extend(extract_ids(row.get("GABA_Neruon_cell_ids")))
    # Optional column in some tables.
    ids.extend(extract_ids(row.get("Non_Neruon_cell_ids")))
    return ids


def read_sample_rows(sample_dir: Path) -> List[Dict[str, str]]:
    rows: List[Dict[str, str]] = []
    for path in sorted(sample_dir.glob("*_merged_regions_table_cell_id.txt")):
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


def class_value(row: Dict[str, str]) -> str:
    return (row.get("class") or row.get("subclass") or "").strip()


def build_target_classes(rows: Sequence[Dict[str, str]]) -> set[str]:
    classes = {class_value(r) for r in rows if class_value(r)}
    glut = {c for c in classes if GLUT_PATTERN.match(c)}
    return glut | GABA_LIST


def process_sample(sample_dir: Path, output_dir: Path) -> tuple[str, int, Path]:
    rows = read_sample_rows(sample_dir)
    target_classes = build_target_classes(rows)
    rows = [r for r in rows if class_value(r) in target_classes]

    groups: Dict[tuple[str, str], List[Dict[str, str]]] = {}
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

            # Keep same behavior as R script (%in% over vector; duplicates count).
            b_set = set(ids_b)
            a_set = set(ids_a)
            overlap_a_in_b = sum(i in b_set for i in ids_a) / n_a
            overlap_b_in_a = sum(i in a_set for i in ids_b) / n_b

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
                    "a_cell_Neuron_type": row_a.get("cell_Neuron_type", ""),
                    "b_cell_Neuron_type": row_b.get("cell_Neuron_type", ""),
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

    output_dir.mkdir(parents=True, exist_ok=True)
    out_file = output_dir / f"{sample_dir.name}_cellid_overlap_summary_filtered.csv"
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

    with out_file.open("w", encoding="utf-8", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(out_rows)

    return sample_dir.name, len(out_rows), out_file


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Compute per-sample overlap summaries for data1.")
    parser.add_argument("--data-dir", type=Path, default=Path("data1"), help="Folder containing sample_* subfolders")
    parser.add_argument("--output-dir", type=Path, default=Path("overlap_results"), help="Output folder")
    parser.add_argument("--workers", type=int, default=20, help="Parallel workers (default: 20)")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    sample_dirs = sorted(p for p in args.data_dir.glob("sample_*") if p.is_dir())
    if not sample_dirs:
        raise SystemExit(f"No sample_* directories found in {args.data_dir}")

    workers = max(1, args.workers)
    print(f"Found {len(sample_dirs)} samples, running with {workers} workers.")

    futures = []
    with ProcessPoolExecutor(max_workers=workers) as pool:
        for sample_dir in sample_dirs:
            futures.append(pool.submit(process_sample, sample_dir, args.output_dir))

        for fut in as_completed(futures):
            sample, n, out_file = fut.result()
            print(f"[{sample}] pairs={n} -> {out_file}")


if __name__ == "__main__":
    main()
