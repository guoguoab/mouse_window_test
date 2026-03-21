#!/usr/bin/env python3
"""按 sample 统计 overlap 分箱相关性，并按 slide+无向配对输出 overlap_a_in_b 记录。"""

from __future__ import annotations

import argparse
import csv
import math
from collections import defaultdict
from datetime import datetime
from multiprocessing import Pool
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Tuple

PAIR_ORDER = ["Gaba-Gaba", "Gaba-Glut", "Glut-Gaba", "Glut-Glut"]
BINS = ["0-20", "20-40", "40-60", "60-80", "80-100"]


def log(msg: str) -> None:
    ts = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    print(f"[{ts}] {msg}", flush=True)


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
    if "gaba" in low:
        return "Gaba"
    if "glut" in low:
        return "Glut"
    if "non" in low and "neuron" in low:
        return "NonNeuron"
    return token


def fill_type(raw_type: str | None, cls: str, mapping: Dict[str, str]) -> str:
    ntype = normalize_type(raw_type)
    if ntype != "Unknown":
        return ntype
    return normalize_type(mapping.get(cls, cls))


def bin_overlap(overlap: float) -> Optional[str]:
    """overlap 按百分比分箱: (0,20], (20,40], ... (80,100]。"""
    if overlap <= 0:
        return None
    if overlap <= 0.2:
        return "0-20"
    if overlap <= 0.4:
        return "20-40"
    if overlap <= 0.6:
        return "40-60"
    if overlap <= 0.8:
        return "60-80"
    return "80-100"


def safe_float(value: str | None) -> Optional[float]:
    if value is None:
        return None
    s = value.strip()
    if not s:
        return None
    try:
        return float(s)
    except ValueError:
        return None


def safe_log2(value: str | None) -> Optional[float]:
    x = safe_float(value)
    if x is None or x <= 0:
        return None
    return math.log2(x)


def pearson_r(xs: List[float], ys: List[float]) -> float:
    n = len(xs)
    if n < 2:
        return float("nan")
    mx = sum(xs) / n
    my = sum(ys) / n
    cov = 0.0
    vx = 0.0
    vy = 0.0
    for x, y in zip(xs, ys):
        dx = x - mx
        dy = y - my
        cov += dx * dy
        vx += dx * dx
        vy += dy * dy
    if vx <= 0 or vy <= 0:
        return float("nan")
    return cov / math.sqrt(vx * vy)


def canonical_class_pair(a_class: str, b_class: str) -> str:
    a = a_class.strip()
    b = b_class.strip()
    return " || ".join(sorted([a, b]))


def summarize_chunk(
    rows: List[Dict[str, str]],
    mapping: Dict[str, str],
) -> Tuple[
    Dict[Tuple[str, str], List[Tuple[Optional[float], Optional[float], Optional[float], Optional[float]]]],
    List[Tuple[str, str, str, float]],
]:
    corr_bucket: Dict[Tuple[str, str], List[Tuple[Optional[float], Optional[float], Optional[float], Optional[float]]]] = defaultdict(list)
    class_pair_records: List[Tuple[str, str, str, float]] = []

    for row in rows:
        a_class = (row.get("a_class") or "").strip()
        b_class = (row.get("b_class") or "").strip()
        a_type = fill_type(row.get("a_cell_Neuron_type"), a_class, mapping)
        b_type = fill_type(row.get("b_cell_Neuron_type"), b_class, mapping)

        pair = f"{a_type}-{b_type}"
        if pair not in PAIR_ORDER:
            continue

        a2b = safe_float(row.get("overlap_a_in_b"))
        b2a = safe_float(row.get("overlap_b_in_a"))
        if a2b is None or b2a is None:
            continue

        lx = safe_log2(row.get("a_cell_id_n"))
        ly = safe_log2(row.get("b_cell_id_n"))
        ex = safe_float(row.get("a_ei"))
        ey = safe_float(row.get("b_ei"))

        # 与原统计方式一致：同类配对使用双向；异类配对按方向归入对应配对。
        if pair in {"Gaba-Gaba", "Glut-Glut"}:
            for ov in (a2b, b2a):
                b = bin_overlap(ov)
                if b is not None:
                    corr_bucket[(pair, b)].append((lx, ly, ex, ey))
        elif pair == "Gaba-Glut":
            b = bin_overlap(a2b)
            if b is not None:
                corr_bucket[("Gaba-Glut", b)].append((lx, ly, ex, ey))
            b = bin_overlap(b2a)
            if b is not None:
                corr_bucket[("Glut-Gaba", b)].append((lx, ly, ex, ey))
        elif pair == "Glut-Gaba":
            b = bin_overlap(a2b)
            if b is not None:
                corr_bucket[("Glut-Gaba", b)].append((lx, ly, ex, ey))
            b = bin_overlap(b2a)
            if b is not None:
                corr_bucket[("Gaba-Glut", b)].append((lx, ly, ex, ey))

        if a2b > 0:
            slide = (row.get("slide") or "").strip()
            class_pair_records.append((slide, a_class, b_class, a2b))

    return corr_bucket, class_pair_records


def merge_corr_bucket(
    total: Dict[Tuple[str, str], List[Tuple[Optional[float], Optional[float], Optional[float], Optional[float]]]],
    part: Dict[Tuple[str, str], List[Tuple[Optional[float], Optional[float], Optional[float], Optional[float]]]],
) -> None:
    for k, vals in part.items():
        total[k].extend(vals)


def chunked_rows(path: Path, chunk_size: int) -> Iterable[List[Dict[str, str]]]:
    with path.open("r", encoding="utf-8", newline="") as fh:
        reader = csv.DictReader(fh)
        batch: List[Dict[str, str]] = []
        for row in reader:
            batch.append(row)
            if len(batch) >= chunk_size:
                yield batch
                batch = []
        if batch:
            yield batch


def process_file(path: Path, sample_name: str, mapping: Dict[str, str], n_workers: int, chunk_size: int):
    merged_corr: Dict[Tuple[str, str], List[Tuple[Optional[float], Optional[float], Optional[float], Optional[float]]]] = defaultdict(list)
    merged_class_pairs: List[Tuple[str, str, str, float]] = []

    chunks = [(chunk, mapping) for chunk in chunked_rows(path, chunk_size)]
    if not chunks:
        log(f"{sample_name}: no rows in {path}")
        return merged_corr, merged_class_pairs

    with Pool(processes=n_workers) as pool:
        for corr_part, pair_part in pool.starmap(summarize_chunk, chunks):
            merge_corr_bucket(merged_corr, corr_part)
            merged_class_pairs.extend(pair_part)

    log(f"finished reading {sample_name}: {path} (chunks={len(chunks)}, class_pair_records={len(merged_class_pairs)})")
    return merged_corr, merged_class_pairs


def compute_corr_rows(sample_name: str, corr_data):
    out_rows = []
    for pair in PAIR_ORDER:
        for b in BINS:
            pts = corr_data.get((pair, b), [])

            logx = [x for x, y, ex, ey in pts if x is not None and y is not None]
            logy = [y for x, y, ex, ey in pts if x is not None and y is not None]
            eix = [ex for x, y, ex, ey in pts if ex is not None and ey is not None]
            eiy = [ey for x, y, ex, ey in pts if ex is not None and ey is not None]

            out_rows.append(
                {
                    "sample": sample_name,
                    "pair_type": pair,
                    "overlap_bin": b,
                    "n_rows": len(pts),
                    "n_log2_pairs": len(logx),
                    "log2_cell_id_pearson_r": pearson_r(logx, logy),
                    "n_ei_pairs": len(eix),
                    "ei_pearson_r": pearson_r(eix, eiy),
                }
            )
    return out_rows


def write_corr_output(rows: List[Dict[str, object]], out_path: Path) -> None:
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fields = [
        "sample",
        "pair_type",
        "overlap_bin",
        "n_rows",
        "n_log2_pairs",
        "log2_cell_id_pearson_r",
        "n_ei_pairs",
        "ei_pearson_r",
    ]
    with out_path.open("w", encoding="utf-8", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fields)
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def write_class_pair_output(rows: List[Dict[str, object]], out_path: Path) -> None:
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fields = [
        "sample",
        "slide",
        "pair_key_unordered",
        "a_class",
        "b_class",
        "overlap_a_in_b",
    ]
    with out_path.open("w", encoding="utf-8", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fields)
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def main() -> None:
    parser = argparse.ArgumentParser(description="按 sample 统计 overlap 分箱相关性（并行 chunk 读取）")
    parser.add_argument(
        "--true-table",
        type=Path,
        default=Path("true_data_overlap_results/true_data_total_overlap_table.csv"),
    )
    parser.add_argument(
        "--sample-glob-dir",
        type=Path,
        default=Path("overlap_results"),
    )
    parser.add_argument(
        "--type-mapping",
        type=Path,
        default=Path("Merfish_brain_cell_type_subclass.txt"),
    )
    parser.add_argument("--workers", type=int, default=20)
    parser.add_argument("--chunk-size", type=int, default=20000)
    parser.add_argument(
        "--corr-output",
        type=Path,
        default=Path("overlap_results/sample_pair_bin_correlations.csv"),
    )
    parser.add_argument(
        "--pair-output",
        type=Path,
        default=Path("overlap_results/sample_slide_class_pair_overlap_a_in_b.csv"),
    )
    args = parser.parse_args()

    mapping = load_type_mapping(args.type_mapping)

    input_files: List[Tuple[str, Path]] = [("true_data", args.true_table)]
    for p in sorted(args.sample_glob_dir.glob("*_cellid_overlap_summary_filtered.csv")):
        sample_name = p.stem.replace("_cellid_overlap_summary_filtered", "")
        input_files.append((sample_name, p))

    all_corr_rows: List[Dict[str, object]] = []
    all_pair_rows: List[Dict[str, object]] = []

    for sample_name, path in input_files:
        if not path.exists():
            log(f"skip missing file: {path}")
            continue
        corr_data, class_pair_data = process_file(
            path=path,
            sample_name=sample_name,
            mapping=mapping,
            n_workers=args.workers,
            chunk_size=args.chunk_size,
        )
        all_corr_rows.extend(compute_corr_rows(sample_name, corr_data))
        for slide, a_class, b_class, overlap_a_in_b in class_pair_data:
            all_pair_rows.append(
                {
                    "sample": sample_name,
                    "slide": slide,
                    "pair_key_unordered": canonical_class_pair(a_class, b_class),
                    "a_class": a_class,
                    "b_class": b_class,
                    "overlap_a_in_b": overlap_a_in_b,
                }
            )

    write_corr_output(all_corr_rows, args.corr_output)
    write_class_pair_output(all_pair_rows, args.pair_output)
    log(f"wrote correlation output: {args.corr_output} (rows={len(all_corr_rows)})")
    log(f"wrote class pair output: {args.pair_output} (rows={len(all_pair_rows)})")


if __name__ == "__main__":
    main()
