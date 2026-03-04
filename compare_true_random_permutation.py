#!/usr/bin/env python3
"""比较 sample_level_summary 与 true_data 的多指标差异，并计算置换检验风格 p 值。"""

from __future__ import annotations

import argparse
import csv
import glob
import math
import os
from collections import Counter, defaultdict
from multiprocessing import Pool
from statistics import mean
from typing import DefaultDict, Dict, Iterable, List, Tuple


Key = Tuple[str, str, str]  # (slide, subclass, layer)
MetricMap = Dict[Key, Dict[str, float]]
RandomMap = Dict[Key, List[Tuple[str, Dict[str, float]]]]


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


def to_float(row: dict, key: str) -> float | None:
    value = row.get(key, "")
    if value is None:
        return None
    value = str(value).strip()
    if value == "":
        return None
    try:
        return float(value)
    except ValueError:
        return None


def permutation_pvalues(true_value: float, random_values: List[float]) -> Dict[str, float]:
    n = len(random_values)
    random_ge = sum(v >= true_value for v in random_values)
    random_le = sum(v <= true_value for v in random_values)
    abs_diff_true = [abs(v - true_value) for v in random_values]
    random_more_extreme = sum(d >= 0 for d in abs_diff_true)  # 恒为 n，使用 centered 版本见下

    center = mean(random_values)
    true_centered = abs(true_value - center)
    random_centered = [abs(v - center) for v in random_values]
    random_more_extreme = sum(v >= true_centered for v in random_centered)

    return {
        "p_one_sided_true_gt": (random_ge + 1.0) / (n + 1.0),
        "p_one_sided_true_lt": (random_le + 1.0) / (n + 1.0),
        "p_two_sided": (random_more_extreme + 1.0) / (n + 1.0),
        "random_ge_count": float(random_ge),
        "random_le_count": float(random_le),
        "random_n": float(n),
    }


def load_random_summary(path: str) -> RandomMap:
    random_map: DefaultDict[Key, List[Tuple[str, Dict[str, float]]]] = defaultdict(list)

    for row in read_tsv(path):
        key = (
            row.get("slide", "").strip(),
            clean_subclass(row.get("subclass", "")),
            normalize_layer(row.get("layer", "")),
        )
        sample = row.get("sample", "").strip() or "NA"

        if not all(key):
            continue

        metric = {
            "region_num": to_float(row, "region_num"),
            "glut_num": to_float(row, "glut_num"),
            "gaba_num": to_float(row, "gaba_num"),
            "ei_ratio": to_float(row, "glut_gaba_ratio"),
        }

        if metric["region_num"] is None:
            continue

        random_map[key].append((sample, metric))

    return random_map


def choose_true_files(true_dir: str) -> List[str]:
    patt = os.path.join(true_dir, "*_merged_regions_table_cell_id.txt")
    files = sorted(glob.glob(patt))
    if files:
        return files

    patt2 = os.path.join(true_dir, "*_merged_rgions_loc.txt")
    return sorted(glob.glob(patt2))


def load_true_metrics(true_dir: str) -> Tuple[MetricMap, Counter]:
    true_files = choose_true_files(true_dir)
    if not true_files:
        raise FileNotFoundError(f"在 {true_dir} 中找不到可用真实文件")

    per_key_region: DefaultDict[Key, List[float]] = defaultdict(list)
    per_key_glut: DefaultDict[Key, List[float]] = defaultdict(list)
    per_key_gaba: DefaultDict[Key, List[float]] = defaultdict(list)
    per_key_ei: DefaultDict[Key, List[float]] = defaultdict(list)
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

        key = (slide, subclass, layer)
        regions = {
            r.get("region", "").strip()
            for r in rows
            if r.get("region", "").strip()
        }
        if not regions:
            skipped["missing_region"] += 1
            continue

        per_key_region[key].append(float(len(regions)))

        glut_vals = [to_float(r, "Glut_Neruon_cell_ids_num") for r in rows]
        gaba_vals = [to_float(r, "GABA_Neruon_cell_ids_num") for r in rows]
        glut_vals = [v for v in glut_vals if v is not None]
        gaba_vals = [v for v in gaba_vals if v is not None]

        if glut_vals:
            per_key_glut[key].append(mean(glut_vals))
        if gaba_vals:
            per_key_gaba[key].append(mean(gaba_vals))

        ratio_vals = []
        for r in rows:
            g1 = to_float(r, "Glut_Neruon_cell_ids_num")
            g2 = to_float(r, "GABA_Neruon_cell_ids_num")
            if g1 is None or g2 is None or g2 == 0:
                continue
            ratio_vals.append(g1 / g2)
        if ratio_vals:
            per_key_ei[key].append(mean(ratio_vals))

    true_map: MetricMap = {}
    all_keys = set(per_key_region.keys())
    for key in all_keys:
        true_map[key] = {
            "region_num": mean(per_key_region[key]),
            "glut_num": mean(per_key_glut[key]) if per_key_glut[key] else math.nan,
            "gaba_num": mean(per_key_gaba[key]) if per_key_gaba[key] else math.nan,
            "ei_ratio": mean(per_key_ei[key]) if per_key_ei[key] else math.nan,
            "true_file_count": float(len(per_key_region[key])),
        }

    return true_map, skipped


def compute_statistics(true_map: MetricMap, random_map: RandomMap) -> Tuple[List[dict], List[dict], Counter]:
    summary_rows: List[dict] = []
    detail_rows: List[dict] = []
    missing = Counter()

    metric_names = ["region_num", "glut_num", "gaba_num", "ei_ratio"]

    for key in sorted(true_map.keys()):
        slide, subclass, layer = key
        true_metrics = true_map[key]

        random_values = random_map.get(key)
        if not random_values:
            missing["no_random_match"] += 1
            continue

        summary_row = {
            "slide": slide,
            "subclass": subclass,
            "layer": layer,
            "true_file_count": int(true_metrics["true_file_count"]),
            "random_sample_n": len(random_values),
        }

        for metric in metric_names:
            true_value = true_metrics.get(metric)
            if true_value is None or math.isnan(true_value):
                summary_row[f"true_{metric}"] = "nan"
                continue

            rand_for_metric = [m[metric] for _, m in random_values if m.get(metric) is not None]
            if not rand_for_metric:
                summary_row[f"true_{metric}"] = f"{true_value:.6f}"
                summary_row[f"{metric}_p_one_sided_true_gt"] = "nan"
                summary_row[f"{metric}_p_one_sided_true_lt"] = "nan"
                summary_row[f"{metric}_p_two_sided"] = "nan"
                continue

            pvals = permutation_pvalues(true_value, rand_for_metric)
            summary_row[f"true_{metric}"] = f"{true_value:.6f}"
            summary_row[f"{metric}_p_one_sided_true_gt"] = f"{pvals['p_one_sided_true_gt']:.6g}"
            summary_row[f"{metric}_p_one_sided_true_lt"] = f"{pvals['p_one_sided_true_lt']:.6g}"
            summary_row[f"{metric}_p_two_sided"] = f"{pvals['p_two_sided']:.6g}"

        summary_rows.append(summary_row)

        for sample, metric_dict in random_values:
            detail_row = {
                "slide": slide,
                "subclass": subclass,
                "layer": layer,
                "sample": sample,
            }
            for metric in metric_names:
                tv = true_metrics.get(metric)
                rv = metric_dict.get(metric)
                detail_row[f"true_{metric}"] = "nan" if tv is None or math.isnan(tv) else f"{tv:.6f}"
                detail_row[f"random_{metric}"] = "nan" if rv is None else f"{rv:.6f}"
                if rv is None or tv is None or math.isnan(tv):
                    detail_row[f"true_gt_random_{metric}"] = "nan"
                else:
                    detail_row[f"true_gt_random_{metric}"] = int(tv > rv)
            detail_rows.append(detail_row)

    return summary_rows, detail_rows, missing


def _compute_chunk(
    chunk_keys: List[Key],
    true_map: MetricMap,
    random_map: RandomMap,
) -> Tuple[List[dict], List[dict], Counter]:
    chunk_true_map = {key: true_map[key] for key in chunk_keys}
    return compute_statistics(chunk_true_map, random_map)


def chunked_compute_statistics(
    true_map: MetricMap,
    random_map: RandomMap,
    workers: int,
    chunk_size: int,
    show_progress: bool = True,
) -> Tuple[List[dict], List[dict], Counter]:
    keys = sorted(true_map.keys())
    if not keys:
        return [], [], Counter()

    workers = max(1, workers)
    if workers == 1 or len(keys) <= chunk_size:
        if show_progress:
            print("[进度] 已完成 1/1 chunks (100.0%)")
        return compute_statistics(true_map, random_map)

    chunks = [keys[i : i + chunk_size] for i in range(0, len(keys), chunk_size)]
    args = [(chunk, true_map, random_map) for chunk in chunks]

    summary_rows: List[dict] = []
    detail_rows: List[dict] = []
    missing = Counter()
    total_chunks = len(chunks)
    done_chunks = 0
    with Pool(processes=workers) as pool:
        for chunk_summary, chunk_detail, chunk_missing in pool.starmap(_compute_chunk, args):
            summary_rows.extend(chunk_summary)
            detail_rows.extend(chunk_detail)
            missing.update(chunk_missing)
            done_chunks += 1
            if show_progress:
                percent = (done_chunks / total_chunks) * 100.0
                print(f"[进度] 已完成 {done_chunks}/{total_chunks} chunks ({percent:.1f}%)", flush=True)

    summary_rows.sort(key=lambda r: (r["slide"], r["subclass"], r["layer"]))
    return summary_rows, detail_rows, missing


def compute_global_totals(true_map: MetricMap, random_map: RandomMap) -> List[dict]:
    true_totals = {
        "region_num": sum(v["region_num"] for v in true_map.values() if not math.isnan(v["region_num"])),
        "glut_num": sum(v["glut_num"] for v in true_map.values() if not math.isnan(v["glut_num"])),
        "gaba_num": sum(v["gaba_num"] for v in true_map.values() if not math.isnan(v["gaba_num"])),
    }

    sample_totals: DefaultDict[str, Dict[str, float]] = defaultdict(lambda: {"region_num": 0.0, "glut_num": 0.0, "gaba_num": 0.0})
    for _, sample_metrics in random_map.items():
        for sample, metric in sample_metrics:
            for metric_name in ["region_num", "glut_num", "gaba_num"]:
                value = metric.get(metric_name)
                if value is not None:
                    sample_totals[sample][metric_name] += value

    rows = []
    for metric_name in ["region_num", "glut_num", "gaba_num"]:
        random_values = [totals[metric_name] for totals in sample_totals.values()]
        pvals = permutation_pvalues(true_totals[metric_name], random_values)
        rows.append(
            {
                "metric": metric_name,
                "true_total": f"{true_totals[metric_name]:.6f}",
                "random_total_mean": f"{mean(random_values):.6f}",
                "random_total_min": f"{min(random_values):.6f}",
                "random_total_max": f"{max(random_values):.6f}",
                "sample_n": len(random_values),
                "p_one_sided_true_gt": f"{pvals['p_one_sided_true_gt']:.6g}",
                "p_one_sided_true_lt": f"{pvals['p_one_sided_true_lt']:.6g}",
                "p_two_sided": f"{pvals['p_two_sided']:.6g}",
            }
        )
    return rows


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
    parser = argparse.ArgumentParser(description="比较真实数据与随机 sample 的 region/glut/gaba/EI 多指标")
    parser.add_argument("--sample-summary", default="sample_level_summary.txt", help="number.R 生成的 sample_level_summary.txt")
    parser.add_argument("--true-dir", default="true_data", help="真实数据目录")
    parser.add_argument("--out-summary", default="true_vs_random_permutation_summary.tsv", help="逐键汇总输出")
    parser.add_argument("--out-detail", default="true_vs_random_permutation_detail.tsv", help="逐sample详情输出")
    parser.add_argument("--out-global", default="true_vs_random_global_totals.tsv", help="总体统计输出")
    parser.add_argument("--workers", type=int, default=20, help="并行进程数")
    parser.add_argument("--chunk-size", type=int, default=0, help="每个并行任务处理的 key 数；0 自动")
    args = parser.parse_args()

    random_map = load_random_summary(args.sample_summary)
    if not random_map:
        raise RuntimeError(f"未从 {args.sample_summary} 读取到有效随机汇总")

    true_map, skipped = load_true_metrics(args.true_dir)
    if not true_map:
        raise RuntimeError(f"未从 {args.true_dir} 读取到有效真实统计")

    chunk_size = args.chunk_size
    if chunk_size <= 0:
        chunk_size = max(1, math.ceil(len(true_map) / max(1, args.workers)))

    summary_rows, detail_rows, missing = chunked_compute_statistics(
        true_map=true_map,
        random_map=random_map,
        workers=args.workers,
        chunk_size=chunk_size,
    )
    global_rows = compute_global_totals(true_map, random_map)

    write_tsv(args.out_summary, summary_rows)
    write_tsv(args.out_detail, detail_rows)
    write_tsv(args.out_global, global_rows)

    print("完成")
    print(f"- 随机键数量: {len(random_map)}")
    print(f"- 真实键数量: {len(true_map)}")
    print(f"- 成功匹配键数量: {len(summary_rows)}")
    print(f"- 逐sample比较记录数: {len(detail_rows)}")
    print(f"- 总体统计指标数: {len(global_rows)}")
    print(f"- 并行进程数: {max(1, args.workers)}")
    print(f"- 分块大小: {chunk_size}")
    if skipped:
        print(f"- 跳过真实文件统计: {dict(skipped)}")
    if missing:
        print(f"- 未匹配统计: {dict(missing)}")
    print(f"- 汇总输出: {args.out_summary}")
    print(f"- 详情输出: {args.out_detail}")
    print(f"- 总体输出: {args.out_global}")


if __name__ == "__main__":
    main()
