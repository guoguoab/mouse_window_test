#!/usr/bin/env python3
"""直接基于 superclass 百分比表绘制 Sample vs True 对比图（每 4 个 pair 一张 SVG）。

输入：
1) sample_pair_overlap_percent_superclass.csv
2) true_data_pair_overlap_percent_superclass.csv（或任意同结构 true 百分比表）
3) sample_vs_true_permutation_pvalues_superclass.csv（可选；如不存在则自动重算）

输出：
- {out-prefix}_part01.svg, {out-prefix}_part02.svg, ...（每张最多 4 个 pair）
- sample_vs_true_permutation_pvalues_superclass.recomputed.csv（当未提供/不存在 pvalue CSV 时）
"""

from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path
from statistics import mean
from typing import Dict, List, Tuple

BINS = ["0–20%", "20–40%", "40–60%", "60–80%", "80–100%"]
COLORS = ["#d9d2ad", "#b7d4ca", "#57a9a5", "#3d80ad", "#294f72"]
SAMPLE_BAR_COLOR = "#7b4ab8"
TRUE_LABEL_COLOR = "#1f1f1f"
SAMPLE_LABEL_COLOR = "#4b237a"


def chunked(items: List[str], n: int) -> List[List[str]]:
    return [items[i : i + n] for i in range(0, len(items), n)]


def load_true_distribution(path: Path) -> Dict[str, Dict[str, float]]:
    out: Dict[str, Dict[str, float]] = {}
    with path.open("r", encoding="utf-8", newline="") as fh:
        reader = csv.DictReader(fh)
        for row in reader:
            pair = (row.get("pair_type") or "").strip()
            b = (row.get("overlap_bin") or "").strip()
            percent = row.get("percent") or row.get("true_percent") or "0"
            if not pair or b not in BINS:
                continue
            out.setdefault(pair, {bb: 0.0 for bb in BINS})
            out[pair][b] = float(percent)
    return out


def load_sample_distribution(path: Path) -> Dict[str, Dict[str, Dict[str, float]]]:
    out: Dict[str, Dict[str, Dict[str, float]]] = {}
    with path.open("r", encoding="utf-8", newline="") as fh:
        reader = csv.DictReader(fh)
        for row in reader:
            sample = (row.get("sample") or "").strip()
            pair = (row.get("pair_type") or "").strip()
            b = (row.get("overlap_bin") or "").strip()
            if not sample or not pair or b not in BINS:
                continue
            out.setdefault(sample, {})
            out[sample].setdefault(pair, {bb: 0.0 for bb in BINS})
            out[sample][pair][b] = float(row.get("percent") or 0.0)
    return out


def infer_pair_order(sample_dists: Dict[str, Dict[str, Dict[str, float]]], true_dist: Dict[str, Dict[str, float]]) -> List[str]:
    sample_pairs = set()
    for sample in sample_dists.values():
        sample_pairs.update(sample.keys())
    true_pairs = set(true_dist.keys())

    def sort_key(pair: str) -> Tuple[str, str]:
        # pair 可能含多个 '-'，仅用于稳定排序
        left = pair.split("-", 1)[0]
        return left, pair

    all_pairs = sorted(sample_pairs | true_pairs, key=sort_key)
    return all_pairs


def permutation_direct_pvalue(random_vals: List[float], true_val: float) -> float:
    n = len(random_vals)
    if n == 0:
        return float("nan")
    extreme = sum(1 for x in random_vals if x >= true_val)
    return (extreme + 1.0) / (n + 1.0)


def stars_for_p(p: float) -> str:
    if math.isnan(p):
        return "NA"
    if p < 0.001:
        return "***"
    if p < 0.01:
        return "**"
    if p < 0.05:
        return "*"
    return "ns"


def compute_stats(
    pair_order: List[str],
    true_dist: Dict[str, Dict[str, float]],
    sample_dists: Dict[str, Dict[str, Dict[str, float]]],
) -> Dict[Tuple[str, str], Dict[str, float]]:
    stats: Dict[Tuple[str, str], Dict[str, float]] = {}
    sample_names = sorted(sample_dists)
    for pair in pair_order:
        for b in BINS:
            vals = [sample_dists[s].get(pair, {}).get(b, 0.0) for s in sample_names]
            tv = true_dist.get(pair, {}).get(b, 0.0)
            m = mean(vals) if vals else 0.0
            sd = math.sqrt(sum((x - m) ** 2 for x in vals) / (len(vals) - 1)) if len(vals) > 1 else 0.0
            pv = permutation_direct_pvalue(vals, tv)
            stats[(pair, b)] = {"true": tv, "mean": m, "sd": sd, "p": pv, "sig": stars_for_p(pv)}
    return stats


def load_stats_from_pvalue_csv(
    path: Path,
    pair_order: List[str],
    true_dist: Dict[str, Dict[str, float]],
    sample_dists: Dict[str, Dict[str, Dict[str, float]]],
) -> Dict[Tuple[str, str], Dict[str, float]]:
    stats = compute_stats(pair_order, true_dist, sample_dists)
    with path.open("r", encoding="utf-8", newline="") as fh:
        reader = csv.DictReader(fh)
        for row in reader:
            pair = (row.get("pair_type") or "").strip()
            b = (row.get("overlap_bin") or "").strip()
            if (pair, b) not in stats:
                continue
            stats[(pair, b)]["p"] = float(row.get("p_value") or "nan")
            stats[(pair, b)]["sig"] = (row.get("significance") or "").strip() or stars_for_p(stats[(pair, b)]["p"])
    return stats


def write_recomputed_pvalue_csv(
    pair_order: List[str],
    stats: Dict[Tuple[str, str], Dict[str, float]],
    out_csv: Path,
    sample_n: int,
) -> None:
    out_csv.parent.mkdir(parents=True, exist_ok=True)
    with out_csv.open("w", encoding="utf-8", newline="") as fh:
        writer = csv.DictWriter(
            fh,
            fieldnames=["pair_type", "overlap_bin", "true_percent", "sample_n", "sample_mean", "sample_sd", "p_value", "significance"],
        )
        writer.writeheader()
        for pair in pair_order:
            for b in BINS:
                st = stats[(pair, b)]
                writer.writerow(
                    {
                        "pair_type": pair,
                        "overlap_bin": b,
                        "true_percent": f"{st['true']:.6f}",
                        "sample_n": sample_n,
                        "sample_mean": f"{st['mean']:.6f}",
                        "sample_sd": f"{st['sd']:.6f}",
                        "p_value": "nan" if math.isnan(st["p"]) else f"{st['p']:.6g}",
                        "significance": st["sig"],
                    }
                )


def draw_plot_for_pairs(
    pair_subset: List[str],
    true_dist: Dict[str, Dict[str, float]],
    sample_dists: Dict[str, Dict[str, Dict[str, float]]],
    stats: Dict[Tuple[str, str], Dict[str, float]],
    out_svg: Path,
) -> None:
    width, height = 1500, 760
    left, right, top, bottom = 90, 40, 120, 120
    chart_w = width - left - right
    chart_h = height - top - bottom

    max_true = max((true_dist.get(p, {}).get(b, 0.0) for p in pair_subset for b in BINS), default=0.0)
    max_sample = max((sample_dists[s].get(p, {}).get(b, 0.0) for s in sample_dists for p in pair_subset for b in BINS), default=0.0)
    max_mean_sd = max((stats[(p, b)]["mean"] + stats[(p, b)]["sd"] for p in pair_subset for b in BINS), default=0.0)
    y_max = max(60.0, max_true, max_sample, max_mean_sd) + 15.0

    def y_to_px(v: float) -> float:
        return top + chart_h - (v / y_max) * chart_h

    group_w = chart_w / len(pair_subset)
    bar_w = 18
    pair_gap = 3
    bin_gap = 8

    svg: List[str] = []
    svg.append(f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}">')
    svg.append('<rect width="100%" height="100%" fill="white"/>')
    svg.append('<text x="750" y="40" text-anchor="middle" font-size="24" font-weight="700">Sample vs True Superclass Overlap Percent (Permutation Test)</text>')
    svg.append(f'<line x1="{left}" y1="{top}" x2="{left}" y2="{top + chart_h}" stroke="#222" stroke-width="3"/>')
    svg.append(f'<line x1="{left}" y1="{top + chart_h}" x2="{left + chart_w}" y2="{top + chart_h}" stroke="#222" stroke-width="3"/>')

    for tick in range(0, int(y_max) + 1, 10):
        y = y_to_px(float(tick))
        svg.append(f'<line x1="{left - 6}" y1="{y:.2f}" x2="{left}" y2="{y:.2f}" stroke="#222" stroke-width="2"/>')
        svg.append(f'<text x="{left - 12}" y="{y + 4:.2f}" text-anchor="end" font-size="14">{tick}%</text>')

    lx, ly = 140, 72
    for i, b in enumerate(BINS):
        x = lx + i * 145
        svg.append(f'<rect x="{x}" y="{ly - 12}" width="24" height="14" fill="{COLORS[i]}" stroke="#222" stroke-width="1.6"/>')
        svg.append(f'<text x="{x + 32}" y="{ly}" font-size="15">{b}</text>')

    for g, pair in enumerate(pair_subset):
        center = left + group_w * (g + 0.5)
        cluster_w = len(BINS) * (2 * bar_w + pair_gap) + (len(BINS) - 1) * bin_gap
        start_x = center - cluster_w / 2

        for i, b in enumerate(BINS):
            bin_start = start_x + i * (2 * bar_w + pair_gap + bin_gap)
            true_x = bin_start
            mean_x = bin_start + bar_w + pair_gap
            tv = true_dist.get(pair, {}).get(b, 0.0)
            y_true = y_to_px(tv)
            h_true = top + chart_h - y_true
            svg.append(f'<rect x="{true_x:.2f}" y="{y_true:.2f}" width="{bar_w}" height="{h_true:.2f}" fill="{COLORS[i]}" fill-opacity="0.35" stroke="#222" stroke-width="1.8"/>')
            svg.append(f'<text x="{true_x + bar_w / 2:.2f}" y="{max(12, y_true - 4):.2f}" text-anchor="middle" font-size="9" font-weight="600" fill="{TRUE_LABEL_COLOR}">{tv:.1f}%</text>')

            st = stats[(pair, b)]
            m = st["mean"]
            sd = st["sd"]
            y_mean = y_to_px(m)
            h_mean = top + chart_h - y_mean
            svg.append(f'<rect x="{mean_x:.2f}" y="{y_mean:.2f}" width="{bar_w}" height="{h_mean:.2f}" fill="{SAMPLE_BAR_COLOR}" fill-opacity="0.9" stroke="#222" stroke-width="1.8"/>')
            svg.append(f'<text x="{mean_x + bar_w / 2:.2f}" y="{max(12, y_mean - 4):.2f}" text-anchor="middle" font-size="9" font-weight="600" fill="{SAMPLE_LABEL_COLOR}">{m:.1f}%</text>')

            ex = mean_x + bar_w / 2
            y1 = y_to_px(max(0.0, m - sd))
            y2 = y_to_px(m + sd)
            svg.append(f'<line x1="{ex:.2f}" y1="{y1:.2f}" x2="{ex:.2f}" y2="{y2:.2f}" stroke="#111" stroke-width="2.2"/>')
            svg.append(f'<line x1="{ex - 4:.2f}" y1="{y1:.2f}" x2="{ex + 4:.2f}" y2="{y1:.2f}" stroke="#111" stroke-width="2.2"/>')
            svg.append(f'<line x1="{ex - 4:.2f}" y1="{y2:.2f}" x2="{ex + 4:.2f}" y2="{y2:.2f}" stroke="#111" stroke-width="2.2"/>')

            sig_y = y_to_px(max(tv, m + sd) + 5)
            sig_x = (true_x + mean_x + bar_w) / 2
            svg.append(f'<text x="{sig_x:.2f}" y="{max(16, sig_y):.2f}" text-anchor="middle" font-size="13" font-weight="700">{st["sig"]}</text>')

        label = pair if len(pair) <= 24 else pair[:24] + "…"
        svg.append(f'<text x="{center:.2f}" y="{top + chart_h + 45}" text-anchor="middle" font-size="17">{label}</text>')

    svg.append('</svg>')
    out_svg.parent.mkdir(parents=True, exist_ok=True)
    out_svg.write_text("\n".join(svg), encoding="utf-8")


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Directly plot superclass sample/true overlap from precomputed CSV files.")
    p.add_argument("--sample-csv", type=Path, default=Path("overlap_results/sample_pair_overlap_percent_superclass.csv"))
    p.add_argument("--true-csv", type=Path, default=Path("true_data_overlap_results/true_data_pair_overlap_percent_superclass.csv"))
    p.add_argument("--pvalue-csv", type=Path, default=Path("overlap_results/sample_vs_true_permutation_pvalues_superclass.csv"))
    p.add_argument("--out-dir", type=Path, default=Path("overlap_results"))
    p.add_argument("--out-prefix", type=str, default="sample_vs_true_superclass_overlap_plot")
    p.add_argument("--pairs-per-figure", type=int, default=4, help="每张图展示的 pair 数（默认 4）")
    return p.parse_args()


def main() -> None:
    args = parse_args()
    for p in [args.sample_csv, args.true_csv]:
        if not p.exists():
            raise SystemExit(f"Missing required input: {p}")

    sample_dists = load_sample_distribution(args.sample_csv)
    if not sample_dists:
        raise SystemExit(f"No sample data loaded from: {args.sample_csv}")

    true_dist = load_true_distribution(args.true_csv)
    pair_order = infer_pair_order(sample_dists, true_dist)
    if not pair_order:
        raise SystemExit("No pair_type found in sample/true csv.")

    if args.pvalue_csv.exists():
        stats = load_stats_from_pvalue_csv(args.pvalue_csv, pair_order, true_dist, sample_dists)
    else:
        stats = compute_stats(pair_order, true_dist, sample_dists)
        out_csv = args.out_dir / "sample_vs_true_permutation_pvalues_superclass.recomputed.csv"
        write_recomputed_pvalue_csv(pair_order, stats, out_csv, sample_n=len(sample_dists))
        print(f"Recomputed p-value table -> {out_csv}")

    groups = chunked(pair_order, max(1, args.pairs_per_figure))
    generated: List[Path] = []
    for idx, pair_subset in enumerate(groups, start=1):
        out_svg = args.out_dir / f"{args.out_prefix}_part{idx:02d}.svg"
        draw_plot_for_pairs(pair_subset, true_dist, sample_dists, stats, out_svg)
        generated.append(out_svg)

    print(f"Samples: {len(sample_dists)}")
    print(f"Total pair types: {len(pair_order)}")
    print(f"Figures generated: {len(generated)}")
    for p in generated:
        print(f"Plot -> {p}")


if __name__ == "__main__":
    main()
