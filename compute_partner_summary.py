#!/usr/bin/env python3
import csv
from pathlib import Path
from typing import Dict, List, Tuple

TRUE_DATA_FILE = Path('true_data_overlap_results/true_data_total_overlap_table.csv')
SAMPLE_GLOB_DIR = Path('overlap_results')
SAMPLE_GLOB_PATTERN = 'sample_*_cellid_overlap_summary_filtered.csv'
OUTPUT_FILE = Path('overlap_results/sample_true_partner_summary.csv')


def parse_num(value: str) -> float:
    try:
        return float(value)
    except (TypeError, ValueError):
        return 0.0


def load_rows(path: Path) -> List[dict]:
    with path.open('r', newline='', encoding='utf-8') as f:
        reader = csv.DictReader(f)
        return list(reader)


def summarize_rows(rows: List[dict], sample_name: str) -> List[dict]:
    unique_keys: List[Tuple[str, str, str]] = []
    seen = set()
    for row in rows:
        key = (row['slide'], row['a_class'], row['a_merge_region'])
        if key not in seen:
            seen.add(key)
            unique_keys.append(key)

    enrich_sum_by_class: Dict[str, float] = {}
    for slide, a_class, a_merge_region in unique_keys:
        # use the unique key row's a_enrich_class_cell_ids_num value
        for row in rows:
            if row['slide'] == slide and row['a_class'] == a_class and row['a_merge_region'] == a_merge_region:
                enrich_sum_by_class[a_class] = enrich_sum_by_class.get(a_class, 0.0) + parse_num(
                    row.get('a_enrich_class_cell_ids_num', '')
                )
                break

    results: List[dict] = []
    for slide, a_class, a_merge_region in unique_keys:
        partner_number = 0
        for row in rows:
            if row['slide'] == slide and row['a_class'] == a_class and row['a_merge_region'] == a_merge_region:
                if parse_num(row.get('overlap_a_in_b', '')) > 0:
                    partner_number += 1

        results.append(
            {
                'sample_name': sample_name,
                'a_class': a_class,
                'enrich_cell_num_sum': enrich_sum_by_class.get(a_class, 0.0),
                'a_merge_region': a_merge_region,
                'partner_number': partner_number,
            }
        )

    return results


def infer_sample_name(file_path: Path) -> str:
    # sample_1_cellid_overlap_summary_filtered.csv -> sample_1
    name = file_path.name
    return name.replace('_cellid_overlap_summary_filtered.csv', '')


def main() -> None:
    all_results: List[dict] = []

    true_rows = load_rows(TRUE_DATA_FILE)
    all_results.extend(summarize_rows(true_rows, 'true_data'))

    sample_files = sorted(SAMPLE_GLOB_DIR.glob(SAMPLE_GLOB_PATTERN))
    for sample_file in sample_files:
        sample_rows = load_rows(sample_file)
        sample_name = infer_sample_name(sample_file)
        all_results.extend(summarize_rows(sample_rows, sample_name))

    OUTPUT_FILE.parent.mkdir(parents=True, exist_ok=True)
    with OUTPUT_FILE.open('w', newline='', encoding='utf-8') as f:
        fieldnames = [
            'sample_name',
            'a_class',
            'enrich_cell_num_sum',
            'a_merge_region',
            'partner_number',
        ]
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(all_results)

    print(f'Wrote {len(all_results)} rows to {OUTPUT_FILE}')


if __name__ == '__main__':
    main()
