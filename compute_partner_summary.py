#!/usr/bin/env python3
import csv
import math
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
from typing import Dict, Iterable, List, Tuple

TRUE_DATA_FILE = Path('true_data_overlap_results/true_data_total_overlap_table.csv')
SAMPLE_GLOB_DIR = Path('overlap_results')
SAMPLE_GLOB_PATTERN = 'sample_*_cellid_overlap_summary_filtered.csv'
OUTPUT_FILE = Path('overlap_results/sample_true_partner_summary.csv')
MAX_WORKERS = 20

# Worker globals for process pool
_WORK_KEY_TO_ROWS: Dict[Tuple[str, str, str], List[dict]] = {}
_WORK_ENRICH_SUM_BY_CLASS: Dict[str, float] = {}
_WORK_SAMPLE_NAME = ''


def parse_num(value: str) -> float:
    try:
        return float(value)
    except (TypeError, ValueError):
        return 0.0


def load_rows(path: Path) -> List[dict]:
    with path.open('r', newline='', encoding='utf-8') as f:
        reader = csv.DictReader(f)
        return list(reader)


def split_chunks(items: List[Tuple[str, str, str]], chunk_size: int) -> Iterable[List[Tuple[str, str, str]]]:
    for i in range(0, len(items), chunk_size):
        yield items[i:i + chunk_size]


def build_key_to_rows(rows: List[dict]) -> Dict[Tuple[str, str, str], List[dict]]:
    key_to_rows: Dict[Tuple[str, str, str], List[dict]] = {}
    for row in rows:
        key = (row['slide'], row['a_class'], row['a_merge_region'])
        key_to_rows.setdefault(key, []).append(row)
    return key_to_rows



def init_worker(
    key_to_rows: Dict[Tuple[str, str, str], List[dict]],
    enrich_sum_by_class: Dict[str, float],
    sample_name: str,
) -> None:
    global _WORK_KEY_TO_ROWS, _WORK_ENRICH_SUM_BY_CLASS, _WORK_SAMPLE_NAME
    _WORK_KEY_TO_ROWS = key_to_rows
    _WORK_ENRICH_SUM_BY_CLASS = enrich_sum_by_class
    _WORK_SAMPLE_NAME = sample_name


def process_chunk(chunk_keys: List[Tuple[str, str, str]]) -> List[dict]:
    output: List[dict] = []
    for slide, a_class, a_merge_region in chunk_keys:
        rows = _WORK_KEY_TO_ROWS[(slide, a_class, a_merge_region)]
        partner_number = sum(1 for row in rows if parse_num(row.get('overlap_a_in_b', '')) > 0)
        output.append(
            {
                'sample_name': _WORK_SAMPLE_NAME,
                'a_class': a_class,
                'enrich_cell_num_sum': _WORK_ENRICH_SUM_BY_CLASS.get(a_class, 0.0),
                'a_merge_region': a_merge_region,
                'partner_number': partner_number,
            }
        )
    return output


def summarize_rows_serial(rows: List[dict], sample_name: str) -> List[dict]:
    key_to_rows = build_key_to_rows(rows)
    unique_keys = list(key_to_rows.keys())

    enrich_sum_by_class: Dict[str, float] = {}
    for key in unique_keys:
        a_class = key[1]
        row0 = key_to_rows[key][0]
        enrich_sum_by_class[a_class] = enrich_sum_by_class.get(a_class, 0.0) + parse_num(
            row0.get('a_enrich_class_cell_ids_num', '')
        )

    results: List[dict] = []
    for slide, a_class, a_merge_region in unique_keys:
        partner_number = sum(
            1
            for row in key_to_rows[(slide, a_class, a_merge_region)]
            if parse_num(row.get('overlap_a_in_b', '')) > 0
        )
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


def summarize_rows_parallel(rows: List[dict], sample_name: str, max_workers: int = MAX_WORKERS) -> List[dict]:
    key_to_rows = build_key_to_rows(rows)
    unique_keys = list(key_to_rows.keys())

    enrich_sum_by_class: Dict[str, float] = {}
    for key in unique_keys:
        a_class = key[1]
        row0 = key_to_rows[key][0]
        enrich_sum_by_class[a_class] = enrich_sum_by_class.get(a_class, 0.0) + parse_num(
            row0.get('a_enrich_class_cell_ids_num', '')
        )

    if not unique_keys:
        return []

    workers = max(1, min(max_workers, len(unique_keys)))
    chunk_size = max(1, math.ceil(len(unique_keys) / (workers * 4)))
    chunks = list(split_chunks(unique_keys, chunk_size))

    all_results: List[dict] = []
    with ProcessPoolExecutor(
        max_workers=workers,
        initializer=init_worker,
        initargs=(key_to_rows, enrich_sum_by_class, sample_name),
    ) as executor:
        for chunk_result in executor.map(process_chunk, chunks):
            all_results.extend(chunk_result)

    return all_results


def infer_sample_name(file_path: Path) -> str:
    # sample_1_cellid_overlap_summary_filtered.csv -> sample_1
    name = file_path.name
    return name.replace('_cellid_overlap_summary_filtered.csv', '')


def main() -> None:
    all_results: List[dict] = []

    true_rows = load_rows(TRUE_DATA_FILE)
    all_results.extend(summarize_rows_serial(true_rows, 'true_data'))

    sample_files = sorted(SAMPLE_GLOB_DIR.glob(SAMPLE_GLOB_PATTERN))
    for idx, sample_file in enumerate(sample_files, start=1):
        sample_rows = load_rows(sample_file)
        sample_name = infer_sample_name(sample_file)
        print(f'[sample {idx}/{len(sample_files)}] start: {sample_name}, input_rows={len(sample_rows)}', flush=True)
        sample_results = summarize_rows_parallel(sample_rows, sample_name, max_workers=MAX_WORKERS)
        all_results.extend(sample_results)
        print(
            f'[sample {idx}/{len(sample_files)}] done: {sample_name}, output_rows={len(sample_results)}',
            flush=True,
        )

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
