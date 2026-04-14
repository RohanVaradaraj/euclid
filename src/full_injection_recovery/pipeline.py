#!/usr/bin/env python3

"""
Queue-driven multi-band injection pipeline.

For each wave of image sets:
1) launch up to `processing.batch_size` multi-band set jobs,
2) wait for that wave to complete,
3) optionally clean temporary products,
4) move on until `source_injection.n_images` total sets are processed.
"""

from __future__ import annotations

from pathlib import Path
import glob
import json
import os
import time
import numpy as np
from astropy.table import Table

from utils import load_config


def _safe_job_name(text: str) -> str:
    return ''.join(ch if ch.isalnum() or ch in ['_', '-'] else '_' for ch in text)


def _resolve_images_lis_path(config: dict) -> Path:
    field_name = config['field']['name']
    processing_cfg = config.get('processing', {})

    if processing_cfg.get('images_lis'):
        return Path(processing_cfg['images_lis'])

    # Default location used elsewhere in this project.
    return Path.cwd().parents[3] / 'data' / field_name / 'images.lis'


def _read_images_table(config: dict) -> Table:
    images_lis = _resolve_images_lis_path(config)
    if not images_lis.is_file():
        raise FileNotFoundError(f"Could not find images.lis at: {images_lis}")

    table = Table.read(images_lis, format='ascii.commented_header')
    required = ['Name', 'Image', 'Weight', 'directory']
    missing = [col for col in required if col not in table.colnames]
    if missing:
        raise ValueError(f"images.lis missing required columns: {missing}")

    processing_cfg = config.get('processing', {})
    selected_names = processing_cfg.get('selected_names', [])
    if selected_names:
        keep = [str(row['Name']) in selected_names for row in table]
        table = table[keep]

    return table


def _chunk_indices(total_count: int, batch_size: int):
    for start in range(0, total_count, batch_size):
        stop = min(start + batch_size, total_count)
        yield range(start, stop)


def _write_worker_shell(script_path: Path, task_file: Path, done_file: Path):
    script_path.parent.mkdir(parents=True, exist_ok=True)

    cmd = (
        f"{Path.cwd() / 'queued_worker.py'} "
        f"--task-file {task_file} "
        f"--done-file {done_file}"
    )

    with open(script_path, 'w') as f:
        f.write('#!/bin/bash\n')
        f.write(f'cd {Path.cwd()}\n')
        f.write(cmd + '\n')

    os.system(f'chmod u+x {script_path}')


def _cleanup_transient_images():
    cutout_dir = Path.cwd() / 'images' / 'cutouts'
    cutout_wht_dir = cutout_dir / 'weights'
    injected_dir = Path.cwd() / 'images' / 'injected'

    for pattern in [
        cutout_dir / '*.fits',
        cutout_wht_dir / '*.fits',
        injected_dir / '*.fits',
        (injected_dir / 'weights') / '*.fits',
    ]:
        for fp in glob.glob(str(pattern)):
            try:
                os.remove(fp)
            except FileNotFoundError:
                pass


def _cleanup_runtime_state():
    for directory in [
        Path.cwd() / 'runtime_tasks',
        Path.cwd() / 'runtime_done',
        Path.cwd() / 'shell_scripts',
    ]:
        directory.mkdir(exist_ok=True)
        for fp in directory.glob('*'):
            if fp.is_file():
                fp.unlink()


def _match_truth_to_output(truth_table: Table, output_table: Table, match_tolerance_pixels: float) -> np.ndarray:
    truth_x = np.asarray(truth_table['x'], dtype=float)
    truth_y = np.asarray(truth_table['y'], dtype=float)
    output_x = np.asarray(output_table['X_IMAGE'], dtype=float)
    output_y = np.asarray(output_table['Y_IMAGE'], dtype=float)

    matched_indices = np.full(len(truth_table), -1, dtype=int)
    available_output = set(range(len(output_table)))
    max_dist2 = match_tolerance_pixels ** 2

    for truth_idx in range(len(truth_table)):
        if not available_output:
            break

        candidate_indices = np.fromiter(available_output, dtype=int)
        dx = output_x[candidate_indices] - truth_x[truth_idx]
        dy = output_y[candidate_indices] - truth_y[truth_idx]
        dist2 = dx ** 2 + dy ** 2

        nearest = int(np.argmin(dist2))
        if dist2[nearest] <= max_dist2:
            output_idx = int(candidate_indices[nearest])
            matched_indices[truth_idx] = output_idx
            available_output.remove(output_idx)

    return matched_indices


def _make_unmatched_column(template_column, n_rows: int):
    dtype_kind = template_column.dtype.kind
    if dtype_kind in ['i', 'u', 'f']:
        return np.full(n_rows, -99, dtype=template_column.dtype)
    if dtype_kind == 'b':
        return np.full(n_rows, False, dtype=template_column.dtype)
    return np.full(n_rows, 'NONE', dtype=template_column.dtype)


def _safe_column_suffix(text: str) -> str:
    return ''.join(ch if ch.isalnum() else '_' for ch in str(text))


def _append_matched_filter_columns(
    matched_table: Table,
    truth_table: Table,
    output_table: Table,
    filter_name: str,
    match_tolerance_pixels: float,
):
    matched_indices = _match_truth_to_output(truth_table, output_table, match_tolerance_pixels)
    recovered_mask = matched_indices >= 0
    suffix = _safe_column_suffix(filter_name)
    matched_table[f'recovered_{suffix}'] = recovered_mask

    for column_name in output_table.colnames:
        fill_values = _make_unmatched_column(output_table[column_name], len(truth_table))
        if np.any(recovered_mask):
            fill_values[recovered_mask] = np.asarray(output_table[column_name][matched_indices[recovered_mask]])
        matched_table[f'{column_name}_{suffix}'] = fill_values


def _match_and_cleanup_output_catalogues(match_tolerance_pixels: float):
    input_dir = Path.cwd() / 'catalogues' / 'input'
    output_dir = Path.cwd() / 'catalogues' / 'output'
    matched_dir = Path.cwd() / 'catalogues' / 'matched'
    matched_dir.mkdir(parents=True, exist_ok=True)

    for summary_file in sorted(input_dir.glob('set_*_summary.json')):
        with open(summary_file, 'r') as f:
            summary = json.load(f)

        truth_catalog = input_dir / summary['truth_catalog']
        if not truth_catalog.is_file():
            continue

        truth_table = Table.read(truth_catalog)
        cutouts = summary.get('cutouts', [])
        if not cutouts:
            continue

        detection_cutout = summary.get('detection_cutout')
        if detection_cutout is None:
            continue

        detection_stem = detection_cutout.replace('.fits', '')
        matched_table = truth_table.copy(copy_data=True)
        output_catalogues_to_delete = []

        for row in cutouts:
            measurement_stem = row['cutout_name'].replace('.fits', '')
            output_name = f'{measurement_stem}__det_{detection_stem}_cat.fits'
            output_catalog = output_dir / output_name
            if not output_catalog.is_file():
                continue

            output_table = Table.read(output_catalog)
            _append_matched_filter_columns(
                matched_table=matched_table,
                truth_table=truth_table,
                output_table=output_table,
                filter_name=row['filter_name'],
                match_tolerance_pixels=match_tolerance_pixels,
            )
            output_catalogues_to_delete.append(output_catalog)

        matched_name = f"{summary['set_id']}_matched_catalog.fits"
        matched_table.write(matched_dir / matched_name, overwrite=True)

        for output_catalog in output_catalogues_to_delete:
            output_catalog.unlink()

        truth_catalog.unlink()
        summary_file.unlink()


def RunFullInjectionRecoveryPipeline(overwrite=True):
    config = load_config('config.yaml')
    se_cfg = config.get('source_extraction', {})
    injection_cfg = config.get('source_injection', {})
    processing_cfg = config.get('processing', {})

    queue = processing_cfg.get('queue', 'cmb')
    queue_memory_gb = int(processing_cfg.get('queue_memory_gb', 2))
    batch_size = int(processing_cfg.get('batch_size', 50))
    check_interval = int(processing_cfg.get('check_interval', 15))
    cleanup_batch_images = bool(processing_cfg.get('cleanup_batch_images', True))
    match_tolerance_pixels = float(processing_cfg.get('match_tolerance_pixels', 3.0))

    rows = list(_read_images_table(config))
    if len(rows) == 0:
        raise RuntimeError('No usable rows found in images.lis')

    configured_filters = list(injection_cfg.get('filters', []))
    row_by_name = {str(row['Name']): row for row in rows}
    selected_rows = [row_by_name[flt] for flt in configured_filters if flt in row_by_name]
    missing_filters = [flt for flt in configured_filters if flt not in row_by_name]

    if missing_filters:
        raise RuntimeError(
            f"Configured filters missing from images.lis selection: {missing_filters}"
        )
    if len(selected_rows) == 0:
        raise RuntimeError('No selected filter rows available for multi-band set creation.')

    n_image_sets = int(injection_cfg.get('n_images', 0))
    if n_image_sets <= 0:
        raise RuntimeError('source_injection.n_images must be a positive integer.')

    output_cat_dir = Path.cwd() / 'catalogues' / 'output'
    matched_cat_dir = Path.cwd() / 'catalogues' / 'matched'
    if overwrite and bool(processing_cfg.get('overwrite_output_catalogues', False)):
        for file in glob.glob(str(output_cat_dir / '*.fits')):
            os.remove(file)
        for file in glob.glob(str(matched_cat_dir / '*.fits')):
            os.remove(file)

    if cleanup_batch_images:
        _cleanup_transient_images()

    _cleanup_runtime_state()

    task_dir = Path.cwd() / 'runtime_tasks'
    shell_dir = Path.cwd() / 'shell_scripts'
    done_dir = Path.cwd() / 'runtime_done'

    task_dir.mkdir(exist_ok=True)
    shell_dir.mkdir(exist_ok=True)
    done_dir.mkdir(exist_ok=True)

    all_batches = list(_chunk_indices(n_image_sets, batch_size))
    print(
        f'Processing {n_image_sets} multi-band image sets in {len(all_batches)} waves of up to {batch_size}'
    )
    print(f'Each set uses {len(selected_rows)} filters: {[str(row["Name"]) for row in selected_rows]}')

    rows_payload = [
        {
            'name': str(row['Name']),
            'image': str(row['Image']),
            'weight': str(row['Weight']),
            'directory': str(row['directory']),
        }
        for row in selected_rows
    ]

    for bidx, batch_indices in enumerate(all_batches, start=1):
        batch_indices = list(batch_indices)
        print(
            f'\n=== Starting wave {bidx}/{len(all_batches)} with {len(batch_indices)} image sets '
            f'({batch_indices[0] + 1}-{batch_indices[-1] + 1} of {n_image_sets}) ==='
        )

        done_files = []
        fail_files = []
        wave_shell_files = []

        for set_idx in batch_indices:
            set_id = f'set_{set_idx + 1:05d}'
            task_id = _safe_job_name(set_id)
            task_file = task_dir / f'{task_id}.json'
            done_file = done_dir / f'{task_id}.done'
            fail_file = done_dir / f'{task_id}.failed'
            shell_file = shell_dir / f'run_worker_{task_id}.sh'

            if done_file.exists():
                done_file.unlink()
            if fail_file.exists():
                fail_file.unlink()

            with open(task_file, 'w') as f:
                json.dump(
                    {
                        'config_file': str(Path.cwd() / 'config.yaml'),
                        'set_id': set_id,
                        'rows': rows_payload,
                        'se_aperture_diameter_arcsec': float(se_cfg.get('aperture_diameter_arcsec', 2.0)),
                        'se_overwrite': bool(se_cfg.get('overwrite', True)),
                    },
                    f,
                    indent=2,
                )

            _write_worker_shell(shell_file, task_file, done_file)
            submit_shell = Path('shell_scripts') / shell_file.name
            os.system(f"addqueue -c {task_id} -m {queue_memory_gb} -q {queue} ./{submit_shell}")

            done_files.append(done_file)
            fail_files.append(fail_file)
            wave_shell_files.append(shell_file)

        print(f'Waiting for wave {bidx} jobs to complete...')
        while True:
            n_done = sum(int(fp.exists()) for fp in done_files)
            n_fail = sum(int(fp.exists()) for fp in fail_files)

            if n_fail > 0:
                raise RuntimeError(
                    f'Wave {bidx}: {n_fail} worker(s) failed. See *.failed files in {done_dir}'
                )

            if n_done == len(done_files):
                break

            print(f'Wave {bidx}: {n_done}/{len(done_files)} sets complete. Re-checking in {check_interval}s...')
            time.sleep(check_interval)

        print(f'Wave {bidx} complete.')
        _match_and_cleanup_output_catalogues(match_tolerance_pixels)
        for shell_file in wave_shell_files:
            if shell_file.exists():
                shell_file.unlink()



if __name__ == '__main__':
    RunFullInjectionRecoveryPipeline(overwrite=True)
