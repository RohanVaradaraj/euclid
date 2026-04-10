#!/usr/bin/env python3

"""
Batch-wise queue-driven injection/recovery pipeline.

For each batch of rows from images.lis:
1) launch cutout+inject+recover worker jobs to the queue,
2) wait for batch completion,
3) clean temporary images,
4) move to next batch.
"""

from __future__ import annotations

from astropy.table import Table
from pathlib import Path
import glob
import json
import os
import time

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


def _chunk_rows(rows, batch_size: int):
    for i in range(0, len(rows), batch_size):
        yield rows[i:i + batch_size]


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


def RunFullInjectionRecoveryPipeline(overwrite=True):
    config = load_config('config.yaml')
    se_cfg = config.get('source_extraction', {})
    processing_cfg = config.get('processing', {})

    queue = processing_cfg.get('queue', 'cmb')
    queue_memory_gb = int(processing_cfg.get('queue_memory_gb', 2))
    batch_size = int(processing_cfg.get('batch_size', 50))
    check_interval = int(processing_cfg.get('check_interval', 15))
    cleanup_batch_images = bool(processing_cfg.get('cleanup_batch_images', True))

    rows = list(_read_images_table(config))
    if len(rows) == 0:
        raise RuntimeError('No usable rows found in images.lis')

    output_cat_dir = Path.cwd() / 'catalogues' / 'output'
    if overwrite and bool(processing_cfg.get('overwrite_output_catalogues', False)):
        for file in glob.glob(str(output_cat_dir / '*.fits')):
            os.remove(file)

    task_dir = Path.cwd() / 'runtime_tasks'
    shell_dir = Path.cwd() / 'shell_scripts'
    done_dir = Path.cwd() / 'runtime_done'

    task_dir.mkdir(exist_ok=True)
    shell_dir.mkdir(exist_ok=True)
    done_dir.mkdir(exist_ok=True)

    all_batches = list(_chunk_rows(rows, batch_size))
    print(f'Processing {len(rows)} images in {len(all_batches)} batches of up to {batch_size}')

    for bidx, batch in enumerate(all_batches, start=1):
        print(f'\n=== Starting batch {bidx}/{len(all_batches)} with {len(batch)} images ===')

        done_files = []
        fail_files = []

        for ridx, row in enumerate(batch, start=1):
            name = str(row['Name'])
            task_id = _safe_job_name(f'b{bidx:03d}_{ridx:03d}_{name}')

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
                        'name': str(row['Name']),
                        'image': str(row['Image']),
                        'weight': str(row['Weight']),
                        'directory': str(row['directory']),
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

        print(f'Waiting for batch {bidx} jobs to complete...')
        while True:
            n_done = sum(1 for fp in done_files if fp.exists())
            n_fail = sum(1 for fp in fail_files if fp.exists())

            if n_fail > 0:
                raise RuntimeError(
                    f'Batch {bidx}: {n_fail} worker(s) failed. See *.failed files in {done_dir}'
                )

            if n_done == len(done_files):
                break

            print(f'Batch {bidx}: {n_done}/{len(done_files)} complete. Re-checking in {check_interval}s...')
            time.sleep(check_interval)

        print(f'Batch {bidx} complete.')

        if cleanup_batch_images:
            _cleanup_transient_images()
            print(f'Batch {bidx}: cleaned temporary cutout/injected images.')


if __name__ == '__main__':
    RunFullInjectionRecoveryPipeline(overwrite=True)
