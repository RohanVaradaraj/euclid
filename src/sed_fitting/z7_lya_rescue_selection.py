#!/usr/bin/env python3

"""
z7_lya_rescue_selection.py

Build a final z7 directory directly from visually accepted high-z SEDs.

This skips the dusty and brown-dwarf filtering steps. It copies sources from
det_*_best_highz_good into det_*_z7 if either:

1. The normal high-z fit has 6.5 < zphot < 7.5.
2. The normal fit has zphot < 6.5, but the matching Lya fit has zphot > 6.5
   and a lower chi2 than the normal fit.
"""

from pathlib import Path
import glob
import sys
import json
import shutil

SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = SCRIPT_DIR.parents[1]
sys.path.append(str(SCRIPT_DIR))

from sed_fitting_codes import parse_spec_file


# Defaults for running this script directly from the command line.
filters = {
    'Y+J': {'type': 'stacked-detection', 'value': 5},
    'HSC-G_DR3': {'type': 'non-detection', 'value': 2},
    'HSC-R_DR3': {'type': 'non-detection', 'value': 2},
    'HSC-I_DR3': {'type': 'non-detection', 'value': 2},
}
run_type = 'with_euclid'
field_name = 'COSMOS'

overwrite = True
dry_run = '--dry-run' in sys.argv

z7_min = 6.5
z7_max = 7.5
require_lya_better_chi2 = True
apply_z7_max_to_lya = False

names_param = [
    'Type', 'Nline', 'Model', 'Library', 'Nband', 'Zphot', 'Zinf', 'Zsup',
    'Chi2', 'PDF', 'Extlaw', 'EB-V', 'Lir', 'Age', 'Mass', 'SFR', 'SSFR'
]


args = [arg for arg in sys.argv[1:] if arg != '--dry-run']
if len(args) > 0:
    filters = json.loads(args[0])
    # args[1] and args[2] are bools/all_filters in the pipeline interface.
    run_type = json.loads(args[3])
    field_name = json.loads(args[4])


def make_base_det(filters, run_type):
    det_list = [f for f, t in filters.items() if t['type'] == 'detection']
    if det_list:
        base_det = 'det_' + '_'.join(det_list)
    else:
        stack_list = [f for f, t in filters.items() if t['type'] == 'stacked-detection']
        if not stack_list:
            raise ValueError('No detection or stacked-detection filter found.')
        base_det = 'det_' + '_'.join(stack_list[0].split('+'))

    if run_type:
        base_det += f'_{run_type}'

    return base_det


def spec_sort_key(spec_file):
    return int(Path(spec_file).stem.split('Id')[-1].lstrip('0'))


def read_model_params(spec_file):
    params = parse_spec_file(spec_file).get('model')
    params.rename_columns(params.colnames, names_param)
    return params


def copy_spec(source, destination_dir):
    if dry_run:
        return
    shutil.copy2(source, destination_dir)


base_det = make_base_det(filters, run_type)

zphot_root = PROJECT_ROOT / 'data' / 'sed_fitting' / 'zphot' / field_name
best_fits_dir = zphot_root / 'best_fits'

good_dir = best_fits_dir / f'{base_det}_best_highz_good'
lya_dir = zphot_root / f'{base_det}_lya'
z7_dir = best_fits_dir / f'{base_det}_z7'

print(f'Taking visually accepted sources from: {good_dir}')
print(f'Taking Lya fits from: {lya_dir}')
print(f'Writing final z7 selection to: {z7_dir}')

if not good_dir.exists():
    raise FileNotFoundError(f'Input good directory does not exist: {good_dir}')
if not lya_dir.exists():
    raise FileNotFoundError(f'Lya directory does not exist: {lya_dir}')

if not dry_run:
    z7_dir.mkdir(parents=True, exist_ok=True)

if overwrite and z7_dir.exists() and not dry_run:
    for file in z7_dir.glob('*.spec'):
        file.unlink()

spec_files = glob.glob(str(good_dir / '*.spec'))
spec_files = sorted(spec_files, key=spec_sort_key)

number_normal_z7 = 0
number_lya_rescued = 0
number_not_z7 = 0
number_missing_lya = 0

for spec_file in spec_files:
    spec_file = Path(spec_file)
    params = read_model_params(spec_file)

    zphot_normal = params['Zphot'][0]
    chi2_normal = params['Chi2'][0]

    is_normal_z7 = (zphot_normal > z7_min) & (zphot_normal < z7_max)
    if is_normal_z7:
        number_normal_z7 += 1
        copy_spec(spec_file, z7_dir)
        continue

    if zphot_normal >= z7_min:
        number_not_z7 += 1
        continue

    lya_file = lya_dir / spec_file.name
    if not lya_file.exists():
        number_missing_lya += 1
        number_not_z7 += 1
        continue

    lya_params = read_model_params(lya_file)
    zphot_lya = lya_params['Zphot'][0]
    chi2_lya = lya_params['Chi2'][0]

    is_lya_z7 = zphot_lya > z7_min
    if apply_z7_max_to_lya:
        is_lya_z7 = is_lya_z7 & (zphot_lya < z7_max)

    lya_is_better = True
    if require_lya_better_chi2:
        lya_is_better = chi2_lya < chi2_normal

    if is_lya_z7 and lya_is_better:
        number_lya_rescued += 1
        copy_spec(lya_file, z7_dir)
    else:
        number_not_z7 += 1

print(f'Normal 6.5<z<7.5 sources copied to _z7: {number_normal_z7}')
print(f'Lya-rescued sources copied to _z7: {number_lya_rescued}')
print(f'Total sources copied to _z7: {number_normal_z7 + number_lya_rescued}')
print(f'Sources not copied to _z7: {number_not_z7}')
print(f'Missing Lya files for normal z<{z7_min} sources: {number_missing_lya}')

if dry_run:
    print('Dry run only: no files were copied or deleted.')
