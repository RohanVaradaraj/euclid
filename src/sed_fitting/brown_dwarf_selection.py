"""
dusty_selection.py

Run the brown dwarf cut: chi2_star > 10.

Created: Monday 28th October 2024.
"""

from pathlib import Path
import glob
import sys
import json
import shutil
from astropy.io import fits, ascii

def stellar_type(model):
    stellar_dict = {
        1: 'M4',
        2: 'M5',
        3: 'M6',
        4: 'M7',
        5: 'M8',
        6: 'M9',
        7: 'L0',
        8: 'L1',
        9: 'L2',
        10: 'L3',
        11: 'L4',
        12: 'L5',
        13: 'L6',
        14: 'L7',
        15: 'L8',
        16: 'L9',
        17: 'T0',
        18: 'T1',
        19: 'T2',
        20: 'T3',
        21: 'T4',
        22: 'T5',
        23: 'T6',
        24: 'T7',
        25: 'T8',
    }
    return stellar_dict[model]

if len(sys.argv) > 1:
    filters_json = sys.argv[1]
    filters = json.loads(filters_json)
    bools_json = sys.argv[2]
    bools = json.loads(bools_json)
    all_filters_json = sys.argv[3]
    all_filters = json.loads(all_filters_json)

#! Output PDF name setup from detection filters
det_list = [f for f, t in filters.items() if t['type'] == 'detection']
if len(det_list) != 0:
    base_det = 'det_' + '_'.join(det_list)
else:
    stack_list = [f for f, t in filters.items() if t['type'] == 'stacked-detection']
    stack_filters = stack_list[0].split('+')
    base_det = 'det_' + '_'.join(stack_filters)
    det_list = stack_filters

# Define start and end points of the SED fitting parameter section in the LePhare .spec file
names_param = ['Type', 'Nline', 'Model', 'Library', 'Nband', 'Zphot', 'Zinf', 'Zsup', 'Chi2', 'PDF', 'Extlaw', 'EB-V', 'Lir', 'Age', 'Mass', 'SFR', 'SSFR']
ds_param = 3
de_param = 9

# Directory setup
zphot_folder = base_det + '_notDustyInterlopers_bdBase'
zphot_dir = Path.cwd().parents[1] / 'data' / 'sed_fitting' / 'zphot' / 'best_fits' / zphot_folder

# Make zphot
if not zphot_dir.exists():
    zphot_dir.mkdir(parents=True)

# Good/maybe files, filenames from the not-dusty low-z interloper directory
not_dusty_dir = zphot_dir.parents[0] / (base_det + '_notDustyInterlopers')
good_files = glob.glob(str(not_dusty_dir / '*.spec'))

# Use the above names to copy files from BD directory
original_bd_dir = zphot_dir.parents[1] / (base_det + '_bd')

# Get the same .spec files in original_dusty_dir
good_files = [f.split('/')[-1] for f in good_files]

# If the file exists and is not already copied over, copy it
for file in good_files:
    if (original_bd_dir / file).exists() and not (zphot_dir / file).exists():
        shutil.copy(original_bd_dir / file, zphot_dir)
    elif not (original_bd_dir / file).exists():
        print(f'File {file} does not exist in {original_bd_dir}')
    

# Make the brown dwarf/not brown dwarf directories
not_bd_dir = zphot_dir.parents[0] / (base_det + '_notBD')
if not not_bd_dir.exists():
    not_bd_dir.mkdir(parents=True)

bd_dir = zphot_dir.parents[0] / (base_det + '_BD')
if not bd_dir.exists():
    bd_dir.mkdir(parents=True)

# Go through files in visual selection directory
spec_files = glob.glob(str(zphot_dir / '*.spec'))
spec_files = sorted(spec_files, key=lambda x: int(x.split('/')[-1].split('Id')[-1].lstrip('0').split('.spec')[0]))

number_BD = 0
number_not_BD = 0

for i, spec_file in enumerate(spec_files):

    print(f'Object {i + 1} of {len(spec_files)}')
    ID = spec_file.split('/')[-1].split('Id')[-1].lstrip('0').split('.spec')[0]
    print('ID:', ID)

    params = ascii.read(spec_file, format='basic', data_start=ds_param, data_end=de_param, delimiter=' ', names=names_param)

    star = params['Type'][5]
    chi2_star = params['Chi2'][5]
    # Check if the zphot_primary is larger than zphot_secondary
    solution_is_BD = chi2_star < 10.

    if solution_is_BD:
        number_BD += 1
        # Copy file to bd_dir
        shutil.copy2(spec_file, bd_dir)
    else:
        number_not_BD += 1
        # Copy file to dusty_dir
        shutil.copy2(spec_file, not_bd_dir)

    print('BD PREFERRED:', solution_is_BD)

    zphot_primary = params['Zphot'][0]
    chi2_primary = params['Chi2'][0]

    stellar_model = params['Model'][-1]
    star_type = stellar_type(stellar_model)

    print('Primary:', zphot_primary, chi2_primary)
    print('Stellar:', star_type, chi2_star)
    print('\n')

print(f'Number of BDs: {number_BD}')
print(f'Number of not BDs: {number_not_BD}')



