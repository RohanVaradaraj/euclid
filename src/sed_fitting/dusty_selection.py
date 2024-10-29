"""
dusty_selection.py

Run the dusty low-z interloper cut: take the SEDs with preferred high-z when IRAC included.
Then go through the rejected sources and make sure the IRAC photometry is unconfused.

Created: Monday 28th October 2024.
"""

from pathlib import Path
import glob
import sys
import json
import shutil
from astropy.io import fits, ascii

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
zphot_folder = base_det + '_visualSelection'
zphot_dir = Path.cwd().parents[1] / 'data' / 'sed_fitting' / 'zphot' / 'best_fits' / zphot_folder

# Make the base visual selection directory
if not zphot_dir.exists():
    zphot_dir.mkdir(parents=True)

# Get the filenames from good_dir and maybe_dir, then copy from the original_dusty_dir into zphot_dir
good_dir = zphot_dir.parents[0] / (base_det + f'_best_highz_good')
maybe_dir = zphot_dir.parents[0] / (base_det + f'_best_highz_maybe')
original_dusty_dir = zphot_dir.parents[1] / (base_det + '_dusty')

# Good/maybe files
good_files = glob.glob(str(good_dir / '*.spec'))
maybe_files = glob.glob(str(maybe_dir / '*.spec'))

# Get the same .spec files in original_dusty_dir
good_files = [f.split('/')[-1] for f in good_files]
maybe_files = [f.split('/')[-1] for f in maybe_files]

# If the file exists and is not already copied over, copy it
for file in good_files + maybe_files:
    if (original_dusty_dir / file).exists() and not (zphot_dir / file).exists():
        shutil.copy(original_dusty_dir / file, zphot_dir)
    elif not (original_dusty_dir / file).exists():
        raise FileNotFoundError(f'File {file} not found in original_dusty_dir.')
    elif (zphot_dir / file).exists():
        print(f'File {file} already exists in zphot_dir.')



# Make the not dusty directory
not_dusty_dir = zphot_dir.parents[0] / (base_det + '_notDustyInterlopers')
if not not_dusty_dir.exists():
    not_dusty_dir.mkdir(parents=True)

# Make the dusty directory
dusty_dir = zphot_dir.parents[0] / (base_det + '_dustyInterlopers')
if not dusty_dir.exists():
    dusty_dir.mkdir(parents=True)

# Go through files in visual selection directory
spec_files = glob.glob(str(zphot_dir / '*.spec'))
spec_files = sorted(spec_files, key=lambda x: int(x.split('/')[-1].split('Id')[-1].lstrip('0').split('.spec')[0]))

number_low_z = 0
number_high_z = 0

for i, spec_file in enumerate(spec_files):

    print(f'Object {i + 1} of {len(spec_files)}')
    ID = spec_file.split('/')[-1].split('Id')[-1].lstrip('0').split('.spec')[0]
    print('ID:', ID)

    params = ascii.read(spec_file, format='basic', data_start=ds_param, data_end=de_param, delimiter=' ', names=names_param)

    zphot_primary = round(params['Zphot'][0], 2)
    chi2_primary = round(params['Chi2'][0], 1)

    zphot_secondary = round(params['Zphot'][1], 2)
    chi2_secondary = round(params['Chi2'][1], 1)

    # Check if the zphot_primary is larger than zphot_secondary
    solution_is_highz = zphot_primary > zphot_secondary

    if solution_is_highz:
        number_high_z += 1
        # Copy file to not_dusty_dir
        shutil.copy2(spec_file, not_dusty_dir)
    else:
        number_low_z += 1
        # Copy file to dusty_dir
        shutil.copy2(spec_file, dusty_dir)

    print('HIGH-Z PREFERRED:', solution_is_highz)

    print('Primary:', zphot_primary, chi2_primary)
    print('Secondary:', zphot_secondary, chi2_secondary)
    print('\n')

print(f'Number of high-z preferred: {number_high_z}')
print(f'Number of low-z preferred: {number_low_z}')

