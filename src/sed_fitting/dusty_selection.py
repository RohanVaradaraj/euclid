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
import sys
sed_path = Path.cwd().parents[0] / 'sed_fitting'
sys.path.append(str(sed_path))
from sed_fitting_codes import parse_spec_file

verbose = False

overwrite = True

if len(sys.argv) > 1:
    filters_json = sys.argv[1]
    filters = json.loads(filters_json)
    bools_json = sys.argv[2]
    bools = json.loads(bools_json)
    all_filters_json = sys.argv[3]
    all_filters = json.loads(all_filters_json)
    run_type_json = sys.argv[4]
    run_type = json.loads(run_type_json)
    field_name_json = sys.argv[5]
    field_name = json.loads(field_name_json)

#! Output PDF name setup from detection filters
det_list = [f for f, t in filters.items() if t['type'] == 'detection']
if det_list:
    base_det = 'det_' + '_'.join(det_list)
else:
    stack_filters = next((f.split('+') for f, t in filters.items() if t['type'] == 'stacked-detection'), [])
    base_det = 'det_' + '_'.join(stack_filters)
    det_list = stack_filters

# Append run type if provided
base_det += f'_{run_type}' if run_type else ''

# Define start and end points of the SED fitting parameter section in the LePhare .spec file
names_param = ['Type', 'Nline', 'Model', 'Library', 'Nband', 'Zphot', 'Zinf', 'Zsup', 'Chi2', 'PDF', 'Extlaw', 'EB-V', 'Lir', 'Age', 'Mass', 'SFR', 'SSFR']

# Directory setup
zphot_folder = base_det + '_visualSelection'
zphot_dir = Path.cwd().parents[1] / 'data' / 'sed_fitting' / 'zphot' / field_name / 'best_fits' / zphot_folder

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
        if verbose:
            print(f'File {file} already exists in zphot_dir.')

# Make the not dusty directory
not_dusty_dir = zphot_dir.parents[0] / (base_det + '_notDustyInterlopers')
if not not_dusty_dir.exists():
    not_dusty_dir.mkdir(parents=True)

# Make the dusty directory
dusty_dir = zphot_dir.parents[0] / (base_det + '_dustyInterlopers')
if not dusty_dir.exists():
    dusty_dir.mkdir(parents=True)

# Print directories
if verbose:
    print('Visual selection directory:', zphot_dir)
    print('Not dusty directory:', not_dusty_dir)
    print('Dusty directory:', dusty_dir)

# If overwrite is True, delete all previous files in the above
if overwrite:
    for directory in [not_dusty_dir, dusty_dir]:
        for file in directory.glob('*.spec'):
            file.unlink()
            
# Go through files in visual selection directory
# spec_files = glob.glob(str(zphot_dir / '*.spec'))
# spec_files = sorted(spec_files, key=lambda x: int(x.split('/')[-1].split('Id')[-1].lstrip('0').split('.spec')[0]))

spec_files = glob.glob(str(good_dir / '*.spec'))
spec_files = sorted(spec_files, key=lambda x: int(x.split('/')[-1].split('Id')[-1].lstrip('0').split('.spec')[0]))

number_low_z = 0
number_high_z = 0

print(len(spec_files), 'files to process.')

for i, spec_file in enumerate(spec_files):

    ID = spec_file.split('/')[-1].split('Id')[-1].lstrip('0').split('.spec')[0]

    if verbose:
        print(f'Object {i + 1} of {len(spec_files)}')
        print('ID:', ID)

    # Read spec file from dusty directory. Name is Id and 9 digits with leading zeros.
    spec_file = original_dusty_dir / spec_file.split('/')[-1]
    #print(i, spec_file)

    params = parse_spec_file(spec_file).get('model')
    params.rename_columns(params.colnames, names_param)

    zphot_primary = round(params['Zphot'][0], 2)
    chi2_primary = round(params['Chi2'][0], 1)

    zphot_secondary = round(params['Zphot'][1], 2)
    chi2_secondary = round(params['Chi2'][1], 1)

    #! --------------------------------------------------------------------------
    #! Selection step: Check if the zphot_primary is larger than zphot_secondary
    #! --------------------------------------------------------------------------
    solution_is_highz = (zphot_primary > zphot_secondary) & (zphot_primary > 5) & (zphot_secondary < 5) & (chi2_secondary - chi2_primary > 4)

    if solution_is_highz:
        number_high_z += 1
        # Copy file to not_dusty_dir
        shutil.copy2(spec_file, not_dusty_dir)
    else:
        number_low_z += 1
        # Copy file to dusty_dir
        shutil.copy2(spec_file, dusty_dir)

    if verbose:
        print('HIGH-Z PREFERRED:', solution_is_highz)
        print('Primary:', zphot_primary, chi2_primary)
        print('Secondary:', zphot_secondary, chi2_secondary)
        print('\n')

print(f'Number of high-z preferred: {number_high_z}')
print(f'Number of low-z preferred: {number_low_z}')

