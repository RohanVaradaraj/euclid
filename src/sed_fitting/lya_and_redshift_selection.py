"""
lya_and_redshift_selection.py

The final step of the galaxy selection.

Take galaxies with 6.5<zphot<7.5, and galaxies at z<6.5 but with a preferred Lya solution that boosts it to z>6.5.

Created: Tuesday 29th October 2024.
"""

from pathlib import Path
import glob
import sys
import json
import shutil
from astropy.io import ascii

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
zphot_folder = base_det + '_lyaBase'
zphot_dir = Path.cwd().parents[1] / 'data' / 'sed_fitting' / 'zphot' / 'best_fits' / zphot_folder

# Make zphot
if not zphot_dir.exists():
    zphot_dir.mkdir(parents=True)

# Good/maybe files
not_BD_dir = zphot_dir.parents[0] / (base_det + '_notBD')
good_files = glob.glob(str(not_BD_dir / '*.spec'))

# Use the above names to copy files from original lya directory
original_lya_dir = zphot_dir.parents[1] / (base_det + '_lya')

# Get the same .spec files in original_dusty_dir
good_files = [f.split('/')[-1] for f in good_files]

# If the file exists and is not already copied over, copy it
for file in good_files:
    if (original_lya_dir / file).exists() and not (zphot_dir / file).exists():
        shutil.copy(original_lya_dir / file, zphot_dir)
    elif not (original_lya_dir / file).exists():
        print(f'File {file} does not exist in {original_lya_dir}')

# Make the lya/not lya directories. 
not_lya_dir = zphot_dir.parents[0] / (base_det + '_not_lya')
if not not_lya_dir.exists():
    not_lya_dir.mkdir(parents=True)

lya_dir = zphot_dir.parents[0] / (base_det + '_lya')
if not lya_dir.exists():
    lya_dir.mkdir(parents=True)

# Make the z7 directory
z7_dir = zphot_dir.parents[0] / (base_det + '_z7')
if not z7_dir.exists():
    z7_dir.mkdir(parents=True)

not_z7_dir = zphot_dir.parents[0] / (base_det + '_not_z7')
if not not_z7_dir.exists():
    not_z7_dir.mkdir(parents=True)

#! NOTE: _lya means the source has original z below z=6.5 but has a preferred Lya solution that boosts it to z>6.5.
#! So need to compare to the original SED fitting.
normal_sed_dir = zphot_dir.parents[1] / base_det

# Go through files in visual selection directory
spec_files = glob.glob(str(zphot_dir / '*.spec'))
spec_files = sorted(spec_files, key=lambda x: int(x.split('/')[-1].split('Id')[-1].lstrip('0').split('.spec')[0]))

number_z7 = 0
number_not_z7 = 0

number_lya = 0
number_not_lya = 0

for i, spec_file in enumerate(spec_files):

    print(f'Object {i + 1} of {len(spec_files)}')
    ID = spec_file.split('/')[-1].split('Id')[-1].lstrip('0').split('.spec')[0]
    print('ID:', ID)

    params = ascii.read(spec_file, format='basic', data_start=ds_param, data_end=de_param, delimiter=' ', names=names_param)

    # Open original SED. 
    # spec file name is made of 9 digits, with left filled with zeros
    original_sed_file = normal_sed_dir / f'Id{ID.zfill(9)}.spec'
    original_params = ascii.read(original_sed_file, format='basic', data_start=ds_param, data_end=de_param, delimiter=' ', names=names_param)

    # First check if the original SED has z>6.5
    zphot_primary = original_params['Zphot'][0]
    chi2_primary = original_params['Chi2'][0]
    is_z7 = (zphot_primary > 6.5) & (zphot_primary < 7.5)

    # If at z=7, copy over to z7_dir, else copy to not_z7_dir
    if is_z7:
        number_z7 += 1
        shutil.copy2(original_sed_file, z7_dir)
    else:
        number_not_z7 += 1
        shutil.copy2(original_sed_file, not_z7_dir)

    # Now, if it is at z<6.5, check if the Lya solution boosts it to z>6.5
    if not is_z7:

        print('Checking Lya solution...')

        zphot_lya = params['Zphot'][0]
        chi2_lya = params['Chi2'][0]
        is_lya = (zphot_lya > 6.5) & (chi2_lya < chi2_primary)

        if is_lya:
            number_lya += 1
            shutil.copy2(spec_file, lya_dir)
        else:
            number_not_lya += 1
            shutil.copy2(spec_file, not_lya_dir)
    else:
        continue

print(f'Number of z=7 galaxies: {number_z7}')
print(f'Number of not z=7 galaxies: {number_not_z7}')

print(f'Number of Lya galaxies: {number_lya}')
print(f'Number of not Lya galaxies: {number_not_lya}')
