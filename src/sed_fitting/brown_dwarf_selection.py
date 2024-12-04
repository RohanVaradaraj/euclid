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
import sys
sed_path = Path.cwd().parents[0] / 'sed_fitting'
sys.path.append(str(sed_path))
from sed_fitting_codes import parse_spec_file

verbose = False 

overwrite = True

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
    run_type_json = sys.argv[4]
    run_type = json.loads(run_type_json)

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
zphot_folder = base_det + '_notDustyInterlopers_bdBase'
zphot_dir = Path.cwd().parents[1] / 'data' / 'sed_fitting' / 'zphot' / 'best_fits' / zphot_folder

# Make zphot
if not zphot_dir.exists():
    zphot_dir.mkdir(parents=True)

# If overwrite, clear the directory
if overwrite:
    for file in zphot_dir.glob('*.spec'):
        file.unlink()

# Normal SED fitting dir
lbg_dir = Path.cwd().parents[1] / 'data' / 'sed_fitting' / 'zphot' / 'best_fits' / (base_det + '_best_highz')

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

#* Also make a directory for strong brown dwarfs when running the strong/weak BD selection. Copy weak ones into normal BD directory.
strong_bd_dir = zphot_dir.parents[0] / (base_det + '_strongBD')
if not strong_bd_dir.exists():
    strong_bd_dir.mkdir(parents=True)

# If overwrite is True, delete all previous files in the above
if overwrite:
    for directory in [not_bd_dir, bd_dir]:
        for file in directory.glob('*.spec'):
            file.unlink()

# Go through files in visual selection directory
spec_files = glob.glob(str(zphot_dir / '*.spec'))
print('Taking files from:', zphot_dir)
spec_files = sorted(spec_files, key=lambda x: int(x.split('/')[-1].split('Id')[-1].lstrip('0').split('.spec')[0]))

number_BD = 0
number_not_BD = 0

number_weak_BD = 0
number_strong_BD = 0

for i, spec_file in enumerate(spec_files):

    ID = spec_file.split('/')[-1].split('Id')[-1].lstrip('0').split('.spec')[0]

    if verbose:
        print(f'Object {i + 1} of {len(spec_files)}')
        print('ID:', ID)

    params = parse_spec_file(spec_file).get('model')
    params.rename_columns(params.colnames, names_param)

    # Open corresponding .spec file in lbg dir
    lbg_file = lbg_dir / f'Id{ID.zfill(9)}.spec'
    lbg_params = parse_spec_file(lbg_file).get('model')
    chi2_highz = lbg_params['Chi2'][0]

    star = params['Type'][5]
    chi2_star = params['Chi2'][5]

    #! ---------------------------------------------------------------------------------------------------------
    #! Selection step: VISTA chi2<10, otherwise chi2_star > chi2_highz and then split into strong and weak BDs
    #! ---------------------------------------------------------------------------------------------------------

    #! Check if the solution is a Brown Dwarf (BD) or not
    if run_type == '':
        #? For VISTA samples, use the chi2=10 cut from Bowler+15
        solution_is_BD = chi2_star < 10
    else:
        #? For other run types, compare chi2_highz and chi2_star
        solution_is_BD = chi2_highz > chi2_star

    #* 1) Solution is a brown dwarf
    if solution_is_BD:
        number_BD += 1  # Increment BD counter

        #? If run type is provided, split into strong and weak BDs
        if run_type != '':
            # Split into strong and weak BD based on delta chi2
            delta_chi2 = chi2_highz - chi2_star
            print('For this object:', delta_chi2)
            if delta_chi2 > 4:
                solution_strength = 'strong'
                number_strong_BD += 1 # Increment strong BD counter
                shutil.copy2(spec_file, strong_bd_dir)  # Copy to strong BD directory
            else:
                solution_strength = 'weak'
                number_weak_BD += 1  # Increment weak BD counter
                shutil.copy2(spec_file, bd_dir)  # Copy to normal BD directory

        #? Otherwise for VISTA samples, copy to BD directory
        else:
            shutil.copy2(spec_file, bd_dir)  # Copy to BD directory for VISTA samples
    
    #* 2) Solution is not a brown dwarf
    else:
        number_not_BD += 1  # Increment non-BD counter
        shutil.copy2(spec_file, not_bd_dir)  # Copy to not BD directory

    zphot_primary = params['Zphot'][0]
    chi2_primary = params['Chi2'][0]

    stellar_model = params['Model'][-1]
    star_type = stellar_type(stellar_model)

    if verbose:
        print('BD PREFERRED:', solution_is_BD)
        print('Primary:', zphot_primary, chi2_primary)
        print('Stellar:', star_type, chi2_star)
        print('\n')

print(f'Number of BDs: {number_BD}')
if run_type != '':
    print(f'Number of strong BDs: {number_strong_BD}')
    print(f'Number of weak BDs: {number_weak_BD}')
print(f'Number of not BDs: {number_not_BD}')



