#!/usr/bin/env python3

"""
Look at the number of all brown dwarfs and compare to Rebecca's models.

Created: Tuesday 1st July 2025.
"""

from pathlib import Path
import sys
sed_fitting_path = Path.cwd().parents[0] / 'sed_fitting'
sys.path.append(str(sed_fitting_path))

from sed_fitting_codes import parse_spec_file
import matplotlib.pyplot as plt
import glob
import numpy as np
from astropy.table import Table

plt.rcParams['axes.linewidth'] = 2.5
plt.rcParams.update({'font.size': 15})
plt.rcParams['figure.dpi'] = 100

field_name = 'COSMOS'

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
        -1: 'No fit'
    }
    return stellar_dict[model]

#! Directories of brown dwarfs
sed_dir = Path.cwd().parents[1] / 'data' / 'sed_fitting' / 'zphot' / field_name / 'det_HSC-Z_DR3_bd'

# Get all spec files
spec_files = glob.glob(str(sed_dir / '*.spec'))

# Make a new astropy Table to store the results
bd_table = Table(
    names=('ID', 'mJ', 'stellar_type', 'chi2_star'),
    dtype=('i8', 'f4', 'U10', 'f4')  # U10 means Unicode string of max length 10
)

results = []

for i, spec_file in enumerate(spec_files):
    if i % 10 == 0:
        print(f'Processing file {i+1}/{len(spec_files)}: {spec_file.split("/")[-1]}')

    file = parse_spec_file(spec_file)
    params = file.get('model')
    phot = file.get('phot')

    chi2_star = params[params['ID'] == 'STAR']['Chi2'][0]
    chi2_lbg = params[params['ID'] == 'GAL-1']['Chi2'][0]
    chi2_dusty = params[params['ID'] == 'GAL-2']['Chi2'][0]

    if chi2_star < chi2_lbg and chi2_star < chi2_dusty:
        ID = spec_file.split('/')[-1].split('.')[0].split('Id')[-1].lstrip('0')
        model = int(params[params['ID'] == 'STAR']['Model'][0])
        stellar_model = stellar_type(model)
        mJ = phot['col1'][5]

        # Store row as tuple
        results.append((int(ID), mJ, stellar_model, chi2_star))

# Only now, convert to Astropy Table
bd_table = Table(rows=results, names=('ID', 'mJ', 'stellar_type', 'chi2_star'),
                 dtype=('i8', 'f4', 'U10', 'f4'))

# Save the table to a file
output_file = f'brown_dwarfs_{field_name}.fits'
bd_table.write(output_file, format='fits', overwrite=True)
