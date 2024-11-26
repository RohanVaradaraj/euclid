#!/usr/bin/env python3

"""
chi2_sed_cuts.py

Script for applying chi2 cuts to SED fits from LePhare.
Uses the specific runs for each scenario - e.g. fewer bands for brown dwarfs, higher Av for dusty galaxies.
Dumps the best-fitting high-z, low-z and BD solutions into new folders.
These are det_{filter}_best_highz, det_{filter}_best_lowz and det_{filter}_best_bd.

Created: Wednesday 9th October 2022
"""

from pathlib import Path
import os
from astropy.io import fits, ascii
from astropy.table import Table
import matplotlib.pyplot as plt
import numpy as np
import sys
import json
from sed_fitting_codes import remove_items

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

run_type_str = '_' + run_type if run_type != '' else ''

#! Set up the directories

# Get the detection filters
detection_filters = [filter_name for filter_name, threshold in filters.items() if threshold['type'] == 'detection']
if len(detection_filters) == 0:
    stack_filters = [f for f, t in filters.items() if t['type'] == 'stacked-detection']
    detection_filters = stack_filters[0].split('+')

# Base directory for SED outputs
base_dir = Path.cwd().parents[1] / 'data' / 'sed_fitting' / 'zphot'

# Determine the base folder name
lbg_dir = 'det_' + '_'.join(detection_filters)

# Set up the other directories
bd_dir = lbg_dir + run_type_str + '_bd'
dusty_dir = lbg_dir + run_type_str + '_dusty'
lya_dir = lbg_dir + run_type_str + '_lya'

# Need the outside_footprint directory to exclude objects outside the Euclid footprint
outside_footprint_dir = base_dir / (lbg_dir + '_outside_footprint')
print(outside_footprint_dir)

# Set up the output directories where good SEDs will be stored
output_dir = Path.cwd().parents[1] / 'data' / 'sed_fitting' / 'zphot' / 'best_fits'

lbg_out_dir = lbg_dir + run_type_str + '_best_highz' # Has chi2 < 2sigma threshold
bd_out_dir = lbg_dir + run_type_str + '_best_bd' # Has chi2 < 10
dusty_out_dir = lbg_dir + run_type_str + '_best_dusty' # Has chi2_lowz < chi2_highz
lya_out_dir = lbg_dir + run_type_str + '_best_lya' # Has chi2_lya < chi2_lbg
pristine_out_dir = lbg_dir + run_type_str + '_pristine' # Satisfies chi2 < 2sigma threshold, chi2_bd < 10, chi2_highz + 4 < chi2_lowz

# Finally modify the lbg dir in place
lbg_dir += run_type_str

# Create the output directories if they don't exist
create_dir = lambda dir_name: (output_dir / dir_name).mkdir(parents=True) if not (output_dir / dir_name).exists() else None

create_dir(lbg_out_dir)
create_dir(bd_out_dir)
create_dir(dusty_out_dir)
create_dir(lya_out_dir)
create_dir(pristine_out_dir)

if overwrite:
    for dir_name in [lbg_out_dir, bd_out_dir, dusty_out_dir, lya_out_dir, pristine_out_dir]:
        for file in (output_dir / dir_name).glob('*.spec'):
            file.unlink()
        print(f'Deleted all previous .spec files in {output_dir / dir_name}')

#! Determine the number of degrees of freedom from the total number of filters

# But first, we need to know if filters are removed for brown dwarfs/non-dusty things
if bools[0] == True:
    filters_to_remove = ['CFHT-u', 'CFHT-g', 'CFHT-r', 'HSC-G_DR3', 'HSC-R_DR3', 'f277w', 'f444w', 'ch1cds', 'ch2cds']
    all_filters = remove_items(all_filters, filters_to_remove)

if bools[1] == False:
    filters_to_remove = ['f444w', 'ch1cds', 'ch2cds']
    all_filters = remove_items(all_filters, filters_to_remove)

print('All filters:', all_filters)

n = len(all_filters)
dof = n - 6 # 6 degrees of freedom for the 6 parameters in the SED fitting

# Compute the 2sigma threshold, which is the chi2 threshold for a 95% confidence interval
two_sigma_thresh = dof + 2*np.sqrt(2*dof)

print('SED fits must have chi2 < ', two_sigma_thresh, ' to be considered good fits')

#! Get the lephare output tables

# Read in the ascii tables for the SED fits
cat_dir = Path.home() / 'lephare' / 'lephare_dev' / 'test'

# Generate the output catalogue names

out_names = {
    'lbg': 'det_' + '_'.join(detection_filters) + run_type_str + '.out',
    'bd': 'det_' + '_'.join(detection_filters) + run_type_str + '_bd.out',
    'dusty': 'det_' + '_'.join(detection_filters) + run_type_str + '_dusty.out',
    'lya': 'det_' + '_'.join(detection_filters) + run_type_str + '_lya.out'
}

# Read the tables
tables = {}

for key, out_name in out_names.items():
    print(f'Reading {out_name}...')
    tables[key] = ascii.read(cat_dir / out_name, format='no_header')

lbg_table = tables['lbg']
bd_table = tables['bd']
dusty_table = tables['dusty']
lya_table = tables['lya']

#! Apply the chi2 cuts

# LBGs
condition_1 = lbg_table[:]['col2'] > 6                       # LBG solution at z>6
condition_2 = lbg_table[:]['col6'] < two_sigma_thresh        # High-z solution has chi2 < 2sigma threshold
condition_3 = lbg_table[:]['col6'] + 4 < lbg_table[:]['col15'] # delta-chi2 between high-z and low-z is more than 4
good_lbg = lbg_table[:][condition_1 & condition_2 & condition_3]

# Brown dwarfs
good_bd = bd_table[:][bd_table['col21'] < 10]              # BD solution has chi2_star < chi2_LBG

# Dusty galaxies 
condition_1 = dusty_table[:]['col2'] < 6                      # Best solution at z<6
condition_2 = dusty_table[:]['col6'] + 4 < dusty_table[:]['col15']   # Dusty solution has significantly lower chi2 than LBG solution
condition_3 = dusty_table[:]['col6'] < two_sigma_thresh       # High-z solution has chi2 < 2sigma threshold

good_dusty = dusty_table[:][condition_1 & condition_2 & condition_3]

# Lya emitters
condition_1 = lya_table[:]['col6'] + 4 < lbg_table[:]['col6']     # Lya solution has significantly lower chi2 than LBG solution
condition_2 = lya_table['col2'] > 6                           # Lya solution at z>6
condition_3 = lya_table['col6'] < two_sigma_thresh             # High-z solution has chi2 < 2sigma threshold

good_lya = lya_table[:][condition_1 & condition_2 & condition_3]

# Pristine
condition_1 = lbg_table[:]['col2'] > 6                        # LBG solution at z>6
condition_2 = lbg_table[:]['col6'] < two_sigma_thresh         # High-z solution has chi2 < 2sigma threshold
condition_3 = bd_table[:]['col21'] < 10                       # BD solution has chi2 < 10
condition_4 = dusty_table[:]['col2'] > 6                      # Dusty solution at z>6
condition_5 = dusty_table[:]['col6'] + 4 < dusty_table[:]['col15'] # Dusty solution has chi2_highz + 4 < chi2_highz

pristine = lbg_table[:][condition_1 & condition_2 & condition_3 & condition_4 & condition_5]

#! Copy the best fits to the output directories

# Get the IDs
lbg_ids = good_lbg['col1']
bd_ids = good_bd['col1']
dusty_ids = good_dusty['col1']
lya_ids = good_lya['col1']
pristine_ids = pristine['col1']

# Dictionary to map categories to their ID lists and output directories
categories = {
    'LBGs': (lbg_ids, lbg_out_dir),
    'BDs': (bd_ids, bd_out_dir),
    'Dusty galaxies': (dusty_ids, dusty_out_dir),
    'Lya emitters': (lya_ids, lya_out_dir),
    'Pristine objects': (pristine_ids, pristine_out_dir)
}

print(len(lbg_ids), 'LBGs')

orig_dir = {
    'LBGs': lbg_dir,
    'BDs': bd_dir,
    'Dusty galaxies': dusty_dir,
    'Lya emitters': lya_dir,
    'Pristine objects': lbg_dir
}

# Print the out_dirs of the categories
for category, (ids, out_dir) in categories.items():
    print(f'{category} will be copied to {output_dir / out_dir}')

# Print where files will be copied from
for category, (ids, out_dir) in categories.items():
    print(f'Files will be copied from {base_dir / orig_dir[category]}')

for category, (ids, out_dir) in categories.items():
    for ID in ids:

        file_name = f'Id{str(ID).zfill(9)}.spec'
        
        # Check if the file exists
        if not (base_dir / orig_dir[category] / file_name).exists():
            print(f'File {file_name} does not exist in {base_dir / orig_dir[category]}')
            continue

        # If the spec file is in the outside footprint directory, ignore it
        if (outside_footprint_dir / file_name).exists():
            continue

        os.system(f'cp {str(base_dir / orig_dir[category] / file_name)} {output_dir / out_dir}')
    print(f'Copied {len(ids)} {category} to {output_dir / out_dir}')
