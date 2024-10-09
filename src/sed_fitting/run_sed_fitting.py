#!/bin/bash
"""
Run the SED fitting.

Created: Friday 12th July 2024.
"""

from pathlib import Path
import subprocess
import sys
import os
from sed_fitting_codes import *
import json

#! SETUP AND SWITCHES

# Run brown dwarf fits only with filters that cover the template wavelength range
run_brown_dwarfs = False

# Run dusty galaxies by including long-wavelength filters and with higher Av in the models
run_dusty = False

# Run Lyman-alpha emitter fits by using modified BC03 templates
run_lya = False

# BOOL ARRAY
#[RUN_BROWN_DWARFS, RUN_LYA, RUN_DUSTY]

#! Filter sets

#! All filters
all_filters = ['HSC-G_DR3', 'HSC-R_DR3', 'HSC-I_DR3', 'HSC-NB0816_DR3', 'HSC-Z_DR3', 'HSC-NB0921_DR3', 'HSC-Y_DR3', 
                'Y', 'J', 'H', 'Ks', 
                'f115w', 'f150w', 'f277w' , 'f444w', 
                'VIS', 'Ye', 'Je', 'He']  


#! Ye
filters = {
    'Ye': {'type': 'detection', 'value': 5},
    'Je': {'type': 'detection', 'value': 2},
    'He': {'type': 'detection', 'value': 2},
    'HSC-G_DR3': {'type': 'non-detection', 'value': 2},
    'HSC-R_DR3': {'type': 'non-detection', 'value': 2},
    'HSC-I_DR3': {'type': 'non-detection', 'value': 2},
    'HSC-Z_DR3': {'type': 'non-detection', 'value': 2},
    #'VIS': {'type': 'non-detection', 'value': 2}, # Might interfere with e.g. z=6.9 candidates
}

#! Je
# filters = {
#     'Je': {'type': 'detection', 'value': 5},
#     'HSC-G_DR3': {'type': 'non-detection', 'value': 2},
#     'HSC-R_DR3': {'type': 'non-detection', 'value': 2},
#     'HSC-I_DR3': {'type': 'non-detection', 'value': 2},
#     'HSC-Z_DR3': {'type': 'non-detection', 'value': 2},
#     'HSC-Y_DR3': {'type': 'non-detection', 'value': 2},
#     'Y': {'type': 'non-detection', 'value': 2},
#     'VIS': {'type': 'non-detection', 'value': 2},
# }

#! He
# filters = {
#     'He': {'type': 'detection', 'value': 5},
#     'HSC-G_DR3': {'type': 'non-detection', 'value': 2},
#     'HSC-R_DR3': {'type': 'non-detection', 'value': 2},
#     'HSC-I_DR3': {'type': 'non-detection', 'value': 2},
#     'HSC-Z_DR3': {'type': 'non-detection', 'value': 2},
#     'HSC-Y_DR3': {'type': 'non-detection', 'value': 2},
#     'Y': {'type': 'non-detection', 'value': 2},
#     'J': {'type': 'non-detection', 'value': 2},
#     'VIS': {'type': 'non-detection', 'value': 2},
#     'f115w': {'type': 'non-detection', 'value': 2},
# }

def run_sed_fitting():

    # Serialize filters and bool switches to a JSON string
    filters_json = json.dumps(filters)
    bools_json = json.dumps([run_brown_dwarfs, run_dusty, run_lya])
    all_filters_json = json.dumps(all_filters)

    selection_script = Path.cwd() / 'selection.py'
    convert_script = Path.cwd() / 'convert_fits_txt.py'

    #! Selection
    print("Running selection.py...")
    subprocess.run(['python3', str(selection_script), filters_json], check=True)

    #! Convert catalogue to lephare format
    print("Running convert_fits_txt.py...")
    subprocess.run(['python3', str(convert_script), filters_json, bools_json, all_filters_json], check=True)

    #! Generate lephare config file
    det_filters = [f for f, t in filters.items() if t['type'] == 'detection']
    nondet_filters = [f for f, t in filters.items() if t['type'] == 'non-detection']

    GenerateLePhareConfig(all_filters, run_brown_dwarfs, run_dusty, run_lya)

    #! Run LePhare

    #! Set up the correct folders for outputting .spec files
    if run_brown_dwarfs:
        det_folder = 'det_' + '_'.join(det_filters) + '_bd'
    if run_lya:
        det_folder = 'det_' + '_'.join(det_filters) + '_lya'
    if run_dusty:
        det_folder = 'det_' + '_'.join(det_filters) + '_dusty'
    else:
        det_folder = 'det_' + '_'.join(det_filters)

        #! Create directories if necessary
        zphot_dir = Path.cwd().parents[1] / 'data' / 'sed_fitting' / 'zphot' / det_folder
        if not zphot_dir.exists():
            zphot_dir.mkdir(parents=True)

        #! Build the LePhare libraries
        buildLePhareLibrary(parameter_file='euclid.para', build_libs=True, build_filters=True, build_mags=True)

        #! Run the photometric redshifts
        runPhotometricRedshifts(parameter_file='euclid.para', zphot_dir=zphot_dir, overwrite=True)

    return None

if __name__ == "__main__":
    run_sed_fitting()