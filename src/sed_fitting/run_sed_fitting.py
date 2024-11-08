#!/usr/bin/env python3

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

# Overwrite existing .spec files
overwrite = False

#! BOOL ARRAY
#![RUN_BROWN_DWARFS, RUN_DUSTY, RUN_LYA]

#! Filter sets

#! All filters
all_filters = ['HSC-G_DR3', 'HSC-R_DR3', 'HSC-I_DR3', 'HSC-NB0816_DR3', 'HSC-Z_DR3', 'HSC-NB0921_DR3', 'HSC-Y_DR3', 
               'Y', 'J', 'H', 'Ks', 'ch1cds', 'ch2cds']


#! VISTA-Y
filters = {
     'Y+J': {'type': 'stacked-detection', 'value': 5},
     'HSC-G_DR3': {'type': 'non-detection', 'value': 2},
     'HSC-R_DR3': {'type': 'non-detection', 'value': 2},
     'HSC-I_DR3': {'type': 'non-detection', 'value': 2},
}

#! Ye
# filters = {
#     'Ye': {'type': 'detection', 'value': 5},
#     'Je': {'type': 'detection', 'value': 2},
#     'He': {'type': 'detection', 'value': 2},``
#     # 'HSC-G_DR3': {'type': 'non-detection', 'value': 2},
#     # 'HSC-R_DR3': {'type': 'non-detection', 'value': 2},
#     # 'HSC-I_DR3': {'type': 'non-detection', 'value': 2},
#     # 'HSC-Z_DR3': {'type': 'non-detection', 'value': 2},
#     'VIS': {'type': 'non-detection', 'value': 2}, # Might interfere with e.g. z=6.9 candidates
# }

# #! Je
# filters = {
#     'Je': {'type': 'detection', 'value': 5},
#     'He'  : {'type': 'detection', 'value': 2},
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

    #! ############## Run preparation scripts ##############

    # #! Selection
    # print("Running selection.py...")
    # selection_script = Path.cwd() / 'selection.py'
    # subprocess.run(['python3', str(selection_script), filters_json], check=True)

    # #! Convert catalogue to lephare format
    # print("Running convert_fits_txt.py...")
    # convert_script = Path.cwd() / 'convert_fits_txt.py'
    # subprocess.run(['python3', str(convert_script), filters_json, bools_json, all_filters_json], check=True)

    # #! Generate lephare config file
    # det_filters = [f for f, t in filters.items() if t['type'] == 'detection']
    # if len(det_filters) == 0:
    #     stack_filters = [f for f, t in filters.items() if t['type'] == 'stacked-detection']
    #     det_filters = stack_filters[0].split('+')
    # nondet_filters = [f for f, t in filters.items() if t['type'] == 'non-detection']

    # GenerateLePhareConfig(all_filters, det_filters, run_brown_dwarfs, run_dusty, run_lya)

    # ############! Run LePhare ############

    # #! Set up the correct folders for outputting .spec files
    # if run_brown_dwarfs:
    #     det_folder = 'det_' + '_'.join(det_filters) + '_bd'
    # if run_lya:
    #     det_folder = 'det_' + '_'.join(det_filters) + '_lya'
    # if run_dusty:
    #     det_folder = 'det_' + '_'.join(det_filters) + '_dusty'
    # if not run_brown_dwarfs and not run_lya and not run_dusty:
    #     det_folder = 'det_' + '_'.join(det_filters)

    # #! Create directories if necessary
    # zphot_dir = Path.cwd().parents[1] / 'data' / 'sed_fitting' / 'zphot' / det_folder
    # if not zphot_dir.exists():
    #     zphot_dir.mkdir(parents=True)

    # #! Build the LePhare libraries
    # buildLePhareLibrary(parameter_file='euclid.para', build_libs=True, build_filters=True, build_mags=True)

    # #! Run the photometric redshifts
    # if overwrite:
    #     for file in zphot_dir.glob('*.spec'):
    #         file.unlink()
    #     print('Deleted all previous .spec files.')

    # runPhotometricRedshifts(parameter_file='euclid.para', zphot_dir=zphot_dir)

    #! Extract good SED fits once all of the above have run.
    # good_seds_script = Path.cwd() / 'chi2_sed_cuts.py'
    # subprocess.run(['python3', str(good_seds_script), filters_json, bools_json, all_filters_json], check=True)

    #! Run the plotting code
    plot_script = Path.cwd() / 'plot_SEDs.py'
    subprocess.run(['python3', str(plot_script), filters_json, bools_json, all_filters_json], check=True)

    #! Run the visual selection code
    # visual_script = Path.cwd() / 'visual_selection.py'
    # subprocess.run(['python3', str(visual_script), filters_json, bools_json, all_filters_json], check=True)

    #! Run the dusty selection code
    # dusty_script = Path.cwd() / 'dusty_selection.py'
    # subprocess.run(['python3', str(dusty_script), filters_json, bools_json, all_filters_json], check=True)

    #! Run the brown dwarf selection code
    # brown_dwarf_script = Path.cwd() / 'brown_dwarf_selection.py'
    # subprocess.run(['python3', str(brown_dwarf_script), filters_json, bools_json, all_filters_json], check=True)

    # #! Run the Lyman-alpha and redshift selection code
    # lya_script = Path.cwd() / 'lya_and_redshift_selection.py'
    # subprocess.run(['python3', str(lya_script), filters_json, bools_json, all_filters_json], check=True)


    return None

if __name__ == "__main__":

    # #! First run the normal sed fitting
    # run_brown_dwarfs = False
    # run_dusty = False
    # run_lya = False

    # run_sed_fitting()

    # #! Then run brown dwarfs
    # run_brown_dwarfs = True
    # run_dusty = False
    # run_lya = False

    # run_sed_fitting()

    #! Then run dusty galaxies
    # run_brown_dwarfs = False
    # run_dusty = True
    # run_lya = False

    # run_sed_fitting()

    #! Then run Lyman-alpha emitters
    # run_brown_dwarfs = False
    # run_dusty = False
    # run_lya = True

    run_sed_fitting()

    # ! Extract good SED fits once all of the above have run.
    # good_seds_script = Path.cwd() / 'chi2_sed_cuts.py'
    # subprocess.run(['python3', str(good_seds_script), filters_json, bools_json, all_filters_json], check=True)

    # #! Run the plotting code
    # plot_script = Path.cwd() / 'plot_SEDs.py'
    # subprocess.run(['python3', str(plot_script), filters_json, bools_json, all_filters_json], check=True)