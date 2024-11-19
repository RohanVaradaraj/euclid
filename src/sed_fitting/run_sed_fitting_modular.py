#!/usr/bin/env python3

"""
Run the SED fitting with various configurations.

Created: Friday 12th July 2024.
"""

from pathlib import Path
import subprocess
import sys
import os
from sed_fitting_codes import *
import json

# Configuration flags
config = {
    "run_type": '',  # Options: '', 'with_euclid', 'just_euclid', 'CDS', 'all_filters'
    "overwrite": False,
    "steps": {
        "selection": False,     # Initial dropout selection
        "lephare": False,       # Run LePhare. Converts the fits file into text, and builds the LePhare config file too.
        "extract_seds": False,  # Take all the good SEDs from the LePhare fitting.
        "plotting": False,      # Plot the SEDs
        "visual_selection": True, # Visual selection of SEDs
        "final_selection": False  # Final selection of SEDs with BD, dusty, lya and z>6.5 cuts.
    }
}

# Base filter sets
base_filters = {
    '': ['HSC-G_DR3', 'HSC-R_DR3', 'HSC-I_DR3', 'HSC-NB0816_DR3', 'HSC-Z_DR3', 'HSC-NB0921_DR3', 'HSC-Y_DR3', 'Y', 'J', 'H', 'Ks', 'ch1cds', 'ch2cds'],
    'with_euclid': ['HSC-G_DR3', 'HSC-R_DR3', 'HSC-I_DR3', 'HSC-NB0816_DR3', 'HSC-Z_DR3', 'HSC-NB0921_DR3', 'HSC-Y_DR3', 'VIS', 'Ye', 'Je', 'He', 'Y', 'J', 'H', 'Ks', 'ch1cds', 'ch2cds'],
    'just_euclid': ['VIS', 'Ye', 'Je', 'He'],
    'CDS': ['HSC-G_DR3', 'HSC-R_DR3', 'HSC-I_DR3', 'HSC-Z_DR3', 'VIS', 'Ye', 'Je', 'He', 'ch1cds', 'ch2cds'],
    'all_filters': ['HSC-G_DR3', 'HSC-R_DR3', 'HSC-I_DR3', 'HSC-NB0816_DR3', 'HSC-Z_DR3', 'HSC-NB0921_DR3', 'HSC-Y_DR3', 'VIS', 'Ye', 'Je', 'He', 'f115w', 'f150w', 'f277w', 'f444w', 'Y', 'J', 'H', 'Ks', 'ch1cds', 'ch2cds']
}

filters = {
    'Y+J': {'type': 'stacked-detection', 'value': 5},
    'HSC-G_DR3': {'type': 'non-detection', 'value': 2},
    'HSC-R_DR3': {'type': 'non-detection', 'value': 2},
    'HSC-I_DR3': {'type': 'non-detection', 'value': 2},
}

def run_sed_fitting(run_type, run_brown_dwarfs, run_dusty, run_lya, config):

    # Serialize filters and bool switches to a JSON string
    filters_json = json.dumps(filters)
    bools_json = json.dumps([run_brown_dwarfs, run_dusty, run_lya])
    all_filters_json = json.dumps(base_filters[run_type])
    run_type_json = json.dumps(run_type)

    print(f'Running with config: {run_type}, Brown Dwarfs: {run_brown_dwarfs}, Dusty: {run_dusty}, Lya: {run_lya}')

    # Step 1: Run selection if required
    if config["steps"]["selection"]:
        print("Running selection step...")
        selection_script = Path.cwd() / 'selection.py'
        subprocess.run(['python3', str(selection_script), filters_json], check=True)

    # Step 2: Run LePhare if required
    if config["steps"]["lephare"]:
        print("Running LePhare step...")
        # Convert catalogue to lephare format
        convert_script = Path.cwd() / 'convert_fits_txt.py'
        subprocess.run(['python3', str(convert_script), filters_json, bools_json, all_filters_json, run_type_json], check=True)

        # Build LePhare libraries
        buildLePhareLibrary(parameter_file='euclid.para', build_libs=True, build_filters=True, build_mags=True)
        
        # Run the photometric redshifts
        zphot_dir = Path.cwd().parents[1] / 'data' / 'sed_fitting' / 'zphot' / run_type
        if config['overwrite']:
            for file in zphot_dir.glob('*.spec'):
                file.unlink()
            print('Deleted all previous .spec files.')
        runPhotometricRedshifts(parameter_file='euclid.para', zphot_dir=zphot_dir)

    # Step 3: Extract good SEDs if required
    if config["steps"]["extract_seds"]:
        print("Running extract SEDs step...")
        good_seds_script = Path.cwd() / 'chi2_sed_cuts.py'
        subprocess.run(['python3', str(good_seds_script), filters_json, bools_json, all_filters_json, run_type_json], check=True)

    # Step 4: Run plotting if required
    if config["steps"]["plotting"]:
        print("Running plotting step...")
        plot_script = Path.cwd() / 'plot_SEDs.py'
        subprocess.run(['python3', str(plot_script), filters_json, bools_json, all_filters_json, run_type_json], check=True)

    # Step 5: Run visual selection if required
    if config["steps"]["visual_selection"]:
        print("Running visual selection step...")
        visual_script = Path.cwd() / 'visual_selection.py'
        subprocess.run(['python3', str(visual_script), filters_json, bools_json, all_filters_json, run_type_json], check=True)

    # Step 6: Run final selection if required
    if config["steps"]["final_selection"]:
        print("Running final selection step...")
        visual_script = Path.cwd() / 'visual_selection.py'
        subprocess.run(['python3', str(visual_script), filters_json, bools_json, all_filters_json, run_type_json], check=True)


if __name__ == "__main__":
    # List of run_types to iterate through
    run_types = [''] #, 'with_euclid', 'just_euclid', 'CDS', 'all_filters']  # Modify as needed

    # Specific combinations of flags
    flag_combinations = [
        (False, False, False),  # All False = Normal SED fitting
        #(True, False, False),   # Only run_brown_dwarfs = True
        #(False, True, False),   # Only run_dusty = True
        #(False, False, True)    # Only run_lya = True
    ]

    # Loop through all combinations of run_types and flag combinations
    for run_type in run_types:
        for run_brown_dwarfs, run_dusty, run_lya in flag_combinations:
            run_sed_fitting(run_type, run_brown_dwarfs, run_dusty, run_lya, config)
