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

#! Configuration flags. Best to run steps one at a time.
config = {
    "run_type": '',                 #? Options: '', 'with_euclid', 'just_euclid', 'CDS', 'all_filters'
    "overwrite": True,
    "steps": {
        "selection": False,         #? Initial dropout selection
        "lephare": False,            #? Run LePhare. Converts the fits file into text, and builds the LePhare config file too.
        "extract_seds": False,      #? Take all the good SEDs from the LePhare fitting.
        "plotting": True,          #? Plot the SEDs
        "visual_selection": False,  #? Visual selection of SEDs
        "final_selection": False   #? Final selection of SEDs with BD, dusty, lya and z>6.5 cuts.
    }
}

#! Different run types 
run_types = ['', 'with_euclid', 'just_euclid', 'CDS', 'all_filters']

#! Whether to run the masking of the data to the euclid footprint
mask_euclid = False

#! Specific combinations of flags.
flag_combinations = [
    (False, False, False),  #? All False = Normal SED fitting
    #(True, False, False),   #? Only run_brown_dwarfs = True
    #(False, True, False),   #? Only run_dusty = True
    #(False, False, True)    #? Only run_lya = True
]

#! IF LOOPING RUN TYPES
# ?If we want to run more than one run_type, we can loop through the run_types list
loop_run_types = False

#! IF PLOTTING:
#? Define the type of object to plot, which goes into the SED code to name the PDF and find the correct folder
#? E.g. in rohan/euclid/data/sed_fitting/zphot/best_fits/, if your desired folder is det_Y_J_with_euclid_z7, below is 'z7'
plot_object_type = 'z7' # 'BD_PLUS_EUCLID_PHOT' #'best_highz' # 'best_bd'

#! Base filter sets
base_filters = {
    '': ['HSC-G_DR3', 'HSC-R_DR3', 'HSC-I_DR3', 'HSC-NB0816_DR3', 'HSC-Z_DR3', 'HSC-NB0921_DR3', 'HSC-Y_DR3', 'Y', 'J', 'H', 'Ks', 'ch1cds', 'ch2cds'],
    'with_euclid': ['HSC-G_DR3', 'HSC-R_DR3', 'HSC-I_DR3', 'HSC-NB0816_DR3', 'HSC-Z_DR3', 'HSC-NB0921_DR3', 'HSC-Y_DR3', 'VIS', 'Ye', 'Je', 'He', 'Y', 'J', 'H', 'Ks', 'ch1cds', 'ch2cds'],
    'just_euclid': ['VIS', 'Ye', 'Je', 'He'],
    'CDS': ['HSC-G_DR3', 'HSC-R_DR3', 'HSC-I_DR3', 'HSC-Z_DR3', 'VIS', 'Ye', 'Je', 'He', 'ch1cds', 'ch2cds'],
    'all_filters': ['HSC-G_DR3', 'HSC-R_DR3', 'HSC-I_DR3', 'HSC-NB0816_DR3', 'HSC-Z_DR3', 'HSC-NB0921_DR3', 'HSC-Y_DR3', 'VIS', 'Ye', 'Je', 'He', 'f115w', 'f150w', 'f277w', 'f444w', 'Y', 'J', 'H', 'Ks', 'ch1cds', 'ch2cds']
}

#! Dropout se;ection filters
filters = {
    'Y+J': {'type': 'stacked-detection', 'value': 5},
    'HSC-G_DR3': {'type': 'non-detection', 'value': 2},
    'HSC-R_DR3': {'type': 'non-detection', 'value': 2},
    'HSC-I_DR3': {'type': 'non-detection', 'value': 2},
}

def run_sed_fitting(run_type, run_brown_dwarfs, run_dusty, run_lya, config):

    """
    Run the SED fitting with the given configuration.

    Parameters
    ----------

    run_type : str
        The type of run to perform. Options are '', 'with_euclid', 'just_euclid', 'CDS', 'all_filters'.
        '' is the default run with the HSC + VISTA filters.
        'with_euclid' is the run with HSC+VISTA+Euclid.
        'just_euclid' is the run with just Euclid filters.
        'CDS' is the run with HSC+Euclid, simulating the Euclid Deep Field North.
        'all_filters' is the run with all filters, including JWST.

    run_brown_dwarfs : bool
        Whether to run the brown dwarf selection. If so, the bluest filters are removed since there is no information for the BD SED.

    run_dusty : bool
        Whether to run the dusty selection. If so, Spitzer/IRAC and the reddest JWST bands are included.

    run_lya : bool
        Whether to run the Lya selection. If so, runs modified BC03 models with Lya emission.

    config : dict
        The configuration dictionary outlining the steps to run. Defined at top of this code.

    Returns
    -------
    None

    """

    # Serialize filters and bool switches to a JSON string
    filters_json = json.dumps(filters)
    bools_json = json.dumps([run_brown_dwarfs, run_dusty, run_lya])
    all_filters_json = json.dumps(base_filters[run_type])
    run_type_json = json.dumps(run_type)
    object_type_json = json.dumps(plot_object_type)

    run_flag_dict = {
        (False, False, False): 'normal SED fitting',
        (True, False, False): 'brown dwarf fitting',
        (False, True, False): 'low-redshift dusty galaxy fitting',
        (False, False, True): 'lyman-alpha emitter fitting'
    }

    print('------------------------------------------------------------------------------------------')
    print(f'Running with config: {run_type}, Brown Dwarfs: {run_brown_dwarfs}, Dusty: {run_dusty}, Lya: {run_lya}. So running {run_flag_dict[(run_brown_dwarfs, run_dusty, run_lya)]}')
    print('------------------------------------------------------------------------------------------')



    #! -----------------------------
    #! Step 1: Run dropout selection 
    #! -----------------------------
    if config["steps"]["selection"]:

        print("Running selection step...")

        selection_script = Path.cwd() / 'selection.py'
        subprocess.run(['python3', str(selection_script), filters_json], check=True)



    #! -----------------------------
    #! Step 2: Run LePhare 
    #! -----------------------------
    if config["steps"]["lephare"]:

        print("Running LePhare step...")

        #? Convert catalogue to lephare format
        convert_script = Path.cwd() / 'convert_fits_txt.py'
        subprocess.run(['python3', str(convert_script), filters_json, bools_json, all_filters_json, run_type_json], check=True)

        #? Set up the output file names and folders
        #* 1) Get the detection and non-detection filters
        det_filters = [f for f, t in filters.items() if t['type'] in {'detection', 'stacked-detection'}]

        if not det_filters:
            det_filters = filters[next(f for f, t in filters.items() if t['type'] == 'stacked-detection')].split('+')

        nondet_filters = [f for f, t in filters.items() if t['type'] == 'non-detection']

        # Replace '+' with '_' in detection filters
        formatted_det_filters = [f.replace('+', '_') for f in det_filters]

        #* 2) Set up the correct folders for outputting .spec files
        tags = [
            '_'.join(formatted_det_filters) if formatted_det_filters else '',
            run_type,
            'bd' if run_brown_dwarfs else '',
            'lya' if run_lya else '',
            'dusty' if run_dusty else ''
        ]

        det_folder = 'det_' + '_'.join(filter(None, tags))

        #? Generate LePhare parameter file
        GenerateLePhareConfig(base_filters[run_type], det_filters, run_type=run_type, run_brown_dwarfs=run_brown_dwarfs, run_dusty=run_dusty, run_lya=run_lya)

        #? Build LePhare libraries
        buildLePhareLibrary(parameter_file='euclid.para', build_libs=True, build_filters=True, build_mags=True)
        
        #? Run the photometric redshifts
        zphot_dir = Path.cwd().parents[1] / 'data' / 'sed_fitting' / 'zphot' / det_folder

        # With an overwrite step
        if config['overwrite']:
            for file in zphot_dir.glob('*.spec'):
                file.unlink()
            print('Deleted all previous .spec files.')

        runPhotometricRedshifts(parameter_file='euclid.para', zphot_dir=zphot_dir)



    #! -----------------------------
    #! Step 3: Extract good SEDs 
    #! -----------------------------
    if config["steps"]["extract_seds"]:

        #? For VISTA+Euclid, get the good SEDs within the Euclid footprint. Don't need to run chi2_sed_cuts.py on top.
        #if run_type == 'with_euclid':

        if mask_euclid:
            print('Extracting good SEDs that are within the Euclid footprint...')

            mask_script = Path.cwd() / 'mask_euclid_pointing.py'
            subprocess.run(['python3', str(mask_script), filters_json, bools_json, all_filters_json, run_type_json], check=True)

        #? Otherwise, Extract the SEDs from the chi2 etc.
        else:

            print("Running extract SEDs step...")

            good_seds_script = Path.cwd() / 'chi2_sed_cuts.py'
            subprocess.run(['python3', str(good_seds_script), filters_json, bools_json, all_filters_json, run_type_json], check=True)



    #! -----------------------------
    #! Step 4: Run plotting 
    #! -----------------------------
    if config["steps"]["plotting"]:

        print("Running plotting step...")

        plot_script = Path.cwd() / 'plot_SEDs.py'
        subprocess.run(['python3', str(plot_script), filters_json, bools_json, all_filters_json, run_type_json, object_type_json], check=True)



    #! -----------------------------
    #! Step 5: Run visual selection
    #! -----------------------------
    if config["steps"]["visual_selection"]:

        print("Running visual selection step...")

        visual_script = Path.cwd() / 'visual_selection.py'
        subprocess.run(['python3', str(visual_script), filters_json, bools_json, all_filters_json, run_type_json], check=True)



    #! -----------------------------
    #! Step 6: Run final selection 
    #! -----------------------------
    if config["steps"]["final_selection"]:
        print("Running final selection step...")

        dusty_script = Path.cwd() / 'dusty_selection.py'
        bd_script = Path.cwd() / 'brown_dwarf_selection.py'
        lya_script = Path.cwd() / 'lya_and_redshift_selection.py'

        print('Running dusty selection...')
        subprocess.run(['python3', str(dusty_script), filters_json, bools_json, all_filters_json, run_type_json], check=True)
        print('Running brown dwarf selection...')
        subprocess.run(['python3', str(bd_script), filters_json, bools_json, all_filters_json, run_type_json], check=True)
        print('Running Lya and z>6.5 selection...')
        subprocess.run(['python3', str(lya_script), filters_json, bools_json, all_filters_json, run_type_json], check=True)

    return None



if __name__ == "__main__":

    # List of run_types to iterate through
    if not loop_run_types:
        run_types = [config["run_type"]]

    # Loop through all combinations of run_types and flag combinations
    for run_type in run_types:
        for run_brown_dwarfs, run_dusty, run_lya in flag_combinations:
            run_sed_fitting(run_type, run_brown_dwarfs, run_dusty, run_lya, config)
