#!/usr/bin/env python3

"""
convert_fits.txt.py

Script for converting fits format catalogues into text files for LePhare to read.

Created: Tuesday 19th April 2022

"""

from astropy.table import Table, Column
from astropy.io import fits, ascii
import numpy as np
from pathlib import Path
import os
import sys
import json
from selection import generate_selection_name, generate_input_name
from sed_fitting_codes import remove_items

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
    run_all_json = sys.argv[6]
    run_all_objects = json.loads(run_all_json)
    custom_name = sys.argv[7] if len(sys.argv) > 7 else None
    custom_name = json.loads(custom_name) if custom_name else None

#! BOOL ARRAY
#![RUN_BROWN_DWARFS, RUN_DUSTY, RUN_LYA]

#????? PAPER CORRECTION: CHECK WHEN SPITZER ERROR FLOOR IS 5% OF FLUX
spitzer_five_percent = False

#???? PAPER CORRECTION: CUSTOM CATALOGUES. WANT TO RERUN U+E SAMPLE WITH ONLY EUCLID PHOTOMETRY
custom_cat = False
#custom_cat_name = 'Euclid_UltraVISTA_z7_sample.fits'
custom_cat_name = 'COSMOS_5sig_Y_J_nonDet_HSC_G_nonDet_HSC_R_nonDet_HSC_I_with_irac_handbook.fits'

# Example usage
if __name__ == "__main__":

    # Generate catalogue name from input filters
    cat_name = generate_selection_name(field_name, filters)

    '''SETUP'''
    cat_dir = Path.cwd().parents[1] / 'data' / 'catalogues'

    if field_name == 'COSMOS':
        image_dir = Path.cwd().parents[3] / 'data' / field_name
    else:
        #? We need image_dir to get images.lis to get filter names. So just use one of the XMM/CDFS tiles.
        image_dir = Path.cwd().parents[3] / 'data' / (field_name+'1')

    # Directory
    out_dir = Path.home().parents[1] / 'hoy' / 'temporaryFilesROHAN' / 'lephare' / 'inputs' / 'euclid'

    # Base name for output
    base_output_name = 'euclid'

    # Get input name
    out_name = generate_input_name(filters, run_type, *bools)
    if run_all_objects:
        out_name = custom_name + '.in'
    print(out_name)

    # Delete blue filters if we are running fitting for brown dwarfs
    if bools[0] == True:
        filters_to_remove = ['HSC-G_DR3', 'HSC-R_DR3', 'f277w', 'f444w', 'ch1cds', 'ch2cds']
        all_filters = remove_items(all_filters, filters_to_remove)
        if field_name == 'CDFS':
            filters_to_remove = ['HSC-G', 'HSC-R', 'u', 'g', 'r']
            all_filters = remove_items(all_filters, filters_to_remove)
        print('Running brown dwarfs: blue filters and long-wavelength filters removed in input catalogue.')

    # Delete reddest filters if we are NOT running dusty galaxies
    if bools[1] == False:
        if field_name != 'XMM':
            filters_to_remove = ['f444w', 'ch1cds', 'ch2cds']
        if field_name == 'XMM':
            filters_to_remove = ['f444w', 'ch1servs', 'ch2servs']
        all_filters = remove_items(all_filters, filters_to_remove)
        print('Not running dusty galaxies: reddest filters removed in input catalogue.')

    # Read input filters
    inputs = Table.read(image_dir / 'images.lis', format='ascii.commented_header')
    availFilters = np.array(inputs['Name'])

    # Read data table
    if not run_all_objects:
        #????? SETTING SPITZER ERROR FLOOR TO 5% OF FLUX
        if spitzer_five_percent:
            cat_name = cat_name.replace('.fits', '_5percent_IRACfloor.fits')
        if custom_cat:
            cat_name = custom_cat_name
            cat_dir = cat_dir #/ 'candidates'

        print(f'Reading catalogue: {cat_dir / cat_name}')
        t = Table.read(cat_dir / cat_name, format='fits', hdu=1)
    if run_all_objects:
        cat_names = {
            'XMM': Path.cwd().parents[3] / 'data' / 'catalogues' / 'XMMFULL' / 'XMMFULL_DR3_MASKVISTADET_HSC-Z_DR3_2.0as_IRAC2.8as_2024_01_18.fits',
            'COSMOS': Path.cwd().parents[3] / 'data' / 'catalogues' / 'finalCOSMOS' / 'other' / 'COSMOSFULL_DR3_MASKVISTADET_HSC-Z_DR3_2025_06_05_1.8as_IRAC_2.8as_ALL.fits',
        }
        t = Table.read(cat_names[field_name], format='fits', hdu=1)

    # Limit to where JWST errors are positive
    #t = t[t['err_f115w'] > 0]

    '''CREATE TABLE'''

    # Create arrays for flux and error column names

    # Create base arrays
    fluxes = ['ID']
    remainder = []

    # Loop and add desired filters to start of list
    for filterName in all_filters:
        fluxes.append(f'flux_{filterName}')
        fluxes.append(f'err_{filterName}')

    # Find unused filters
    remainder = list(set(availFilters) - set(all_filters))

    # Uncomment to append unused filters to the end
    '''
    for filterName in remainder:
        fluxes.append(f'flux_{filterName}')
        fluxes.append(f'err_{filterName}')
    '''

    # Find the context based on filters used, and their order.
    n = len(all_filters)

    # Begin sum at zero
    sum = 0

    # Add the powers
    for i in range(0, n):
        sum = sum + 2**i

    # Define the result as the context
    context = sum

    t['Context'] = context
    t['ID'] = t['ID'].astype(int)

    # Add spec_z
    fluxes.append('Context')

    # Rearrange table accordingly
    t = t[fluxes]

    # Write to the desired format: no header names, .in file.
    ascii.write(t, out_dir / out_name, format='commented_header', overwrite=True)
    print(t)
    print(f'Context: {context}')
