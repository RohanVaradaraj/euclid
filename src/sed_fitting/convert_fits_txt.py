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

#! BOOL ARRAY
#![RUN_BROWN_DWARFS, RUN_DUSTY, RUN_LYA]

# Example usage
if __name__ == "__main__":

    # Generate catalogue name from input filters
    cat_name = generate_selection_name('COSMOS', filters)


    '''SETUP'''
    cat_dir = Path.cwd().parents[1] / 'data' / 'catalogues'
    image_dir = Path.cwd().parents[3] / 'data' / 'COSMOS'

    # Directory
    out_dir = Path.home().parents[1] / 'hoy' / 'temporaryFilesROHAN' / 'lephare' / 'inputs' / 'euclid'

    # Base name for output
    base_output_name = 'euclid'

    # Get input name
    out_name = generate_input_name(filters, run_type, *bools)
    print(out_name)

    # Delete blue filters if we are running fitting for brown dwarfs
    if bools[0] == True:
        filters_to_remove = ['CFHT-u', 'CFHT-g', 'CFHT-r', 'HSC-G_DR3', 'HSC-R_DR3', 'f277w', 'f444w', 'ch1cds', 'ch2cds']
        all_filters = remove_items(all_filters, filters_to_remove)
        print('Running brown dwarfs: blue filters and long-wavelength filters removed in input catalogue.')

    # Delete reddest filters if we are NOT running dusty galaxies
    if bools[1] == False:
        filters_to_remove = ['f444w', 'ch1cds', 'ch2cds']
        all_filters = remove_items(all_filters, filters_to_remove)
        print('Not running dusty galaxies: reddest filters removed in input catalogue.')

    # Read input filters
    inputs = Table.read(image_dir / 'images.lis', format='ascii.commented_header')
    availFilters = np.array(inputs['Name'])

    # Read data table
    t = Table.read(cat_dir / cat_name, format='fits', hdu=1)

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
