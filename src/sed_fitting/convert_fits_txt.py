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

# Example usage
if __name__ == "__main__":

    # Get catalogue name from environment variable
    #cat_name = 'COSMOS_5sig_Je_3sig_J_nonDet_HSC_G_nonDet_HSC_R_nonDet_HSC_I_nonDet_HSC_Z_nonDet_HSC_Y_nonDet_Y.fits'
    cat_name = 'COSMOS_5sig_Ye_2sig_Y_nonDet_HSC_G_nonDet_HSC_R_nonDet_HSC_I.fits'


    '''SETUP'''
    cat_dir = Path.cwd().parents[1] / 'data' / 'catalogues'
    image_dir = Path.cwd().parents[3] / 'data' / 'COSMOS'

    # Directory
    out_dir = Path.home().parents[1] / 'hoy' / 'temporaryFilesROHAN' / 'lephare' / 'inputs' / 'euclid'
    #out_dir = Path.home() / 'lephare' / 'lephare_dev' / 'test'

    # Base name for output
    base_output_name = 'euclid'

    # Example filters dictionary from the first script
    filters = {
        'Ye': {'type': 'detection', 'value': 5},
        'HSC-Y_DR3': {'type': 'detection', 'value': 2},
        'HSC-G_DR3': {'type': 'non-detection', 'value': 2},
        'HSC-R_DR3': {'type': 'non-detection', 'value': 2},
        'HSC-I_DR3': {'type': 'non-detection', 'value': 2},
    }

    # Get input name
    out_name = 'det_Ye_Y.in'
    # out_name = 'hsc_nb.in'
    #out_name = 'det_NB0921_Ye.in'
    
    all_filters = ['CFHT-u', 'CFHT-g', 'CFHT-r', 'CFHT-z', 'HSC-G_DR3', 'HSC-R_DR3', 'HSC-I_DR3', 'HSC-NB0816_DR3', 'HSC-Z_DR3', 'HSC-NB0921_DR3', 'HSC-Y_DR3', 'Y', 'J', 'H', 'Ks', 'VIS', 'Ye', 'Je', 'He']

    if 'no_euclid' in out_name:
        all_filters = ['CFHT-u', 'CFHT-g', 'CFHT-r', 'CFHT-z', 'HSC-G_DR3', 'HSC-R_DR3', 'HSC-I_DR3', 'HSC-NB0816_DR3', 'HSC-Z_DR3', 'HSC-NB0921_DR3', 'HSC-Y_DR3', 'Y', 'J', 'H', 'Ks']

    print(all_filters)

    # Read input filters
    inputs = Table.read(image_dir / 'images.lis', format='ascii.commented_header')
    availFilters = np.array(inputs['Name'])

    # Read data table
    t = Table.read(cat_dir / cat_name, format='fits', hdu=1)

    '''CREATE TABLE'''

    print(t.columns)
    print(t)
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
    print(context)

    t['Context'] = context
    t['ID'] = t['ID'].astype(int)

    # Add spec_z
    fluxes.append('Context')

    # Rearrange table accordingly
    t = t[fluxes]

    # Write to the desired format: no header names, .in file.
    ascii.write(t, out_dir / out_name, format='commented_header', overwrite=True)
    print(t)
