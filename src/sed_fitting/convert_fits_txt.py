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

'''SETUP'''

cat_dir = Path.cwd().parents[1] /'data' / 'catalogues'
cat_name = 'COSMOS_5sig_Ye_2sig_VISTA_Y_nonDet_HSC_G_nonDet_HSC_R_nonDet_HSC_I.fits'

image_dir = Path.cwd().parents[3] / 'data' / 'COSMOS'

# Driec
out_dir = Path.home().parents[1] / 'hoy' / 'temporaryFilesROHAN' / 'lephare' / 'inputs' / 'euclid'
out_name = 'euclid_test.in'

filters = ['CFHT-u', 'CFHT-g', 'CFHT-r', 'CFHT-iy', 'CFHT-z', 'HSC-G_DR3', 'HSC-R_DR3', 'HSC-I_DR3', 'HSC-NB0816_DR3', 'HSC-Z_DR3', 'HSC-NB0921_DR3', 'HSC-Y_DR3', 'Y', 'J', 'H', 'Ks', 'VIS', 'Ye', 'Je', 'He'] 

print(filters)

# Read input filters
inputs = Table.read(image_dir / 'images.lis', format='ascii.commented_header')
availFilters = np.array(inputs['Name'])

# Read data table
t = Table.read(cat_dir / cat_name, format='fits', hdu=1) # SIMULATION

'''CREATE TABLE'''

print(t.columns)
print(t)
# Create arrays for flux and error column names

# Create base arrays
fluxes = ['ID']
remainder = []

# Loop and add desired filters to start of list
for i, filterName in enumerate(filters):

    fluxes = fluxes + ['flux_{0}'.format(filterName)]

    fluxes = fluxes + ['err_{0}'.format(filterName)]

# Find unused filters
remainder = list(set(availFilters) - set(filters))

# Append unused filters to the end
'''
for i, filterName in enumerate(remainder):

    fluxes = fluxes + ['flux_{0}'.format(filterName)]

    fluxes = fluxes + ['err_{0}'.format(filterName)]
'''
# Find the context based on filters used, and their order.

# How many filters
n = len(filters)

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
fluxes = fluxes + ['Context']

# Rearrange table accordingly
t = t[fluxes] #, z_spec]

# Write to the desired format: no header names, .in file.
ascii.write(t, out_dir / out_name, format='commented_header', overwrite = True)
print(t)
