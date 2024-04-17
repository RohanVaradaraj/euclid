#!/usr/bin/env python3

"""
convert_fits.txt.py

Script for converting fits format catalogues into text files for LePhare to read.

Created: Tuesday 19th April 2022

"""

from astropy.table import Table, Column
from astropy.io import fits, ascii
import numpy as np

'''SETUP'''

#lephareParamsStr = 'lephare_params_DR3.txt' # XMM DR3 Y+J
#lephareParamsStr = 'lephare_params.txt' # CDFS Y+J
#lephareParamsStr = 'COS_YJ_lephare_params.txt' # COSMOS Y+J

lephareParamsStr = 'jwst/beta_params.txt' # FITTING BETA SLOPE

# SED SIMULATION STUFF
lephareParamsStr = 'lephare_params_sim.txt'


# Read catlaogue name in from lephare input txt file
lephareParams = Table.read(lephareParamsStr, format='ascii.commented_header')
cat = lephareParams['catalogue']
field = str(lephareParams['field'][0])

# Catalogue directory
catDir = '/mnt/vardy/vardygroupshare/data/catalogues/' # NORMAL
catDir = '/mnt/vardy/vardygroupshare/HSC_SSP_DR3/ref_catalogues/nathan/' # SOME SIDE JWST STUFF
catDir = f'/mnt/vardy/vardygroupshare/HSC_SSP_DR3/simulations/mockCats/{field}/' # SIMULATION
catDir = f'/mnt/vardy/vardygroupshare/HSC_SSP_DR3/simulations/mockCats/{field}/UDEEP/' # SIMULATION in UDEEP

# Inputs directory from images.lis
if field == 'XMM' or 'CDFS':
    inputDir = '/mnt/vardy/vardygroupshare/data/{0}1/'.format(field)
if field == 'COSMOS':
    inputDir = '/mnt/vardy/vardygroupshare/data/{0}/'.format(field)

# Output directory
outDir = '/mnt/hoy/temporaryFilesROHAN/lephare/inputs/'

cat = str(cat[0])
print(cat)
print(field)

# Filters to be used in the SED fitting. For my z=7 work. These will be placed first in the text file.
if field == 'XMM':
    #filters = ['HSC-G_DR3', 'HSC-R_DR3', 'HSC-I_DR3', 'HSC-NB0816_DR3', 'HSC-Z_DR3', 'HSC-NB0921_DR3', 'HSC-Y_DR3', 'Y', 'J', 'H', 'Ks', 'ch1', 'ch2']
    #filters = ['HSC-G', 'HSC-R_DR3', 'HSC-I', 'HSC-NB0816_DR3', 'HSC-Z_DR3', 'HSC-NB0921_DR3', 'HSC-Y', 'Y', 'J', 'H', 'Ks'] #, 'ch1', 'ch2'] #, 'ch1servs', 'ch2servs']
    filters = ['HSC-I_DR3', 'HSC-NB0816_DR3', 'HSC-Z_DR3', 'HSC-NB0921_DR3', 'HSC-Y_DR3', 'Y', 'J', 'H', 'Ks'] #, 'ch1servs', 'ch2servs'] # Without GR for BD
if field == 'CDFS':
    #filters = ['u', 'g', 'r', 'i', 'HSC-G', 'HSC-R', 'HSC-I', 'HSC-Z', 'Y', 'J', 'H', 'Ks', 'ch1cds', 'ch2cds']
    filters = ['i', 'HSC-I', 'HSC-Z', 'Y', 'J', 'H', 'Ks'] # Without ugr for BD
if field == 'COSMOS':
#    filters = ['HSC-G_DR3', 'HSC-R_DR3', 'HSC-I_DR3', 'HSC-NB0816_DR3', 'HSC-Z_DR3', 'HSC-NB0921_DR3', 'HSC-Y_DR3', 'Y', 'J', 'H', 'Ks', 'ch1cds', 'ch2cds']
#    filters = ['u', 'HSC-G', 'HSC-R', 'HSC-I', 'HSC-Z', 'Y', 'J', 'H', 'Ks'] # SOME SIDE JWST STUFF

    filters = ['Y', 'J', 'H'] # BETA SLOPES

## SLIGHTLY DIFFERENT NAMING CONVENTION FOR SIMULATIONS
filters = [name.replace('-', '_').split('_DR3')[0] for name in filters] # SIMULATION


print(filters)

# Read input filters
inputs = Table.read(inputDir + 'images.lis', format='ascii.commented_header')
availFilters = np.array(inputs['Name'])


# Read data
if field == 'XMM' or field == 'CDFS':
    #t = Table.read(catDir + '{0}FULL/cutting/'.format(field) + cat, format='fits', hdu=1) # Normal data
    t = Table.read(catDir + cat, format='fits', hdu=1) # SIMULATION

if field == 'COSMOS':
#    t = Table.read(catDir + 'final{0}/cutting/'.format(field) + cat, format='fits', hdu=1) # NORMAL
    t = Table.read(catDir + cat, format='fits', hdu=1) # SOME SIDE JWST STUFF


'''CREATE TABLE'''

print(t.columns)
print(t)
# Create arrays for flux and error column names

# Create base arrays
fluxes = ['ID']
#fluxes = ['UID'] # SOME SIDE JWST STUFF
remainder = []

# Get the photometric redshifts
#zphots = t['Z_Phot']
zphots = t['ZPhot']

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

# Set spec_z to -99 if it is nan
#for i in range(0, len(t['z_spec'])):
#    if np.isnan(t['z_spec'][i]):
#        t['z_spec'][i] = -99. # LePhare ignores -99

# Spectroscopic redshifts
#t['z_spec'] = zphots


# Write to the desired format: no header names, .in file.
ascii.write(t, outDir+str(lephareParams['inputName'][0]), format='commented_header', overwrite = True)
print(lephareParams['inputName'][0])
print(t)
