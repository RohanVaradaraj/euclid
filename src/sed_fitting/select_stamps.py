#!/usr/bin/env python3

"""
select_stamps.py

Chopped down version of cut_stamps.py

Lists objects from a catalogue and asks the user if they want to keep each.

Created: Tuesday 19th July 2022

"""

''' Import modules '''

import matplotlib as mpl
mpl.use('Agg')

import os
from os.path import exists
import numpy as np
from astropy.io import fits, ascii
from astropy.table import Table, Column
from astropy.coordinates import SkyCoord
import astropy.units as u
from new_catalogue_codes import cut_out, which_field, label_ct
from astropy.wcs import WCS
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import cv2
import scipy.ndimage as ndimage

''' Set-up '''

# Read catlaogue name and other things in from lephare input txt file
lephareParams = Table.read('lephare_params.txt', format='ascii.commented_header')

# Catalogue
catName = str(lephareParams['catalogue'][0])

# Field
field = str(lephareParams['field'][0])

# Get name for file from the source directory
sourceDir = str(lephareParams['dir_source'][0])
targetDir = str(lephareParams['dir_target'][0])
title = sourceDir.split('/')[-2]

# Plot dir from this dir too.
plotDir = sourceDir + '../../plots/'

# Catalogue directory
catDir = '/mnt/vardy/vardygroupshare/data/catalogues/{0}FULL/cutting/'.format(field)

# Maybe object directory
#maybeDir = '/mnt/hoy/temporaryFilesROHAN/lephare/vis_checked/maybe_{0}/'.format(field)
maybeDir = sourceDir + '../2maybe/'

# Image directory
imDir = '/mnt/vardy/vardygroupshare/data/'

# Where to save stamps
outDir = '/mnt/hoy/temporaryFilesROHAN/stamps/'

# Looping through directory
loop = True

# Switch to do visual check
interactive = True

# Desired filters
if field == 'XMM':
    reqFilters = ['YJ', 'GRI', 'HSC-G', 'HSC-R_DR3', 'HSC-I', 'HSC-NB0816_DR3', 'HSC-Z_DR3', 'HSC-NB0921_DR3', 'HSC-Y', 'Y', 'J', 'H', 'Ks', 'ch1servs', 'ch2servs']
if field == 'CDFS':
    #reqFilters = ['u', 'g', 'r', 'i', 'HSC-G', 'HSC-R', 'HSC-I', 'HSC-Z', 'Y', 'J', 'H', 'Ks', 'ch1cds', 'ch2cds']
    reqFilters = ['YJ', 'griGRI', 'u', 'g', 'HSC-G', 'r', 'i', 'HSC-I', 'HSC-Z', 'Y', 'J', 'H', 'Ks', 'ch1cds', 'ch2cds']

# Stamp size
stampSize = 10.0 # Arcseconds on each side

''' EXTRACT IMAGES '''

# Flag crosstalk artefacts. Only need to run once
#label_ct(catDir+catName, 'XMMFULL')


# Loop through objects
if loop:

    # Read in
    hdu = fits.open(catDir+catName)
    header = hdu[0].header
    cat = hdu[1].data

    for i, obj in enumerate(cat['ID'][347:]):

        print('OBJECT NUMBER: {0}'.format(i+1+347))
        print('OBJECT ID: {0}'.format(obj))

        # Interactive removal of objects
        if interactive:

            # Ask user
            keep = input("Keep this object? y, n, or m for maybe. Type stop to end the code. ")

            # If we keep the object, move to vischeck folder
            if keep == 'y':
                strID = 'Id' + str(obj).zfill(9) + '.spec'

                # Visual check after IRAC and z>6.5
#                os.system('cp /mnt/hoy/temporaryFilesROHAN/lephare/secondFitting/{0}_irac_z6.5/{1} /mnt/hoy/temporaryFilesROHAN/lephare/secondFitting/{2}_vis/'.format(field, strID, field))
                os.system('cp ' + sourceDir + strID + ' ' + targetDir)

                #cat['vis_check'][i] = 1.0
                print('Object kept')

            # Maybe objects go in a separate folder
            if keep == 'm':
                strID = 'Id' + str(obj).zfill(9) + '.spec'

                # Visual check after IRAC and z>6.5
#                os.system('cp /mnt/hoy/temporaryFilesROHAN/lephare/secondFitting/{0}_irac_z6.5/{1} /mnt/hoy/temporaryFilesROHAN/lephare/secondFitting/{2}_vis/'.format(field, strID, field))
                os.system('cp ' + sourceDir + strID + ' ' + maybeDir)

                #cat['vis_check'][i] = 2.0
                print('Object put into maybe file')

            # If we don't keep it, do nothing.
            if keep == 'n':
                print('Object not kept')

            # Exit clause
            if keep == 'stop':
                print('Stopped at object number {0}, ID {1}'.format(i+1, obj))
                exit()

    # Save final catalogue, with vis_check flag
    #if field == 'CDFS':
    #    hdu.writeto('/mnt/vardy/vardygroupshare/data/catalogues/CDFSFULL/cutting/CDFSFULL_DR2_MASKVISTADET_YJ_1.8as_IRAC2.8as_2022_04_2_vischeck_2.fits', overwrite=True)
    #if field == 'XMM':
    #    hdu.writeto('/mnt/vardy/vardygroupshare/data/catalogues/XMMFULL/XMMFULL_DR2_MASKVISTADET_YJ_1.8as_IRAC2.8as_2022_03_24_vischeck_2.fits', overwrite=True)
