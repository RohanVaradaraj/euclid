#!/usr/bin/env python3

"""
stack_images.py

Create an image stack in Euclid.

Created: Tueday 11th June 2024.

"""

# Import modules ###############################################

import numpy as np
import os
import astropy.io.fits as fits
import gc
from pathlib import Path

# Setup ########################################################

# Image dir for UltraVISTA
imageDir = Path.cwd().parents[3] / 'data' / 'euclid' / 'images'

stackDir = Path.home() / 'euclid'

filter1 = 'Y' # Filters to be stacked
filter2 = 'J'
filter3 = 'H' # Set filter3 to 'None' if you only want to stack 2.

filters = [filter1, filter2, filter3]

fields = ['COSMOS']

print('Stacking {0} in {1}'.format(filters, fields))

# Loop through fields ##########################################

for i, fieldName in enumerate(fields):

    stackDir = stackDir / fieldName.upper()

    print('Running loop. Stacking {0} images in {1}'.format(filters, fieldName))

    # Set up images
    image1 = imageDir / f'{fieldName}_{filter1}_MOSAIC.fits'
    image2 = imageDir / f'{fieldName}_{filter2}_MOSAIC.fits'
    if filter3 != 'None':
        image3 = imageDir / f'{fieldName}_{filter3}_MOSAIC.fits'

    wht1 = imageDir / f'{fieldName}_{filter1}_MOSAIC_WHT.fits'
    wht2 = imageDir / f'{fieldName}_{filter2}_MOSAIC_WHT.fits'
    if filter3 != 'None':
        wht3 = imageDir / f'{fieldName}_{filter3}_MOSAIC_WHT.fits'

    print('First image: ', image1)
    print('With weight: ', wht1)
    print('Second image: ', image2)
    print('With weight: ', wht2)
    if filter3 != 'None':
        print('Third image: ', image3)
        print('With weight: ', wht3)


    # Load images in a memory safe manner
    with fits.open(image1) as imHDU1:
        image1 = imHDU1[0].data
        header1 = imHDU1[0].header
    with fits.open(image2) as imHDU2:
        image2 = imHDU2[0].data
        header2 = imHDU2[0].header
    with fits.open(wht1) as whtHDU1:
        wht1 = whtHDU1[0].data
        wheader1 = whtHDU1[0].header
    with fits.open(wht2) as whtHDU2:
        wht2 = whtHDU2[0].data
        wheader2 = whtHDU2[0].header

    print('##########Check header of image 1!#############')
    print(header1)

    # Extract third image if stacking
    if filter3 != 'None':
        with fits.open(image3) as imHDU3:
            imHDU3 = fits.open(image3, memmap=True)
            image3 = imHDU3[0].data
            header3 = imHDU3[0].header

        with fits.open(wht3) as whtHDU3:
            whtHDU3 = fits.open(wht3, memmap=True)
            wht3 = whtHDU3[0].data
            wheader3 = whtHDU3[0].header

    if filter3 != 'None':
        print('Computing stack of {0}, {1}, {2}'.format(filter1, filter2, filter3))

        print('Computing sum of the weights...')
        weightSum = wht1 + wht2 + wht3

        print('Computing weighted sum...')
        w1_1 = wht1 * image1
        w2_2 = wht2 * image2
        w3_3 = wht3 * image3

        finalImage = (w1_1 + w2_2 + w3_3) / weightSum

    if filter3 == 'None':
        print('Computing stack of {0}, {1}'.format(filter1, filter2))

        print('Computing sum of the weights...')
        weightSum = wht1 + wht2

        print('Computing weighted sum...')
        w1_1 = wht1 * image1

    #! Save image
    print('########Saving image.#########')

    if filter3 == 'None':
        fits.writeto(stackDir / fieldName.upper() / '/{0}_{1}{2}_STACK.fits'.format(fieldName.upper(), name1, name2), finalImage, header1, overwrite=True)
    if filter3 != 'None':
        fits.writeto(stackDir / fieldName.upper() / '/{0}_{1}{2}{3}_STACK.fits'.format(fieldName.upper(), filter1, filter2, filter3), finalImage, header1, overwrite=True)

    print('###########Saved image to ', stackDir / fieldName.upper())

    #! Save weight
    print('########Saving weight#########')

    if filter3 == 'None':
        fits.writeto(stackDir / fieldName.upper() / '/{0}_{1}{2}_STACK_WHT.fits'.format(fieldName.upper(), name1, name2), weightSum, wheader1, overwrite=True)
    if filter3 != 'None':
        fits.writeto(stackDir / fieldName.upper() / '/{0}_{1}{2}{3}_STACK_WHT.fits'.format(fieldName.upper(), filter1, filter2, filter3), weightSum, wheader1, overwrite=True)

    print('###########Saved weight to ', stackDir / fieldName.upper())

