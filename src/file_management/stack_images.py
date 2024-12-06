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

# Image dir for Euclid
#imageDir = Path.cwd().parents[3] / 'data' / 'euclid' / 'images'

# Stack dir for regular COSMOS directory
imageDir = Path.cwd().parents[3] / 'data' / 'COSMOS'

stackDir = imageDir

filter1 = 'Y' # Filters to be stacked
filter2 = 'J'
filter3 = 'None' # Set filter3 to 'None' if you only want to stack 2.
filter4 = 'None'

filters = [filter1, filter2, filter3, filter4]

fields = ['COSMOS']

print('Stacking {0} in {1}'.format(filters, fields))

# Loop through fields ##########################################

for i, fieldName in enumerate(fields):

    print('Running loop. Stacking {0} images in {1}'.format(filters, fieldName))

    # Set up images
    #image1 = imageDir / f'{fieldName}_{filter1}_MOSAIC.fits'
    #image2 = imageDir / f'{fieldName}_{filter2}_MOSAIC.fits'
    # image1 = imageDir / f'Euclid_{filter1}_vista_matched.fits'
    #image2 = imageDir / f'Euclid_{filter2}_vista_matched.fits'
    #image1 = imageDir / f'{fieldName}_{filter1}_resamp.fits'
    #image2 = imageDir / f'{fieldName}_{filter2}_resamp.fits'
    

    #! VISTA
    image1 = imageDir / f'UVISTA_{filter1}_DR6.fits'
    image2 = imageDir / f'UVISTA_{filter2}_DR6.fits'

    if filter3 != 'None':
        #image3 = imageDir / f'{fieldName}_{filter3}_MOSAIC.fits'
        #image3 = imageDir / f'Euclid_{filter3}_vista_matched.fits'
        #image3 = imageDir / f'UVISTA_{filter3}_DR6.fits'
        image3 = imageDir / f'{fieldName}_{filter3}_resamp.fits'

    if filter4 != 'None':
        image4 = imageDir / f'UVISTA_{filter4}_DR6.fits'

    #wht1 = imageDir / f'{fieldName}_{filter1}_MOSAIC_WHT.fits'
    #wht2 = imageDir / f'{fieldName}_{filter2}_MOSAIC_WHT.fits'
    #wht1 = imageDir / f'Euclid_{filter1}_vista_matched_WHT.fits'
    #wht2 = imageDir / f'Euclid_{filter2}_vista_matched_WHT.fits'

    #! VISTA
    wht1 = imageDir / f'UVISTA_{filter1}_DR6_wht.fits'
    wht2 = imageDir / f'UVISTA_{filter2}_DR6_wht.fits'
    # wht1 = imageDir / f'{fieldName}_{filter1}_resamp_wht.fits'
    # wht2 = imageDir / f'{fieldName}_{filter2}_resamp_wht.fits'


    if filter3 != 'None':
        #wht3 = imageDir / f'{fieldName}_{filter3}_MOSAIC_WHT.fits'
        #wht3 = imageDir / f'Euclid_{filter3}_vista_matched_WHT.fits'
        #wht3 = imageDir / f'UVISTA_{filter3}_DR6_wht.fits'
        wht3 = imageDir / f'{fieldName}_{filter3}_resamp_wht.fits'
    
    if filter4 != 'None':
        wht4 = imageDir / f'UVISTA_{filter4}_DR6_wht.fits'
        
    print('First image: ', image1)
    print('With weight: ', wht1)
    print('Second image: ', image2)
    print('With weight: ', wht2)
    if filter3 != 'None':
        print('Third image: ', image3)
        print('With weight: ', wht3)
    if filter4 != 'None':
        print('Fourth image: ', image4)
        print('With weight: ', wht4)


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

    #! Extract third image if stacking
    if filter3 != 'None':
        with fits.open(image3) as imHDU3:
            imHDU3 = fits.open(image3, memmap=True)
            image3 = imHDU3[0].data
            header3 = imHDU3[0].header

        with fits.open(wht3) as whtHDU3:
            whtHDU3 = fits.open(wht3, memmap=True)
            wht3 = whtHDU3[0].data
            wheader3 = whtHDU3[0].header

    #! Extract fourth image if stacking
    if filter4 != 'None':
        with fits.open(image4) as imHDU4:
            imHDU4 = fits.open(image4, memmap=True)
            image4 = imHDU4[0].data
            header4 = imHDU4[0].header
        
        with fits.open(wht4) as whtHDU4:
            whtHDU4 = fits.open(wht4, memmap=True)
            wht4 = whtHDU4[0].data
            wheader4 = whtHDU4[0].header




    #! STACKING PRECISELY TWO IMAGES
    if filter3 == 'None':
        print('Computing stack of {0}, {1}'.format(filter1, filter2))

        print('Computing sum of the weights...')
        weightSum = wht1 + wht2

        # Replace zeros with ones
        print('Replacing zeros with ones...')
        weightSum[weightSum == 0] = 1

        print('Computing weighted sum...')
        w1_1 = wht1 * image1
        w2_2 = wht2 * image2

        finalImage = (w1_1 + w2_2) / weightSum



    #! STACKING PRECISELY THREE IMAGES
    if filter3 != 'None' and filter4 == 'None':
        print('Computing stack of {0}, {1}, {2}'.format(filter1, filter2, filter3))

        print('Computing sum of the weights...')
        weightSum = wht1 + wht2 + wht3

        # Replace zeros with ones
        print('Replacing zeros with ones...')
        weightSum[weightSum == 0] = 1

        print('Computing weighted sum...')
        w1_1 = wht1 * image1
        w2_2 = wht2 * image2
        w3_3 = wht3 * image3

        finalImage = (w1_1 + w2_2 + w3_3) / weightSum




    #! STACKING PRECISELY FOUR IMAGES
    if filter4 != 'None':
        print('Computing stack of {0}, {1}, {2}, {3}'.format(filter1, filter2, filter3, filter4))

        print('Computing sum of the weights...')
        weightSum = wht1 + wht2 + wht3 + wht4

        # Replace zeros with ones
        print('Replacing zeros with ones...')
        weightSum[weightSum == 0] = 1

        print('Computing weighted sum...')
        w1_1 = wht1 * image1
        w2_2 = wht2 * image2
        w3_3 = wht3 * image3
        w4_4 = wht4 * image4

        finalImage = (w1_1 + w2_2 + w3_3 + w4_4) / weightSum

    #! Save image
    print('########Saving image.#########')

    if filter3 == 'None':
        #fits.writeto(stackDir / '{0}_{1}{2}_STACK.fits'.format(fieldName.upper(), name1, name2), finalImage, header1, overwrite=True)
        #fits.writeto(stackDir / f'Euclid_{filter1}{filter2}_vista_matched.fits', finalImage, header1, overwrite=True)
        fits.writeto(stackDir / f'UVISTA_{filter1}{filter2}_DR6.fits', finalImage, header1, overwrite=True)
    if filter3 != 'None' and filter4 == 'None':
        #fits.writeto(stackDir / '{0}_{1}{2}{3}_STACK.fits'.format(fieldName.upper(), filter1, filter2, filter3), finalImage, header1, overwrite=True)
        #fits.writeto(stackDir / f'Euclid_{filter1}{filter2}{filter3}_vista_matched.fits', finalImage, header1, overwrite=True)
        fits.writeto(stackDir / f'COSMOS_{filter1}{filter2}{filter3}_resamp.fits', finalImage, header1, overwrite=True)

    if filter4 != 'None':
        fits.writeto(stackDir / f'UVISTA_{filter1}{filter2}{filter3}{filter4}_DR6.fits', finalImage, header1, overwrite=True)

    print('###########Saved image to ', stackDir / fieldName.upper())

    #! Save weight
    print('########Saving weight#########')

    if filter3 == 'None' and filter4 == 'None':
        #fits.writeto(stackDir / '{0}_{1}{2}_STACK_WHT.fits'.format(fieldName.upper(), name1, name2), weightSum, wheader1, overwrite=True)
        #fits.writeto(stackDir / f'Euclid_{filter1}{filter2}_vista_matched_WHT.fits', weightSum, wheader1, overwrite=True)
        fits.writeto(stackDir / f'UVISTA_{filter1}{filter2}_DR6_wht.fits', weightSum, header1, overwrite=True)
    if filter3 != 'None':
        #fits.writeto(stackDir / '{0}_{1}{2}{3}_STACK_WHT.fits'.format(fieldName.upper(), filter1, filter2, filter3), weightSum, wheader1, overwrite=True)
        #fits.writeto(stackDir / f'Euclid_{filter1}{filter2}{filter3}_vista_matched_WHT.fits', weightSum, wheader1, overwrite=True)
        fits.writeto(stackDir / f'COSMOS_{filter1}{filter2}{filter3}_resamp_wht.fits', weightSum, header1, overwrite=True)
    if filter4 != 'None':
        fits.writeto(stackDir / f'UVISTA_{filter1}{filter2}{filter3}{filter4}_DR6_wht.fits', weightSum, wheader1, overwrite=True)

    print('###########Saved weight to ', stackDir / fieldName.upper())



