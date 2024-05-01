#!/usr/bin/env python3

# TEST

"""

area_calculation.py

Calculate the area of the Euclid image. Modified from area_mask.py

Originally created: Monday 3rd October 2022.
Modified: Tuesday 30th April 2024.

"""

# --- Import libraries ---

import numpy as np
from numpy import linalg as LA
import os
from astropy.io import fits
from astropy.table import Table
import matplotlib.pyplot as plt
from astropy.coordinates import SkyCoord
from astropy import units as u
from astropy.wcs import WCS
from regions import Regions
from shapely.geometry import Point, Polygon
from pathlib import Path

# Run the binary image creation?
make_binary = False

# Calculate the total area?
area_calc = True

dataDir = Path.home() / 'euclid' / 'COSMOS'

maskDir = Path.cwd().parent.parent.parent.parent / 'data' / 'masks'

if make_binary:
    mosaic = dataDir / 'COSMOS_Y_MOSAIC.fits'
else:
    mosaic = dataDir / 'junk' / 'COSMOS_BINARY_MOSAIC.fits'






#! ----- Read in image ---------

hdu = fits.open(mosaic, memmap=True)
header = hdu[0].header
image = hdu[0].data
print('Image loaded')

imageSize = image.shape # Size of mosaic
print('Image size: ', imageSize)

# Get WCS and pixel size
pixSizeDeg = abs(header['CD1_1'])
wsci = WCS(header)
print('pixel size in degrees: ', pixSizeDeg)


#! ------ APPLY EUCLID FOOTPRINT -------

if make_binary:

    # --- Read in the Euclid footprint ---
    euclid_mask = Regions.read(maskDir / 'Euclid_square_COSMOS.reg', format='ds9')

    print(euclid_mask)

    reg = euclid_mask[0]

    # Convert region to string to find out if it is circle or rectangle.
    regStr = str(reg)
    shape = regStr[8:14]

    # Polygons
    if shape == 'Polygo':

        print('This region is a polygon')

        # Get vertex RA
        RAs = reg.vertices.ra #.degree
        DECs = reg.vertices.dec #.degree

        x, y = wsci.all_world2pix(RAs, DECs, 1)

        # Points are labelled from bottom left, clockwise.
        # Want to compute gradients of the lines.
        m1 = (y[1]-y[0]) / (x[1]-x[0])
        m2 = (y[2]-y[1]) / (x[2]-x[1])
        m3 = (y[3]-y[2]) / (x[3]-x[2])
        m4 = (y[3]-y[0]) / (x[3]-x[0])


        print('Gradients are :', m1, m2, m3, m4)

        c1 = y[0] - m1*x[0]
        c2 = y[1] - m2*x[1]
        c3 = y[2] - m3*x[2]
        c4 = y[3] - m4*x[3]

        print('Intercepts are: ', c1, c2, c3, c4)

        Y, X = np.ogrid[:imageSize[0], :imageSize[1]]  # Make a mask to then apply to the image.
        print('Open grid created')

        # --- COMPUTE THE LINES ---
        line1 = (Y-c1)/m1 # Left line
        line2 = m2*X + c2 # Top line
        line3 = (Y-c3)/m3 # Right line
        line4 = m4*X + c4 # Bottom line
        print('Lines created')


        # Bound in lines
        cond = (X < line3) & (X > line1) & (Y < line2) & (Y > line4)
        print('Condition created')

        image[cond] = 1
        image[~cond] = 0
        print('Condition applied')

        del X, Y, cond
    
    ##################################
    #! Read in the mask for the image
    ##################################

    mask = Regions.read(maskDir / 'COSMOS' / 'Y_euclid.reg', format='ds9')

    for ii, reg in enumerate(mask):

        print(reg)

        shape = str(reg)[8:11]

        if shape == 'Pol':

            print('This region is a polygon')

            # Get vertex RA
            RAs = reg.vertices.ra #.degree
            DECs = reg.vertices.dec #.degree

            x, y = wsci.wcs_world2pix(RAs, DECs, 1)

            # Points are labelled from bottom left, clockwise.
            # Want to compute gradients of the lines.
            m1 = (y[1]-y[0]) / (x[1]-x[0])
            m2 = (y[2]-y[1]) / (x[2]-x[1])
            m3 = (y[3]-y[2]) / (x[3]-x[2])
            m4 = (y[3]-y[0]) / (x[3]-x[0])


            print('Gradients are :', m1, m2, m3, m4)

            c1 = y[0] - m1*x[0]
            c2 = y[1] - m2*x[1]
            c3 = y[2] - m3*x[2]
            c4 = y[3] - m4*x[3]

            print('Intercepts are: ', c1, c2, c3, c4)

            Y, X = np.ogrid[:imageSize[0], :imageSize[1]]  # Make a mask to then apply to the image.
            print('Open grid created')

            # --- COMPUTE THE LINES ---
            line1 = (Y-c1)/m1 # Left line
            line2 = m2*X + c2 # Top line
            line3 = (Y-c3)/m3 # Right line
            line4 = m4*X + c4 # Bottom line
            print('Lines created')


            # Bound in lines
            cond = (X < line3) & (X > line1) & (Y < line2) & (Y > line4)
            print('Condition created')

#                if (np.max(x) - np.min(x) > 200):
            image[cond] = 0
            print('Condition applied')

            del X, Y, cond

        if shape == 'Rec':

            print('This region is a rectangle')
            #continue

            # RA DEC of centre
            x = reg.center.ra
            y = reg.center.dec

            print('X and Y of region: ', x, y)

            Hw = reg.width.value/2

            # Convert from arcsec to pixels
            Hw = Hw / (pixSizeDeg * 3600)
            print('half width: ', Hw)

            # Half height
            Hh = reg.height.value/2

            # Convert from arcsec to pixels
            Hh = Hh / (pixSizeDeg * 3600)
            print('half height: ', Hh)

            # Convert to pixels in the mosaic
            xc, yc = wsci.wcs_world2pix(x, y, 1)

            # Points starting from bottom left, anticlockwise.
            p1 = [xc-Hw, yc-Hh]
            p2 = [xc+Hw, yc-Hh]
            p3 = [xc+Hw, yc+Hh]
            p4 = [xc-Hw, yc+Hh]

            Y, X = np.ogrid[:imageSize[0], :imageSize[1]]  # Make a mask to then apply to the image.
            print('Open grid created')

            # Compute region in box
            box = (X > p1[0]) & (X < p2[0]) & (Y < p3[1]) & (Y > p2[1])

            image[box] = 0

            del X, Y, box

    hdu.writeto(dataDir / 'junk' / 'COSMOS_BINARY_MOSAIC.fits', overwrite=True)


if area_calc:
    #! ---------- COMPUTE FINAL AREA!!! -------------
    Npix = np.sum(image)

    area = pixSizeDeg * pixSizeDeg * Npix

    print('AREA IS {0} SQUARE DEGREES'.format(area))
    exit()

