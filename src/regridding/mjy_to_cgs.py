#!/usr/bin/env python3

"""
This module contains a function that converts flux densities from mJy/steradian to cgs units.
"""

from astropy.io import fits
import numpy as np
from pathlib import Path

image_dir = Path.cwd().parents[3] / 'data' / 'COSMOS'

filter_names = ['f115w', 'f150w', 'f277w', 'f444w']

for filter_name in filter_names:

    print(f'Converting {filter_name} image to cgs units')
    image_name = f'CWEB_{filter_name}_vista_matched.fits'

    with fits.open(image_dir / image_name) as hdul:
        image_data = hdul[0].data
        image_header = hdul[0].header
    print('Opened image')

    #! Convert from MJy/steradian to cgs

    # solid angle of a pixel
    solid_angle = ( 0.15 * (np.pi/180) * (1/3600) )**2

    # Find area of pixel in steradians
    pix_area = solid_angle

    # Zeropoint in mAB depends on pixel area when the units are MJy/sr. 
    # See https://jwst-docs.stsci.edu/jwst-near-infrared-camera/nircam-performance/nircam-absolute-flux-calibration-and-zeropoints#gsc.tab=0
    zpt = -6.10 - 2.5*np.log10(solid_angle)
    print(f'The AB zeropoint for a pixel scale of 0.15 arcsec/pix is {zpt}')

    zpt_factor = (48.6 + 28.08) / 2.5

    # conversion factor
    conversion_factor = 1e-23 * 1e6 * solid_angle * 10**zpt_factor
    print(f'The conversion factor is {conversion_factor}')

    print('Converting to cgs')
    image_data_cgs = image_data / conversion_factor
    print('Converted to cgs')
    print(np.nanmean(image_data_cgs))
    print(np.nanmean(image_data))

    #! Save the image
    print('Saving image')
    image_name_cgs = f'CWEB_{filter_name}_vista_matched_cgs.fits'

    hdu = fits.PrimaryHDU(image_data_cgs, header=image_header)
    hdu.writeto(image_dir / image_name_cgs, overwrite=True)
    print('Saved image')



