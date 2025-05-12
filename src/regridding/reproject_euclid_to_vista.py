#!/usr/bin/env python3

"""
reproject the euclid images to VIDEO. Modified from reproject_image.py. More general - automatically find image centres, dimensions, etc.

Created: Friday 14th June 2024.

Modified: Monday 12th May 2025.
"""

from pathlib import Path
import os
from astropy.io import fits
from astropy.wcs import WCS
import glob
import numpy as np

fields = ['CDFS1', 'CDFS2', 'CDFS3']

filter_names = ['VIS', 'Y', 'J', 'H'] # Euclid
#filter_names = ['YJH'] # Euclid stack

data_dir = Path.cwd().parents[3] / 'data'
calib_dir = data_dir / 'bertin_config'

input_swarp_filename = 'regrid_euclid_to_vista.swarp'
input_swarp = Path.cwd() / input_swarp_filename

for field in fields:

    reference_dir = data_dir / f'{field}'
    reference_file = reference_dir / 'HSC-G.fits'

    # Open the reference file
    with fits.open(reference_file) as ref_hdu:
        ref_wcs = WCS(ref_hdu[0].header)
        ref_header = ref_hdu[0].header

        # Image size
        naxis1 = ref_header['NAXIS1']
        naxis2 = ref_header['NAXIS2']

        # Pixel scale
        pixel_scale = np.abs(ref_header['CD1_1']) * 3600.0 # Convert to arcseconds

        # Find centre pixel coordinate
        centre_x = naxis1 / 2
        centre_y = naxis2 / 2

        # Find the RA and DEC of the centre pixel
        centre_ra, centre_dec = ref_wcs.all_pix2world(centre_x, centre_y, 0) 
        print(f'Centre RA: {centre_ra}, Centre DEC: {centre_dec}')
        continue

    euclid_dir = data_dir / 'euclid' / 'images'

    for filter_name in filter_names:

        if filter_name != 'YJH':
            euclid_image = euclid_dir / f'{field}_{filter_name}_MOSAIC.fits' #! Original images
            weight_image = euclid_dir / f'{field}_{filter_name}_MOSAIC_RMS.fits'

        if filter_name == 'YJH':
            euclid_image = euclid_dir / f'{field}_{filter_name}_STACK.fits' #! YJH stack
            weight_image = euclid_dir / f'{field}_{filter_name}_STACK_WHT.fits'

        save_image = euclid_dir / f'{field}_{filter_name}_resamp.fits'
        save_weight = euclid_dir / f'{field}_{filter_name}_resamp_wht.fits'

        #? FIRST THE SCI IMAGE
        # Extra memory keywords
        keywords = f'-COMBINE_BUFSIZE 2048 -COMBINE_TYPE WEIGHTED -VMEM_MAX 16384 -VMEM_DIR . \
                    -MEM_MAX 2048 -COMBINE_BUFSIZE 2048 -IMAGE_SIZE {naxis1},{naxis2}  \
                    -PIXEL_SCALE {pixel_scale} -CENTER {centre_x},{centre_y} '

        # Now reproject the image
        os.system(f'~/swarp/bin/swarp {euclid_image} -WEIGHT_IMAGE {weight_image} -c {input_swarp_filename} -IMAGEOUT_NAME {str(save_image)} {keywords}')

        #? NOW THE RMS IMAGE
        # Extra memory keywords
        keywords = '-WEIGHT_TYPE NONE -COMBINE_BUFSIZE 2048 -COMBINE_TYPE AVERAGE -VMEM_MAX 16384 -VMEM_DIR . \
                    -MEM_MAX 2048 -COMBINE_BUFSIZE 2048 -IMAGE_SIZE {naxis1},{naxis2}  \
                    -PIXEL_SCALE {pixel_scale} -CENTER {centre_x},{centre_y} '

        # Now reproject the weight image
        os.system(f'~/swarp/bin/swarp {weight_image} -c regrid.swarp -IMAGEOUT_NAME {str(save_weight)} {keywords}')
