#!/usr/bin/env python3

"""
reproject the euclid images to the VISTA images

Created: Friday 14th June 2024.
"""

from pathlib import Path
import os

fields = ['COSMOS']
filter_names = ['VIS'] #, 'Y', 'J', 'H']

data_dir = Path.cwd().parents[3] / 'data'
calib_dir = data_dir / 'bertin_config'

input_sex = calib_dir / 'rohan_scamp_sex.sex'
input_scamp = calib_dir / 'tweak_default.scamp'
input_swarp = Path.cwd() / 'regrid.swarp'


for field in fields:

    reference_dir = data_dir / f'{field}'
    reference_file = reference_dir / 'UVISTA_Y_DR6_cropped.fits'

    euclid_dir = data_dir / 'euclid' / 'images'

    for filter_name in filter_names:

        euclid_image = euclid_dir / f'{field}_{filter_name}_MOSAIC.fits'
        weight_image = euclid_dir / f'{field}_{filter_name}_MOSAIC_WHT.fits'

        

        #os.system(f'fitsheader {str(reference_file)} > ref_header.txt')
        os.system(f'~/swarp/bin/swarp {str(euclid_image)} -WEIGHT_IMAGE {str(weight_image)} -c regrid.swarp -IMAGEOUT_NAME {field}_{filter_name}_resamp.fits -WEIGHTOUT_NAME {field}_{filter_name}_resamp_wht.fits')
