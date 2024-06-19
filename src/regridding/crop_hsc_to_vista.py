#!/usr/bin/env python3

"""
DR6 images are larger than DR5 etc, so crop down.

Created: Wednesday 19th June 2024.
"""

from pathlib import Path
import os
from astropy.io import fits

fields = ['COSMOS']

filter_names = ['Y'] #, 'J', 'H', 'K']
#filter_names = ['YJH']

ref_image = Path.cwd().parents[3] / 'data' / 'COSMOS' / 'HSC-Y_DR3.fits'


data_dir = Path.cwd().parents[3] / 'data'
input_swarp = Path.cwd() / 'crop.swarp'

for field in fields:

    output_dir = data_dir / 'COSMOS'

    for filter_name in filter_names:

        print('Cropping: ', filter_name)

        euclid_image = output_dir / f'UVISTA_p2m-{filter_name}_full_hr_v10.fits'
        weight_image = output_dir / f'UVISTA_p2m-{filter_name}_full_hr_v10_weight.fits'

        save_image = output_dir / f'UVISTA_{filter_name}_DR6.fits'
        save_weight = output_dir / f'UVISTA_{filter_name}_DR6_WHT.fits'

        # Now open output image and update the header CRPIX1 and CRPIX2 values of the reference image
        with fits.open(euclid_image, mode='update') as hdul:
            hdul[0].header['CRPIX1'] = fits.open(ref_image)[0].header['CRPIX1']
            hdul[0].header['CRPIX2'] = fits.open(ref_image)[0].header['CRPIX2']
            hdul.flush()

        with fits.open(euclid_image, mode='update') as hdul:
            hdul[0].header['CRPIX1'] = fits.open(ref_image)[0].header['CRPIX1']
            hdul[0].header['CRPIX2'] = fits.open(ref_image)[0].header['CRPIX2']
            hdul.flush()

        #os.system(f'fitsheader {str(reference_file)} > ref_header.txt')
        os.system(f'~/swarp/bin/swarp {str(euclid_image)} -WEIGHT_IMAGE {str(weight_image)} -c crop.swarp -IMAGEOUT_NAME {str(save_image)} -WEIGHTOUT_NAME {str(save_weight)}')

