#!/usr/bin/env python3

"""
Not quite a crop - just need to adjust the centre of the final Euclid frames.

Created: Wednesday 19th June 2024.
"""

from pathlib import Path
import os

fields = ['COSMOS']

filter_names = ['VIS', 'Y', 'J', 'H']
#filter_names = ['YJH']


data_dir = Path.cwd().parents[3] / 'data'
input_swarp = Path.cwd() / 'crop.swarp'

for field in fields:

    euclid_dir = data_dir / 'COSMOS'

    for filter_name in filter_names:

        euclid_image = euclid_dir / f'Euclid_{filter_name}_vista_matched.fits'
        weight_image = euclid_dir / f'Euclid_{filter_name}_vista_matched_WHT.fits'

        save_image = f'Euclid_{filter_name}_cropped.fits'
        save_weight = f'Euclid_{filter_name}_cropped_WHT.fits'

        #os.system(f'fitsheader {str(reference_file)} > ref_header.txt')
        os.system(f'~/swarp/bin/swarp {str(euclid_image)} -WEIGHT_IMAGE {str(weight_image)} -c crop.swarp -IMAGEOUT_NAME {str(save_image)} -WEIGHTOUT_NAME {str(save_weight)}')


    