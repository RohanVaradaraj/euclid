#!/usr/bin/env python3

"""
convert_rms_to_weight.py

Create inverse variance images from the RMS images. Want to see if this helps with PSFEx etc.

Created: Tuesday 9th April 2024
"""

from astropy.io import fits
from pathlib import Path
import glob

filter_names = ['VIS', 'Y', 'J', 'H']

euclid_dir = Path.home() / 'euclid'

for filter_name in filter_names:
    print(f'Running in {filter_name}')

    image_dir = euclid_dir / filter_name / 'COSMOS'

    # Get all files containing "RMS" in the name
    rms_files = glob.glob(str(image_dir / '*RMS*'))

    for rms_file in rms_files:

        # Open the RMS image
        with fits.open(Path(rms_file)) as hdul:
            data = hdul[0].data

        # Create the inverse variance image
        inverse_variance = 1 / data**2

        # Save the inverse variance image
        inverse_variance_file = rms_file.replace('RMS', 'MAP_WEIGHT')
        print(inverse_variance_file)
        with fits.open(rms_file) as hdul:
            hdul[0].data = inverse_variance
            hdul.writeto(inverse_variance_file, overwrite=True)