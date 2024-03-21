#!/usr/bin/env python3

"""
swarp.py

Swarp together the Euclid tiles to make a full mosaic.

Created: Thursday 21st March 2024.
"""

import os
from astropy.io import fits
import glob
from pathlib import Path

euclid_dir = Path.home() / 'euclid'
field_name = 'COSMOS'
filter_names = ['VIS'] #, 'J', 'H', 'VIS']

# Define the sorting function
def sort_by_number(filename):
    return int(filename.split('_')[-1].split('.')[0])



for filter_name in filter_names:

    # Get all the fits files for the filter.
    fits_files = glob.glob(str(euclid_dir / filter_name / field_name / '*BGSUB*.fits'))
    rms_files = glob.glob(str(euclid_dir / filter_name / field_name / '*RMS*.fits'))

    # Sort fits_files and rms_files based on the numeric part of the filenames
    fits_files = sorted(fits_files, key=sort_by_number)
    rms_files = sorted(rms_files, key=sort_by_number)

    # Create a list of the images to swarp together.
    images = ' '.join(fits_files)

    # Create a list of the rms images
    rms = '-WEIGHT_IMAGE ' + ','.join(rms_files)

    # Output directory.
    output_dir = euclid_dir / field_name

    # Output file.
    output_file = output_dir / f'{field_name}_{filter_name}_MOSAIC.fits'

    outrms_file = output_dir / f'{field_name}_{filter_name}_MOSAIC_RMS.fits'

    # Key words
    config_string = '-c ./euclid.swarp'

    # Swarp command.
    swarp_command = f'swarp {images} {rms} {config_string} -IMAGEOUT_NAME {output_file} -WEIGHTOUT_NAME {outrms_file}'
    print(swarp_command)

    # Run the swarp command.
    #os.system(swarp_command)
