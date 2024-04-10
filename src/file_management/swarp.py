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
filter_names = ['Y'] #, 'J', 'H', 'VIS']

# Define the sorting function
def sort_by_number(filename):
    return int(filename.split('_')[-1].split('.')[0])



for filter_name in filter_names:

    # Get all the fits files for the filter.
    fits_files = glob.glob(str(euclid_dir / filter_name / field_name / '*BGSUB*.fits'))
    rms_files = glob.glob(str(euclid_dir / filter_name / field_name / '*MAP_WEIGHT*.fits'))

    # Sort fits_files and rms_files based on the numeric part of the filenames
    fits_files = sorted(fits_files, key=sort_by_number)
    rms_files = sorted(rms_files, key=sort_by_number)

    # Create a list of the images to swarp together.
    images = fits_files[0]
    for i in range(1, len(fits_files)):
        images = images + ',' + fits_files[i]

    # Same for rms images
    rms = rms_files[0]
    for i in range(1, len(rms_files)):
        rms = rms + ',' + rms_files[i]


    # Output directory.
    output_dir = euclid_dir / field_name

    # Output file.
    output_file = output_dir / f'{field_name}_{filter_name}_MOSAIC.fits'

    outrms_file = output_dir / f'{field_name}_{filter_name}_MOSAIC_WHT.fits'

    # Key words
    config_string = '-c ./euclid.swarp'

    #! Run swarp on Euclid image
    keywords = '-WEIGHT_TYPE MAP_WEIGHT -COMBINE_BUFSIZE 2048 -COMBINE_TYPE WEIGHTED -VMEM_MAX 16384 -VMEM_DIR . \
                -MEM_MAX 2048 -COMBINE_BUFSIZE 2048 -PIXELSCALE_TYPE MEDIAN -PIXEL_SCALE 0.0 -IMAGE_SIZE 0'
    
    swarp_command = f'~/swarp/bin/swarp {images} -WEIGHT_IMAGE {rms} {config_string} -IMAGEOUT_NAME {output_file} {keywords}' 
    print(swarp_command)
    #os.system(swarp_command)

    #! And again on the weight image
    keywords = ' -WEIGHT_TYPE ' + 'NONE' + \
               " -PIXELSCALE_TYPE MEDIAN " + "-PIXEL_SCALE 0.0 -IMAGE_SIZE 0 "+ \
               " -VMEM_MAX 16384  -VMEM_DIR . " \
               + " -MEM_MAX 2048 -COMBINE_BUFSIZE 2048 -COMBINE_TYPE AVERAGE "
    
    swarp_command_rms = f'~/swarp/bin/swarp {rms} {config_string} -IMAGEOUT_NAME {outrms_file} {keywords}'
    print(swarp_command_rms)
    os.system(swarp_command_rms)
