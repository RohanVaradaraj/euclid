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
import pickle
import numpy as np

field_name = 'CDFS'
filter_names = ['Y', 'J', 'H', 'VIS']

video_tiles = ['CDFS1', 'CDFS2', 'CDFS3']

euclid_dir = Path.cwd().parents[3] / 'data' / 'euclid' / 'euclid_deep_field_fornax'

# For CDFS, load the tile names that lie in the VIDEO footprints.
euclid_in_video = pickle.load(open('euclid_within_video.pkl', 'rb'))
print(euclid_in_video['CDFS1'])

# Loop over VIDEO tiles
for video_tile in video_tiles:
    # Loop over Euclid filters
    for filter_name in filter_names:

        # Get all the fits files for the filter, which lie in this video tile
        fits_files = glob.glob(str(euclid_dir / filter_name / '*BGSUB*.fits'))
        rms_files = glob.glob(str(euclid_dir / filter_name / '*RMS*.fits'))

        # Get the tile numbers
        all_tile_numbers = [file_name.split('/')[-1].split('-')[3].split('TILE')[1] for file_name in fits_files] 
        all_tile_numbers_rms = [file_name.split('/')[-1].split('-')[3].split('TILE')[1] for file_name in rms_files]
        
        #! Limit fits_files to those that lie in the VIDEO tile
        fits_files = [fits_file for fits_file, tile_number in zip(fits_files, all_tile_numbers) if tile_number in euclid_in_video[video_tile]]
        rms_files = [rms_file for rms_file, tile_number in zip(rms_files, all_tile_numbers_rms) if tile_number in euclid_in_video[video_tile]]

        # Sort files and rms to have same order, based on the tile numbers
        fits_files = sorted(fits_files, key=lambda x: int(x.split('/')[-1].split('-')[3].split('TILE')[1])) 
        rms_files = sorted(rms_files, key=lambda x: int(x.split('/')[-1].split('-')[3].split('TILE')[1]))

        # Create a list of the images to swarp together.
        images = fits_files[0]
        for i in range(1, len(fits_files)):
            images = images + ',' + fits_files[i]

        # Same for rms images
        rms = rms_files[0]
        for i in range(1, len(rms_files)):
            rms = rms + ',' + rms_files[i]

        # Output directory.
        output_dir = euclid_dir

        # Output file.
        output_file = output_dir / f'{video_tile}_{filter_name}_MOSAIC.fits'

        outrms_file = output_dir / f'{video_tile}_{filter_name}_MOSAIC_RMS.fits'

        # Key words
        config_string = '-c ./euclid.swarp'

        #! Run swarp on Euclid image
        keywords = '-WEIGHT_TYPE MAP_RMS -COMBINE_BUFSIZE 2048 -COMBINE_TYPE WEIGHTED -VMEM_MAX 16384 -VMEM_DIR . \
                    -MEM_MAX 2048 -COMBINE_BUFSIZE 2048 -PIXELSCALE_TYPE MEDIAN -PIXEL_SCALE 0.0 -IMAGE_SIZE 0'
        
        swarp_command = f'~/swarp/bin/swarp {images} -WEIGHT_IMAGE {rms} {config_string} -IMAGEOUT_NAME {output_file} {keywords}' 
        print(swarp_command)
        os.system(swarp_command)

        #! And again on the weight image
        keywords = ' -WEIGHT_TYPE ' + 'NONE' + \
                " -PIXELSCALE_TYPE MEDIAN " + "-PIXEL_SCALE 0.0 -IMAGE_SIZE 0 "+ \
                " -VMEM_MAX 16384  -VMEM_DIR . " \
                + " -MEM_MAX 2048 -COMBINE_BUFSIZE 2048 -COMBINE_TYPE AVERAGE "
        
        swarp_command_rms = f'~/swarp/bin/swarp {rms} {config_string} -IMAGEOUT_NAME {outrms_file} {keywords}'
        print(swarp_command_rms)
        os.system(swarp_command_rms)
