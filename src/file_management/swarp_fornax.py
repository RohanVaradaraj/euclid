#!/usr/bin/env python3

"""
swarp.py

Swarp together the Euclid tiles to make a full mosaic, resampled to VISTA and on same plate scale.

Created: Thursday 21st March 2025.
"""

import os
from astropy.io import fits
import glob
from pathlib import Path
import pickle
import numpy as np
from astropy.wcs import WCS
import re

def get_tile_id(filename):
    match = re.search(r'TILE(\d+)', filename)
    return match.group(1) if match else None

#! Switches/variables
filter_names = ['VIS'] #? Euclid filter(s) we want to resample/swarp
video_tiles = ['CDFS1', 'CDFS2', 'CDFS3'] #? VIDEO tile(s) we want to match the Euclid mosaics to.
test = False #? If true, run on three tiles only

# Path to the Euclid tiles
#euclid_dir = Path.cwd().parents[3] / 'data' / 'euclid' / 'euclid_deep_field_fornax' #! Q1
euclid_dir = Path.home() / 'euclid' / 'CDFS' #! DR1

# Output directory
#output_dir = euclid_dir / 'tmp #! Q1
output_dir = euclid_dir.parent #! DR1

#! Load the Euclid tile names that lie in the VIDEO footprints.
# euclid_in_video = pickle.load(open('euclid_within_video.pkl', 'rb')) #! Q1
euclid_in_video = pickle.load(open('euclid_DR1_within_video.pkl', 'rb')) #! DR1

# Loop over VIDEO tiles
for video_tile in video_tiles:

    # update output dir
    final_dir = output_dir / video_tile
    print(f'Final mosaic will be in {final_dir}')

    # Open the VIDEO tile
    video_dir = Path.home().parents[1] / 'hoy' / 'VIDEO_FINAL'
    video_name = f'{video_tile.lower()}_Y.fits'
    with fits.open(video_dir / video_name) as hdul:

        # Get pixel scale
        pixel_scale = np.abs(hdul[0].header['CD1_1']) * 3600
        print(f'Pixel scale for {video_tile} is {pixel_scale} arcsec/pixel')

        # Get the image size
        image_size = [hdul[0].header['NAXIS1'], hdul[0].header['NAXIS2']]
        print(f'Image size for {video_tile} is {image_size[0]},{image_size[1]} pixels')

        # Find the centre of the image in ra dec
        wcs = WCS(hdul[0].header)
        centre = wcs.pixel_to_world(image_size[0] / 2, image_size[1] / 2)
        print(f'Centre of {video_tile} is {centre.ra.deg},{centre.dec.deg}')

    # Loop over Euclid filters
    for filter_name in filter_names:

        # Get all the fits files for the filter
        # fits_files = glob.glob(str(euclid_dir / filter_name / '*BGSUB*.fits'))    #! Q1
        # rms_files = glob.glob(str(euclid_dir / filter_name / '*MAP_WEIGHT*.fits'))
        fits_files = glob.glob(str(euclid_dir / '*BGSUB*.fits'))      #! DR1
        rms_files = glob.glob(str(euclid_dir / '*MAP_WEIGHT*.fits'))

        # Get the tile numbers from the file names
        all_tile_numbers = [get_tile_id(file_name) for file_name in fits_files]
        all_tile_numbers_rms = [get_tile_id(file_name) for file_name in rms_files]
        
        #! Limit fits_files to those that lie in this specific VIDEO tile
        fits_files = [fits_file for fits_file, tile_number in zip(fits_files, all_tile_numbers) if tile_number in euclid_in_video[video_tile]]
        rms_files = [rms_file for rms_file, tile_number in zip(rms_files, all_tile_numbers_rms) if tile_number in euclid_in_video[video_tile]]

        #! Also limit the tile numbers to those that lie in this specific VIDEO tile
        tile_numbers_here = [tile_number for tile_number in all_tile_numbers if tile_number in euclid_in_video[video_tile]]

        # Build a dictionary of {tile_id: filename} for each
        sci_dict = {get_tile_id(f): f for f in fits_files}
        rms_dict = {get_tile_id(f): f for f in rms_files}

        # Get intersection of available tile IDs in both lists
        common_tile_ids = sorted(set(sci_dict) & set(rms_dict), key=int)

        # Build aligned, ordered lists
        fits_files = [sci_dict[tile_id] for tile_id in common_tile_ids]
        rms_files = [rms_dict[tile_id] for tile_id in common_tile_ids]
        
        # Also sort the tile numbers!
        tile_numbers_here = sorted(tile_numbers_here, key=lambda x: int(x))

        if test:
            # Limit to the first three tiles for testing
            fits_files = fits_files[0:2]
            rms_files = rms_files[0:2]
            tile_numbers_here = tile_numbers_here[0:2]

        #! ############################ OPTION B: SWARP ALL EUCLID TILES TOGETHER ############################
        # Output file names
        output_file = output_dir / f'{video_tile}_{filter_name}_MOSAIC.fits'
        outrms_file = output_dir / f'{video_tile}_{filter_name}_MOSAIC_WHT.fits'

        # Check if the output file and RMS file already exist
        if output_file.exists():
            print(f'Skipping {output_file} as it already exists.')
            continue
            
        if outrms_file.exists():
            print(f'Skipping {outrms_file} as it already exists.')
            continue

        #! If the file doesn't exist, go ahead and swarp it to match VIDEO!
        # Key words
        config_string = f'-c ./euclid_{video_tile.lower()}.swarp'
        #config_string = f'-c ./euclid.swarp'

        #! Run swarp on Euclid image
        keywords = f'-WEIGHT_TYPE MAP_WEIGHT -COMBINE_TYPE WEIGHTED -VMEM_MAX 16384 -VMEM_DIR . \
                    -MEM_MAX 2048 -COMBINE_BUFSIZE 2048 ' 

        # Convert fits_files and rms_files to a comma-separated string
        fits_files = ','.join(fits_files)
        rms_files = ','.join(rms_files)
        
        swarp_command = f'~/swarp/bin/swarp {fits_files} -WEIGHT_IMAGE {rms_files} {config_string} -IMAGEOUT_NAME {output_file} {keywords}' 
        print(swarp_command)
        os.system(swarp_command)

        #! And again on the weight image
        keywords = "-WEIGHT_TYPE NONE \
                -VMEM_MAX 16384  -VMEM_DIR . \
                -MEM_MAX 2048 -COMBINE_BUFSIZE 2048 -COMBINE_TYPE AVERAGE "
        
        swarp_command_rms = f'~/swarp/bin/swarp {rms_files} {config_string} -IMAGEOUT_NAME {outrms_file} {keywords}'
        print(swarp_command_rms)
        os.system(swarp_command_rms)