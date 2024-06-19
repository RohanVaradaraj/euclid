#!/usr/bin/env python3

"""
Swarp together the CWEB tiles resampled to the VISTA pixel scale.

Created: Tuesday 18th June 2024.
"""

import os
from astropy.io import fits
import glob
from pathlib import Path

jwst_dir = Path.cwd().parents[2] / 'images' / 'CWEB'

filter_names = ['F115W'] #, 'F150W', 'F277W', 'F444W']

for filter_name in filter_names:

    print('Swarping together the CWEB tiles for filter:', filter_name)

    fits_files = glob.glob(str(jwst_dir / 'mosaic_0*' / f'CWEB-{filter_name}-*resamp.fits'))
    rms_files = glob.glob(str(jwst_dir / 'mosaic_0*' / f'CWEB-{filter_name}-*resamp_rms.fits'))

    fits_files = sorted(fits_files)
    rms_files = sorted(rms_files)

    # Make string of image names to pass to Swarp
    images = fits_files[0]
    for i in range(1, len(fits_files)):
        images = images + ',' + fits_files[i]

    # Same for rms images
    rms = rms_files[0]
    for i in range(1, len(rms_files)):
        rms = rms + ',' + rms_files[i]

    # Output file.
    output_file = jwst_dir / f'CWEB_{filter_name}_resamp.fits'
    outrms_file = jwst_dir / f'CWEB_{filter_name}_resamp_rms.fits'

    # Key words
    config_string = '-c ./euclid.swarp'

    #! Run swarp on Euclid image
    keywords = '-WEIGHT_TYPE MAP_RMS -COMBINE_BUFSIZE 2048 -COMBINE_TYPE AVERAGE -VMEM_MAX 16384 -VMEM_DIR . \
                -MEM_MAX 2048 -COMBINE_BUFSIZE 2048'
    
    # swarp_command = f'~/swarp/bin/swarp {images} -WEIGHT_IMAGE {rms} {config_string} -IMAGEOUT_NAME {output_file} {keywords}' # \
    #                 #-WEIGHTOUT_NAME {outrms_file}'

    swarp_command = f'~/swarp/bin/swarp {images} -WEIGHT_IMAGE {rms} {config_string} -IMAGEOUT_NAME {output_file} {keywords}' # \
                    #-WEIGHTOUT_NAME {outrms_file}' 
    print(swarp_command)
    os.system(swarp_command)

    # #! And again on the weight image
    # keywords = ' -WEIGHT_TYPE ' + 'NONE' + \
    #            " -VMEM_MAX 16384  -VMEM_DIR . " \
    #            + " -MEM_MAX 2048 -COMBINE_BUFSIZE 2048 -COMBINE_TYPE AVERAGE "
    
    # swarp_command_rms = f'~/swarp/bin/swarp {rms} {config_string} -IMAGEOUT_NAME {outrms_file} {keywords}'
    # print(swarp_command_rms)
    # os.system(swarp_command_rms)