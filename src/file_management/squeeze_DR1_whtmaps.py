#!/usr/bin/env python3

"""
Squeeze the final axis of the DR1 WHT NISP CDFS maps.

Modified from convert_rms_to_weight_cdfs_DR1.py
"""

from astropy.io import fits
from pathlib import Path
import glob
import numpy as np

euclid_dir = Path.home().parents[1] / 'extraspace' / 'varadaraj' / 'euclid' / 'DR1'


for field in ['CDFS1', 'CDFS2', 'CDFS3']:
    print(f'Running in {field}')

    image_dir = euclid_dir / field
    print(image_dir)

    # Get all files containing "var"
    rms_files = glob.glob(str(image_dir / '*.inverse_var.*'))
    print('Inverse variance files')
    print(rms_files)

    for i, rms_file in enumerate(rms_files):
        print(rms_file)
        print(f'File no. {i+1} of {len(rms_files)}')

        with fits.open(rms_file) as hdul:
            data = np.squeeze(hdul[0].data)
            # Print the shape of the data
            print(f'Shape of data: {data.shape}')
            hdul[0].data = data
            hdul.writeto(rms_file, overwrite=True)





