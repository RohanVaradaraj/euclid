#!/usr/bin/env python3

"""
convert_rms_to_weight_cdfs_DR1.py

Modified from convert_rms_to_weight.py 
"""

from astropy.io import fits
from pathlib import Path
import glob

euclid_dir = Path.home().parents[1] / 'extraspace' / 'varadaraj' / 'euclid' / 'DR1'


for field in ['CDFS1', 'CDFS2', 'CDFS3']:
    print(f'Running in {field}')

    image_dir = euclid_dir / field
    print(image_dir)

    # Get all files containing "var"
    rms_files = glob.glob(str(image_dir / '*.var.*'))
    print('Variance files')
    print(rms_files)

    for i, rms_file in enumerate(rms_files):
        print(rms_file)
        print(f'File no. {i+1} of {len(rms_files)}')

        # Open the VAR image
        with fits.open(Path(rms_file)) as hdul:
            data = hdul[0].data

        # Create the inverse variance image
        inverse_variance = 1 / data

        # Save the inverse variance image
        inverse_variance_file = rms_file.replace('.var.', '.inverse_var.')
        print(inverse_variance_file)
        with fits.open(rms_file) as hdul:
            hdul[0].data = inverse_variance
            hdul.writeto(inverse_variance_file, overwrite=True)
        
        # Delete the RMS map safely.
        print(f'Deleting {rms_file}')
        Path(rms_file).unlink(missing_ok=True)