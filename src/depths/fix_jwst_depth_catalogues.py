"""
fix_jwst_depth_catalogues.py

The output catalogues from running the depth codes on the JWST tiles returns an incorrect RA DEC.
The resulting catalogues have the incorrect overlap.
This code assumes the pixel coordinates of the objects are correct, reads in the WCS from the image headers and applies the transformation.

Created: Wednesday 2nd October 2024.
"""

from astropy.io import fits
from astropy.wcs import WCS
from astropy.table import Table
from pathlib import Path
import glob
import fnmatch
import re

image_dir = Path.cwd().parents[3] / 'data' / 'CWEB'
cat_dir = Path.cwd().parents[1] / 'data' / 'depths' / 'COSMOS' / 'catalogues'

filter_names = ['f115w', 'f150w', 'f277w', 'f444w']

for filter_name in filter_names:

    tile_pattern = '[0-7][AB]'
    files = glob.glob(str(cat_dir / f'd{filter_name}_{tile_pattern}.fits'))


    print(files)

    for file_name in files:

        tile = file_name.split('_')[1].split('.')[0]

        print(tile)

        mosaic_dir = image_dir / f'mosaic_{tile}'

        with fits.open(mosaic_dir / f'CWEB-{filter_name.upper()}-{tile}_i2dnobg_small.fits') as hdul:
            wcs = WCS(hdul[1].header)

        t = Table.read(file_name)

        ra, dec = wcs.all_pix2world(t['X_IMAGE'], t['Y_IMAGE'], 0)

        t['RA'] = ra
        t['DEC'] = dec

        t.remove_columns(['ALPHA_J2000', 'DELTA_J2000'])
        t.rename_column('RA', 'ALPHA_J2000')
        t.rename_column('DEC', 'DELTA_J2000')

        t.write(file_name, overwrite=True)

