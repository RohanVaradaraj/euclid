"""
Generate the 10x10 arcsec cutout files as {source_name}_{filter_name}.fits for the UltraVISTA+Euclid sample.

Created: Friday 23rd January 2026.
"""

from pathlib import Path
import os
from astropy.io import fits
from astropy.wcs import WCS
import numpy as np
from astropy.table import Table
from astropy.nddata import Cutout2D
from astropy.coordinates import SkyCoord
import astropy.units as u
import matplotlib.pyplot as plt

#! Directories
image_dir = Path.home() / 'euclid' / 'COSMOS'

output_dir = Path.cwd() / 'cutout_files'
os.makedirs(output_dir, exist_ok=True)

#! Load candidate catalogue
cat_filename = Path.cwd() / 'UltraVISTA_plus_Euclid_z7_sample.fits'
t = Table.read(cat_filename)

t.sort('Muv')

print(t)

#! Filters to generate cutouts for
filter_names = ['Y', 'J', 'H', 'VIS']
cutout_size_arcsec = 10.0 

filter_dict = {
    'Y': 'YE',
    'J': 'JE',
    'H': 'HE',
    'VIS': 'IE'
}

#! Loop through sources
for i, source in enumerate(t):

    ID = source['SOURCE_ID']
    ra = source['RA']
    dec = source['DEC']

    name = ID.replace(' ', '')

    print('Source {}/{}:'.format(i+1, len(t)))

    print(f'Generating cutouts for source ID {ID} at RA={ra}, DEC={dec}')
    for filter_name in filter_names:

        with fits.open(image_dir / f'COSMOS_{filter_name}_MOSAIC.fits') as hdul:

            img = hdul[0].data
            hdr = hdul[0].header
            wcs = WCS(hdr)

            position = SkyCoord(ra, dec, unit='deg')
            cutout_size = u.Quantity((cutout_size_arcsec, cutout_size_arcsec), u.arcsec)
            cutout = Cutout2D(img, position, cutout_size, wcs=wcs)

            # Generate output filename
            output_filename = output_dir / f'{name}_{filter_dict[filter_name]}.fits'
            print(f'  Saving cutout to {output_filename}')

            # Save cutout to FITS file
            hdu = fits.PrimaryHDU(data=cutout.data, header=cutout.wcs.to_header())
            hdu.writeto(output_filename, overwrite=True)


            # Plot check
            # plt.imshow(cutout.data, origin='lower', cmap='gray')
            # plt.savefig('check_stamp.pdf')
            # plt.show()
            # exit()




