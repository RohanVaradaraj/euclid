#!/usr/bin/env python3

'''
crop_mosaics.py

The outputs of swarp are much larger than they need to be.

Crop the images down to a smaller, more manageable size.

Created: Friday 22nd March 2024.
'''

from astropy.io import fits
from astropy.wcs import WCS
from astropy.nddata.utils import Cutout2D
from pathlib import Path

mosaic_dir = Path.home() / 'euclid' / 'COSMOS'

filter_names = ['J', 'H', 'VIS']

xmin = 5200
xmax = 44000
ymin = 7600
ymax = 47600

xcen, ycen = (xmin + xmax) / 2, (ymin + ymax) / 2
x_size, y_size = xmax - xmin, ymax - ymin

for filter_name in filter_names:

    print(f'Cropping {filter_name}')

    #! Crop the image
    # Get the mosaic file.
    mosaic_file = mosaic_dir / f'COSMOS_{filter_name}_MOSAIC.fits'
    print('Doing science image')

    with fits.open(mosaic_file, memmap=True) as hdul:
        data = hdul[0].data
        header = hdul[0].header
        print('Image opened')

        # Get the wcs.
        wcs = WCS(header)

        # Cutout.
        cutout = Cutout2D(data, (xcen, ycen), size=(y_size, x_size), wcs=wcs)
        print('Cutout made')

        # Update the header.
        header.update(cutout.wcs.to_header())

        # Save the cutout.
        cutout_file = mosaic_dir / f'COSMOS_{filter_name}_MOSAIC.fits'
        cutout_data = cutout.data
        hdu = fits.PrimaryHDU(cutout_data, header=header)
        hdu.writeto(cutout_file, overwrite=True)
        print(f'Cutout saved to {cutout_file}')

    #! Crop the RMS
    # Get the mosaic file.
    rms_file = mosaic_dir / f'COSMOS_{filter_name}_MOSAIC_RMS.fits'
    print('Doing RMS image')

    with fits.open(rms_file, memmap=True) as hdul:
        data = hdul[0].data
        header = hdul[0].header
        print('Image opened')

        # Get the wcs.
        wcs = WCS(header)

        # Cutout.
        cutout = Cutout2D(data, (xcen, ycen), size=(y_size, x_size), wcs=wcs)
        print('Cutout made')

        # Update the header.
        header.update(cutout.wcs.to_header())

        # Save the cutout.
        cutout_file = mosaic_dir / f'COSMOS_{filter_name}_MOSAIC_RMS.fits'
        cutout_data = cutout.data
        hdu = fits.PrimaryHDU(cutout_data, header=header)
        hdu.writeto(cutout_file, overwrite=True)
        print(f'Cutout saved to {cutout_file}')
