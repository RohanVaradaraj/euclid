#!/usr/bin/env python3

"""
Crop the euclid image to make a smaller test region to test different sextractor configs.

Created: Monday 8th July 2024.
"""

from astropy.io import fits
from pathlib import Path
from astropy.wcs import WCS

image_dir = Path.cwd().parents[3] / 'data' / 'COSMOS'
image_name = 'Euclid_Y_vista_matched.fits'
weight_name = 'Euclid_Y_vista_matched_WHT.fits'


with fits.open(image_dir / image_name) as hdu:
    image = hdu[0].data
    hdr = hdu[0].header

with fits.open(image_dir / weight_name) as hdu:
    weight = hdu[0].data
    hdr_w = hdu[0].header

# Crop the image to central 5000x1000 pixels
x_cen, y_cen = image.shape[0]//2, image.shape[1]//2

image_crop = image[x_cen-2500:x_cen+2500, y_cen-2500:y_cen+2500]
weight_crop = weight[x_cen-2500:x_cen+2500, y_cen-2500:y_cen+2500]

# Update the header
hdr['NAXIS1'] = 5000
hdr['NAXIS2'] = 5000

hdr_w['NAXIS1'] = 5000
hdr_w['NAXIS2'] = 5000

hdu = fits.PrimaryHDU(image_crop, header=hdr)
hdu.writeto(image_dir / 'test_Y.fits', overwrite=True)

hdu = fits.PrimaryHDU(weight_crop, header=hdr_w)
hdu.writeto(image_dir / 'test_Y_wht.fits', overwrite=True)