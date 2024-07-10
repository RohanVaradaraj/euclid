#!/usr/bin/env python3

"""
plot the euclid and HSC footprints and adjust to get them on top of one another.

Edit: this code just helped me realise CRPIX and CRVAL are NOT centre pix and centre val.
Instead, must be coordinate reference or something. I am stupid.

Created: Tuesday 25th June 2024
"""

from pathlib import Path
from astropy.io import fits
import matplotlib.pyplot as plt
from astropy.wcs import WCS
import numpy as np

plt.rcParams['axes.linewidth'] = 2.5
plt.rcParams.update({'font.size': 15})
plt.rcParams['figure.dpi'] = 100

data_dir = Path.cwd().parents[3] / 'data'

hsc_image = data_dir / 'COSMOS' / 'HSC-Y_DR3.fits'

vista_image = data_dir / 'COSMOS' / 'UVISTA_Y_DR6.fits'

with fits.open(vista_image) as hdul:
    euclid_data = hdul[0].data
    euclid_header = hdul[0].header

with fits.open(hsc_image) as hdul:
    hsc_data = hdul[0].data
    hsc_header = hdul[0].header

# euclid_header['CRPIX1'] = 16316.5
# euclid_header['CRPIX2'] = 14968.5

wcs_euclid = WCS(euclid_header)
wcs_hsc = WCS(hsc_header)

# Find centre RA DEC of the HSC image
ra_hsc, dec_hsc = wcs_hsc.all_pix2world(hsc_data.shape[1] / 2, hsc_data.shape[0] / 2, 0)
print(ra_hsc, dec_hsc)

# Find differences in coords in bottom left pixel of each wcs in pixels
ra_euclid, dec_euclid = wcs_euclid.all_pix2world(0, 0, 0)
ra_hsc, dec_hsc = wcs_hsc.all_pix2world(0, 0, 0)
print(ra_euclid, dec_euclid)

# Convert to pixels
ra_euclid_pix, dec_euclid_pix = wcs_euclid.all_world2pix(ra_hsc, dec_hsc, 0)
ra_hsc_pix, dec_hsc_pix = wcs_hsc.all_world2pix(ra_euclid, dec_euclid, 0)
print('this is the ra and dec in pixels', ra_euclid_pix, dec_euclid_pix)

# Find differences in pixels
diff_ra_pix = int(ra_euclid_pix - 0)
diff_dec_pix = int(dec_euclid_pix - 0)
print('this is the difference in pixels', diff_ra_pix, diff_dec_pix)


# Shift Euclid image by these differences
euclid_header['CRPIX1'] -= int(diff_ra_pix)
euclid_header['CRPIX2'] -= int(diff_dec_pix)

# Shift CRVAL by this much
ra_euclid, dec_euclid = wcs_euclid.all_pix2world(0, 0, 0)
print('this is the ra and dec in world', ra_euclid, dec_euclid)
euclid_header['CRVAL1'] = ra_euclid
euclid_header['CRVAL2'] = dec_euclid

wcs_euclid = WCS(euclid_header)

# Modify the data to match this: create a new array of zeros with the same shape as the Euclid data
# Then copy the Euclid data into the new array, shifted by the differences in pixels
print('Shifting the image')
euclid_data_shifted = np.zeros_like(euclid_data)
euclid_data_shifted[diff_dec_pix:, diff_ra_pix:] = euclid_data[:-diff_dec_pix, :-diff_ra_pix]

# Overwrite the Euclid data with the shifted data
euclid_data = euclid_data_shifted

# Overwrite and save the Euclid image
print('Writing image')
euclid_image_shifted = Path.cwd() / 'Euclid_VIS_cropped_shifted.fits'
fits.writeto(euclid_image_shifted, euclid_data, euclid_header, overwrite=True)

# Get footprints
euclid_footprint = wcs_euclid.calc_footprint()
hsc_footprint = wcs_hsc.calc_footprint()

# Make patches for footprints
euclid_patch = plt.Polygon(euclid_footprint, fill=None, edgecolor='blue', lw=2.5, label='VISTA')
hsc_patch = plt.Polygon(hsc_footprint, fill=None, edgecolor='red', lw=2.5, label='HSC Y')


# Plot footprints
plt.figure(figsize=(10, 10))
plt.gca().add_patch(euclid_patch)
plt.gca().add_patch(hsc_patch)

plt.xlabel('RA')
plt.ylabel('Dec')

plt.xlim(149, 151)
plt.ylim(1.45, 2.9)

# Reverse x-axis
plt.gca().invert_xaxis()

# Ensure axes have same scale/aspect ratio
plt.axis('equal')

plt.legend()
plt.show()


