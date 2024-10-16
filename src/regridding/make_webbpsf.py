"""
Use WEBB PSF to generate PSFs for the CWEB imaging.

#!Run in the conda env webbpsf_env

Created: Friday 26th July 2024.
"""

import webbpsf
from astropy.io import fits
import matplotlib.pyplot as plt
from pathlib import Path
from matplotlib.colors import LogNorm
from scipy.ndimage import rotate
import numpy as np

# Load the CWEB image to get the header
image_dir = Path.cwd().parents[2] / 'images' / 'CWEB' / 'mosaic_0A'

psf_dir = Path.cwd().parents[1] / 'data' / 'psf' / 'COSMOS' / 'results'

filename = 'CWEB-F150W-0A.fits'
with fits.open(image_dir / filename) as hdul:
    header = hdul[0].header
    # Assuming the position angle is stored in a keyword like 'PA_V3' (check your actual header keyword)
    pa_v3 = header['PA_V3']
    print(f'PA_V3: {pa_v3}')

# Define filters to process
filters = ['F115W', 'F150W', 'F277W', 'F444W']
psf_data_rotated_dict = {}

# Calculate and rotate PSFs for each filter
for filter_name in filters:
    nc = webbpsf.NIRCam()
    nc.filter = filter_name
    
    # Set the position angle
    nc.pupil_rotation = -pa_v3
    
    # Calculate the PSF
    if filter_name[1] == '1':
        psf = nc.calc_psf(oversample=1, fov_arcsec=6.0)
    else:
        psf = nc.calc_psf(oversample=2, fov_arcsec=6.0)
    
    # Extract PSF data
    psf_data = psf[0].data
    
    # Rotate the PSF data to match the position angle
    angle_degrees = -pa_v3
    psf_data_rotated = rotate(psf_data, angle_degrees, reshape=False, order=1, mode='constant')
    
    psf_data_rotated_dict[filter_name] = psf_data_rotated

    # Save the PSF data to a FITS file
    psf_filename = f'webbpsf_{filter_name}.fits'
    fits.writeto(psf_dir / psf_filename, psf_data_rotated, overwrite=True)


# Plotting the PSFs side by side
fig, axes = plt.subplots(1, 4, figsize=(20, 5), sharex=True, sharey=True)

for ax, filter_name in zip(axes, filters):
    ax.imshow(psf_data_rotated_dict[filter_name], cmap='viridis', norm=LogNorm())
    ax.set_title(f'Filter: {filter_name}')
    ax.set_xlabel('X pixel')
    ax.set_ylabel('Y pixel')

plt.tight_layout()
plt.colorbar(axes[0].images[0], ax=axes, orientation='vertical', fraction=0.02, pad=0.04)
plt.show()
