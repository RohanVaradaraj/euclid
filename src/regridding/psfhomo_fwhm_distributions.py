"""
plot the fwhms of the PSFs in the PSF-homogenised euclid and jwst images.

Created: Monday 8th July 2024.
"""

"""
psf_fwhm_plots.py

Plot the FWHM distributions of the stars in Euclid, compare with Euclid pipeline PSFEx outputs.

Created: Monday 15th April 2024
"""

from astropy.io import ascii
from astropy.table import Table
import matplotlib.pyplot as plt
from pathlib import Path
import numpy as np
from scipy.stats import norm

plt.rcParams['axes.linewidth'] = 2.5
plt.rcParams.update({'font.size': 15})
plt.rcParams['figure.dpi'] = 100

psfex_dir = Path.cwd().parents[3] / 'data' / 'psf' / 'COSMOS' / 'catalogues'
plot_dir = Path.cwd().parent.parent / 'plots' / 'psf'

# Plot pipeline or our PSF FWHMs? True for pipeline, False for our PSFEx
plot_pipe = False

filter_names = ['VIS', 'Ye', 'Je', 'He']

compare_names = ['HSC-I_DR3', 'Y_DR6', 'J_DR6', 'H_DR6']

colors = ['blue', 'green', 'red', 'purple']

pix_scale = 0.15 # arcsec / pix

plt.figure(figsize=(10, 6))

for i, filter_name in enumerate(filter_names):

    psfex_stars = psfex_dir / f'{filter_name}_stars.fits'
    psfex_data = Table.read(psfex_stars)

    # FWHM is in pixels, convert to arcsec with pixel scale
    psfex_data['FWHM_IMAGE'] *= pix_scale
    print(psfex_data['FWHM_IMAGE'])

    # Fit a Gaussian to the pipeline data to measure the peak
    mu, std = norm.fit(psfex_data['FWHM_IMAGE'])

    plt.hist(psfex_data['FWHM_IMAGE'], bins=np.arange(0.6, 1, 0.005), color=colors[filter_names.index(filter_name)], alpha=0.8, label=filter_name + f', peak={mu:.2f}as', density=True, histtype='step')


compare_stars = psfex_dir / f'Y_DR6_stars.fits'
compare_data = Table.read(compare_stars)
compare_data['FWHM_IMAGE'] *= pix_scale
print(compare_data['FWHM_IMAGE'])
plt.hist(compare_data['FWHM_IMAGE'], bins=np.arange(0.6, 1, 0.005), color='gray', alpha=0.2, density=True, histtype='step')

plt.xlabel('FWHM (arcsec)')
plt.ylabel('Normalized count')
plt.legend()
plt.show()