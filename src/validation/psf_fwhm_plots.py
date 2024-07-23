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

stars_dir = Path.cwd().parent.parent / 'data' / 'ref_catalogues' / 'stars'
psfex_dir = Path.cwd().parent.parent / 'data' / 'psf' / 'COSMOS' / 'catalogues'
plot_dir = Path.cwd().parent.parent / 'plots' / 'psf'

# Plot pipeline or our PSF FWHMs? True for pipeline, False for our PSFEx
plot_pipe = False

filter_names = ['VIS', 'Y', 'J', 'H']

colors = ['blue', 'green', 'red', 'purple']

pix_scale = 0.1 # arcsec / pix

plt.figure(figsize=(10, 6))

for filter_name in filter_names:

    #! 1) First load the pipeline PSFs
    if plot_pipe:
        pipeline_stars = stars_dir / f'{filter_name}_euclid_pipeline_psf_fwhms.txt'
        pipeline_data = ascii.read(pipeline_stars)

        #Fit a Gaussian to the pipeline data to measure the peak
        mu, std = norm.fit(pipeline_data['FWHM'])
        
        plt.hist(pipeline_data['FWHM'], bins=np.arange(0.1, 0.6, 0.001), color=colors[filter_names.index(filter_name)], alpha=0.8, label=filter_name + f', peak={round(mu, 2)}as', density=True)

    #! 2) Now load our PSFEx PSFs
    else:
        psfex_stars = psfex_dir / f'{filter_name}_stars.fits'
        psfex_data = Table.read(psfex_stars)

        # FWHM is in pixels, convert to arcsec with pixel scale
        psfex_data['FWHM_IMAGE'] *= pix_scale

        # Fit a Gaussian to the pipeline data to measure the peak
        mu, std = norm.fit(psfex_data['FWHM_IMAGE'])

        plt.hist(psfex_data['FWHM_IMAGE'], bins=np.arange(0.1, 0.6, 0.005), color=colors[filter_names.index(filter_name)], alpha=0.8, label=filter_name + f', peak={mu:.2f}as', density=True)

plt.xlabel('FWHM (arcsec)')
plt.ylabel('Normalized count')
plt.legend()
if plot_pipe:
    plt.savefig(plot_dir / 'pipeline_psf_fwhms.png')
else:
    plt.savefig(plot_dir / 'manual_selection_psf_fwhms.png')

plt.show()