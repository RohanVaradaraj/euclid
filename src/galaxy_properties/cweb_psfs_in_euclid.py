"""
Take PSFs from COSMOS Web. Measure their FWHMs in Euclid as a function of magnitude.

Created: Thursday 9th January 2025.
"""

import numpy as np
from astropy.io import fits
from pathlib import Path
from astropy.table import Table
import matplotlib.pyplot as plt
from sklearn.utils import resample


plt.rcParams['axes.linewidth'] = 2.5
plt.rcParams.update({'font.size': 15})
plt.rcParams['figure.dpi'] = 100

cat_dir = Path.cwd().parents[1] / 'data' / 'ref_catalogues' / 'f444w_stars_faint_PSFs_XMATCHWITH_Euclid_photCat_NONSAT.fits'

t = Table.read(cat_dir)

print(t.colnames)

mag_J = -2.5*np.log10(t['flux_Je']) - 48.6

fwhm = t['FWHM_IMAGE_2']

# Convert to arcsec
fwhm *= 0.1

# Bin FWHM by J mag, and compute median and std

# Bins of width 0.5 mag
bins = np.arange(20, 30, 0.5)

indices = np.digitize(mag_J, bins)
fwhm_mean = [np.mean(fwhm[indices == i]) for i in range(1, len(bins))]
fwhm_median = [np.median(fwhm[indices == i]) for i in range(1, len(bins))]
fwhm_std = [np.std(fwhm[indices == i]) for i in range(1, len(bins))]

plt.figure(figsize=(10, 7))

plt.scatter(mag_J, fwhm, color='black', alpha=0.3, s=13, edgecolor='none')
plt.xlabel(r'$J_{E}$' + ' (mag)')
plt.ylabel('FWHM (arcsec)')

# iqr = [np.percentile(fwhm[indices == i], 75) - 
#        np.percentile(fwhm[indices == i], 25) for i in range(1, len(bins))]

# Get upper and lower ranges containing 68% of the data
yerr_upp = [np.percentile(fwhm[indices == i], 84) - fwhm_median[i-1] for i in range(1, len(bins))]
yerr_low = [fwhm_median[i-1] - np.percentile(fwhm[indices == i], 16) for i in range(1, len(bins))]

# Concatenate the upper and lower errors
yerr = np.array([yerr_low, yerr_upp])

# Save fwhm median and std to npy files
np.save('cweb_psfs_euclid_fwhm_median.npy', fwhm_median)
np.save('cweb_psfs_euclid_fwhm_std.npy', fwhm_std)
np.save('cweb_psfs_euclid_fwhm_1sigma.npy', yerr)


# Plot binned FWHM
#plt.errorbar(bins[:-1], fwhm_mean, yerr=np.array(fwhm_std), fmt='o', color='red', markersize=15, elinewidth=3.5, markeredgecolor='black')
#plt.errorbar(bins[:-1], fwhm_mean, yerr=np.array(iqr)/2, fmt='o', color='red', markersize=15, elinewidth=3.5, markeredgecolor='black')
plt.errorbar(bins[:-1], fwhm_median, yerr=[yerr_low, yerr_upp], fmt='o', color='red', markersize=15, elinewidth=3.5, markeredgecolor='black')

# Horizontal line at 0.51
plt.axhline(0.51, color='deepskyblue', linestyle='--', linewidth=3.5, label=r'$\mathrm{Euclid \ PSF \ FWHM} \ (J_{E})$')

plt.xlim(21.1, 27.4)
plt.ylim(0.25, 1.75)

plt.legend()
plt.tight_layout()
plot_dir = Path.cwd().parents[1] / 'plots' / 'sizes'
if not plot_dir.exists():
    plot_dir.mkdir()
plt.savefig(plot_dir / 'cweb_psfs_euclid_fwhm_vs_Jmag_IQR.pdf')
plt.show()

