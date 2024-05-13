"""
correlate_euclid_chips.py

For stars with astrometric offset larger than the pixel scale but smaller than the PSF FWHM, we want to see if they are aligned with the Euclid chips.
In this code we run some tests for this!

Created: Wednesday 8th May 2024.
"""

from astropy.io import ascii, fits
from astropy.table import Table
import numpy as np
from scipy.spatial import cKDTree
import matplotlib.pyplot as plt
from pathlib import Path
from astroML.correlation import two_point
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
from scipy.stats import kstest, ks_2samp

# Making plots look nice
plt.rcParams['axes.linewidth'] = 2.5
plt.rcParams.update({'font.size': 15})
plt.rcParams['figure.dpi'] = 100

def KS_D_test(sig, n, m):
    val = np.sqrt( - np.log(sig/2) * ((1+m/n)/(2*m)) )
    return val

# Set up directories
stars_dir = Path.cwd().parent.parent / 'data' / 'ref_catalogues' / 'stars'
euclid_dir = Path.home() / 'euclid' / 'COSMOS'

# Open the star coordinates
#t = ascii.read(stars_dir / 'Y_euclid_gaia_coords.ascii')
t = Table.read(stars_dir / 'Y_all_gaia_euclid_outside_pixScale.fits')

# Create a KDTree from the star coordinates
#stars = np.array([t['RA_euclid'], t['DEC_euclid']]).T
stars = np.array([t['RA_1'], t['DEC_1']]).T


#########################################################!
#######! TWO POINT CORRELATION FUNCTION    ##############
#########################################################!

#! Random points
np.random.seed(0)
X = np.random.random((len(t), 2))
bins = np.linspace(0, 0.73, 30)
corr = two_point(X, bins)
print(np.allclose(corr, 0, atol=0.02))

#! Stars
star_corr = two_point(stars, bins)
print(np.allclose(star_corr, 0, atol=0.02))

# Plot
plt.figure(figsize=(10, 8))
plt.plot(bins[:-1], star_corr, marker='s', label='Stars')
plt.plot(bins[:-1], corr, marker='o', label='Random points')
plt.legend()
#plt.yscale('log')
plt.xlabel(r'$\theta$ (degrees)')
plt.ylabel('Two-point correlation function')
plt.tight_layout()
plt.show()
exit()

#########################################################!
#######! GETTING WEIGHT VALUES AT THE STARS   ############
#########################################################!

hdu = fits.open(euclid_dir / 'COSMOS_Y_MOSAIC_WHT.fits', memmap=True)

header = hdu[0].header
data = hdu[0].data

wcs = WCS(header)

# Get the pixel coordinates of the stars
star_pix_coords = wcs.all_world2pix(stars, 0)

# Don't want to measure right on top of stars - go 20 pixels up and right
star_pix_coords = np.round(star_pix_coords + 20).astype(int)

# Get the weight values at the star coordinates
star_weights = data[star_pix_coords[:, 1].astype(int), star_pix_coords[:, 0].astype(int)]

# Given the shape of the data, generate random pixel coordinates to get the weight values
np.random.seed(42)
rand_pix_coords = np.random.randint(0, data.shape[1], (len(stars), 2))
rand_weights = data[rand_pix_coords[:, 1], rand_pix_coords[:, 0]]

# Remove zeros from both weight arrays
star_weights = star_weights[star_weights != 0.]
rand_weights = rand_weights[rand_weights != 0.]

bins=np.linspace(0.1, 5.5, 15)

#! Extracting the weight values
# Plot the weight values
plt.figure(figsize=(10, 8))
plt.hist(star_weights, bins=bins, histtype='step', lw=2.5, label='Stars', density=False)
plt.hist(rand_weights, bins=bins, histtype='step', lw=2.5, label='Random points', density=False, linestyle='--')
plt.legend()
plt.xlabel('Weight map value')
plt.ylabel('N')
#plt.yscale('log')
plt.show()

#! Run a Kolmogorov-Smirnov test on the weight values
ks_stat, ks_pval = ks_2samp(star_weights, rand_weights)
#ks_stat, ks_pval = ks_2samp(rand_weights, rand_weights)

print(f'KS statistic: {ks_stat}')
print(f'KS p-value: {ks_pval}')

d_value = KS_D_test(0.05, len(star_weights), len(rand_weights))

print('KS D test value:', d_value)
print('Distributions are significantly different: ', d_value > ks_stat)
