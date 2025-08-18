"""
Main feedback from Adam Muzzin: since some galaxies are extended, aperture fluxes might be underestimated.
So try Kron photometry.

I never ran the original catalogues with MAG_AUTO, like an idiot.
So I reran them and have nice new catalogues. MAG_AUTO is exactly a Kron aperture measurement.
In TOPCAT, I have crossmatched this catalogue with our candidates.
This script will then compare the aperture fluxes to Kron to see the difference.

Aim a (magnitude dependent?) correction for Muv.

Created: Friday 15th August 2025.
"""

from astropy.table import Table
from pathlib import Path
import numpy as np
from astropy.coordinates import SkyCoord
import astropy.units as u
import matplotlib.pyplot as plt
from astropy import cosmology as cosmo
import pickle

# Plotting configuration
plt.rcParams.update({
    'axes.linewidth': 2.5,
    'font.size': 15,
    'figure.dpi': 100
})

plt.rcParams.update({
    # Ticks on all sides, pointing inwards
    'xtick.top': True, 'xtick.bottom': True,
    # 'ytick.left': True, 'ytick.right': True,
    'xtick.direction': 'in', 'ytick.direction': 'in',

    # Major tick size and width
    'xtick.major.size': 6.5, 'ytick.major.size': 6.5,
    'xtick.major.width': 2, 'ytick.major.width': 2,

    # Minor tick size and width
    'xtick.minor.size': 3, 'ytick.minor.size': 3,
    'xtick.minor.width': 1.5, 'ytick.minor.width': 1.5,
})

def flux_to_mag(flux):
    return -2.5 * np.log10(flux) - 48.6

def mag_to_flux(mag):
    return 10 ** ((-48.6 - mag) / 2.5)

def crossmatch_and_get_fwhms(candidates, depth_cat, filter_name, max_distance=1 * u.arcsec):
    """Crossmatch candidates with the depth catalog and return FWHMs and magnitudes."""
    coords_candidates = SkyCoord(ra=candidates['RA'], dec=candidates['DEC'], frame='icrs')
    coords_depth = SkyCoord(ra=depth_cat['ALPHA_J2000'], dec=depth_cat['DELTA_J2000'], frame='icrs')
    idx, d2d, _ = coords_candidates.match_to_catalog_sky(coords_depth)
    matched_idx = idx[d2d < max_distance]
    fwhms = depth_cat['FWHM_IMAGE'][matched_idx] * 0.2  # Convert to arcsec
    flux = candidates[f'flux_{filter_name}e'][d2d < max_distance]
    mag = flux_to_mag(flux)
    return fwhms, mag, candidates[d2d < max_distance]


# paths
cat_dir = Path.cwd().parents[1] / 'data' / 'catalogues' / 'candidates'
euclid_depth_dir = Path.cwd().parents[1] / 'data' / 'depths' / 'COSMOS' / 'catalogues'

# Sample cat
cat_name = 'Euclid_UltraVISTA_z7_sample_MAG_AUTO.fits' #? U+E
cat_name = 'UltraVISTA_only_z7_sample_MAG_AUTO.fits' #? U-only
#cat_name = 'COSMOS_5sig_Y_J_nonDet_HSC_G_nonDet_HSC_R_nonDet_HSC_I_candidates_2025_02_14_with_euclid_WITH_MAG_AUTO.fits'
t_orig = Table.read(cat_dir / cat_name)
#t = t[t['Muv'] < 0]



t = t_orig

# Euclid depth cat with FWHMs
# euclid_depth_cat = Table.read(euclid_depth_dir / f'dJ.fits')

# euclid_lbg_fwhms, euclid_lbg_mag, t = crossmatch_and_get_fwhms(t_orig, euclid_depth_cat, 'J')
# print(len(euclid_lbg_fwhms), len(euclid_lbg_mag))
# extended = np.where(euclid_lbg_fwhms > 1.1)


# Take extended sources where euclid_lbg_fwhms * 0.1 > 1.2
#t = t[extended]

# Create a new column Muv_kron
t['Muv_kron'] = t['Muv']
t['Muv_linear_fit'] = t['Muv']

filter_name = 'J'

aper_mag = flux_to_mag(t[f'flux_{filter_name}'])
kron_mag = t[f'MAG_AUTO']

#! Fluxes
aper_flux = t[f'flux_{filter_name}']
kron_count = t['FLUX_AUTO']
value = -(48.6 + 30)/2.5
kron_flux = (10**value)*kron_count

ratio = kron_flux / aper_flux

# Mask out ratio > 2 and < 0.5
mask = (ratio < 2) & (ratio > 0.5)
#mask = (ratio < 10) & (ratio > 0.1)
aper_flux = aper_flux[mask]
kron_flux = kron_flux[mask]
t = t[mask]
ratio = ratio[mask]

# Get fwhms


# plt.scatter(aper_flux, ratio, s=100, alpha=0.8, c=euclid_lbg_fwhms[mask], cmap='coolwarm', edgecolor='black',
#             vmin=0.4, vmax=2)
plt.scatter(aper_flux, ratio, s=100, alpha=0.8, c=t['FWHM_IMAGE']*0.15, cmap='coolwarm', edgecolor='black',
vmin=0.6, vmax=3)

# Fit a straight line to the data
coeffs = np.polyfit(np.log10(aper_flux), ratio, 1)
x_fit = np.linspace(np.log10(min(aper_flux)), np.log10(max(aper_flux)), 100)
fit_line = np.polyval(coeffs, x_fit)

# Write the fit coefficients to a pkl file
fit_coeffs = {'slope': coeffs[0], 'intercept': coeffs[1]}
with open('kron_fit_coeffs_U_only.pkl', 'wb') as f:
    pickle.dump(fit_coeffs, f)

plt.plot(10**x_fit, fit_line, color='black', linestyle='--', linewidth=2, label='Linear fit')

# Also plot the same line, but piecewise: where the fit is < 1, plot horizontal line at 1, where the fit is > 1, plot the fit line
x_fit = np.linspace(np.log10(min(aper_flux)), np.log10(max(aper_flux)), 100)
fit_line_piecewise = np.where(fit_line < 1, 1, fit_line)
plt.plot(10**x_fit, fit_line_piecewise, color='red', linestyle='--', linewidth=2, label='Piecewise fit')
plt.legend()

# Also bin in log(flux) and calculate mean and std and plot
# bin_edges = np.logspace(np.log10(min(aper_flux)), np.log10(max(aper_flux)), 5)
# bin_centers = []
# bin_means = []
# bin_stds = []
# bin_counts = []
# for i in range(len(bin_edges) - 1):
#     mask_bin = (aper_flux >= bin_edges[i]) & (aper_flux < bin_edges[i+1])
#     if np.any(mask_bin):
#         bin_centers.append(np.mean([bin_edges[i], bin_edges[i+1]]))
#         bin_means.append(np.mean(ratio[mask_bin]))
#         bin_stds.append(np.std(ratio[mask_bin]))
#         bin_counts.append(np.sum(mask_bin))
# plt.errorbar(bin_centers, bin_means, yerr=bin_stds/np.sqrt(bin_counts),
#              fmt='o', color='red', markersize=13, markeredgecolor='black')



plt.xscale('log')
plt.xlabel('Aperture flux')
plt.ylabel('Kron flux / Aperture flux')
plt.colorbar(label='FWHM (arcsec)', cmap='coolwarm')
plt.show()

#
exit()

# # Calculate the difference
# diff = kron_mag - aper_mag

# # # Restrict table to where diff > -2, i.e. physical
# mask = diff < 2
# t = t[mask]
# aper_mag = aper_mag[mask]
# kron_mag = kron_mag[mask]
# diff = diff[mask]

# t['Muv_kron'] += diff

# print(np.array(t['ID']))
# pairs = zip(np.array(t['ID']), diff)
# print(list(pairs))


# # Choose bin edges (adjust as needed)
# bin_width = 0.25
# bins = np.arange(min(t['Muv']), max(t['Muv']) + bin_width, bin_width)

# # Assign each point to a bin
# bin_idx = np.digitize(t['Muv'], bins)

# # Prepare arrays for results
# bin_centers = []
# bin_means = []
# bin_stds = []
# bin_counts = []

# for b in range(1, len(bins)):
#     mask_bin = bin_idx == b
#     if np.any(mask_bin):
#         bin_centers.append(np.mean([bins[b-1], bins[b]]))
#         bin_means.append(np.mean(diff[mask_bin]))
#         bin_stds.append(np.std(diff[mask_bin]))
#         bin_counts.append(np.sum(mask_bin))

# bin_centers = np.array(bin_centers)
# bin_means = np.array(bin_means)
# bin_stds = np.array(bin_stds)
# bin_counts = np.array(bin_counts)

# # Plot
# plt.figure(figsize=(10, 6))
# plt.scatter(t['Muv'], diff, s=70, alpha=0.5)
# plt.errorbar(bin_centers, bin_means, yerr=bin_stds/np.sqrt(bin_counts), 
#             fmt='o', color='red', markersize=13, markeredgecolor='black')

# # Fit straight line to the binned points
# coeffs = np.polyfit(bin_centers, bin_means, 1)
# fit_line = np.polyval(coeffs, bin_centers)
# plt.plot(bin_centers, fit_line, color='black', linestyle='--', linewidth=2)

# # Add linear fit line from t['Muv_linear_fit'], only if it is in the extended sources
# # for i in range(len(t)):
# #     if extended[i]:
# #         plt.plot([t['Muv'][i], t['Muv'][i]], [t['Muv_linear_fit'][i], t['Muv_kron'][i]], color='green', linestyle='--', linewidth=1.5)

    

# plt.xlabel(r'$M_{\rm{UV}}$', size=15)
# plt.ylabel(r'$\Delta m = \rm{Kron Linear} - \rm{Aper}$', size=15)

# plt.title(filter_name)
# plt.gca().invert_xaxis()
# plt.ylim(-1, 1.3)
# plt.show()
# plt.close()

# #plt.hist(t['Muv_kron'], bins=np.arange(-24, -18, 0.25), histtype='step', lw=2, label='Kron', color='tab:red', alpha=0.7)
# plt.hist(t['Muv'], bins=np.arange(-24, -18, 0.25), histtype='step', lw=2, label='Aperture', color='tab:blue', alpha=0.7)
# plt.hist(t['Muv_kron'], bins=np.arange(-24, -18, 0.25), histtype='step', lw=2, label='Linear Fit', color='tab:green', alpha=0.7)
# plt.legend()
# plt.xlabel(r'$M_{\rm{UV}}$', size=15)
# plt.show()
# plt.close()

# plt.scatter(t['Muv'], t['Muv_linear_fit'], s=70)

# # Plot y=x 
# x = np.linspace(-22.6, -20.3, 100)
# plt.plot(x, x, color='black', linestyle='--', linewidth=2)

# plt.xlabel(r'$M_{\rm{UV}}$', size=15)
# plt.ylabel(r'$M_{\rm{UV, Kron}}$', size=15)
# plt.gca().invert_xaxis()
# plt.gca().invert_yaxis()
# plt.show()

# #t.write(cat_dir / cat_name, overwrite=True)



