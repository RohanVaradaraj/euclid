"""
Crossmatch the LBG and BD candidates from VISTA and Euclid with the depth catalogues to get the FWHMs.

Created: Tuesday 26th November 2024.
"""

from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table
from astropy.coordinates import SkyCoord
import astropy.units as u

plt.rcParams['axes.linewidth'] = 2.5
plt.rcParams.update({'font.size': 15})
plt.rcParams['figure.dpi'] = 100

#! Set up the directories

filter_name = 'J'

# LBG and BD samples
candidate_dir = Path.cwd().parents[1] / 'data' / 'catalogues' / 'candidates'

# Depth directories
euclid_depth_dir = Path.cwd().parents[1] / 'data' / 'depths' / 'COSMOS' / 'catalogues'
vista_depth_dir = Path.home().parent.parent / 'vardy' / 'vardygroupshare' / 'data' / 'depths' / 'COSMOS' / 'catalogues'

#! Read in LBG and BD candidates
# euclid_lbg_catname = 'COSMOS_5sig_Y_J_nonDet_HSC_G_nonDet_HSC_R_nonDet_HSC_I_candidates_2024_11_21_with_euclid.fits'
# euclid_bd_catname = 'COSMOS_5sig_Y_J_nonDet_HSC_G_nonDet_HSC_R_nonDet_HSC_I_best_bd_INTERLOPERS_2024_11_26_with_euclid.fits'

euclid_lbg_catname = 'COSMOS_5sig_Y_J_nonDet_HSC_G_nonDet_HSC_R_nonDet_HSC_I_really_good_LBGs_INTERLOPERS_2024_11_26_with_euclid.fits' #? Visually selected
euclid_bd_catname = 'COSMOS_5sig_Y_J_nonDet_HSC_G_nonDet_HSC_R_nonDet_HSC_I_really_good_BDs_INTERLOPERS_2024_11_26_with_euclid.fits'

vista_lbg_catname = 'COSMOS_5sig_Y_J_nonDet_HSC_G_nonDet_HSC_R_nonDet_HSC_I_candidates_2024_11_20.fits'
vista_bd_catname = 'COSMOS_5sig_Y_J_nonDet_HSC_G_nonDet_HSC_R_nonDet_HSC_I_BD_INTERLOPERS_2024_11_26.fits'

euclid_lbg = Table.read(candidate_dir / euclid_lbg_catname)
euclid_bd = Table.read(candidate_dir / euclid_bd_catname)
print(len(euclid_lbg), len(euclid_bd))

vista_lbg = Table.read(candidate_dir / vista_lbg_catname)
vista_bd = Table.read(candidate_dir / vista_bd_catname)

#! Read in the depth catalogues
euclid_depth_catname = f'd{filter_name}.fits'
vista_depth_catname = f'd{filter_name}.fits'

euclid_depth_cat = Table.read(euclid_depth_dir / euclid_depth_catname)
vista_depth_cat = Table.read(vista_depth_dir / vista_depth_catname)

#! Crossmatch the LBG and BD candidates with the depth catalogues
coords_euclid_lbg = SkyCoord(ra=euclid_lbg['RA'], dec=euclid_lbg['DEC'], frame='icrs')
coords_euclid_bd = SkyCoord(ra=euclid_bd['RA'], dec=euclid_bd['DEC'], frame='icrs')

coords_vista_lbg = SkyCoord(ra=vista_lbg['RA'], dec=vista_lbg['DEC'], frame='icrs')
coords_vista_bd = SkyCoord(ra=vista_bd['RA'], dec=vista_bd['DEC'], frame='icrs')

coords_euclid_depth = SkyCoord(ra=euclid_depth_cat['ALPHA_J2000'], dec=euclid_depth_cat['DELTA_J2000'], frame='icrs')
coords_vista_depth = SkyCoord(ra=vista_depth_cat['ALPHA_J2000'], dec=vista_depth_cat['DELTA_J2000'], frame='icrs')

#! CROSSMATCHING

# Maximum distance for crossmatching
max_distance = 1 * u.arcsec

#? 1) Crossmatch the Euclid LBG candidates with the Euclid depth catalogue
idx, d2d, d3d = coords_euclid_lbg.match_to_catalog_sky(coords_euclid_depth)
matched_idx1 = idx[d2d < max_distance]

# Get the FWHMs
euclid_lbg_fwhms = euclid_depth_cat['FWHM_IMAGE'][matched_idx1] * 0.1 # Convert to arcsec

# Get the J-band fmag from the LBG cat
euclid_lbg_flux = euclid_lbg[f'flux_{filter_name}e'][d2d < max_distance]
euclid_lbg_mag = -2.5 * np.log10(euclid_lbg_flux) - 48.6

#? 2) Crossmatch the Euclid BD candidates with the Euclid depth catalogue
idx, d2d, d3d = coords_euclid_bd.match_to_catalog_sky(coords_euclid_depth)
matched_idx2 = idx[d2d < max_distance]

# Get the FWHMs
euclid_bd_fwhms = euclid_depth_cat['FWHM_IMAGE'][matched_idx2] * 0.1 # Convert to arcsec

# Get the J-band fmag from the BD cat
euclid_bd_flux = euclid_bd[f'flux_{filter_name}e'][d2d < max_distance]
euclid_bd_mag = -2.5 * np.log10(euclid_bd_flux) - 48.6

#? 3) Crossmatch the Vista LBG candidates with the Vista depth catalogue
idx, d2d, d3d = coords_vista_lbg.match_to_catalog_sky(coords_vista_depth)
matched_idx3 = idx[d2d < max_distance]

# Get the FWHMs
vista_lbg_fwhms = vista_depth_cat['FWHM_IMAGE'][matched_idx3] * 0.15 # Convert to arcsec

# Get the J-band fmag from the LBG cat
vista_lbg_flux = vista_lbg[f'flux_{filter_name}e'][d2d < max_distance]
vista_lbg_mag = -2.5 * np.log10(vista_lbg_flux) - 48.6

#? 4) Crossmatch the Vista BD candidates with the Vista depth catalogue
idx, d2d, d3d = coords_vista_bd.match_to_catalog_sky(coords_vista_depth)
matched_idx4 = idx[d2d < max_distance]

# Get the FWHMs
vista_bd_fwhms = vista_depth_cat['FWHM_IMAGE'][matched_idx4] * 0.15 # Convert to arcsec

# Get the J-band fmag from the BD cat
vista_bd_flux = vista_bd[f'flux_{filter_name}e'][d2d < max_distance]
vista_bd_mag = -2.5 * np.log10(vista_bd_flux) - 48.6

#! Plot the mag vs FHWM
plt.figure(figsize=(10, 8))
plt.scatter(euclid_lbg_mag, euclid_lbg_fwhms, label='Euclid LBG', color='blue', marker='o')
plt.scatter(euclid_bd_mag, euclid_bd_fwhms, label='Euclid BD', color='red', marker='o')
plt.xlabel(f'{filter_name} mag')
plt.ylabel('FWHM (arcsec)') 
plt.xlim(23, 28)
plt.show()
plt.close()

plt.figure(figsize=(10, 8))
plt.scatter(vista_lbg_mag, vista_lbg_fwhms, label='Vista LBG', color='blue', marker='o')
plt.scatter(vista_bd_mag, vista_bd_fwhms, label='Vista BD', color='red', marker='o')
plt.xlabel(f'{filter_name} mag')
plt.ylabel('FWHM (arcsec)') 
plt.xlim(23, 28)
plt.show()
plt.close()









