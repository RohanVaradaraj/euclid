from astropy.table import Table
from astropy.coordinates import SkyCoord
import astropy.units as u
from pathlib import Path
from regions import Regions
from astropy.io import fits
from astropy.wcs import WCS
import matplotlib.pyplot as plt
import numpy as np

plt.rcParams['axes.linewidth'] = 2.5
plt.rcParams.update({'font.size': 20})
plt.rcParams['figure.dpi'] = 100

plt.rcParams.update({
    # Ticks on all sides, pointing inwards
    'xtick.top': True, 'xtick.bottom': True,
    'ytick.left': True, 'ytick.right': True,
    'xtick.direction': 'in', 'ytick.direction': 'in',

    # Major tick size and width
    'xtick.major.size': 6.5, 'ytick.major.size': 6.5,
    'xtick.major.width': 3, 'ytick.major.width': 3,

    # Minor tick size and width
    'xtick.minor.size': 3, 'ytick.minor.size': 3,
    'xtick.minor.width': 2, 'ytick.minor.width': 2,
})

# === File paths ===
cat_dir = Path.cwd().parents[1] / 'data' / 'catalogues' / 'candidates'
plot_dir = Path.cwd().parents[1] / 'plots' / 'sample'

ultravista_file = 'COSMOS_5sig_Y_J_nonDet_HSC_G_nonDet_HSC_R_nonDet_HSC_I_candidates_2025_02_14.fits'
euclid_file = 'COSMOS_5sig_Y_J_nonDet_HSC_G_nonDet_HSC_R_nonDet_HSC_I_candidates_2025_02_14_with_euclid.fits'

euclid_dir = Path.cwd().parents[3] / 'data' / 'euclid' / 'images'
r = Regions.read(euclid_dir / 'euclid.reg')
print(r)

#! Read in a wcs
with fits.open(euclid_dir / 'COSMOS_H_resamp.fits') as hdul:
    wcs = WCS(hdul[0].header)

# === Load the catalogs ===
ultravista = Table.read(cat_dir / ultravista_file)
euclid = Table.read(cat_dir / euclid_file)

print(ultravista.colnames)

# === Create SkyCoord objects for crossmatching ===
uv_coords = SkyCoord(ra=ultravista['RA'], dec=ultravista['DEC'])
eu_coords = SkyCoord(ra=euclid['RA'], dec=euclid['DEC'])


#! MASKING UVISTA TO EUCLID FOOTPRINT
# === Identify region types ===
polygon = r[0]
box = r[1]
print(f"Polygon: {polygon}, Box: {box}")

# Convert to shapely geometries

# === Create masks ===
in_polygon = polygon.contains(uv_coords, wcs=wcs)
in_box = box.contains(uv_coords, wcs=wcs)

# === Final mask: inside polygon AND NOT inside box ===
final_mask = in_polygon & ~in_box

# === Apply the mask ===
print(len(ultravista), "sources in UltraVISTA catalog before masking.")
ultravista = ultravista[final_mask]
print(len(ultravista), "sources in UltraVISTA catalog after masking.")

# Mask the coordinates as well
uv_coords = uv_coords[final_mask]

#! CROSSMATCHING
# === Crossmatch: UltraVISTA to Euclid ===
idx_uv2eu, d2d_uv2eu, _ = uv_coords.match_to_catalog_sky(eu_coords)
match_mask_uv = d2d_uv2eu < 1.0 * u.arcsec

# === (a) UltraVISTA-only sources ===
ultravista_only = ultravista[~match_mask_uv]

# === Crossmatch: Euclid to UltraVISTA ===
idx_eu2uv, d2d_eu2uv, _ = eu_coords.match_to_catalog_sky(uv_coords)
match_mask_eu = d2d_eu2uv < 1.0 * u.arcsec

# === (b) Euclid-only sources ===
euclid_only = euclid[~match_mask_eu]

# == (c) in both ==
both_uv = ultravista[match_mask_uv]
both_eu = euclid[idx_uv2eu[match_mask_uv]]

# === Save output ===
# ultravista_only.write('ultravista_only.fits', overwrite=True)
# euclid_only.write('euclid_only.fits', overwrite=True)

print(f"Done! Found {len(ultravista_only)} UltraVISTA-only sources.")
print(f"Found {len(euclid_only)} UltraVISTA+Euclid-only sources.")
print(f"Found {len(both_uv)} sources in both UltraVISTA and Euclid.")

#! COMPARE Muv

# Filter catalogues so that Muv<0
both_uv = both_uv[both_eu['Muv'] < 0]
both_eu = both_eu[both_eu['Muv'] < 0]

# Set up the LF bins
Muv_bins = [-22.4, -22., -21.8, -21.6, -21.4, -21.2, -21.0, -20.8, -20.6, -20.4, -20.2]
Muv_uv = np.array(both_uv['Muv'])
Muv_eu = np.array(both_eu['Muv'])
uv_bin_indices = np.digitize(Muv_uv, Muv_bins) - 1
eu_bin_indices = np.digitize(Muv_eu, Muv_bins) - 1
bin_shift_mask = uv_bin_indices != eu_bin_indices

# compare Muv for sources in both catalogues
# Plot all points in light grey as base
#plt.scatter(Muv_uv, Muv_eu, s=50, alpha=0.6, color='grey')

# Highlight only those that shift bins
plt.figure(figsize=(8, 6))
plt.scatter(Muv_uv[bin_shift_mask], Muv_eu[bin_shift_mask], s=50, alpha=1, color='black', zorder=10)

# Error bars
plt.errorbar(Muv_uv, Muv_eu, 
             xerr=[both_uv['dMuv_inf'], both_uv['dMuv_sup']], 
             yerr=[both_eu['dMuv_inf'], both_eu['dMuv_sup']], 
             fmt='none', alpha=0.5, color='black')

# Extract Muv from the unmatched sources
uvista_only_muv = np.array(ultravista_only['Muv'])
euclid_only_muv = np.array(euclid_only['Muv'])

# Apply Muv < 0 filter (like before)
uvista_only_muv = uvista_only_muv[uvista_only_muv < 0]
euclid_only_muv = euclid_only_muv[euclid_only_muv < 0]

# # === UltraVISTA-only: horizontal markers at y = lower edge ===
# plt.scatter(uvista_only_muv, [-19.95]*len(uvista_only_muv), marker='^', color='dodgerblue', label='UltraVISTA-only', s=50, alpha=0.7)

# # === Euclid-only: vertical markers at x = lower edge ===
# plt.scatter([-19.95]*len(euclid_only_muv), euclid_only_muv, marker='>', color='magenta', label='Euclid-only', s=50, alpha=0.7)

# Bin lines
# for b in Muv_bins:
#     plt.axvline(b, linestyle=':', color='gray', alpha=0.5)
#     plt.axhline(b, linestyle=':', color='gray', alpha=0.5)

# Diagonal line for equality
lims = [-22.7, -19.9]  # Set limits for the diagonal line
plt.plot(lims, lims, color='black', linestyle='--')

plt.minorticks_on()    

plt.xlabel(r'$M_{\rm{UV}}$ (UltraVISTA-only)')
plt.ylabel(r'$M_{\rm{UV}}$ (UltraVISTA+$Euclid$)')
plt.xlim(lims)
plt.ylim(lims)
plt.gca().invert_xaxis()
plt.gca().invert_yaxis()
plt.tight_layout()
plt.savefig(plot_dir / 'compare_Muv_uvista_euclid.pdf', bbox_inches='tight')
plt.show()
exit()

#! COMPARE ZPHOT
plt.figure(figsize=(8, 8))
plt.scatter(both_uv['Zphot'], both_eu['Zphot'], s=50, alpha=1, color='black')

plt.errorbar(both_uv['Zphot'], both_eu['Zphot'], 
            xerr=[both_uv['Zphot']-both_uv['Zinf'], 
                   both_uv['Zsup']-both_uv['Zphot']],
            yerr=[both_eu['Zphot']-both_eu['Zinf'],
                     both_eu['Zsup']-both_eu['Zphot']],
             fmt='none', alpha=0.4, color='black')

plt.xlabel('Zphot (UltraVISTA)')
plt.ylabel('Zphot (Euclid)')



# Add straight line
plt.plot([6.3, 7.6], [6.3, 7.6], color='black', linestyle='--')

plt.xlim(6.3, 7.6)
plt.ylim(6.3, 7.6)

plt.show()
