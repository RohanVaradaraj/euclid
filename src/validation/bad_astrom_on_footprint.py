"""
bad_astrom_on_footprint.py

Plot the locations of stars which have large astrometric offsets between VISTA and Euclid.

Created: Wednesday 17th April 2024.
"""

from astropy.io import fits, ascii
from astropy.table import Table
import matplotlib.pyplot as plt
from pathlib import Path
import glob
from astropy.wcs import WCS
from matplotlib.patches import Polygon
import numpy as np

# Making plots look nice
plt.rcParams['axes.linewidth'] = 2.5
plt.rcParams.update({'font.size': 15})
plt.rcParams['figure.dpi'] = 100

euclid_dir = Path.home() / 'euclid' / 'Y' / 'COSMOS'
stars_dir = Path.cwd().parent.parent / 'data' / 'ref_catalogues' / 'stars'
detector_dir = Path.home() / 'euclid' / 'DETECTOR_LAYERING'


images = glob.glob(str(euclid_dir / '*BGSUB-*'))

# Compare against VISTA...
vista = False

# Or JWST
jwst = False

# Or gaia
gaia = True

# List to store centroid coordinates and corresponding image filenames
centroid_list = []

# Create a dictionary to store the mapping of tile ID to tile index
tile_index_map = {}

for i, image in enumerate(images):

    with fits.open(image) as hdul:

        data = hdul[0].data
        header = hdul[0].header

        # Get wcs from header
        wcs = WCS(header)

        footprint = wcs.calc_footprint()

        # Calculate centroid
        centroid_x = np.mean(footprint[:, 0])
        centroid_y = np.mean(footprint[:, 1])

        # Append centroid coordinates and image filename to the list
        centroid_list.append((centroid_x, centroid_y, image))

# Sort the list based on centroid coordinates
centroid_list.sort(key=lambda x: (-x[1], -x[0]), reverse=True)

plt.figure(figsize=(10, 8))
ax = plt.subplot()

# # Plot the images in sorted order
# for i, (centroid_x, centroid_y, image) in enumerate(centroid_list):
#     with fits.open(image) as hdul:
#         tile = image.split('TILE_')[-1].split('.fits')[0]
#         tile_index_map[tile] = i + 1
#         header = hdul[0].header
#         wcs = WCS(header)
#         footprint = wcs.calc_footprint()
#         r = Polygon(np.array(footprint), closed=True, edgecolor='r', facecolor='none', lw=2.5, alpha=0.8)
#         ax.add_patch(r)
#         ax.text(centroid_x, centroid_y, str(i + 1), color='black', ha='center', va='center')


# Save the tile_index_map dictionary
np.save(Path.cwd().parent.parent / 'data' / 'mosaic' / 'tile_index_map.npy', tile_index_map)

#! Also add the UltraVISTA tile and PRIMER tile
uvista_Y = Path.home().parent.parent / 'vardy' / 'vardygroupshare' / 'data' / 'COSMOS' / 'UVISTA_Y_DR6_cropped.fits'
with fits.open(uvista_Y) as hdu_uvista:
    wcs_uvista = WCS(hdu_uvista[0].header)
    uvista_footprint = wcs_uvista.calc_footprint()

print(uvista_footprint)

r = Polygon(np.array(uvista_footprint), closed=True, edgecolor='b', facecolor='none', lw=2.5, label='UltraVISTA', alpha=0.8)
ax.add_patch(r)


#! Add the Euclid mask which is roughly the footprint
e = Polygon(np.array([[150.6599884,2.6947483], [149.8463557,2.9082861], [149.6063068,2.0034576], [150.4217294,1.7836002]]), closed=True, edgecolor='r', facecolor='none', lw=2.5, label='Euclid', alpha=0.8)
ax.add_patch(e)

#! Plot positions of large offset stars
filter_name = 'Y'
if vista:
    stars = ascii.read(stars_dir / f'{filter_name}_outside_pixscale_vista_euclid_coords.ascii')
    stars_ra = stars['RA_vista']
    stars_dec = stars['DEC_vista']
if jwst:
    stars = ascii.read(stars_dir / f'{filter_name}_outside_fwhm_jwst_euclid_coords.ascii')
    stars_ra = stars['RA_jwst']
    stars_dec = stars['DEC_jwst']
if gaia:
    stars = ascii.read(stars_dir / f'{filter_name}_outside_fwhm_euclid_gaia_coords.ascii') #! Outside fwhm
    stars_ra = stars['RA_Gaia']
    stars_dec = stars['DEC_Gaia']
    stars_pmra = stars['pmra']
    stars_pmdec = stars['pmdec']
    stars_pm = stars['pm']
    stars = Table.read(stars_dir / 'Y_all_gaia_euclid_outside_pixScale.fits') #! Outside pixel scale

    stars_ra_pix = stars['RA_1']
    stars_dec_pix = stars['DEC_1']
    stars_pmra_pix = stars['pmra']
    stars_pmdec_pix = stars['pmdec']
    stars_pm_pix = np.array([np.sqrt(x**2 + y**2) for x, y in zip(stars_pmra_pix, stars_pmdec_pix)])


    # for pm values, replace '--' with 0
    stars_pmra = np.array([0 if x == '--' else x for x in stars_pmra])
    stars_pmdec = np.array([0 if x == '--' else x for x in stars_pmdec])
    stars_pm = np.array([0 if x == '--' else x for x in stars_pm])

    # Convert the values in the arrays to floats
    stars_pmra = stars_pmra.astype(float)
    stars_pmdec = stars_pmdec.astype(float)
    stars_pm = stars_pm.astype(float)

    # stars_pix = ascii.read(stars_dir / f'{filter_name}_outside_pixscale_euclid_gaia_coords.ascii') #! Outside pixel scale
    # stars_ra_pix = stars_pix['RA_Gaia']
    # stars_dec_pix = stars_pix['DEC_Gaia']
    # stars_pmra_pix = stars_pix['pmra']
    # stars_pmdec_pix = stars_pix['pmdec']
    # stars_pm_pix = stars_pix['pm']

    # # for pm values, replace '--' with 0
    # stars_pmra_pix = np.array([0 if x == '--' else x for x in stars_pmra_pix])
    # stars_pmdec_pix = np.array([0 if x == '--' else x for x in stars_pmdec_pix])
    # stars_pm_pix = np.array([0 if x == '--' else x for x in stars_pm_pix])

    # # Convert the values in the arrays to floats
    # stars_pmra_pix = stars_pmra_pix.astype(float)
    # stars_pmdec_pix = stars_pmdec_pix.astype(float)
    # stars_pm_pix = stars_pm_pix.astype(float)

#! Plot the ultra-deep stripes
strip1 = Polygon(np.array([[150.43, 2.76], [150.57, 2.76], [150.57, 1.66], [150.43, 1.66]]), closed=True, edgecolor='orange', facecolor='none', lw=2.5, alpha=0.8)
strip2 = Polygon(np.array([[150.06, 2.76], [150.20, 2.76], [150.20, 1.66], [150.06, 1.66]]), closed=True, edgecolor='orange', facecolor='none', lw=2.5, alpha=0.8)
strip3 = Polygon(np.array([[149.70, 2.76], [149.83, 2.76], [149.83, 1.66], [149.70, 1.66]]), closed=True, edgecolor='orange', facecolor='none', lw=2.5, alpha=0.8)
strip4 = Polygon(np.array([[149.33, 2.76], [149.46, 2.76], [149.46, 1.66], [149.33, 1.66]]), closed=True, edgecolor='orange', facecolor='none', lw=2.5, alpha=0.8, label='Ultra-deep stripes')

ax.add_patch(strip1)
ax.add_patch(strip2)
ax.add_patch(strip3)
ax.add_patch(strip4)
                 

#! Add euclid detector layering

files = glob.glob(str(detector_dir / '*fits'))


counter = 0

for i, file_name in enumerate(files):

    print(file_name)

    t = Table.read(detector_dir / file_name)
    counter += len(t)

    for j, row in enumerate(t):

        # Get polygon from 'POLYGON' column
        polygon = row['POLYGON']

        # Convert from radians to degrees
        polygon = np.degrees(polygon)

        # Now it is in form [x1, y1, x2, y2, x3, y3, x4, y4]. Plot as polygon
        p = Polygon(np.array(polygon).reshape(4, 2), closed=True, edgecolor='gray', facecolor='none', lw=1, alpha=0.002)
        ax.add_patch(p)

print('Number of exposures: ', counter)


# Plot the bad stars
#ax.scatter(stars_ra, stars_dec, marker='x', color='black', s=10, label=f'astrom. offset > {filter_name} PSF FWHM')

if gaia:
    ax.scatter(stars_ra_pix, stars_dec_pix, marker='x', color='red', s=10, label=f'astrom. offset > {filter_name} pix scale (0.1\")', zorder=10)
    factor = 1e-3
    ax.quiver(stars_ra_pix, stars_dec_pix, stars_pmra_pix*factor, stars_pmdec_pix*factor, color='black', scale=1, width=0.002, headwidth=3, headlength=3, headaxislength=3)

# Set axis limits and labels
ax.set_xlim(150.9, 149.1)
ax.set_ylim(1.5, 3.2)
ax.set_xlabel('RA (deg)')
ax.set_ylabel('DEC (deg)')
ax.legend(loc='upper right')
plt.gca().set_aspect('equal', adjustable='box')
plt.tight_layout()
# if vista:
#     plt.savefig(Path.cwd().parent.parent / 'plots' / 'mosaic' / f'{filter_name}_bad_astrom_on_footprint_vista.png')
# if jwst:
#     plt.savefig(Path.cwd().parent.parent / 'plots' / 'mosaic' / f'{filter_name}_bad_astrom_on_footprint_jwst.png')
#if gaia:
#    plt.savefig(Path.cwd().parent.parent / 'plots' / 'mosaic' / f'{filter_name}_bad_astrom_on_footprint_gaia.png')
plt.show()



