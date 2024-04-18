"""
bad_astrom_on_footprint.py

Plot the locations of stars which have large astrometric offsets between VISTA and Euclid.

Created: Wednesday 17th April 2024.
"""

from astropy.io import fits, ascii
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

images = glob.glob(str(euclid_dir / '*BGSUB-*'))

# Compare against VISTA or JWST? If true, compare against VISTA
vista = False

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
else:
    stars = ascii.read(stars_dir / f'{filter_name}_outside_pixscale_jwst_euclid_coords.ascii')
    stars_ra = stars['RA_jwst']
    stars_dec = stars['DEC_jwst']

#! Plot the ultra-deep stripes
strip1 = Polygon(np.array([[150.43, 2.76], [150.57, 2.76], [150.57, 1.66], [150.43, 1.66]]), closed=True, edgecolor='orange', facecolor='none', lw=2.5, alpha=0.8)
strip2 = Polygon(np.array([[150.06, 2.76], [150.20, 2.76], [150.20, 1.66], [150.06, 1.66]]), closed=True, edgecolor='orange', facecolor='none', lw=2.5, alpha=0.8)
strip3 = Polygon(np.array([[149.70, 2.76], [149.83, 2.76], [149.83, 1.66], [149.70, 1.66]]), closed=True, edgecolor='orange', facecolor='none', lw=2.5, alpha=0.8)
strip4 = Polygon(np.array([[149.33, 2.76], [149.46, 2.76], [149.46, 1.66], [149.33, 1.66]]), closed=True, edgecolor='orange', facecolor='none', lw=2.5, alpha=0.8, label='Ultra-deep stripes')

ax.add_patch(strip1)
ax.add_patch(strip2)
ax.add_patch(strip3)
ax.add_patch(strip4)
                 

# Plot the bad stars
ax.scatter(stars_ra, stars_dec, marker='x', color='black', s=10, label=f'astrom. offset > {filter_name} PSF FWHM')

# Set axis limits and labels
ax.set_xlim(150.9, 149.1)
ax.set_ylim(1.5, 3.2)
ax.set_xlabel('RA (deg)')
ax.set_ylabel('DEC (deg)')
ax.legend(loc='upper right')
plt.gca().set_aspect('equal', adjustable='box')
plt.tight_layout()
plt.savefig(Path.cwd().parent.parent / 'plots' / 'mosaic' / f'{filter_name}_bad_astrom_on_footprint.png')
plt.show()



