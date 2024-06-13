"""
plot_mosaic.py

The Euclid pointing in COSMOS is split into 22 mosaic tiles. This code will plot the tiles and assign each tile a number 1 through 22. 
The tiles are numbered from left to right, top to bottom, i.e. in increasing order of RA and decreasing order of DEC.

Created: Wednesday 20th March 2024.
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

images = glob.glob(str(euclid_dir / '*BGSUB-*'))

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

# Modify shape so that one degree is the same on both axes
#plt.gca().set_aspect('equal', adjustable='box')
#ax.set_aspect(aspect='equal', adjustable='box')
plt.axis('equal')

# # Plot the images in sorted order
# for i, (centroid_x, centroid_y, image) in enumerate(centroid_list):
#     with fits.open(image) as hdul:
        
#         tile = image.split('TILE_')[-1].split('.fits')[0]
#         tile_index_map[tile] = i + 1

#         header = hdul[0].header

#         wcs = WCS(header)
#         footprint = wcs.calc_footprint()

#         ax = plt.subplot()

#         r = Polygon(np.array(footprint), closed=True, edgecolor='black', facecolor='none', lw=2.5, alpha=0.2)
#         ax.add_patch(r)

# Euclid mosaic tiles
#dummy = Polygon([[0,0],[0,0],[0,0],[0,0]], closed=True, edgecolor='gray', facecolor='none', lw=2.5, alpha=0.1, linestyle='dotted', label='OTF mosaic tiles')
#ax.add_patch(dummy)

#! Add the Euclid mask which is roughly the footprint
e = Polygon(np.array([[150.6599884,2.6947483], [149.8463557,2.9082861], [149.6063068,2.0034576], [150.4217294,1.7836002]]), closed=True, edgecolor='blue', facecolor='none', lw=2.5, label='Euclid', alpha=0.8)
ax.add_patch(e)

# Save the tile_index_map dictionary
np.save(Path.cwd().parent.parent / 'data' / 'mosaic' / 'tile_index_map.npy', tile_index_map)

#! Also add the UltraVISTA tile and PRIMER tile
uvista_Y = Path.home().parent.parent / 'vardy' / 'vardygroupshare' / 'data' / 'COSMOS' / 'UVISTA_Y_DR6_cropped.fits'
uvista_Y = Path.home().parent.parent / 'vardy' / 'vardygroupshare' / 'data' / 'COSMOS' / 'UVISTA_Y_dr5_rc1.fits'
with fits.open(uvista_Y) as hdu_uvista:
    wcs_uvista = WCS(hdu_uvista[0].header)
    uvista_footprint = wcs_uvista.calc_footprint()

primer = Path.home() / 'JWST' / 'primer_cosmos_nircam_v0.5_f444w_30mas_sci.fits'
with fits.open(primer) as hdu_primer:
    wcs_primer = WCS(hdu_primer[0].header)
    primer_footprint = wcs_primer.calc_footprint()

r = Polygon(np.array(uvista_footprint), closed=True, edgecolor='none', facecolor='darkgray', lw=2.5, label='UltraVISTA', alpha=0.8, zorder=-1) # Gray for Rebecca
#r = Polygon(np.array(uvista_footprint), closed=True, edgecolor='b', facecolor='none', lw=2.5, label='UltraVISTA', alpha=0.8)

ax.add_patch(r)
#p = Polygon(np.array(primer_footprint), closed=True, edgecolor='g', facecolor='none', lw=2.5, label='PRIMER', alpha=0.8)
#ax.add_patch(p)

#! Add the COSMOS-Web tiles
jwst_dir = Path.cwd().parents[3] / 'data' / 'CWEB'
jwst_files = glob.glob(str(jwst_dir / '*' / '*F444*'))

for jwst_file in jwst_files:

    # Get tile name from directory the image is in (one level up)
    tile = jwst_file.split('/')[-2].split('mosaic_')[-1]

    with fits.open(jwst_file) as hdu_jwst:
        wcs_jwst = WCS(hdu_jwst[1].header)
        jwst_footprint = wcs_jwst.calc_footprint()
        r = Polygon(np.array(jwst_footprint), closed=True, edgecolor='red', facecolor='none', lw=2.5, alpha=0.8)
        ax.add_patch(r)

        # Add the tile number at centre of polygon
        centroid_x = np.mean(jwst_footprint[:, 0])
        centroid_y = np.mean(jwst_footprint[:, 1])
        ax.text(centroid_x, centroid_y, tile, fontsize=8, color='black', ha='center', va='center')

dummy = Polygon([[0,0], [1,1], [1,0], [0,1]], closed=True, edgecolor='red', facecolor='none', lw=2.5, label='CWEB')
ax.add_patch(dummy)

#! Add the Hubble 3D-DASH footprint
dash = Path.home().parent.parent / 'hoy' / 'DASH' / 'hlsp_3d-dash_hst_wfc3_combined-cosmos_f160w_v1.0_drz-sci.fits'
with fits.open(dash) as hdu_dash:
    wcs_dash = WCS(hdu_dash[0].header)
    hubble_footprint = wcs_dash.calc_footprint()

    #r = Polygon(np.array(hubble_footprint), closed=True, edgecolor='deepskyblue', facecolor='none', lw=2.5, label='HST 3D-DASH', alpha=0.6, linestyle='dashed')
    #ax.add_patch(r)


#! Plot positions of REBELS galaxies
# rebels = ascii.read(Path.cwd().parent.parent / 'data' / 'mosaic' / 'REBELS.csv', format='csv')
# rebels_ra = rebels['RA']
# rebels_dec = rebels['Dec']

# # Plot the REBELS galaxies
# ax.scatter(rebels_ra, rebels_dec, marker='x', color='black', s=10, label='REBELS')


# Set axis limits and labels
ax.set_xlim(150.9, 149.1)
ax.set_ylim(1.5, 3.0)

# Squish ra axis by np.cos(dec) to account for declination


ax.set_xlabel('RA (deg)')
ax.set_ylabel('DEC (deg)')
ax.legend(loc='upper right')

plt.tight_layout()
#plt.savefig(Path.cwd().parent.parent / 'plots' / 'mosaic' / 'mosaic_gray_fullCWEB.png')
plt.show()



