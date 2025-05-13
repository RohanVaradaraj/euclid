"""
Reads in catalogues and determines which VIDEO tile they lie in.

Created: Tuesday 13th May 2025.
"""

from astropy.io import fits
from astropy.table import Table
from pathlib import Path
from astropy.wcs import WCS
from shapely.geometry import Polygon, Point
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.colors as mcolors

plt.rcParams['axes.linewidth'] = 2.5
plt.rcParams.update({'font.size': 20})
plt.rcParams['figure.dpi'] = 100

#! Field names
field_name = 'XMM'
video_tiles = [field_name+x for x in ['1', '2', '3']]

#! Read in catalogue to label
cat_dir = Path.cwd().parents[1] / 'data' / 'catalogues' /'candidates'
cat_name = 'XMM_5sig_HSC_Z_nonDet_HSC_G_nonDet_HSC_R_candidates_2025_05_13.fits'
t = Table.read(cat_dir / cat_name)

#! Directory containing VIDEO-matched images
video_dir = Path.cwd().parents[3] / 'data'

footprints = []

#! Read the VIDEO footprints
for i, tile in enumerate(video_tiles):

    if field_name == 'XMM':
        with fits.open(video_dir / tile / 'HSC-G_DR3.fits') as hdu:
            wcs = WCS(hdu[0].header)
            tile_footprint = wcs.calc_footprint()
            polygon = Polygon(tile_footprint)
            footprints.append(polygon)

    if field_name == 'CDFS':
        with fits.open(video_dir / tile / 'HSC-G.fits') as hdu:
            wcs = WCS(hdu[0].header)
            tile_footprint = wcs.calc_footprint()
            polygon = Polygon(tile_footprint)
            footprints.append(polygon)


#! Get table coords
ra = t['RA']
dec = t['DEC']

# Optional: column to store the matched tile name
tile_matches = []

# Loop through each object
for ra, dec in zip(ra, dec):
    point = Point(ra, dec)

    matched_tile = None
    for idx, polygon in enumerate(footprints):
        if polygon.contains(point):
            matched_tile = video_tiles[idx]  # Assign the FIRST tile match
            break  # Stop checking after the first match

    tile_matches.append(matched_tile)

# Add tile info to the table
t['video_tile'] = tile_matches

print(t)

# Save the updated table
t.write(cat_dir / cat_name, overwrite=True)

# Unique tile names (excluding None)
unique_tiles = sorted(set(tile for tile in tile_matches if tile is not None))

# Assign colors to tiles
cmap = cm.get_cmap('tab10', len(unique_tiles)) 
tile_to_color = {tile: cmap(i) for i, tile in enumerate(unique_tiles)}

fig, ax = plt.subplots(figsize=(10, 10))

# Plot footprints
for i, polygon in enumerate(footprints):
    x, y = polygon.exterior.xy
    ax.fill(x, y, alpha=0.3, fc='lightgray', ec='black', lw=3)

# Plot sources color-coded by matched tile
for i, (ra_val, dec_val, tile) in enumerate(zip(t['RA'], t['DEC'], tile_matches)):
    if tile is not None:
        ax.plot(ra_val, dec_val, 'o', markersize=5, color=tile_to_color[tile])
    else:
        ax.plot(ra_val, dec_val, 'x', markersize=5, color='red')  # Not matched to any tile

# Legend
for tile, color in tile_to_color.items():
    ax.plot([], [], 'o', color=color, label=tile)
if None in tile_matches:
    ax.plot([], [], 'x', color='red', label='No tile match')

ax.set_xlabel('RA [deg]')
ax.set_ylabel('Dec [deg]')
ax.legend(loc='best')
ax.set_title(f'{field_name} candidates and VIDEO tile footprints')

# Flip x axis
ax.set_xlim(ax.get_xlim()[::-1])

plt.tight_layout()
plt.show()


