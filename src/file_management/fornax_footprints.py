"""
Check the coverage of the Euclid tiles against MIGHTEE and the VIDEO fields.

Created: Tuesday 6th May 2025 (first new code as a postdoc!)
"""

from astropy.io import fits
import matplotlib.pyplot as plt
from astropy.wcs import WCS
from pathlib import Path
import glob
import os
from matplotlib.patches import Circle
from astropy.coordinates import SkyCoord
import astropy.units as u
from shapely.geometry import Polygon
import pickle
import numpy as np

plt.rcParams.update({'font.size': 15})
plt.rcParams['axes.linewidth'] = 4
plt.rcParams['figure.dpi'] = 100

#! Euclid Deep Field Fornax directory
fornax_dir = Path.cwd().parents[3] / 'data' / 'euclid' / 'euclid_deep_field_fornax' / 'Y'
fornax_dir = Path.home() / 'euclid' / 'CDFS'

#! VIDEO directory
vista_dir = Path.cwd().parents[3] / 'data'
tiles = ['CDFS1', 'CDFS2', 'CDFS3']

plt.figure(figsize=(10,10))

#! MIGHTEE pointings
# === Load the data from your text file ===
coordinates = []
with open('MIGHTEE_pointings.txt', 'r') as f:
    for line in f:
        if line.startswith("#") or not line.strip():
            continue  # Skip comments or empty lines
        parts = line.strip().split()
        if len(parts) == 3:
            name, ra, dec = parts
            coordinates.append((name, float(ra), float(dec)))

ax = plt.gca()  # Get current axes

# === Plot circles and points ===
radius = 0.7136  # degrees

#! MIGHTEE pointings
# for name, ra, dec in coordinates:
#     circle = Circle((ra, dec), radius, edgecolor='orange', facecolor='none', lw=2, alpha=0.8)
#     ax.add_patch(circle)
    #plt.plot(ra, dec, 'k.', markersize=2)  # center point
    # Optional: label each point
    # plt.text(ra, dec, name, fontsize=6)


pattern = os.path.join(fornax_dir, f'*BGSUB-MOSAIC-VIS*')

# List all matching files
euclid_files = glob.glob(pattern)

# Get the tile names
#tile_names = [file_name.split('-')[3].split('VIS_TILE')[1] for file_name in euclid_files]
tile_names = [file_name.split('-')[2].split('VIS_TILE')[1] for file_name in euclid_files]

print(tile_names)
footprints = []

for i, file in enumerate(euclid_files[1:2]):
    print(f"file {i} of {len(euclid_files)}")
    print(file.split('/')[-1])

    with fits.open(file, ignore_missing_simple=True) as hdu:
        wcs = WCS(hdu[0].header)

    # Get the footprint
    footprint = wcs.calc_footprint()

    # Append footprint
    footprints.append(footprint)


#! LSST footprints
lsst_dir = Path.home().parents[1] / 'extraspace' / 'varadaraj' / 'lsst' / 'y'
lsst_files = glob.glob(str(lsst_dir / '*SCI.fits'))
lsst_tilenames = [file_name.split('LSST_dp1_')[-1].split('_y_SCI.fits')[0] for file_name in lsst_files]
lsst_footprints = []

for i, file in enumerate(lsst_files):
    print(f"LSST file {i} of {len(lsst_files)}")
    with fits.open(file) as hdu:
        print(hdu.info())
        wcs = WCS(hdu[0].header)
        print(wcs)

    footprint = wcs.calc_footprint()
    lsst_footprints.append(footprint)



def plot_footprints(footprints, close_polygons=True, color='blue', lw=2, alpha=0.8, tilename=False, zorder=1):
    for i, fp in enumerate(footprints):
        fp = fp.tolist()
        if close_polygons:
            fp.append(fp[0])  # Close the polygon
        ra, dec = zip(*fp)
        plt.plot(ra, dec, '-', color=color, alpha=alpha, lw=lw, zorder=zorder)
        if tilename:
            ra_center = sum(ra) / len(ra)
            dec_center = sum(dec) / len(dec)
            plt.text(ra_center, dec_center, tile_names[i], fontsize=6, ha='center', va='center', color=color)
    
    plt.xlabel("RA (deg)")
    plt.ylabel("Dec (deg)")
    plt.title("Euclid Deep Field Fornax")

#! Euclid pointings
# plot_footprints(footprints, color='green', lw=2, alpha=0.6, tilename=False)

# ! VIDEO pointings
vista_footprints = []
for i, tile in enumerate(tiles):
    tile_dir = vista_dir / tile
    with fits.open(tile_dir / 'HSC-G.fits') as hdu:
        wcs = WCS(hdu[0].header)

    footprint = wcs.calc_footprint()
    vista_footprints.append(footprint)

plot_footprints(vista_footprints, color='black', lw=4, alpha=1, zorder=10)

#! LSST pointings
plot_footprints(lsst_footprints, color='deepskyblue', lw=2, alpha=0.8, tilename=False)


#! FIND WHICH EUCLID TILES LIE IN VIDEO
# Convert footprints to Shapely Polygons
video_polygons = [Polygon(fp) for fp in vista_footprints]
euclid_polygons = [Polygon(fp) for fp in footprints]

#Build mapping: VIDEO tile -> list of Euclid tile names
euclid_within_video = {}

for j, video_poly in enumerate(video_polygons):
    video_tile_name = tiles[j]
    euclid_within_video[video_tile_name] = []
    
    for i, euclid_poly in enumerate(euclid_polygons):
        if euclid_poly.intersects(video_poly):
            euclid_within_video[video_tile_name].append(tile_names[i])

for video_tile, euclid_tiles in euclid_within_video.items():
    print(f"\nVIDEO tile: {video_tile}")
    if euclid_tiles:
        print("  Euclid tiles overlapping:")
        for tile in euclid_tiles:
            print(f"    - {tile}")
    else:
        print("  No overlapping Euclid tiles.")

video_colors = {
    'CDFS1': 'red',
    'CDFS2': 'orange',
    'CDFS3': 'yellow'
}

#! FIND WHICH LSST TILES LIE IN VIDEO
lsst_polygons = [Polygon(fp) for fp in lsst_footprints]
lsst_within_video = {}


for j, video_poly in enumerate(video_polygons):
    video_tile_name = tiles[j]
    lsst_within_video[video_tile_name] = []
    
    for i, lsst_poly in enumerate(lsst_polygons):
        if lsst_poly.intersects(video_poly):
            lsst_within_video[video_tile_name].append(lsst_tilenames[i])

with open('lsst_DR1_within_video.pkl', 'wb') as f:
    pickle.dump(lsst_within_video, f)

print(lsst_within_video)
print(euclid_within_video)

#! Plotting the Eucid pointings that lie within the VIDEO tiles, one VIDEO tile at a time
# video_tile = tiles[2]  # or just 'CDFS1'
# color = video_colors[video_tile]

# for tile in euclid_within_video[video_tile]:
#     idx = tile_names.index(tile)
#     fp = footprints[idx].tolist()
#     fp.append(fp[0])  # Close polygon
#     ra, dec = zip(*fp)
#     plt.plot(ra, dec, '-', color=color, lw=2, alpha=0.9)#, label=video_tile)

#! Plotting all the Euclid tiles inside VIDEO, coloured by VIDEO tile
# for video_tile, euclid_tile_list in euclid_within_video.items():
#     color = video_colors[video_tile]
#     for tile in euclid_tile_list:
#         idx = tile_names.index(tile)  # Get index of this tile in original footprints list
#         fp = footprints[idx].tolist()
#         fp.append(fp[0])  # Close polygon
#         ra, dec = zip(*fp)
#         plt.plot(ra, dec, '-', color=color, lw=2.5, alpha=0.9)

#! Plotting all the LSST tiles inside VIDEO, coloured by VIDEO tile
for video_tile, lsst_tile_indices in lsst_within_video.items():
    color = video_colors[video_tile]
    for idx in lsst_tile_indices:
        fp = lsst_footprints[idx].tolist()
        fp.append(fp[0])  # Close polygon
        ra, dec = zip(*fp)
        plt.plot(ra, dec, '-', color=color, lw=1.5, alpha=0.8)

#! Plot a circle centred at 53,-28 with radius of my choosing
# radius = 2 # degrees
# center_ra = 53
# center_dec = -28
# circle = Circle((center_ra, center_dec), radius, edgecolor='purple', facecolor='none', lw=2, alpha=0.8)
# ax.add_patch(circle)

# Dummy labels
plt.plot([], [], color='black', lw=3, alpha=1, label='VIDEO')
plt.plot([], [], color='green', lw=2, alpha=0.6, label='Euclid Q1')
# plt.plot([], [], color='orange', lw=2, alpha=0.8, label='MIGHTEE')

#plt.grid(True)
plt.gca().invert_xaxis()  # Optional, matches sky view
plt.tick_params(which='major', length=10, width=3)
plt.tick_params(axis='both', which='minor', length=5, width=2)
plt.legend()
plt.tight_layout()
# equal aspect ratio
plt.axis('equal')
plt.show()

# Save the dictionary of euclid_within_video to a file
with open('euclid_DR1_within_video.pkl', 'wb') as f:
    pickle.dump(euclid_within_video, f)




