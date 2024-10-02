"""
Select stars in COSMOS-Web following methodology of Zhuang and Shen.

Created: Wednesday 2nd October 2024.
"""

from astropy.table import Table
from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np
from astropy.coordinates import SkyCoord
from matplotlib.path import Path as mplPath

plt.rcParams['axes.linewidth'] = 2.5
plt.rcParams.update({'font.size': 15})
plt.rcParams['figure.dpi'] = 100

# Define the directories
cat_dir = Path.cwd().parents[1] / 'data' / 'depths' / 'COSMOS' / 'catalogues'
stars_dir = Path.cwd().parents[1] / 'data' / 'psf' / 'COSMOS' / 'catalogues'
mosaic_dir = Path.cwd().parents[1] / 'data' / 'mosaic'

# Dictionary of FWHM cuts for each band
fwhm_cuts = {'f115w': (1.8, 2.7), 'f150w': [2, 2.7], 'f277w': [3.7, 4.7], 'f444w': [5, 6]}

# Define the filter names
filter_names = ['f115w', 'f150w', 'f277w', 'f444w']
tiles = ['0A', '0B', '1A', '1B', '2A', '2B', '3A', '3B', '4A', '4B', '5A', '5B', '6A', '6B', '7A', '7B']

# Loop through the filters
for filter_name in filter_names:

    t = Table.read(cat_dir / f'd{filter_name}_locus_stars.fits')
    print(len(t))

    signal = t['FLUX_APER'][:, 0]
    noise = t['FLUXERR_APER'][:, 0]

    # Select stars
    stars = (t['ELONGATION'] < 1.5) & (t['CLASS_STAR'] < 0.99) & (t['CLASS_STAR'] > 0.8) & (t['FWHM_IMAGE'] > fwhm_cuts[filter_name][0]) & (t['FWHM_IMAGE'] < fwhm_cuts[filter_name][1])

    t_stars = t[stars]

    print(len(t_stars))

    t_stars.rename_column('NUMBER', 'ID')
    t_stars.remove_columns(['ALPHA_J2000', 'DELTA_J2000'])

    t_stars.rename_column('RA', 'ALPHA_J2000')
    t_stars.rename_column('DEC', 'DELTA_J2000')

    for tile in tiles:

        print(tile)
        tile_footprint = np.load(mosaic_dir / f'CWEB_footprint_{tile}.npy')        

        # Select stars in the tile. tile-footprint is four corners of the tile in ra,dec
        ra = t_stars['ALPHA_J2000']
        dec = t_stars['DELTA_J2000']

        catalog_coords = SkyCoord(ra, dec, unit='deg')
        polygon = mplPath(tile_footprint)
        points = np.vstack((ra, dec)).T
        mask = polygon.contains_points(points)
        objects_in_polygon = t_stars[mask]

        print(objects_in_polygon)

        # plt.scatter(objects_in_polygon['RA'], objects_in_polygon['DEC'], s=5)
        # plt.axis('equal')
        # plt.xlabel('RA')
        # plt.ylabel('Dec')
        # plt.title(filter_name)
        # plt.show()
        # plt.close()

        objects_in_polygon.write(stars_dir / f'{filter_name}_{tile}_stars.fits', overwrite=True)

        # Also save as ascii file. First remove the columns that are not needed
        objects_in_polygon.remove_columns(['FLUX_APER', 'FLUXERR_APER', 'MAG_APER'])
        objects_in_polygon.write(stars_dir / f'{filter_name}_{tile}_stars.ascii', format='ascii.commented_header', overwrite=True)
