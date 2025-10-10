"""
Get all the high-z galaxies at z>4.5 which overlap with COSMOS-Web/COSMOS-3D.

Created: Tuesday 30th September 2025.
"""

from pathlib import Path
import numpy as np
from astropy.table import Table
from shapely.geometry import Point, Polygon

cat_dir = Path.cwd().parents[1] / 'data' / 'ref_catalogues'
radio_xmatch_catname = 'COSMOS-5p2arcsec-CrossMatch_Catalogue_withSpecz-WithMasked_Final_v1.0.fits'

def isCoordInCWEB(ra: np.ndarray, dec: np.ndarray) -> np.ndarray:
    """
    Check if a coordinate lies in the footprints of CWEB.

    If it doesn't, return zero for the coord, else return the mosaic string.

    Parameters
    ----------
    ra : np.array
        Right ascension of the coordinate in degrees.
    dec : np.array
        Declination of the coordinate in degrees.
    
    Returns
    -------
    np.ndarray
        Array with rows corresponding to the input coordinates and columns corresponding to the footprints.
        The value is '0' if it is not in the footprint.
        If it is in the footprint then it takes a numeric string value of the tile label Nathan has assigned (4A, 5A, 5B, 6A, etc.)

    """

    # Check if ra and dec are arrays
    ra = np.array(ra)
    dec = np.array(dec)

    # Reshape ra and dec to make sure they are at least one-dimensional arrays
    ra = np.atleast_1d(ra)
    dec = np.atleast_1d(dec)

    assert ra.shape == dec.shape, "ra and dec arrays must have the same shape"

    points = [Point(x, y) for x, y in zip(ra, dec)]

    mosaics = ['0A', '0B', '1A', '1B', '2A', '2B', '3A', '3B', '4A', '4B', '5A', '5B', '6A', '6B', '7A', '7B']

    polygons = []

    for mosaic in mosaics:
        footprint = np.load(Path.cwd().parent.parent / 'data' / 'mosaic' / f'CWEB_footprint_{mosaic}.npy')

        # Create shapely polygons
        cweb_polygon = Polygon(footprint)
        polygons.append(cweb_polygon)

    points_in_cweb = []
    for point in points:

        # 0 if not in footprint, else the tile label
        in_cweb = '0'

        #! Check if the point lies in any of the COSMOS-Web polygons
        #! If so, find all of them
        for i, polygon in enumerate(polygons):
            if polygon.contains(point):
                in_cweb = mosaics[i]
                break

        # Append results of search to lists
        points_in_cweb.append(in_cweb)

    # Convert lists to numpy arrays
    points_in_cweb = np.array(points_in_cweb)

    return points_in_cweb

t = Table.read(cat_dir / radio_xmatch_catname)
print('Size of original table:', len(t))
print(t.colnames)

# Sort in descending z order
t.sort('z_best_final', reverse=True)

# First restrict table to z>4.5
z_mask = t['z_best_final'] > 4.

t = t[z_mask]
print('Size of z>4.0 table:', len(t))


# Now restrict to where isCoordInCWEB is not zero
cweb_mask = isCoordInCWEB(t['RA_host'], t['Dec_host']) != '0'
t = t[cweb_mask]
print('Size of z>4.0 & in CWEB table:', len(t))

for i, row in enumerate(t):
    RA = t['RA_host'][i]
    Dec = t['Dec_host'][i]
    z = t['z_best_final'][i]
    id = t['Source_Name'][i]
    print(f'{RA:.6f}, {Dec:.6f}, {z:.3f}, {id}')