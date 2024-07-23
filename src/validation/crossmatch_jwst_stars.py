#!/usr/bin/env python3

"""
crossmatch_jwst_stars.py

Created: Wednesday 17th April 2024.
"""


from astropy.io import fits, ascii
from astropy.table import Table
from pathlib import Path
import astropy.units as u
from astropy.coordinates import SkyCoord

jwst_dir = Path.home().parent.parent / 'vardy' / 'vardygroupshare' / 'data' / 'psf' / 'COSMOS' / 'catalogues'
euclid_dir = Path.cwd().parent.parent / 'data' / 'psf' / 'COSMOS' / 'catalogues'
save_dir = Path.cwd().parent.parent / 'data' / 'ref_catalogues' / 'stars'

filter_names = ['Y', 'J', 'H']
jwst_filters = ['f115w', 'f150w', 'f200w']

for i, filter_name in enumerate(filter_names):

    # Open corresponding jwst stars catalogue
    jwst_stars_file = jwst_dir / f'{jwst_filters[i]}_stars.fits'
    jwst_stars = Table.read(jwst_stars_file)
    print(jwst_stars)

    jwst_ra, jwst_dec = jwst_stars['ALPHA_J2000'], jwst_stars['DELTA_J2000']

    # And the euclid stars
    euclid_stars_file = euclid_dir / f'{filter_name}_stars.ascii'
    euclid_stars = ascii.read(euclid_stars_file)

    euclid_ra, euclid_dec = euclid_stars['RA'], euclid_stars['DEC']

    # Crossmatch between jwst and euclid coords within 2 arcsec
    jwst_coords = SkyCoord(ra=jwst_ra, dec=jwst_dec, unit=(u.deg, u.deg))
    euclid_coords = SkyCoord(ra=euclid_ra, dec=euclid_dec, unit=(u.deg, u.deg))

    idx_jwst, idx_euclid, d2d, d3d = euclid_coords.search_around_sky(jwst_coords, seplimit=2*u.arcsec)

    # Write matches to a new file with both coordinates
    crossmatch_file = save_dir / f'{filter_name}_jwst_euclid_coords.ascii'

    with open(crossmatch_file, 'w') as f:
            
            # Add header for RA_jwst, DEC_jwst, RA_euclid, DEC_euclid
            f.write('# RA_jwst DEC_jwst RA_euclid DEC_euclid\n')
    
            for i, idx in enumerate(idx_jwst):
                f.write(f'{jwst_ra[idx]} {jwst_dec[idx]} {euclid_ra[idx_euclid[i]]} {euclid_dec[idx_euclid[i]]}\n')

    print(f'Filter {filter_name} done.')




        