#!/usr/bin/env python3

"""
crossmatch_vista_stars.py

Created: Wednesday 17th April 2024.
"""


from astropy.io import fits, ascii
from astropy.table import Table
from pathlib import Path
import astropy.units as u
from astropy.coordinates import SkyCoord

vista_dir = Path.home().parent.parent / 'vardy' / 'vardygroupshare' / 'data' / 'psf' / 'COSMOS' / 'catalogues'
euclid_dir = Path.cwd().parent.parent / 'data' / 'psf' / 'COSMOS' / 'catalogues'
save_dir = Path.cwd().parent.parent / 'data' / 'ref_catalogues' / 'stars'

filter_names = ['Y', 'J', 'H']

for filter_name in filter_names:

    # Open corresponding vista stars catalogue
    vista_stars_file = vista_dir / f'{filter_name}_DR6_stars.ascii'
    vista_stars = ascii.read(vista_stars_file)

    vista_ra, vista_dec = vista_stars['RA'], vista_stars['DEC']

    # And the euclid stars
    euclid_stars_file = euclid_dir / f'{filter_name}_stars.ascii'
    euclid_stars = ascii.read(euclid_stars_file)

    euclid_ra, euclid_dec = euclid_stars['RA'], euclid_stars['DEC']

    # Crossmatch between vista and euclid coords within 2 arcsec
    vista_coords = SkyCoord(ra=vista_ra, dec=vista_dec, unit=(u.deg, u.deg))
    euclid_coords = SkyCoord(ra=euclid_ra, dec=euclid_dec, unit=(u.deg, u.deg))

    idx_vista, idx_euclid, d2d, d3d = euclid_coords.search_around_sky(vista_coords, seplimit=2*u.arcsec)

    # Write matches to a new file with both coordinates
    crossmatch_file = save_dir / f'{filter_name}_vista_euclid_coords.ascii'

    with open(crossmatch_file, 'w') as f:
            
            # Add header for RA_vista, DEC_vista, RA_euclid, DEC_euclid
            f.write('# RA_vista DEC_vista RA_euclid DEC_euclid\n')
    
            for i, idx in enumerate(idx_vista):
                f.write(f'{vista_ra[idx]} {vista_dec[idx]} {euclid_ra[idx_euclid[i]]} {euclid_dec[idx_euclid[i]]}\n')

    print(f'Filter {filter_name} done.')




        