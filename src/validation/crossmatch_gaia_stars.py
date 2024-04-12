#!/usr/bin/env python3

"""
crossmatch_gaia_stars.py

Created: Friday 12th April 2024.
"""


from astropy.io import fits, ascii
from astropy.table import Table
from pathlib import Path
from astroquery.gaia import Gaia
import astropy.units as u
from astropy.coordinates import SkyCoord

# Suppress astroquery info
import logging
logging.getLogger('astroquery').setLevel(logging.ERROR)


stars_dir = Path.cwd().parent.parent / 'data' / 'psf' / 'COSMOS' / 'catalogues'
save_dir = Path.cwd().parent.parent / 'data' / 'ref_catalogues'

filter_names = ['VIS', 'Y', 'J', 'H']

for filter_name in filter_names:

    print(f'Processing filter {filter_name}')

    # Open empty file to star the euclid coords and Gaia coords if there is a match
    euclid_gaia_coords = stars_dir / f'{filter_name}_brightEuclid_gaia_coords.ascii'
    with open(euclid_gaia_coords, 'w') as f:

        # Add header for RA_euclid, DEC_euclid, RA_Gaia, DEC_Gaia
        f.write('# RA_euclid DEC_euclid RA_Gaia DEC_Gaia\n')

        # Read the commented header stars.ascii files
        #stars_file = stars_dir / f'{filter_name}_stars.ascii'
        #stars = ascii.read(stars_file)

        #! Alternate file: very bright stars
        stars_file = save_dir / f'{filter_name}_bright_stars.fits'
        stars = Table.read(stars_file)
        # Change ALPHA_J2000 and DELTA_J2000 to RA and DEC
        stars.rename_column('ALPHA_J2000', 'RA')
        stars.rename_column('DELTA_J2000', 'DEC')

        # Keep rows with CLASS_STAR above a threshold
        stars = stars[stars['CLASS_STAR'] > 0.95]

        # Only keep the second and third columns
        stars = stars['RA', 'DEC']


        # Loop through each row
        for i, star in enumerate(stars):
            print(f'{i+1} of {len(stars)}')
            ra = star['RA']
            dec = star['DEC']

            # Query Gaia for the star
            euclid_coord = SkyCoord(ra=ra, dec=dec, unit=(u.degree, u.degree), frame='icrs')
            job = Gaia.cone_search_async(euclid_coord, radius=u.Quantity(2.0, u.arcsec))
            r = job.get_results()

            if len(r) > 0:

                print('Match found!')

                # Write to file
                f.write(f'{ra} {dec} {r["ra"][0]} {r["dec"][0]}\n')

    

        