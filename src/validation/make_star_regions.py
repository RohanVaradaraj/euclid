"""
make_star_regions.py

Generate a list of regions to open in ds9 around our problematic stars.

Can then blink them to check for proper motion.
"""

from astropy.io import ascii
from astropy.table import Table
from pathlib import Path

stars_dir = Path.cwd().parent.parent / 'data' / 'ref_catalogues' / 'stars'

filter_names = ['Y']

for filter_name in filter_names:

    t = Table.read(stars_dir / f'{filter_name}_outside_fwhm_vista_gaia_coords.ascii', format='ascii.commented_header')

    with open(stars_dir / 'regions' / f'{filter_name}_pmFiltered_fwhm_pixscale_vista_gaia_coords.reg', 'w') as f:

        f.write('global color=green dashlist=8 3 width=1 font="helvetica 10 normal roman" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1\n')
        f.write('fk5\n')

        for i, row in enumerate(t):

            #f.write(f'circle({row["RA_euclid"]},{row["DEC_euclid"]},0.003)\n')
            f.write(f'circle({row["RA_Gaia"]},{row["DEC_Gaia"]},0.003) # pm = {row["pm"]}\n')