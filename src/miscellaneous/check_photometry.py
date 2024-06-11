"""
check photometry of specific objects
"""

import astropy.units as u
from astropy.table import Table
from astropy.coordinates import SkyCoord
import matplotlib.pyplot as plt
import numpy as np

t = Table.read('XMM3_3085.fits')

print(t)

# Filter names of interest, and columns are of form flux_{filter}
order = ['HSC-G_DR3', 'HSC-R_DR3', 'HSC-I_DR3', 'HSC-Z_DR3', 'HSC-Y_DR3', 'Y', 'J', 'H', 'Ks', 'ch1servs', 'ch2servs']

for col in t.colnames:
    if 'flux_' in col and 'YJ' not in col:

        name = col.replace('flux_', '')
        if name not in order:
            continue

        else:
            flux = t[col]
            err = t[col.replace('flux_', 'err_')]

            # Check if flux is less than 2sigma
            if flux / err < 2:
                print('Non-detection in band', col)

                # Convert to magnitude 2sigma upper limit
                mag = -2.5 * np.log10(2 * err[0]) - 48.6
                print('2sigma upper limit in band', col, 'is', mag)
            else:
                sn = flux[0] / err[0]
                print('Detection in band', col, 'with S/N', sn)

                # Convert to mag
                mag = -2.5 * np.log10(flux[0]) - 48.6
                print('Magnitude in band', col, 'is', mag)

                # Get mag errors
                mag_err = 1.0857 * err[0] / flux[0]
                print('Magnitude error in band', col, 'is', mag_err)

        



