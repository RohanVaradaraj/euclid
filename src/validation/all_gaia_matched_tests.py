"""
all_gaia_matched_tests.py

Seeing if downloading all COSMOS Gaia objects from the science archive and then crossmatching gives different results.

Created: Wednesday 8th May 2024.

"""

from astropy.table import Table
from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np

plt.rcParams['axes.linewidth'] = 2.5
plt.rcParams.update({'font.size': 15})
plt.rcParams['figure.dpi'] = 100

stars_dir = Path.cwd().parent.parent / 'data' / 'ref_catalogues' / 'stars'

t = Table.read(stars_dir / 'all_gaia_euclid_matched.fits')

pm = np.sqrt(t['pmra']**2 + t['pmdec']**2)

# Remove objects with proper motion greater than 1 mas/yr
mask = pm < 1
t = t[mask]

pm = np.sqrt(t['pmra']**2 + t['pmdec']**2)
print(np.max(pm))
exit()


delta_ra = t['RA_1'] - t['ra_2']
delta_dec = t['DEC_1'] - t['dec_2']

delta_ra *= 3600 # convert to arcseconds
delta_dec *= 3600 # convert to arcseconds

# plt.hist(t['Separation'], bins=100, histtype='step', color='black', lw=2.5)
# plt.show()
# plt.close()

plt.figure(figsize=(8, 8))
plt.scatter(delta_ra, delta_dec, s=1, marker='x', color='black', alpha=0.6)
plt.xlabel(r'$\Delta \mathrm{RA \ (as)}$')
plt.ylabel(r'$\Delta \mathrm{Dec \ (as)}$')

# Plot a circle corresponding to the Euclid pixel scale
circle1 = plt.Circle((0, 0), 0.1/2, color='r', fill=False, linewidth=2.5)
plt.gca().add_artist(circle1)

# And plot a Circle corresponding to the filter FWHM peak
circle2 = plt.Circle((0, 0), 0.48/2, color='black', fill=False, linewidth=2.5)
plt.gca().add_artist(circle2)

# Get the objects outside the FWHM and pixel scale
outside_fhwm = np.where((np.sqrt(delta_ra**2 + delta_dec**2) > 0.48))[0]
outside_pix = np.where((np.sqrt(delta_ra**2 + delta_dec**2) > 0.1/2) & (np.sqrt(delta_ra**2 + delta_dec**2) < 0.48/2))[0]

# Make region files for these

with open(stars_dir / 'regions' / f'Y_all_gaia_euclid_outside_FWHM.reg', 'w') as f:
    f.write('global color=green dashlist=8 3 width=1 font="helvetica 10 normal roman" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1\n')
    f.write('fk5\n')
    for i, row in enumerate(t[outside_fhwm]):
        f.write(f'circle({row["RA_1"]},{row["DEC_1"]},0.003) # pm = {pm[outside_fhwm][i]}\n')

with open(stars_dir / 'regions' / f'Y_all_gaia_euclid_outside_pixScale.reg', 'w') as f:
    f.write('global color=green dashlist=8 3 width=1 font="helvetica 10 normal roman" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1\n')
    f.write('fk5\n')
    for i, row in enumerate(t[outside_pix]):
        f.write(f'circle({row["RA_1"]},{row["DEC_1"]},0.003) # pm = {pm[outside_pix][i]}\n')




# Dummy plots for the legend
plt.plot([-99, -98], [-99, -98], color='r', lw=2.5, label='Euclid pixel scale')
plt.plot([-99, -98], [-99, -98], color='black', lw=2.5, label=f'Y PSF FWHM')
plt.axis('equal')
plt.xlim(-0.6, 0.6)
plt.ylim(-0.6, 0.6)
plt.tight_layout()
plt.legend()
plt.show()




