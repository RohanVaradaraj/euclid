"""
Select stars in COSMOS-Web following methodology of Zhuang and Shen.

Created: Wednesday 2nd October 2024.
"""

from astropy.table import Table
from pathlib import Path
import matplotlib.pyplot as plt

plt.rcParams['axes.linewidth'] = 2.5
plt.rcParams.update({'font.size': 15})
plt.rcParams['figure.dpi'] = 100

# Define the directories
cat_dir = Path.cwd().parents[1] / 'data' / 'depths' / 'COSMOS' / 'catalogues'

stars_dir = Path.cwd().parents[1] / 'data' / 'psf'

# Dictionary of FWHM cuts for each band
fwhm_cuts = {'f115w': (1.8, 2.7), 'f150w': [2, 2.7], 'f277w': [3.7, 4.7], 'f444w': [5, 6]}

# Define the filter names
filter_names = ['f115w', 'f150w', 'f277w', 'f444w']

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

    t_stars.write(cat_dir / f'd{filter_name}_stars.fits', overwrite=True)

    plt.scatter(t_stars['RA'], t_stars['DEC'], s=5)
    plt.axis('equal')
    plt.xlabel('RA')
    plt.ylabel('Dec')
    plt.title(filter_name)
    plt.show()
    plt.close()
