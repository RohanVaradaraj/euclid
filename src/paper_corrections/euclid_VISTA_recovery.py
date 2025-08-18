"""
Karina Caputi asked me to check what fraction of UltraVISTA galaxies are recoverable with Euclid data alone.

My beautiful code allows me to run the just euclid stuff very quickly.

Now checking the .out file to look at the recovery.

Created: Friday 15th August 2025.
"""

from astropy.io import ascii
from astropy.table import Table
from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np

# Plotting configuration
plt.rcParams.update({
    'axes.linewidth': 2.5,
    'font.size': 15,
    'figure.dpi': 100
})

plt.rcParams.update({
    # Ticks on all sides, pointing inwards
    'xtick.top': True, 'xtick.bottom': True,
    # 'ytick.left': True, 'ytick.right': True,
    'xtick.direction': 'in', 'ytick.direction': 'in',

    # Major tick size and width
    'xtick.major.size': 6.5, 'ytick.major.size': 6.5,
    'xtick.major.width': 2, 'ytick.major.width': 2,

    # Minor tick size and width
    'xtick.minor.size': 3, 'ytick.minor.size': 3,
    'xtick.minor.width': 1.5, 'ytick.minor.width': 1.5,
})

# Actual sample
cand_dir = Path.cwd().parents[1] / 'data' / 'catalogues' / 'candidates'
tUE = Table.read(cand_dir / 'Euclid_UltraVISTA_z7_sample.fits')

file_name = '/mnt/zfsusers/varadaraj/lephare/lephare_dev/test/det_Y_J_just_euclid.out'

t = ascii.read(file_name, format='no_header')
print(t)

print(tUE)

z_orig = tUE['Zphot']
z_recovered = t['col2']

plt.scatter(z_orig, z_recovered, s=50)
plt.plot([1, 8.8], [1, 8.8], color='red', linestyle='--', linewidth=2)

# Also plot |zrec−zorig|/(1 + zorig) >0.15 lines
z = np.arange(1, 8.8, 0.01)
plt.plot(z, z + 0.15 * (1 + z), color='blue', linestyle='--', linewidth=2)
plt.plot(z, z - 0.15 * (1 + z), color='blue', linestyle='--', linewidth=2)

# Calculate the fraction of source outside the |zrec−zorig|/(1 + zorig) >0.15 lines
mask = np.abs(z_recovered - z_orig) / (1 + z_orig) > 0.15
fraction_outside = np.sum(mask) / len(z_orig)
print(f"Fraction of sources outside |zrec−zorig|/(1 + zorig) > 0.15: {fraction_outside:.2f}")
print('Number of sources outside: ', np.sum(mask))

plt.ylim(1, 8.8)
plt.xlim(1, 8.8)
plt.xlabel(r'$z_{\rm{U+E}}$')
plt.ylabel(r'$z_{\rm{E-only}}$')
plt.show()