"""
source_counts.py

Plot the number of sources per square degree as a function of MAG_AUTO in Euclid, comparing to VISTA.

Created: Tuesday 16th April 2024.
"""

from astropy.table import Table
from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np

plt.rcParams['axes.linewidth'] = 2.5
plt.rcParams.update({'font.size': 15})
plt.rcParams['figure.dpi'] = 100

euclid_dir = Path.cwd().parent.parent / 'data' / 'depths' / 'COSMOS' / 'catalogues'
ground_dir = Path.home().parent.parent / 'vardy' / 'vardygroupshare' / 'data' / 'depths' / 'COSMOS' / 'catalogues'

filter_names = ['Y', 'J', 'H']

for filter_name in filter_names:

    # Open the euclid catalogue
    t_euclid = Table.read(euclid_dir / f'd{filter_name}.fits')

    # Open the ground catalogue
    t_ground = Table.read(ground_dir / f'd{filter_name}_DR6.fits')

    # Get the MAG_AUTO values
    mags_e = t_euclid['MAG_AUTO']
    mags_g = t_ground['MAG_AUTO']

    # Remove missing data
    mags_e = mags_e[(mags_e != -99.) | (mags_e != 99.)]
    mags_g = mags_g[(mags_g != -99.) | (mags_g != 99.)]

    # Bin into 0.5 mag bins
    bins = np.arange(12, 32, 0.5)

    hist_e, bins = np.histogram(mags_e, bins=bins)
    hist_g, bins = np.histogram(mags_g, bins=bins)

    # Plot
    plt.figure(figsize=(8, 8))
    plt.scatter(bins[:-1], hist_e, label=f'Euclid {filter_name}', color='black', marker='v', s=100)
    plt.plot(bins[:-1], hist_e, color='black', alpha=0.6)
    plt.scatter(bins[:-1], hist_g, label=f'VISTA {filter_name}', color='blue', marker='o', s=70, alpha=0.7)
    plt.plot(bins[:-1], hist_g, color='blue', alpha=0.6)

    # Find the maximum value and draw vertical lines, with labels of the mag value here
    max_e = np.max(hist_e)
    max_g = np.max(hist_g)
    plt.axvline(x=bins[np.argmax(hist_e)], color='black', linestyle='--', label=f'VISTA max = {bins[np.argmax(hist_e)]:.1f}')
    plt.axvline(x=bins[np.argmax(hist_g)], color='blue', linestyle='--', label=f'VISTA max = {bins[np.argmax(hist_g)]:.1f}')

    plt.legend()
    plt.xlabel('MAG_AUTO')
    plt.ylabel('N / 0.5 mag')
    plt.yscale('log')
    plt.tight_layout()

    plt.savefig(Path.cwd().parent.parent / 'plots' / 'depths' / f'{filter_name}_source_counts.png', bbox_inches='tight')
    plt.show()
    plt.close()

