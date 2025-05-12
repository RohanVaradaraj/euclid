"""
quick script for plotting the difference in LF between UltraVISTA-only and UltraVISTA+Euclid.

Created: March 14th 2025.
"""
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

plt.rcParams.update({'font.size': 15})
plt.rcParams['axes.linewidth'] = 4
plt.rcParams['figure.dpi'] = 100

uvista_euclid = [2.12861662e-06, 6.03210043e-06, 9.16337525e-06, 1.59486164e-05, 3.23233992e-05, 4.24464998e-05, 6.05560308e-05] #, 6.61033853e-05, 5.90742485e-05, 2.17809970e-05]
uvista_only = [3.32642337e-06, 7.34422140e-06, 6.27708092e-06, 1.57390423e-05, 2.74080836e-05, 4.12950686e-05, 3.83475493e-05] #, 5.04532349e-05, 3.29772144e-05, 4.49456387e-06]

Muv_bins = [-22.4, -22., -21.8,  -21.6, -21.4, -21.2, -21.0, -20.8, -20.6, -20.4, -20.2]
bin_widths = np.abs(np.diff(Muv_bins))
bin_centres = 0.5 * (np.array(Muv_bins[:-1]) + np.array(Muv_bins[1:]))[:-3]
print(bin_centres)

# Plot difference
fig = plt.figure(figsize=(10, 4))

plt.scatter(bin_centres, np.array(uvista_euclid) - np.array(uvista_only), s=100, color='black', label='UVISTA + Euclid - UVISTA only')
plt.hlines(0, -22.4, -20.2, color='black', linestyle='--', lw=2)
plt.xlim([-23.2, -20.7])
#plt.yscale('log')
plot_dir = Path.cwd().parents[1] / 'plots' / 'LF'
plt.tight_layout()
plt.savefig(plot_dir / 'diff_LF.pdf', bbox_inches='tight')


