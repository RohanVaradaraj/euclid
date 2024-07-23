"""
Testing the YJH-selected catalogues produced by SE++

Created: Tuesday 25th June 2024.
"""

import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table
from pathlib import Path
from scipy.ndimage import uniform_filter1d

plt.rcParams['axes.linewidth'] = 2.5
plt.rcParams.update({'font.size': 15})
plt.rcParams['figure.dpi'] = 100

def running_median(data, window_size):
    """Compute the running median with a given window size, extending to the ends of the data."""
    extended_data = np.pad(data, (window_size//2, window_size//2), mode='edge')
    median_filtered = np.median(
        np.lib.stride_tricks.sliding_window_view(extended_data, window_size),
        axis=1
    )
    return median_filtered

cat_dir = Path.cwd().parents[1] / 'data' / 'ref_catalogues' / 'SEplus_cats'

# catalogue of stars
cat_name = 'det_YJH_stars.fits'

t = Table.read(cat_dir / cat_name)

print(t.colnames)


mag_Y = -2.5*np.log10(t['flux_Y_2.0as']) -48.6
mag_J = -2.5*np.log10(t['flux_J_2.0as']) -48.6
mag_H = -2.5*np.log10(t['flux_H_2.0as']) -48.6

# Restrict to stellar locus
# -0.2 < Y-J < 0.4
# -0.2 < J-H < 0.3
t = t[(mag_Y-mag_J > -0.2) & (mag_Y-mag_J < 0.4) & (mag_J-mag_H > -0.2) & (mag_J-mag_H < 0.3)]

mag_Y = -2.5*np.log10(t['flux_Y_2.0as']) -48.6
mag_J = -2.5*np.log10(t['flux_J_2.0as']) -48.6
mag_H = -2.5*np.log10(t['flux_H_2.0as']) -48.6

# plt.scatter(mag_Y-mag_J, mag_J-mag_H, s=1, c='k')
# plt.xlabel('Y-J')
# plt.ylabel('J-H')
# plt.show()

mag_Ye = -2.5*np.log10(t['flux_Ye_1.2as']) -48.6
mag_Je = -2.5*np.log10(t['flux_Je_1.2as']) -48.6
mag_He = -2.5*np.log10(t['flux_He_1.2as']) -48.6


######################! MAG - MAG VS MAG ###############################
#! Compare Y bands
plt.figure(figsize=(10, 6))
plt.scatter(mag_Y, mag_Y-mag_Ye, s=1, c='k')
plt.xlabel(r'$Y$')
plt.ylabel(r'$Y-Y_{e}$')

# Sort data by mag_Y
sorted_indices = np.argsort(mag_Y)
mag_Y_sorted = mag_Y[sorted_indices]
mag_Ye_diff_sorted = (mag_Y - mag_Ye)[sorted_indices]

# Compute the running median for mag_Ye_diff_sorted
window_size = 50
mag_Ye_diff_median = running_median(mag_Ye_diff_sorted, window_size)
plt.plot(mag_Y_sorted, mag_Ye_diff_median[1:], color='red', label='Running Median', linewidth=2)
plt.plot([14, 26], [0, 0], color='green', label='Zero', linewidth=2)

plt.ylim(-1, 0.5)
plt.xlim(15, 24)


plt.show()
plt.close()

#! Compare J bands
plt.figure(figsize=(10, 6))
plt.scatter(mag_J, mag_J-mag_Je, s=1, c='k')
plt.xlabel(r'$J$')
plt.ylabel(r'$J-J_{e}$')

# Sort data by mag_Y
sorted_indices = np.argsort(mag_J)
mag_J_sorted = mag_J[sorted_indices]
mag_Je_diff_sorted = (mag_J - mag_Je)[sorted_indices]

# Compute the running median for mag_Ye_diff_sorted
window_size = 50
mag_Je_diff_median = running_median(mag_Je_diff_sorted, window_size)
plt.plot(mag_J_sorted, mag_Je_diff_median[1:], color='red', label='Running Median', linewidth=2)
plt.plot([14, 26], [0, 0], color='green', label='Zero', linewidth=2)

plt.ylim(-1, 0.5)
plt.xlim(15, 24)


plt.show()
plt.close()

#! Compare H bands
plt.figure(figsize=(10, 6))
plt.scatter(mag_H, mag_H-mag_He, s=1, c='k')
plt.xlabel(r'$H$')
plt.ylabel(r'$H-H_{e}$')

# Sort data by mag_H
sorted_indices = np.argsort(mag_H)
mag_H_sorted = mag_H[sorted_indices]
mag_He_diff_sorted = (mag_H - mag_He)[sorted_indices]

# Compute the running median for mag_He_diff_sorted
window_size = 50
mag_He_diff_median = running_median(mag_He_diff_sorted, window_size)
plt.plot(mag_H_sorted, mag_He_diff_median[1:], color='red', label='Running Median', linewidth=2)
plt.plot([14, 26], [0, 0], color='green', label='Zero', linewidth=2)

plt.ylim(-1, 0.5)
plt.xlim(15, 24)

plt.show()
plt.close()



#########################! MAG VS MAG ###########################

# Compare Y and Ye
plt.figure(figsize=(10, 10))
plt.scatter(mag_Y, mag_Ye, s=1, c='k')
plt.plot([16, 24], [16, 24], color='red', label='1:1', linewidth=2)

plt.xlabel(r'$Y$')
plt.ylabel(r'$Y_{e}$')

plt.show()


