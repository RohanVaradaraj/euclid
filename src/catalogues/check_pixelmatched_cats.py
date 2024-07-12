"""
Check the catalogues made from regridding and PSF-homogenising the Euclid and CWEB data!

Created: Friday 12th June 2024.
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

cat_dir = Path.cwd().parents[3] / 'data' / 'catalogues' / 'finalCOSMOS' / 'other'

# catalogue of stars
cat_name = 'COSMOS_detYJH_masked_1.8as_Euclid_CWEB_2024_07_12.fits'

t = Table.read(cat_dir / cat_name)

print(t.colnames)


mag_Y = -2.5*np.log10(t['flux_Y']) -48.6
mag_J = -2.5*np.log10(t['flux_J']) -48.6
mag_H = -2.5*np.log10(t['flux_H']) -48.6

# Restrict to stellar locus
# -0.2 < Y-J < 0.4
# -0.2 < J-H < 0.3
#t = t[(mag_Y-mag_J > -0.2) & (mag_Y-mag_J < 0.4) & (mag_J-mag_H > -0.2) & (mag_J-mag_H < 0.3)]


# plt.scatter(mag_Y-mag_J, mag_J-mag_H, s=1, c='k')
# plt.xlabel('Y-J')
# plt.ylabel('J-H')
# plt.show()

mag_Ye = -2.5*np.log10(t['flux_Ye']) -48.6
mag_Je = -2.5*np.log10(t['flux_Je']) -48.6
mag_He = -2.5*np.log10(t['flux_He']) -48.6

mag_f115w = -2.5*np.log10(t['flux_f115w']) -48.6
mag_y = -2.5*np.log10(t['flux_HSC-Y_DR3']) -48.6

mag_VIS = -2.5*np.log10(t['flux_VIS']) -48.6
mag_i = -2.5*np.log10(t['flux_HSC-I_DR3']) -48.6


######################! MAG - MAG VS MAG ###############################

# Compare VIS and HSC-I
plt.figure(figsize=(10, 10))
plt.scatter(mag_VIS, mag_VIS-mag_i, s=1, c='k', alpha=0.6)
plt.xlabel(r'$VIS$')
plt.ylabel(r'$VIS-i_{HSC}$')

# Sort data by mag_VIS
sorted_indices = np.argsort(mag_VIS)
mag_VIS_sorted = mag_VIS[sorted_indices]
mag_VIS_diff_sorted = (mag_VIS - mag_i)[sorted_indices]

# Compute the running median for mag_VIS_diff_sorted
window_size = 200
mag_VIS_diff_median = running_median(mag_VIS_diff_sorted, window_size)
plt.plot(mag_VIS_sorted, mag_VIS_diff_median[1:], color='red', label='Running Median', linewidth=2)
plt.plot([14, 26], [0, 0], color='green', label='Zero', linewidth=2)

plt.ylim(-1, 0.5)
plt.xlim(15, 24)

plt.show()
plt.close()

# Compare HSC-Y and vista-Y
plt.figure(figsize=(10, 10))
plt.scatter(mag_y, mag_y-mag_Y, s=1, c='k', alpha=0.6)
plt.xlabel(r'$y_{HSC}$')
plt.ylabel(r'$y_{HSC}-Y_{VISTA}$')

# Sort data by mag_Y
sorted_indices = np.argsort(mag_y)
mag_y_sorted = mag_y[sorted_indices]
mag_y_diff_sorted = (mag_y - mag_Y)[sorted_indices]

# Compute the running median for mag_Y_diff_sorted
window_size = 200
mag_y_diff_median = running_median(mag_y_diff_sorted, window_size)
plt.plot(mag_y_sorted, mag_y_diff_median[1:], color='red', label='Running Median', linewidth=2)
plt.plot([14, 26], [0, 0], color='green', label='Zero', linewidth=2)

plt.ylim(-1, 0.5)
plt.xlim(15, 24)

plt.show()
plt.close()


#! Compare Y bands
plt.figure(figsize=(10, 6))
plt.scatter(mag_Y, mag_Y-mag_Ye, s=1, c='k', alpha=0.6)
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
plt.scatter(mag_J, mag_J-mag_Je, s=1, c='k', alpha=0.6)
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
plt.scatter(mag_H, mag_H-mag_He, s=1, c='k', alpha=0.6)
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

# Compare Ye and f115w
plt.figure(figsize=(10, 10))
plt.scatter(mag_Ye, mag_Ye-mag_f115w, s=1, c='k', alpha=0.6)
plt.xlabel(r'$Y_{e}$')
plt.ylabel(r'$Y_{e}-\mathrm{f115w}$')

# Sort data by mag_Y
sorted_indices = np.argsort(mag_Ye)
mag_Ye_sorted = mag_Ye[sorted_indices]
mag_Ye_diff_sorted = (mag_Ye - mag_f115w)[sorted_indices]

# Compute the running median for mag_Ye_diff_sorted
# window_size = 50
# mag_Ye_diff_median = running_median(mag_Ye_diff_sorted, window_size)
# plt.plot(mag_Ye_sorted, mag_Ye_diff_median[1:], color='red', label='Running Median', linewidth=2)
plt.plot([14, 26], [0, 0], color='green', label='Zero', linewidth=2)

plt.ylim(-1, 0.5)
plt.xlim(15, 24)

plt.show()
plt.close()

exit()
#########################! MAG VS MAG ###########################

# Compare Y and Ye
plt.figure(figsize=(10, 10))
plt.scatter(mag_Y, mag_Ye, s=1, c='k', alpha=0.6)
plt.plot([16, 24], [16, 24], color='red', label='1:1', linewidth=2)

plt.xlabel(r'$Y$')
plt.ylabel(r'$Y_{e}$')

plt.show()
plt.close()

# Compare J and Je
plt.figure(figsize=(10, 10))
plt.scatter(mag_J, mag_Je, s=1, c='k', alpha=0.6)
plt.plot([16, 24], [16, 24], color='red', label='1:1', linewidth=2)

plt.xlabel(r'$J$')
plt.ylabel(r'$J_{e}$')

plt.show()
plt.close()

# Compare H and He
plt.figure(figsize=(10, 10))
plt.scatter(mag_H, mag_He, s=1, c='k', alpha=0.6)
plt.plot([16, 24], [16, 24], color='red', label='1:1', linewidth=2)

plt.xlabel(r'$H$')
plt.ylabel(r'$H_{e}$')

plt.show()
plt.close()