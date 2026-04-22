"""
Checking flat stars in the DR1 data to assess the zeropoint issue exploring with Holly.

Created: Tuesday 21st April 2026.
"""

from astropy.table import Table
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt


def sigma_clip_xy(x, y, sigma=3.0):
    finite_mask = np.isfinite(x) & np.isfinite(y)
    x = np.asarray(x)[finite_mask]
    y = np.asarray(y)[finite_mask]

    median_y = np.median(y)
    mad_y = np.median(np.abs(y - median_y))
    robust_sigma = 1.4826 * mad_y

    if robust_sigma == 0 or len(y) < 2:
        return x, y, np.ones(len(y), dtype=bool)

    clip_mask = np.abs(y - median_y) <= sigma * robust_sigma
    return x[clip_mask], y[clip_mask], clip_mask

plt.rcParams.update({'font.size': 15})
plt.rcParams['axes.linewidth'] = 2.5
plt.rcParams['figure.dpi'] = 100

plt.rcParams.update({
    # Ticks on all sides, pointing inwards
    'xtick.top': True, 'xtick.bottom': True,
    'ytick.left': True, 'ytick.right': True,
    'xtick.direction': 'in', 'ytick.direction': 'in',

    # Major tick size and width
    'xtick.major.size': 6.5, 'ytick.major.size': 6.5,
    'xtick.major.width': 2, 'ytick.major.width': 2,

    # Minor tick size and width
    'xtick.minor.size': 3, 'ytick.minor.size': 3,
    'xtick.minor.width': 1.5, 'ytick.minor.width': 1.5,
})


cat_dir = Path.cwd().parents[1] / 'data' / 'catalogues'

t = Table.read(cat_dir / 'COSMOS_Euclid_DR1_det_YE_zpt_check.fits')

y = -2.5*np.log10(t['flux_Y']) - 48.6
j = -2.5*np.log10(t['flux_J']) - 48.6
h = -2.5*np.log10(t['flux_H']) - 48.6
# Select flat stars
flat_star_mask = (np.abs(y - j) < 0.05) & (np.abs(j - h) < 0.05) & (t['CLASS_STAR'] > 0.99)
mag_cut_y_mask = (y > 17.) & (y < 20.)
mag_cut_j_mask = (j > 17.) & (j < 20.)
mag_cut_h_mask = (h > 17.) & (h < 20.)


t_flat_stars_y = t[flat_star_mask & mag_cut_y_mask]
t_flat_stars_j = t[flat_star_mask & mag_cut_j_mask]
t_flat_stars_h = t[flat_star_mask & mag_cut_h_mask]


yye = np.log10(t_flat_stars_y['flux_Y']) - np.log10(t_flat_stars_y['flux_YE_DR1hom'])


jje = np.log10(t_flat_stars_j['flux_J']) - np.log10(t_flat_stars_j['flux_JE_DR1hom'])


hhe = np.log10(t_flat_stars_h['flux_H']) - np.log10(t_flat_stars_h['flux_HE_DR1hom'])

print(len(yye))

data_y = [yye, jje, hhe]
data_x = [y[flat_star_mask & mag_cut_y_mask], j[flat_star_mask & mag_cut_j_mask], h[flat_star_mask & mag_cut_h_mask]]

names_y = [r'$Y-Y_{\mathrm{E}}$', r'$J-J_{\mathrm{E}}$', r'$H-H_{\mathrm{E}}$']
names_x = [r'$Y$', r'$J$', r'$H$']
colours = ['tab:blue', 'tab:green', 'tab:red']

for i, vals in enumerate(data_y):
    clipped_x, clipped_y, clip_mask = sigma_clip_xy(data_x[i], vals, sigma=3.0)

    best_horizontal_line = np.median(clipped_y)
    vertical_shift_to_zero = -best_horizontal_line
    shifted_y = clipped_y + vertical_shift_to_zero
    mad_clipped_y = np.median(np.abs(clipped_y - best_horizontal_line))
    robust_sigma = 1.4826 * mad_clipped_y
    scatter_on_median = 1.2533 * robust_sigma / np.sqrt(len(clipped_y))
    n_removed = len(vals) - len(clipped_y)

    print(
        f"{names_y[i]}: best-fit horizontal line (median) = {best_horizontal_line:.6f}, "
        f"shift to move line to y=0 = {vertical_shift_to_zero:.6f}, "
        f"standard error on median = {scatter_on_median:.6f}, "
        f"removed {n_removed} points with 3-MAD clipping"
    )

    plt.figure(figsize=(10,6))

    plt.scatter(clipped_x, clipped_y, color=colours[i], label='Original')
    plt.scatter(clipped_x, shifted_y, color='black', alpha=0.7, label='Shifted to y=0')
    plt.axhline(best_horizontal_line, color=colours[i], linewidth=2.5,
                label=fr'Best-fit = {best_horizontal_line:.4f}')

    plt.ylabel(names_y[i])
    plt.xlabel(names_x[i])

    plt.axhline(0, color='black')
    plt.axhline(0.05, linestyle='--', color='black')
    plt.axhline(-0.05, linestyle='--', color='black')

    plt.ylim(-0.3, 0.3)
    plt.legend(frameon=False)

    plt.show()
