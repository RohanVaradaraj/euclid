"""
Compare my results to sims, at request of steve wilkins.

Created: Tuesday 19th August 2025.
"""
from pathlib import Path
from astropy.table import Table
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.cm as cm
import cmasher as cmr

plt.rcParams.update({'font.size': 15})
plt.rcParams['axes.linewidth'] = 4
plt.rcParams['figure.dpi'] = 100

plt.rcParams.update({
    # Ticks on all sides, pointing inwards
    'xtick.top': True, 'xtick.bottom': True,
    'ytick.left': True, 'ytick.right': True,
    'xtick.direction': 'in', 'ytick.direction': 'in',

    # Major tick size and width
    'xtick.major.size': 8, 'ytick.major.size': 8,
    'xtick.major.width': 2, 'ytick.major.width': 2,

    # Minor tick size and width
    'xtick.minor.size': 4, 'ytick.minor.size': 4,
    'xtick.minor.width': 1.5, 'ytick.minor.width': 1.5,
})

fontsize=28
labelsize=24

def dpl_fit(M, phiStar, Mstar, alpha, beta):
    """ For fitting the LFs"""
    return np.log(10) * phiStar / (2.5 * (10 ** (0.4 * (alpha + 1) * (M - Mstar)) + 10 ** (0.4 * (beta + 1) * (M - Mstar))))

def sample_asym_gaussian(mean, err_minus, err_plus, size=1, rng=None):
    """Sample from an asymmetric Gaussian (split normal)."""
    if rng is None:
        rng = np.random.default_rng()
    # Choose side
    u = rng.random(size)
    samples = np.empty(size)
    mask = u < (err_minus / (err_minus + err_plus))
    samples[mask]  = rng.normal(mean, err_minus, mask.sum())
    samples[~mask] = rng.normal(mean, err_plus, (~mask).sum())
    return samples

M_star_best = -21.1431
phi_star_best = 9.0598e-05
alpha_best = -2.1029
beta_best = -4.6266

M_star_err_minus, M_star_err_plus   = 0.25054261288004653, 0.27714872059357987
phi_star_err_minus, phi_star_err_plus = 3.806543369065727e-05, 6.704259904493972e-05
alpha_err_minus, alpha_err_plus     = 0.17300217202221013, 0.21467644551429
beta_err_minus, beta_err_plus       = 0.385521158264428, 0.3357791319020107

data_dir = Path.cwd() / 'flags_data' / 'flags_data' / 'data' / 'DistributionFunctions' / 'LUV' / 'models' / 'binned'

# Get all files in data dir
files = list(data_dir.glob('*.ecsv'))

# Read all as tables
tables = [Table.read(file, format='ascii.ecsv') for file in files]



# Restrict all tables to where z=7.0
tables = [table[table['z'] == 7.0] for table in tables]

# Remove empty tables
tables = [table for table in tables if len(table) > 0]

# Restrict tables to 'M' column and last column
for table in tables:
    # If colname is phi, take log and rename to log10phi. If already log10phi, leave it
    if 'log10phi' in table.colnames:
        table['phi'] = 10**table['log10phi']
        table.remove_column('log10phi')
    if 'log10L' in table.colnames:
        table['M'] = -2.5* table['log10L'] +51.6
        table.remove_column('log10L') 
    print(table.colnames, table.meta['name'])

# My LF
Muv_bins = [-22.4, -22., -21.8,  -21.6, -21.4, -21.2, -21.0, -20.8, -20.6, -20.4, -20.2]
bin_widths = np.abs(np.diff(Muv_bins))
bin_centres = 0.5 * (np.array(Muv_bins[:-1]) + np.array(Muv_bins[1:]))

LF_values = [2.12861662e-06, 6.03210043e-06, 9.11791980e-06, 1.59486164e-05,
 3.23233992e-05, 4.41687494e-05, 5.84964344e-05, 6.61033853e-05,
 6.20456212e-05, 3.48707182e-05]
LF_errors = [1.28489889e-06, 3.18765338e-06, 4.01277818e-06, 5.63136911e-06,
 8.82750291e-06, 1.08894452e-05, 1.42224172e-05, 1.62241218e-05,
 1.75968840e-05, 1.64927900e-05]

# Plot all tables
#plt.figure(figsize=(10, 10))

# Colormap for distinct colors
colors = colors = cmr.take_cmap_colors('cmr.guppy', N=len(tables))

# Cycle through line styles
linestyles = ['-', '--', ':']

# Varadaraj+23 LF points
v23x = [-22.175, -22.925, -23.675]#, -24.425])
v23y = [2.70e-6, 2.81e-7, 2.37e-8]#, np.nan])  
v23dy = [0.66e-6, 1.54e-7, 2.50e-8]#, np.nan])  
v23x, v23y, v23dy = np.array(v23x), np.array(v23y), np.array(v23dy)
#plt.errorbar(v23x+0.05, v23y, yerr=v23dy, fmt='o', color='black', 
#            ecolor='black', elinewidth=2, label='Varadaraj+23', markersize=15, markeredgecolor='black', zorder=4)  

# Finkelstein+15
f15x = [-22.0, -21.5, -21.0, -20.5, -20.0, -19.5, -19.0, -18.5, -18.0]
f15y = [0.0046e-3, 0.0187e-3, 0.0690e-3, 0.1301e-3, 0.2742e-3, 0.3848e-3, 0.5699e-3, 2.5650e-3, 3.0780e-3]
f15y_up = [0.0049e-3, 0.0085e-3, 0.0156e-3, 0.0239e-3, 0.0379e-3, 0.0633e-3, 0.2229e-3, 0.8735e-3, 1.0837e-3]
f15y_lo = [0.0028e-3, 0.0067e-3, 0.0144e-3, 0.0200e-3, 0.0329e-3, 0.0586e-3, 0.1817e-3, 0.7161e-3, 0.8845e-3]
f15x, f15y, f15y_up, f15y_lo = np.array(f15x), np.array(f15y), np.array(f15y_up), np.array(f15y_lo)

fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 12), sharex=True,
                               gridspec_kw={'height_ratios': [3, 1]})

# ---- Top panel: LFs ----
for i, (table, c) in enumerate(zip(tables, colors)):
    ls = linestyles[i % len(linestyles)]
    ax1.plot(table['M'], table['phi'],
             color=c, linestyle=ls, linewidth=3,
             label=table.meta.get('name', f'Sim {i+1}'),
             alpha=0.9)

M_fit = np.linspace(-25, -18, 200)
LF_fit = dpl_fit(M_fit, phi_star_best, M_star_best, alpha_best, beta_best)
ax1.plot(M_fit, LF_fit, color='tab:red', linewidth=7,
         label="Best-fit DPL", alpha=0.9)

# Your LF points
ax1.errorbar(bin_centres[:-3], LF_values[:-3], yerr=LF_errors[:-3],
             fmt='o', color='tab:red', markersize=15, zorder=10,
             markeredgecolor='black', label='This work', elinewidth=3,)

# Varadaraj+23
ax1.errorbar(v23x+0.05, v23y, yerr=v23dy, fmt='o', color='black',
             ecolor='black', elinewidth=3, markersize=15,
             markeredgecolor='black', zorder=4, label='Varadaraj+23')

# Finkelstein+15
ax1.errorbar(f15x, f15y, yerr=[f15y_lo, f15y_up],
             fmt='p', color='tab:purple', alpha=0.9,
             markersize=12, zorder=4, label='Finkelstein+15')

ax1.set_ylim([1e-9, 5e-3])
ax1.set_xlim(-24.5, -19)
ax1.set_yscale('log')
ax1.set_ylabel(r'$\log_{10}(\phi\,/\,\rm{mag}^{-1}\,\rm{Mpc}^{-3})$', fontsize=fontsize)

# ---- Bottom panel: residuals ----
for i, (table, c) in enumerate(zip(tables, colors)):
    ls = linestyles[i % len(linestyles)]
    # Interpolate the best-fit LF onto the same M grid
    LF_fit_interp = dpl_fit(table['M'], phi_star_best, M_star_best, alpha_best, beta_best)
    ratio = table['phi'] / LF_fit_interp
    ax2.plot(table['M'], ratio, color=c, linestyle=ls, linewidth=3, alpha=0.9)

ax2.plot(M_fit, np.ones_like(M_fit), color='tab:red', lw=2, alpha=1, zorder=-1,)

#! Shaded region for DPL uncertainties
n_samples = 3000
rng = np.random.default_rng(42)

M_star_samples   = sample_asym_gaussian(M_star_best,   M_star_err_minus,   M_star_err_plus,   n_samples, rng)
phi_star_samples = sample_asym_gaussian(phi_star_best, phi_star_err_minus, phi_star_err_plus, n_samples, rng)
alpha_samples    = sample_asym_gaussian(alpha_best,    alpha_err_minus,    alpha_err_plus,    n_samples, rng)
beta_samples     = sample_asym_gaussian(beta_best,     beta_err_minus,     beta_err_plus,     n_samples, rng)

LF_samples = []
for i in range(n_samples):
    LF_samples.append(dpl_fit(M_fit,
                              phi_star_samples[i],
                              M_star_samples[i],
                              alpha_samples[i],
                              beta_samples[i]))
LF_samples = np.array(LF_samples)

# Compute median and 68% intervals
LF_median = np.median(LF_samples, axis=0)
LF_lo = np.percentile(LF_samples, 16, axis=0)
LF_hi = np.percentile(LF_samples, 84, axis=0)

# Plot best fit + shaded envelope
#ax1.plot(M_fit, LF_median, color='tab:red', lw=7, label="Best-fit DPL", alpha=0.9)
ax1.fill_between(M_fit, LF_lo, LF_hi, color='tab:red', alpha=0.25, zorder=-1, edgecolor='none')

# ax2.axhline(0, color='black', lw=2)  # reference line
ax2.set_xlim(-24.5, -19)
ax2.set_ylabel(r'$\phi_{\rm sim} / \phi_{\rm DPL}$', fontsize=fontsize)
ax2.set_xlabel(r'$M_{\rm UV}$', fontsize=fontsize)

# Draw the LF_lo and LF_hi ratios on the bottom panel
ax2.fill_between(M_fit, LF_lo / LF_median, LF_hi / LF_median, color='tab:red', alpha=0.25, zorder=-1, edgecolor='none')

plt.tight_layout()
# plt.show()

#plt.errorbar(f15x, f15y, yerr=[f15y_lo, f15y_up], fmt='p', color='tab:purple', alpha=0.9, label='Finkelstein+15', markersize=12, zorder=4)

ax1.set_ylim([1e-9, 2e-3])
ax1.set_xlim(-24., -19)

ax2.set_ylim([0, 5])

ax1.set_yscale('log')

yticks = [1e-9, 1e-8, 1e-7, 1e-6, 1e-5, 1e-4, 1e-3]
ax1.set_yticks(yticks)
ax1.set_yticklabels([r"$-9\,$", r"$-8\,$", r"$-7\,$", r"$-6\,$", r"$-5\,$", r"$-4\,$", r"$-3\,$"])
# Make axis labels larger
ax1.tick_params(axis='both', which='major', labelsize=labelsize)

# Add minor ticks on ax1
ax1.minorticks_on()

# Make ax2 label size larger
ax2.tick_params(axis='both', which='major', labelsize=labelsize)
# Set specific ax2 y axis labels, like -1, 0, 1, 2, 3
ax2.set_yticks([0, 1,2, 3, 4, 5])
# Minor ticks on y axis of ax2
ax2.minorticks_on()


ax1.legend(loc='lower right', fontsize=18, frameon=False, ncols=2)
plot_dir = Path.cwd().parents[1] / 'plots' / 'LF'
plt.savefig(plot_dir / 'LF_comparison_to_sims.pdf', bbox_inches='tight')
plt.show()



