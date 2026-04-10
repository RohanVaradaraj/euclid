"""
Some plots showing the distribution of z=6 sources in the VISTA fields.

Created: Friday 27th March 2026.
"""

import matplotlib.pyplot as plt
from astropy.table import Table
import numpy as np
from pathlib import Path

plt.rcParams.update({'font.size': 15})
plt.rcParams['axes.linewidth'] = 3
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

table_dir = Path.cwd().parents[1] / 'data' / 'catalogues' / 'candidates'

xmm_catname = 'XMM_5sig_HSC_Z_nonDet_HSC_G_nonDet_HSC_R_candidates_2025_05_14.fits'
cosmos_catname = 'COSMOS_5sig_HSC_Z_nonDet_HSC_G_nonDet_HSC_R_candidates_2025_06_06.fits'

xmm_table = Table.read(table_dir / xmm_catname)
cosmos_table = Table.read(table_dir / cosmos_catname)

print(cosmos_table.colnames)

plot_dir = Path.cwd().parents[1] / 'plots' / 'z6_sample_plots'

#! Plot histograms of Muv
# plt.figure(figsize=(8,6))

# bins=np.arange(-23, -19, 0.1)

# plt.hist(xmm_table['Muv'], bins=bins, alpha=0.8, label='XMM', color='tab:blue', histtype='step', linewidth=2)
# plt.hist(cosmos_table['Muv'], bins=bins, alpha=0.8, label='COSMOS', color='tab:orange', histtype='step', linewidth=2)
# plt.xlabel(r'$M_{\mathrm{UV}}$')
# plt.ylabel('Count')
# plt.legend()
# plt.savefig(plot_dir / 'Muv_histogram.pdf', bbox_inches='tight')
# plt.show()
# plt.close()

#! Plot histograms of mAB from flux_HSC-Z_DR3
# plt.figure(figsize=(8,6))
# bins=np.arange(24, 27, 0.1)
# mag_cos = -2.5*np.log10(cosmos_table['flux_HSC-Z_DR3']) -48.6
# mag_xmm = -2.5*np.log10(xmm_table['flux_HSC-Z_DR3']) -48.6

# plt.hist(mag_xmm, bins=bins, alpha=0.8, label='XMM', color='tab:blue', histtype='step', linewidth=2)
# plt.hist(mag_cos, bins=bins, alpha=0.8, label='COSMOS', color='tab:orange', histtype='step', linewidth=2)
# plt.xlabel(r'$m_{\mathrm{AB}}$')
# plt.ylabel('Count')
# plt.legend()
# plt.savefig(plot_dir / 'mAB_histogram.pdf', bbox_inches='tight')
# plt.show()
# plt.close()

#! Plot histogram of Zphot
# plt.figure(figsize=(8,6))
# bins=np.arange(5.5, 6.5, 0.05)
# plt.hist(xmm_table['Zphot'], bins=bins, alpha=0.8, label='XMM', color='tab:blue', histtype='step', linewidth=2)
# plt.hist(cosmos_table['Zphot'], bins=bins, alpha=0.8, label='COSMOS', color='tab:orange', histtype='step', linewidth=2)
# plt.xlabel(r'$z_{\mathrm{phot}}$')
# plt.ylabel('Count')
# plt.legend()
# plt.savefig(plot_dir / 'Zphot_histogram.pdf', bbox_inches='tight')
# plt.show()
# plt.close()

#! Scatter plot of Muv vs Zphot, with errors from dMuv_inf, dMuv_sup, Zinf, Zsup
plt.figure(figsize=(9,6))

xmm_zphot_err = [np.abs(xmm_table['Zphot'] - xmm_table['Zinf']), np.abs(xmm_table['Zsup'] - xmm_table['Zphot'])]
cos_zphot_err = [np.abs(cosmos_table['Zphot'] - cosmos_table['Zinf']), np.abs(cosmos_table['Zsup'] - cosmos_table['Zphot'])]

plt.errorbar(xmm_table['Zphot'], xmm_table['Muv'], xerr=xmm_zphot_err, yerr=[xmm_table['dMuv_inf'], xmm_table['dMuv_sup']], 
fmt='o', label='XMM', color='tab:blue', ecolor='lightblue', elinewidth=2, capsize=4, alpha=0.5)

plt.errorbar(cosmos_table['Zphot'], cosmos_table['Muv'], xerr=cos_zphot_err, yerr=[cosmos_table['dMuv_inf'], cosmos_table['dMuv_sup']],
fmt='o', label='COSMOS', color='tab:orange', ecolor='navajowhite', elinewidth=2, capsize=4, alpha=0.5)

# Flip y axis
plt.gca().invert_yaxis()
plt.xlabel(r'$z_{\mathrm{phot}}$')
plt.ylabel(r'$M_{\mathrm{UV}}$')
plt.legend()
plt.savefig(plot_dir / 'Muv_vs_Zphot.pdf', bbox_inches='tight')
plt.show()

