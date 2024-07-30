"""
Compare the SED fitting with and without Euclid filters.

Created: Tuesday 30th July 2024.
"""

import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path
from astropy.table import Table
import seaborn as sns

sns.set_palette("colorblind")
plt.rcParams['axes.linewidth'] = 2.5
plt.rcParams.update({'font.size': 15})
plt.rcParams['figure.dpi'] = 100

# Load the SED fitting results
table_dir = Path.cwd().parents[1] / 'data' / 'sed_fitting' / 'catalogues' / 'comparison'
table_name = 'det_Ye_Y_comparison.fits'
t = Table.read(table_dir / table_name)

#! Redshift distribution of GAL-1
# plt.hist(t['zphot_pri_1'], bins=np.arange(0, 8, 0.1), alpha=0.5, label='No Euclid')
# plt.hist(t['zphot_pri_2'], bins=np.arange(0, 8, 0.1), alpha=0.5, label='Euclid')
# plt.legend()
# plt.xlabel('GAL-1 zphot')
# plt.show()
# plt.close()

#! Redshift distribution of GAL-1, accounting for where BD fit is preferred.
# zphots_1 = t['zphot_pri_1']
# zphots_2 = t['zphot_pri_2']
# zphots_1[t['chi2_star_1'] < t['chi2_pri_1']] = 0
# zphots_2[t['chi2_star_2'] < t['chi2_pri_2']] = 0
# plt.hist(zphots_1, bins=np.arange(0, 8, 0.1), alpha=0.5, label='No Euclid')
# plt.hist(zphots_2, bins=np.arange(0, 8, 0.1), alpha=0.5, label='Euclid')
# plt.legend()
# plt.xlabel('GAL-1 zphot')
# plt.show()
# plt.close()

#! Magnitude of error on GAL-1
# t = t[t['zphot_pri_2'] > 4.5]
# dz_1 = t['dz_sup_1'] - t['dz_inf_1']
# dz_2 = t['dz_sup_2'] - t['dz_inf_2']
# # plt.scatter(t['zphot_pri_1'], dz_1, alpha=0.5, label='No Euclid')
# # plt.scatter(t['zphot_pri_2'], dz_2, alpha=0.5, label='Euclid')
# plt.hist(dz_1, bins=np.arange(0, 1, 0.05), alpha=0.5, label='No Euclid', density=True)
# plt.hist(dz_2, bins=np.arange(0, 1, 0.05), alpha=0.5, label='Euclid', density=True)
# plt.legend()
# plt.title(r'$z > 4.5$')
# #plt.xlabel('GAL-1 zphot')
# plt.xlabel(r'$\delta z$')
# plt.show()

#! Difference between primary chi2 and stellar chi2
# chi2_diff_1 = t['chi2_star_1'] - t['chi2_pri_1']
# chi2_diff_2 = t['chi2_star_2'] - t['chi2_pri_2']
# plt.scatter(t['zphot_pri_1'], chi2_diff_1, alpha=0.5, label='No Euclid')
# plt.scatter(t['zphot_pri_2'], chi2_diff_2, alpha=0.5, label='Euclid')
# plt.xlabel('GAL-1 zphot')
# plt.ylabel(r'$\Delta \chi^{2}$' + ' (Primary - Stellar)')
# plt.legend()
# plt.show()

#! Primary redshift with and without euclid
# zphots_1 = t['zphot_pri_1']
# zphots_2 = t['zphot_pri_2']
# zphots_1[t['chi2_star_1'] < t['chi2_pri_1']] = 0
# zphots_2[t['chi2_star_2'] < t['chi2_pri_2']] = 0
# plt.scatter(zphots_1, zphots_2, alpha=0.5)
# plt.plot([0, 8], [0, 8], 'k--')
# plt.xlabel('z (No Euclid)')
# plt.ylabel('z (Euclid)')
# plt.show()

#! Difference in primary redshift
# zphots_1 = t['zphot_pri_1']
# zphots_2 = t['zphot_pri_2']
# zphots_1[t['chi2_star_1'] < t['chi2_pri_1']] = 0
# zphots_2[t['chi2_star_2'] < t['chi2_pri_2']] = 0
# plt.scatter(zphots_2, zphots_1-zphots_2, alpha=0.5)
# plt.ylabel(r'$\Delta z$')
# plt.xlabel('z (Euclid)')
# plt.show()

#! Reduced chi2 before and after
# chi2_1 = t['chi2_pri_1']
# chi2_2 = t['chi2_pri_2']
# # If stellar chi2 better, replace with stellar chi2
# chi2_1[t['chi2_star_1'] < t['chi2_pri_1']] = t['chi2_star_1'][t['chi2_star_1'] < t['chi2_pri_1']]
# chi2_2[t['chi2_star_2'] < t['chi2_pri_2']] = t['chi2_star_2'][t['chi2_star_2'] < t['chi2_pri_2']]
# # Divide by degrees of freedom.
# chi2_1 /= (15 - 6)
# chi2_2 /= (19 - 6)
# plt.scatter(chi2_1, chi2_2, alpha=0.5)
# plt.xlabel(r'$\chi^{2}_{\mathrm{reduced}}$ (No Euclid)')
# plt.ylabel(r'$\chi^{2}_{\mathrm{reduced}}$ (Euclid)')
# plt.show()

#! Difference between primary chi2 and secondary chi2
chi2_diff_1 = t['chi2_sec_1'] - t['chi2_pri_1']
chi2_diff_2 = t['chi2_sec_2'] - t['chi2_pri_2']
plt.scatter(t['zphot_pri_1'], chi2_diff_1, alpha=0.5, label='No Euclid')
plt.scatter(t['zphot_pri_2'], chi2_diff_2, alpha=0.5, label='Euclid')
plt.xlabel('GAL-1 zphot')
tmp = 3
plt.ylabel(rf'$\Delta \chi^{2}$ {tmp}' + ' (Primary - Secondary)')
plt.legend()
plt.show()

