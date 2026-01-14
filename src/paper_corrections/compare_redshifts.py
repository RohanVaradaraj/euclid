"""
Compare U-only vs U+E photometric redshifts.

Created: Wednesday 7th January 2025.
"""

from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np
from astropy.table import Table

plt.rcParams['axes.linewidth'] = 2.5
plt.rcParams.update({'font.size': 20})
plt.rcParams['figure.dpi'] = 100

plt.rcParams.update({
    # Ticks on all sides, pointing inwards
    'xtick.top': True, 'xtick.bottom': True,
    'ytick.left': True, 'ytick.right': True,
    'xtick.direction': 'in', 'ytick.direction': 'in',

    # Major tick size and width
    'xtick.major.size': 6.5, 'ytick.major.size': 6.5,
    'xtick.major.width': 3, 'ytick.major.width': 3,

    # Minor tick size and width
    'xtick.minor.size': 3, 'ytick.minor.size': 3,
    'xtick.minor.width': 2, 'ytick.minor.width': 2,
})

data_dir = Path.cwd().parent / 'galaxy_properties'
cat_dir = Path.cwd().parents[1] / 'data' / 'catalogues' / 'candidates'

# Load npy files
# ue = np.load(data_dir / 'z_Muv_sample_with_euclid.npy') 
# uo = np.load(data_dir / 'z_Muv_sample_.npy')

ue_cat = Table.read(cat_dir / 'Euclid_UltraVISTA_z7_sample_kron_piecewise.fits')
uo_cat = Table.read(cat_dir / 'UltraVISTA_only_z7_sample_MAG_AUTO.fits')

# Get ID arrays
uo_ids = uo_cat['ID']
ue_ids = ue_cat['ID']

# IDs common to both
common_ids = np.intersect1d(uo_ids, ue_ids)

# Rows in uo_cat that are also in ue_cat
uo_in_ue = uo_cat[np.isin(uo_ids, common_ids)]
print(uo_in_ue)

# Sort by IDs
uo_in_ue.sort('ID')
# Get zphot
z_uo = uo_in_ue['Zphot']

# We also have Zinf and Zsup, from these get errors
z_uo_err_low = z_uo - uo_in_ue['Zinf']
z_uo_err_up = uo_in_ue['Zsup'] - z_uo


# Restrict ue_cat to only those in common
ue_in_uo = ue_cat[np.isin(ue_ids, common_ids)]
print(ue_in_uo)
# Sort by IDs
ue_in_uo.sort('ID')
z_ue = ue_in_uo['Zphot']
z_ue_err_low = z_ue - ue_in_uo['Zinf']
z_ue_err_up = ue_in_uo['Zsup'] - z_ue

# Ploy z_U vs z_UE
plt.figure(figsize=(8, 6))
plt.xlim(6.25, 7.65)
plt.ylim(6.25, 7.65)

plt.errorbar(z_uo, z_ue, xerr=[z_uo_err_low, z_uo_err_up], yerr=[z_ue_err_low, z_ue_err_up], fmt='none', alpha=0.5, color='black')
plt.scatter(z_uo, z_ue, s=50, alpha=1, color='black', zorder=10)
plt.plot([6, 8], [6, 8], color='black', linestyle='--')
plt.xlabel(r'$z_{\rm{phot}}$ (UltraVISTA-only)')
plt.ylabel(r'$z_{\rm{phot}}$ (UltraVISTA+$Euclid$)')
plt.minorticks_on()  
plt.tight_layout()
plt.savefig('U_vs_UE_zphot_comparison.pdf', bbox_inches='tight', dpi=100)
plt.close()
# plt.show()

# Now plot the histograms of the size of the error as a function of Muv
plt.figure(figsize=(8, 6))
bins = np.linspace(-0.2, 0.2, 30)
plt.hist(z_uo_err_up-z_uo_err_low, bins=bins, alpha=0.5, color='dodgerblue', label='UltraVISTA-only')
plt.hist(z_uo_err_up-z_ue_err_low, bins=bins, alpha=0.5, color='magenta', label='UltraVISTA+Euclid')
plt.xlabel(r'$\Delta z$')
plt.ylabel('N')
plt.savefig('zerr_comparison.pdf', bbox_inches='tight', dpi=100)

