"""
Takes the things cut as BDs in the VISTA selection and compares to the SED fitting with Euclid.

Then plot the fraction of sources misidentified as BD as function of apparent magnitude.
"""

import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table
from scipy.stats import binned_statistic
from astropy.io import ascii
import sys
from pathlib import Path
import glob
import scipy.stats.distributions as dist
from astropy.cosmology import FlatLambdaCDM
import matplotlib.ticker as mticker

sed_fiting_path = Path.cwd().parent / 'sed_fitting'
sys.path.append(str(sed_fiting_path))
from sed_fitting_codes import *

plt.rcParams['axes.linewidth'] = 2.5
plt.rcParams.update({'font.size': 15})
plt.rcParams['figure.dpi'] = 100
names_param = ['Type', 'Nline', 'Model', 'Library', 'Nband', 'Zphot', 'Zinf', 'Zsup', 'Chi2', 'PDF', 'Extlaw', 'EB-V', 'Lir', 'Age', 'Mass', 'SFR', 'SSFR']

# Cosmology
H = 70
omegaM = 0.3
omegaV = 0.7
cosmo = FlatLambdaCDM(H0=H, Om0=omegaM)

def apparent_to_absolute(mag):
    DL = cosmo.luminosity_distance(7).value * 10 ** 6
    M1500 = mag - 5*np.log10(DL/10) + 2.5*np.log10(1+7)
    return M1500

def absolute_to_apparent(M1500):
    DL = cosmo.luminosity_distance(7).value * 10 ** 6
    mag = M1500 + 5*np.log10(DL/10) - 2.5*np.log10(1+7)
    return mag

#! Set up directories
base_sed_dir = Path.cwd().parents[1] / 'data' / 'sed_fitting' / 'zphot'

#! VISTA BD fitting
vista_bd_dir = base_sed_dir / 'best_fits' / 'det_Y_J_BD'

#! Euclid BD fitting
euclid_bd_dir = base_sed_dir / 'best_fits' / 'det_Y_J_with_euclid_BD_PLUS_EUCLID_PHOT'

#! Catalogue with photometry
catalogue_path = Path.cwd().parents[1] / 'data' / 'catalogues' / 'COSMOS_5sig_Y_J_nonDet_HSC_G_nonDet_HSC_R_nonDet_HSC_I.fits'
t = Table.read(catalogue_path)

vista_files = glob.glob(str(vista_bd_dir / '*.spec'))

# Properties to compare
IDs = []
VISTA_Y_mags = []
VISTA_J_mags = []
Euclid_Y_mags = []
Euclid_J_mags = []

is_VISTA_BD = []
is_Euclid_BD = []

for i, spec_file in enumerate(vista_files):

    spec_file_name = spec_file.split('/')[-1]

    ID = spec_file.split('/')[-1].split('Id')[-1].lstrip('0').split('.spec')[0]

    vista_file = parse_spec_file(spec_file)
    euclid_file = parse_spec_file(euclid_bd_dir / spec_file_name)

    vista_params = vista_file.get('model')
    vista_params.rename_columns(vista_params.colnames, names_param)
    vista_chi2_gal = vista_params['Chi2'][0]
    vista_chi2_bd = vista_params['Chi2'][-1]

    euclid_params = euclid_file.get('model')
    euclid_params.rename_columns(euclid_params.colnames, names_param)
    euclid_chi2_gal = euclid_params['Chi2'][0]
    euclid_chi2_bd = euclid_params['Chi2'][-1]

    #! Check if the source is a BD in VISTA
    is_VISTA_BD.append(vista_chi2_bd < 10)

    #! Check if the source is a BD in Euclid
    is_Euclid_BD.append(euclid_chi2_bd < euclid_chi2_gal)

    #! Get magnitudes
    # Define the lambda function for magnitude calculation
    calc_magnitude = lambda flux: -2.5 * np.log10(flux) - 48.6

    # Get magnitudes
    VISTA_Y = calc_magnitude(t[t['ID'] == int(ID)]['flux_Y'])[0]
    VISTA_J = calc_magnitude(t[t['ID'] == int(ID)]['flux_J'])[0]
    Euclid_Y = calc_magnitude(t[t['ID'] == int(ID)]['flux_Ye'])[0]
    Euclid_J = calc_magnitude(t[t['ID'] == int(ID)]['flux_Je'])[0]


    IDs.append(ID)
    VISTA_Y_mags.append(VISTA_Y)
    VISTA_J_mags.append(VISTA_J)
    Euclid_Y_mags.append(Euclid_Y)
    Euclid_J_mags.append(Euclid_J)


# Convert to numpy arrays
IDs = np.array(IDs)
VISTA_Y_mags = np.array(VISTA_Y_mags)
VISTA_J_mags = np.array(VISTA_J_mags)
Euclid_Y_mags = np.array(Euclid_Y_mags)
Euclid_J_mags = np.array(Euclid_J_mags)
is_VISTA_BD = np.array(is_VISTA_BD)
is_Euclid_BD = np.array(is_Euclid_BD)

# # Get the fraction of BDs in Euclid over BDs in VISTA
fraction_bd = []
yerr_lower = []
yerr_upper = []

# Choose which magnitude to use
mag_to_use = VISTA_J_mags

# Plot fraction of is_bd in euclid over is_bd in vista in bins of magnitude
bins = np.arange(np.min(mag_to_use), np.max(mag_to_use), 0.3)

c = 0.68 # Confidence level

for i in range(len(bins) - 1):
    mask = (mag_to_use > bins[i]) & (mag_to_use< bins[i + 1])
    
    # Calculate the fraction of misidentified BDs
    fraction = (np.sum(is_Euclid_BD[mask]) / np.sum(is_VISTA_BD[mask]))
    fraction_bd.append(fraction)

    n = np.sum(is_VISTA_BD[mask])
    k = np.sum(is_Euclid_BD[mask])

    # Cameron+10, https://arxiv.org/pdf/1012.0566
    p_lower = dist.beta.ppf((1 - c) / 2., k + 1, n - k + 1)
    p_upper = dist.beta.ppf(1 - (1 - c) / 2., k + 1, n - k + 1)

    dy_low = fraction - p_lower
    dy_up = p_upper - fraction

    if dy_low < 0:
        dy_low = 0
    if dy_up < 0:
        dy_up = 0

    yerr_lower.append(dy_low)
    yerr_upper.append(dy_up)


fraction_bd = np.array(fraction_bd)
yerr_lower = np.array(yerr_lower)
yerr_upper = np.array(yerr_upper)

# Plot with error bars and secondary x-axis
fig, ax1 = plt.subplots(figsize=(8, 6))

# Plot the primary y-axis (fraction of BDs)
ax1.errorbar(bins[:-1], fraction_bd, yerr=[yerr_lower, yerr_upper], 
             fmt='o-', color='black', capsize=3, capthick=2, linewidth=3, markersize=10)
ax1.set_xlabel(r'$J$ (AB mag)')
ax1.set_ylabel('VISTA BD identification success rate')

# Create a secondary x-axis for absolute UV magnitude at z=7
secax = ax1.secondary_xaxis('top', functions=(apparent_to_absolute, absolute_to_apparent))

secax.set_xlabel(r'$M_{\mathrm{UV}}$')
secax.yaxis.set_major_formatter(mticker.ScalarFormatter())
plot_dir = Path.cwd().parents[1] / 'plots' / 'brown_dwarfs'
plt.savefig(plot_dir / 'BD_identification_rate_VISTA.pdf')
plt.show()



