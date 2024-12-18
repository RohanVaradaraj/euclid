"""
Measure UV mags for galaxy candidates.

Created: Friday 8th November 2024.
"""

from astropy.table import Table, Column
from pathlib import Path
import numpy as np
import glob
import os
import matplotlib.pyplot as plt
from astropy.io import ascii, fits
from scipy.integrate import simps
from astropy.cosmology import FlatLambdaCDM
from scipy import stats
import matplotlib.pyplot as plt
import math

plt.rcParams['axes.linewidth'] = 2.5
plt.rcParams.update({'font.size': 15})
plt.rcParams['figure.dpi'] = 100

import sys
sed_path = Path.cwd().parents[0] / 'sed_fitting'
sys.path.append(str(sed_path))
from sed_fitting_codes import parse_spec_file

run_type = 'with_euclid'

def mag_to_flux(m):
    '''Convert mags to flux count'''
    flux = 10**(-0.4*(m+48.6))
    return flux

def flux_to_mag(f):
    mag = -2.5 * np.log10(f) - 48.6
    return mag

# Define the cosmology
H = 70
omegaM = 0.3
omegaV = 0.7
cosmo = FlatLambdaCDM(H0=H, Om0=omegaM)

#! SED Fitting folder
# Name of the directory we want to use to make the catalogue
folder = f'det_Y_J_{run_type}_z7' if run_type != '' else 'det_Y_J_z7'

# Get the list of objects that made it through the SED fitting
obj_dir = Path.cwd().parents[1] / 'data' / 'sed_fitting' / 'zphot' / 'best_fits'
obj_list = glob.glob(str(obj_dir / folder / '*.spec'))

# Get the IDs
IDs = [spec_file.split('/')[-1].split('Id')[-1].lstrip('0').split('.spec')[0] for spec_file in obj_list]
IDs = [int(ID) for ID in IDs]

#! Catalogue of above objects
# Parent catalogue from which to get fluxes
cat_name = 'COSMOS_5sig_Y_J_nonDet_HSC_G_nonDet_HSC_R_nonDet_HSC_I_candidates_2024_11_28_with_euclid.fits'

# Read in the parent catalogue
cat_dir = Path.cwd().parents[1] / 'data' / 'catalogues' / 'candidates'
t = Table.read(cat_dir / cat_name)

#! Add a column for UV mags if it doesnt exist already
if 'UV_mag' not in t.colnames:
    t['Muv'] = Column(np.zeros(len(t)), name='UV_mag', dtype=float)
    t['dMuv_inf'] = Column(np.zeros(len(t)), name='UV_mag_err', dtype=float)
    t['dMuv_sup'] = Column(np.zeros(len(t)), name='UV_mag_err', dtype=float)

# Loop through the objects in the table and get its SED properties
for i, ID in enumerate(IDs):

    # Find the row index in the table with this ID
    row_index = np.where(t['ID'] == ID)[0]
    # Open the SED solution file
    file_name = obj_list[i]
    spec_data = parse_spec_file(file_name)

    # Get SED table, high-z solution is the first one
    sed = spec_data.get('sed')[0]


    # Get Zphot
    zphot = t['Zphot'][row_index][0]
    zinf = t['Zinf'][row_index][0]
    zsup = t['Zsup'][row_index][0]

    if zinf > zphot:
        print(f'Zinf is less than Zphot for ID {ID}.')
        print(f'Zinf: {zinf}, Zphot: {zphot}')
        t['Zinf'][row_index] = math.floor(zinf * 100) / 100
        zinf = math.floor(zinf * 100) / 100

    # Luminosity distance
    DL = cosmo.luminosity_distance(zphot).value * 10 ** 6 # put into pc

    # And the uncertainty on DL
    DL_sup = cosmo.luminosity_distance(zsup).value * 10 ** 6
    DL_inf = cosmo.luminosity_distance(zinf).value * 10 ** 6

    # Extract SED lambda and flux
    wlen = sed['lambda']
    sed = sed['flux']
    wlen = np.array([float(w) for w in wlen])
    sed = np.array([float(s) for s in sed])
    sed = mag_to_flux(sed)

    # Place tophat filter, rest 1500A, width 100A, in the observed frame
    filter_obs = np.zeros(len(wlen))
    for i, lam in enumerate(wlen):
        if (lam > 1450.0*(1+zphot)) and (lam < 1550.0*(1+zphot)):
            filter_obs[i] = 1.0

    # Convolve SED with tophat filter at rest frame 1500. Gets flux density in T1500(1+z).
    convolved_flux = np.sum(sed * filter_obs) / np.sum(filter_obs)

    # Impose minimum error of 5%
    conv_sup = convolved_flux + (convolved_flux* 0.05)
    conv_inf = convolved_flux - (convolved_flux * 0.05)

    # Compute apparent magnitude
    m1500 = -2.5*np.log10(convolved_flux)-48.6
    m1500_sup = -2.5*np.log10(conv_sup)-48.6
    m1500_inf = -2.5*np.log10(conv_inf)-48.6

   # Compute abs mag using distance modulus and the redshift correction
    M1500 = m1500 - 5*np.log10(DL/10) + 2.5*np.log10(1+zphot)

    # Get the supremum and infimum
    K1_sup = 2.5*np.log10(1+zsup)
    K1_inf = 2.5*np.log10(1+zinf)

    # uncertainties
    M1500_sup = m1500_sup - 5*np.log10(DL_sup/10) + K1_sup
    M1500_inf = m1500_inf - 5*np.log10(DL_inf/10) + K1_inf

    plus = M1500 - M1500_sup
    minus = M1500_inf - M1500

    # Add to the table
    t['Muv'][row_index] = M1500
    t['dMuv_inf'][row_index] = minus
    t['dMuv_sup'][row_index] = plus

# Overwrite table
t.write(cat_dir / cat_name, format='fits', overwrite=True)

# Split table into where EW is non-zero
t_lya = t[t['Lyman_alpha_EW'] > 0]

plt.errorbar(
    t['Zphot'], t['Muv'],
    yerr=[t['dMuv_inf'], t['dMuv_sup']],
    xerr=[np.abs(t['Zphot']-t['Zinf']), np.abs(t['Zsup']-t['Zphot'])],
    fmt='o', color='black', markersize=8,
    label='COSMOS',
    alpha=0.5
)

plt.errorbar(
    t_lya['Zphot'], t_lya['Muv'],
    yerr=[t_lya['dMuv_inf'], t_lya['dMuv_sup']],
    xerr=[np.abs(t_lya['Zphot']-t_lya['Zinf']), np.abs(t_lya['Zsup']-t_lya['Zphot'])],
    fmt='D', color='orange', markersize=10,
    label=r'Lyman-$\alpha$ emitters',
)

# Load paper 1 galaxies
cdfs_cands = Table.read('CDFS_cands.fits')
xmm_cands = Table.read('XMM_cands.fits')

xmm_zphot = np.array(xmm_cands['zphot'])
xmm_Muv = np.array(xmm_cands['Muv'])
xmm_z_inf = np.array(xmm_cands['z_inf'])
xmm_z_sup = np.array(xmm_cands['s_sup'])

cdfs_zphot = np.array(cdfs_cands['zphot'])
cdfs_Muv = np.array(cdfs_cands['Muv'])
cdfs_z_inf = np.array(cdfs_cands['z_inf'])
cdfs_z_sup = np.array(cdfs_cands['s_sup'])

xmm_zphot = xmm_zphot.astype(float)
xmm_Muv = xmm_Muv.astype(float)
xmm_z_inf = xmm_z_inf.astype(float)
xmm_z_sup = xmm_z_sup.astype(float)

cdfs_zphot = cdfs_zphot.astype(float)
cdfs_Muv = cdfs_Muv.astype(float)
cdfs_z_inf = cdfs_z_inf.astype(float)
cdfs_z_sup = cdfs_z_sup.astype(float)

# Red circles for XMM
plt.errorbar(
    xmm_zphot, xmm_Muv,
    xerr=[xmm_z_inf, xmm_z_sup],
    fmt='o', color='red', markersize=10,
    label='XMM'
)

# Blue squares for CDFS
plt.errorbar(
    cdfs_zphot, cdfs_Muv,
    xerr=[cdfs_z_inf, cdfs_z_sup],
    fmt='s', color='blue', markersize=10,
    label='CDFS'
)


# Reverse y axis
plt.gca().invert_yaxis()
plt.xlabel(r'$z_{\rm phot}$')
plt.ylabel(r'$M_{\rm UV}$')
plt.legend(loc='upper left')
plt.tight_layout()
plt.xlim(6.0, 7.6)
plt.ylim(-19.8, -25.4)
plot_dir = Path.cwd().parents[1] /'plots' / 'LF'
plt.savefig(plot_dir / 'z_Muv_sample.pdf')
plt.show()


 