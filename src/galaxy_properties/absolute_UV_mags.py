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
import matplotlib.gridspec as gridspec

plt.rcParams['axes.linewidth'] = 2.5
plt.rcParams.update({'font.size': 27})
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

import sys
sed_path = Path.cwd().parents[0] / 'sed_fitting'
sys.path.append(str(sed_path))
from sed_fitting_codes import parse_spec_file

#! Field name
field_name = 'COSMOS'

#! Det/non-det filters
# filters = {
#     'HSC-Z_DR3': {'type': 'detection', 'value': 5},
#     'HSC-G_DR3': {'type': 'non-detection', 'value': 2},
#     'HSC-R_DR3': {'type': 'non-detection', 'value': 2},
# }

filters = {
    'Y+J': {'type': 'stacked-detection', 'value': 5},
    'HSC-G_DR3': {'type': 'non-detection', 'value': 2},
    'HSC-R_DR3': {'type': 'non-detection', 'value': 2},
    'HSC-I_DR3': {'type': 'non-detection', 'value': 2},
}

#! Run type
run_type = ''
#run_type = ''

# Switch to stop computing Muv if we've already done it.
compute_Muv = False

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

# Load completeness
completeness_dir = Path.cwd().parent / 'injection_recovery'
completeness_name = 'completeness_matrix_2.npy'
completeness_matrix = np.load(completeness_dir / completeness_name)

#! SED Fitting folder
# Generate name of the directory we want to use to make the catalogue
det_filters = [f for f, t in filters.items() if t['type'] in ['detection', 'stacked-detection']]
det_filter_str = '_'.join(det_filters)
if run_type != '':
    folder = f'det_{det_filter_str}_{run_type}_z7'
else:
    folder = f'det_{det_filter_str}_z7'

print(f'Folder name: {folder}')

# Get the list of objects that made it through the SED fitting
obj_dir = Path.cwd().parents[1] / 'data' / 'sed_fitting' / 'zphot' / field_name / 'best_fits'
obj_list = glob.glob(str(obj_dir / folder / '*.spec'))

# Get the IDs
IDs = [spec_file.split('/')[-1].split('Id')[-1].lstrip('0').split('.spec')[0] for spec_file in obj_list]
IDs = [int(ID) for ID in IDs]

#! Catalogue of above candidate objects

#? COSMOS
if field_name == 'COSMOS':
    if run_type == 'with_euclid':
        cat_name = 'COSMOS_5sig_Y_J_nonDet_HSC_G_nonDet_HSC_R_nonDet_HSC_I_candidates_2025_02_14_with_euclid.fits' # with euclid
        ref_cat = 'COSMOS_5sig_Y_J_nonDet_HSC_G_nonDet_HSC_R_nonDet_HSC_I_candidates_2025_02_14.fits'
    if run_type == '':
        cat_name = 'COSMOS_5sig_Y_J_nonDet_HSC_G_nonDet_HSC_R_nonDet_HSC_I_candidates_2025_02_14.fits' # just vista
        #cat_name = 'COSMOS_5sig_HSC_Z_nonDet_HSC_G_nonDet_HSC_R_candidates_2025_06_06.fits' # z=6 sample

#? XMM
if field_name == 'XMM':
    #cat_name = 'XMM_5sig_HSC_Z_nonDet_HSC_G_nonDet_HSC_R_candidates_2025_05_13.fits' # Initial 5sigma selection, following bowler+15, 838 galaxies!
    cat_name = 'XMM_5sig_HSC_Z_nonDet_HSC_G_nonDet_HSC_R_candidates_2025_05_14.fits' # chi2_star < chi2_gal to remove BDs, instead of chi2_star < 10


# Read in the parent catalogue
cat_dir = Path.cwd().parents[1] / 'data' / 'catalogues' / 'candidates'
t = Table.read(cat_dir / cat_name)

if run_type == 'with_euclid':
    t_ref = Table.read(cat_dir / ref_cat)

#! Add a column for UV mags if it doesnt exist already
if 'Muv' not in t.colnames:
    t['Muv'] = Column(np.zeros(len(t)), name='UV_mag', dtype=float)
    t['dMuv_inf'] = Column(np.zeros(len(t)), name='UV_mag_err', dtype=float)
    t['dMuv_sup'] = Column(np.zeros(len(t)), name='UV_mag_err', dtype=float)

if compute_Muv:
    # Loop through the objects in the table and get its SED properties
    for i, ID in enumerate(IDs):
        print(i)

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

# Save the redshift and Muv values to a numpy file
# output_file = f'z_Muv_sample_{run_type}.npy'
# np.save(output_file, np.array([t['Zphot'], t['Muv']]))

# Load the numpy files
z_Muv_ = np.load('z_Muv_sample_.npy')
z_Muv_with_euclid = np.load('z_Muv_sample_with_euclid.npy')

#Split table into where EW is non-zero
t_lya = t[t['Lyman_alpha_EW'] > 0]
print(len(t_lya))
print(t['Muv'])

# Plot Muv vs redshift
plt.figure(figsize=(12, 8))

zlo = t['Zphot'] - t['Zinf']
zhi = t['Zsup'] - t['Zphot']
zerr = np.array([zlo, zhi])

plt.scatter(t['Zphot'], t['Muv'], s=100, color='dodgerblue', alpha=0.9, edgecolor='none', label='XMM')
plt.errorbar(t['Zphot'], t['Muv'], yerr=[t['dMuv_inf'], t['dMuv_sup']], xerr=zerr, fmt='o', color='dodgerblue', markersize=16,)

# Reverse y axis
plt.gca().invert_yaxis()

plt.xlabel(r'$z_{\rm phot}$')
plt.ylabel(r'$M_{\rm UV}$')

#plt.show()
plt.close()


#! paper 1 galaxies
paper1_dir = Path.cwd().parents[3] / 'HSC_SSP_DR3' / 'codes'
dataCDFS = Table.read(paper1_dir / 'vmax_CDFS_comp.txt', format='ascii.commented_header')
dataXMM = Table.read(paper1_dir / 'vmax_XMM_comp.txt', format='ascii.commented_header')

errorXMM = Table.read(paper1_dir / 'errorsMin_primary_XMM.txt', format='ascii.commented_header')
errorCDFS = Table.read(paper1_dir / 'errorsMin_primary_CDFS.txt', format='ascii.commented_header')


# plt.errorbar(dataXMM['z'], dataXMM['Muv'], xerr=(errorXMM['zinf'], errorXMM['zsup']), yerr=(np.abs(errorXMM['Muv_inf']), errorXMM['Muv_sup']), 
#     color='black', alpha=1, linestyle='none', marker='^', markersize=12, elinewidth=2.5,
#     label='Varadaraj+23, XMM')


# plt.errorbar(dataCDFS['z'], dataCDFS['Muv'], xerr=(errorCDFS['zinf'], errorCDFS['zsup']), yerr=(np.abs(errorCDFS['Muv_inf']), errorCDFS['Muv_sup']), 
#     color='black', marker='s', alpha=1, linestyle='none', markersize=12, elinewidth=2.5,
#     label='Varadaraj+23, CDFS')

# # Imshow completeness as grayscale in the background
# completeness_z = np.arange(6.5, 7.5, 0.05)
# completeness_Muv = np.arange(-23, -20, 0.1)
# # plt.imshow(completeness_matrix, aspect='auto', cmap='Greys', extent=[6.5, 7.5, -23, -20], alpha=0.5, zorder=-1)

# # Flip matrix in Muv direction
# completeness_matrix = np.flip(completeness_matrix, axis=0)

# # And add contour for each redshift where the completeness crosses 0.3
# contour = plt.contour(completeness_z[1:], completeness_Muv[1:], completeness_matrix, levels=[0.3], colors='darkgray', linewidths=5, zorder=5, linestyles='--')
# contour = plt.contour(completeness_z[1:], completeness_Muv[1:], completeness_matrix, levels=[0.4], colors='gray', linewidths=5, zorder=5, linestyles='--')
# contour = plt.contour(completeness_z[1:], completeness_Muv[1:], completeness_matrix, levels=[0.5], colors='dimgray', linewidths=5, zorder=5, linestyles='--')

# plt.text(7.469, -20.4, '30%', fontsize=20, color='darkgray')
# plt.text(7.469, -20.65, '40%', fontsize=20, color='gray')
# plt.text(7.469, -20.9, '50%', fontsize=20, color='dimgray')

# #! UltraVISTA

# t = t[t['Zphot'] > 6.5]

# if run_type == 'with_euclid':

#     # split into things that are and arent in the full vista sample
#     mask_in = np.isin(t['ID'], t_ref['ID'])
#     mask_out = np.isin(t['ID'], t_ref['ID'], invert=True)

#     # Create a new table with only the unmatched rows
#     t_unmatched = t[mask_out]
#     t_matched = t[mask_in]

#     plt.errorbar(
#         t_matched['Zphot'], t_matched['Muv'],
#         yerr=[t_matched['dMuv_inf'], t_matched['dMuv_sup']],
#         xerr=[np.abs(t_matched['Zphot']-t_matched['Zinf']), np.abs(t_matched['Zsup']-t_matched['Zphot'])],
#         fmt='o', color='dodgerblue', markersize=16,
#         alpha=0.6, markeredgecolor='none', elinewidth=2.5,
#     )

#     plt.errorbar(
#         t_unmatched['Zphot'], t_unmatched['Muv'],
#         yerr=[t_unmatched['dMuv_inf'], t_unmatched['dMuv_sup']],
#         xerr=[np.abs(t_unmatched['Zphot']-t_unmatched['Zinf']), np.abs(t_unmatched['Zsup']-t_unmatched['Zphot'])],
#         fmt='s', color='purple', markersize=16,
#         alpha=0.9, markeredgecolor='none', elinewidth=2.5, zorder=3,
#     )



# else:
#     plt.errorbar(
#         t['Zphot'], t['Muv'],
#         yerr=[t['dMuv_inf'], t['dMuv_sup']],
#         xerr=[np.abs(t['Zphot']-t['Zinf']), np.abs(t['Zsup']-t['Zphot'])],
#         fmt='o', color='dodgerblue', markersize=16,
#         alpha=0.6, markeredgecolor='none', elinewidth=2.5,
#     )

# plt.errorbar(
#     t_lya['Zphot'], t_lya['Muv'],
#     yerr=[t_lya['dMuv_inf'], t_lya['dMuv_sup']],
#     xerr=[np.abs(t_lya['Zphot']-t_lya['Zinf']), np.abs(t_lya['Zsup']-t_lya['Zphot'])],
#     fmt='o', color='dodgerblue', markersize=16,
#     alpha=0.6, markeredgecolor='none', elinewidth=2.5,
#     # fmt='D', color='red', markersize=16, alpha=1., markeredgecolor='none', elinewidth=2.5,
#     # zorder=4,
# )

# Bouwens 2021

dataBouwens = Table.read(paper1_dir.parent / 'ref_catalogues' / 'bouwens21_z7.dat', format='ascii.commented_header')
# plt.scatter(dataBouwens['col11'], dataBouwens['Muv'], label='Bouwens+21', s=40, color='gray', alpha=0.4, zorder=-1, edgecolor='none')

# # # Add 50% completeness limit
# # plt.plot(z_complimit, Muv_complimit, color='black', linestyle='--', lw=2.5)

# # Make ticks thicker
# plt.tick_params(which='major', length=10, width=3)

# # Add minor ticks
# plt.tick_params(axis='both', which='minor', length=5, width=2)

# # Reverse y axis
# plt.gca().invert_yaxis()
# plt.xlabel(r'$z_{\rm phot}$')
# plt.ylabel(r'$M_{\rm UV}$')
# plt.legend(loc='upper left')
# plt.tight_layout()
# plt.xlim(6.0, 7.6)
# plt.ylim(-19.8, -24.1)
# plot_dir = Path.cwd().parents[1] /'plots' / 'LF'

# # Check for 'with_euclid' in run_type
# if run_type == 'with_euclid':
#     plt.text(6.01, -19.9, r'UltraVISTA + $Euclid$', size=25)
#     plt.savefig(plot_dir / 'z_Muv_sample_with_euclid.pdf')
# if run_type == '':
#     plt.text(6.01, -19.9, 'UltraVISTA only', size=25)
#     plt.savefig(plot_dir / 'z_Muv_sample.pdf')
# #plt.show()
# exit()
##############! ##############! ##############! ##############! ##############! ##############! ##############! ##############! 

# Define figure with GridSpec layout
fig = plt.figure(figsize=(12, 8))
gs = gridspec.GridSpec(4, 4, wspace=0, hspace=0, width_ratios=[6, 1, 0.2, 0.2], height_ratios=[1, 4, 0.2, 0.2])

# Main scatter plot
ax_main = plt.subplot(gs[1, 0])

ax_main.minorticks_on()

# Top histogram (Zphot distribution)
ax_histx = plt.subplot(gs[0, 0], sharex=ax_main)

# Right histogram (Muv distribution)
ax_histy = plt.subplot(gs[1, 1], sharey=ax_main)

ax_histx.minorticks_on()
ax_histy.minorticks_on()

ax_main.errorbar(dataXMM['z'], dataXMM['Muv'], xerr=(errorXMM['zinf'], errorXMM['zsup']), yerr=(np.abs(errorXMM['Muv_inf']), errorXMM['Muv_sup']), 
    color='gray', alpha=0.7, linestyle='none', marker='^', markersize=12, elinewidth=2.5,
    label='Varadaraj+23, XMM', zorder=-1)


ax_main.errorbar(dataCDFS['z'], dataCDFS['Muv'], xerr=(errorCDFS['zinf'], errorCDFS['zsup']), yerr=(np.abs(errorCDFS['Muv_inf']), errorCDFS['Muv_sup']), 
    color='gray', marker='s', alpha=0.7, linestyle='none', markersize=12, elinewidth=2.5,
    label='Varadaraj+23, CDFS', zorder=-1)

# Imshow completeness as grayscale in the background
completeness_z = np.arange(6.5, 7.5, 0.05)
completeness_Muv = np.arange(-23, -20, 0.1)
# plt.imshow(completeness_matrix, aspect='auto', cmap='Greys', extent=[6.5, 7.5, -23, -20], alpha=0.5, zorder=-1)

# Flip matrix in Muv direction
completeness_matrix = np.flip(completeness_matrix, axis=0)

# And add contour for each redshift where the completeness crosses 0.3
contour = ax_main.contour(completeness_z[1:], completeness_Muv[1:], completeness_matrix, levels=[0.3], colors='darkgray', linewidths=5, zorder=5, linestyles='--')
contour = ax_main.contour(completeness_z[1:], completeness_Muv[1:], completeness_matrix, levels=[0.4], colors='gray', linewidths=5, zorder=5, linestyles='--')
contour = ax_main.contour(completeness_z[1:], completeness_Muv[1:], completeness_matrix, levels=[0.5], colors='dimgray', linewidths=5, zorder=5, linestyles='--')

ax_main.text(7.469, -20.4, '30%', fontsize=22, color='darkgray')
ax_main.text(7.469, -20.65, '40%', fontsize=22, color='gray')
ax_main.text(7.469, -20.9, '50%', fontsize=22, color='dimgray')

# Plot main scatter plot

if run_type == 'with_euclid':

    # split into things that are and arent in the full vista sample
    mask_in = np.isin(t['ID'], t_ref['ID'])
    mask_out = np.isin(t['ID'], t_ref['ID'], invert=True)

    # Create a new table with only the unmatched rows
    t_unmatched = t[mask_out]
    t_matched = t[mask_in]

    ax_main.errorbar(
        t_matched['Zphot'], t_matched['Muv'],
        # yerr=[t_matched['dMuv_inf'], t_matched['dMuv_sup']],
        # xerr=[np.abs(t_matched['Zphot']-t_matched['Zinf']), np.abs(t_matched['Zsup']-t_matched['Zphot'])],
        yerr=0,
        xerr=0,
        fmt='o', color='tab:red', markersize=12,
        alpha=0.8, markeredgecolor='none', elinewidth=2.5,
    )

    ax_main.errorbar(
        t_unmatched['Zphot'], t_unmatched['Muv'],
        # yerr=[t_unmatched['dMuv_inf'], t_unmatched['dMuv_sup']],
        # xerr=[np.abs(t_unmatched['Zphot']-t_unmatched['Zinf']), np.abs(t_unmatched['Zsup']-t_unmatched['Zphot'])],
        yerr=0,
        xerr=0,
        fmt='s', color='tab:red', markersize=14,
        alpha=0.8, markeredgecolor='black', elinewidth=2.5, zorder=3,
        markeredgewidth=2
    )

    ax_main.errorbar(
            t_lya['Zphot'], t_lya['Muv'],
            # yerr=[t_lya['dMuv_inf'], t_lya['dMuv_sup']],
            # xerr=[np.abs(t_lya['Zphot']-t_lya['Zinf']), np.abs(t_lya['Zsup']-t_lya['Zphot'])],
            yerr=0,
            xerr=0,
            fmt='o', color='tab:red', markersize=12,
            alpha=0.8, markeredgecolor='none', elinewidth=2.5,
            # fmt='D', color='red', markersize=16, alpha=1., markeredgecolor='none', elinewidth=2.5,
            # zorder=4,
        )

    # Draw typical error
    mean_error_inf = np.mean(np.abs(t_matched['Zphot']-t_matched['Zinf']))
    mean_error_sup = np.mean(np.abs(t_matched['Zsup']-t_matched['Zphot']))
    mean_error_Muv_inf = np.mean(t_matched['dMuv_inf'])
    mean_error_Muv_sup = np.mean(t_matched['dMuv_sup'])
    ax_main.errorbar(6.23, -20.5, xerr=mean_error_sup, yerr=mean_error_Muv_sup, fmt='o', color='tab:red', 
                     markersize=12, alpha=1, markeredgecolor='none', elinewidth=3, markerfacecolor='none')


else:
    ax_main.errorbar(
        t['Zphot'], t['Muv'],
        # yerr=[t['dMuv_inf'], t['dMuv_sup']],
        # xerr=[np.abs(t['Zphot']-t['Zinf']), np.abs(t['Zsup']-t['Zphot'])],
        yerr=0,
        xerr=0,
        fmt='o', color='dodgerblue', markersize=12,
        alpha=0.8, markeredgecolor='none', elinewidth=2.5#, markerfacecolor='none',
    )

    ax_main.errorbar(
        t_lya['Zphot'], t_lya['Muv'],
        # yerr=[t_lya['dMuv_inf'], t_lya['dMuv_sup']],
        # xerr=[np.abs(t_lya['Zphot']-t_lya['Zinf']), np.abs(t_lya['Zsup']-t_lya['Zphot'])],
        yerr=0,
        xerr=0,
        fmt='o', color='dodgerblue', markersize=12,
        alpha=0.8, markeredgecolor='none', elinewidth=2.5#, markerfacecolor='none',
        # fmt='D', color='red', markersize=16, alpha=1., markeredgecolor='none', elinewidth=2.5,
        # zorder=4,
    )


    # Draw typical error
    mean_error_inf = np.mean(np.abs(t['Zphot']-t['Zinf']))
    mean_error_sup = np.mean(np.abs(t['Zsup']-t['Zphot']))
    mean_error_Muv_inf = np.mean(t['dMuv_inf'])
    mean_error_Muv_sup = np.mean(t['dMuv_sup'])
    ax_main.errorbar(6.23, -20.5, xerr=mean_error_sup, yerr=mean_error_Muv_sup, fmt='o', color='dodgerblue', 
                     markersize=12, alpha=1, markeredgecolor='none', elinewidth=3, markerfacecolor='none')


# bouwens+21
dataBouwens = Table.read(paper1_dir.parent / 'ref_catalogues' / 'bouwens21_z7.dat', format='ascii.commented_header')
ax_main.scatter(dataBouwens['col11'], dataBouwens['Muv'], label='Bouwens+21', s=70, color='gray', alpha=0.4, zorder=-1, edgecolor='none')


# Plot histograms
density = True
if run_type == 'with_euclid':
    comparison_z = z_Muv_[0]
    comparison_M = z_Muv_[1]
    ax_histx.hist(comparison_z, bins=np.arange(6.45, 7.55, 0.05), histtype='step', color='dodgerblue', density=density, linewidth=2, alpha=0.7)
    ax_histy.hist(comparison_M, bins=np.arange(-22.6, -20.2, 0.1), orientation='horizontal', histtype='step', color='dodgerblue', density=density, linewidth=2, alpha=0.7)
if run_type == '':
    comparison_z = z_Muv_with_euclid[0]
    comparison_M = z_Muv_with_euclid[1]
    ax_histx.hist(comparison_z, bins=np.arange(6.45, 7.55, 0.05), histtype='step', color='tab:red', density=density, linewidth=2, alpha=0.7)
    ax_histy.hist(comparison_M, bins=np.arange(-22.6, -20.2, 0.1), orientation='horizontal', histtype='step', color='tab:red', density=density, linewidth=2, alpha=0.7)


if run_type=='with_euclid':
    color='tab:red'
else:
    color='dodgerblue'
ax_histx.hist(t['Zphot'], bins=np.arange(6.45, 7.55, 0.05), histtype='step', color=color, density=density, linewidth=3)
ax_histy.hist(t['Muv'], bins=np.arange(-22.6, -20.2, 0.1), histtype='step', color=color, orientation='horizontal', density=density, linewidth=3)



# Remove tick labels for clean aesthetics
ax_histx.set_yticklabels([])
ax_histx.set_yticklabels([])
ax_histy.set_xticklabels([])
ax_histy.set_xticklabels([])

ax_histx.tick_params(axis='x', which='both', labelbottom=False)
ax_histy.tick_params(axis='y', which='both', labelleft=False)

# Remove the y ticks of the top and x ticks of the right
ax_histx.tick_params(axis='y', which='both', left=False, right=False)
ax_histy.tick_params(axis='x', which='both', bottom=False, top=False)

# Make tick size and width same as main plot
ax_histx.tick_params(which='major', length=10, width=3)
ax_histy.tick_params(which='major', length=10, width=3)


#ax_main.set_xticklabels([6.0, 6.2, 6.4, 6.6, 6.8, 7.0, 7.2, 7.4, 7.6])

ax_main.tick_params(which='major', length=10, width=3)

# # Add minor ticks
ax_main.tick_params(axis='both', which='minor', length=5, width=2)

# Axis labels
ax_main.set_xlabel(r'$z_{\rm phot}$', fontsize=30)
ax_main.set_ylabel(r'$M_{\rm UV}$', fontsize=30)
plt.xticks(fontsize=25)
plt.yticks(fontsize=25)

# Reverse y-axis for Muv (brighter = lower values)
ax_main.invert_yaxis()

# Set limits
ax_main.set_xlim(6.0, 7.65)
ax_main.set_ylim(-19.8, -24.2)

ax_main.legend(loc='upper left', fontsize=22, ncol=1, frameon=False)

# Add text label
if run_type == 'with_euclid':
    ax_main.text(6.05, -19.95, r'UltraVISTA + $Euclid$', size=29)
else:
    ax_main.text(6.05, -19.95, 'UltraVISTA only', size=29)

plt.tight_layout()

# Save figure
plot_dir = Path.cwd().parents[1] / 'plots' / 'LF'
if run_type == 'with_euclid':
    plt.savefig(plot_dir / f'z_Muv_sample_{run_type}_marg.pdf', bbox_inches='tight')
else:
    plt.savefig(plot_dir / f'z_Muv_sample_marg.pdf', bbox_inches='tight')
#plt.savefig(plot_dir / f'z_Muv_sample_{run_type}_marg.pdf', bbox_inches='tight')

plt.show()



 