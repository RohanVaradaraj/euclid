"""
The code that computes the LF! Beautiful.

Modified from the z=7 version for z=6 to experiment with different zbin cuts, bit cleaner too.
"""

import numpy as np
from astropy.table import Table
import matplotlib.pyplot as plt
from pathlib import Path

plt.rcParams.update({'font.size': 15})
plt.rcParams['axes.linewidth'] = 4
plt.rcParams['figure.dpi'] = 100

plot_dir = Path.cwd().parents[1] / 'plots' / 'LF'

#! Switches

field_name = 'CDFS'

do_completeness = False

#! ############### FUNCTIONS ####################
def dpl(M, Mstar, phiStar, alpha, beta):

    ''' Produces DPL given parameters alpha, beta, normalisation phi*, char. mag. M* and mag array M. '''

    numerator = np.log(10) * phiStar / 2.5

    denomA = 10 ** (0.4 * (alpha + 1) * (M - Mstar))
    denomB = 10 ** (0.4 * (beta + 1)  * (M - Mstar))

    denominator = denomA + denomB

    phi = numerator/denominator

    return phi



def schechter(M, Mstar, phiStar, alpha):

    ''' Produces Schechter function given parameters: normalisation phi*, faint slope alpha, mag array M and char. mag M*. '''

    coeff = np.log(10) / 2.5

    faint = (10 ** (0.4 * (Mstar - M))) ** (alpha+1)


    bright_exponent = -10 ** (0.4 * (Mstar - M))
    bright = np.exp(bright_exponent)

    phi = coeff * phiStar * faint * bright

    return phi



def dpl_fit(M, Mstar, phiStar, alpha, beta):
    """ For fitting the LFs"""
    return np.log(10) * phiStar / (2.5 * (10 ** (0.4 * (alpha + 1) * (M - Mstar)) + 10 ** (0.4 * (beta + 1) * (M - Mstar))))

def schechter_fit(M, phiStar, Mstar, alpha):
    """ For fitting the LFs"""
    coeff = np.log(10) / 2.5
    faint = (10 ** (0.4 * (Mstar - M))) ** (alpha+1)
    bright_exponent = -10 ** (0.4 * (Mstar - M))
    bright = np.exp(bright_exponent)
    phi = coeff * phiStar * faint * bright
    return phi



def concatenate_arrays(*arrays):
    """
    Concatenates multiple sets of (x, y, dy) arrays, to combine LF data I pass in.
    """
    if len(arrays) < 2:
        raise ValueError("At least two sets of (x, y, dy) arrays are required for concatenation.")

    # Unpack and concatenate each component separately
    x_concat = np.concatenate([arr[0] for arr in arrays])
    y_concat = np.concatenate([arr[1] for arr in arrays])
    dy_concat = np.concatenate([arr[2] for arr in arrays])

    return x_concat, y_concat, dy_concat




#? Read catalogue
cat_dir = Path.cwd().parents[1] / 'data' / 'catalogues' / 'candidates'
if field_name == 'XMM':
    cat_name = 'XMM_5sig_HSC_Z_nonDet_HSC_G_nonDet_HSC_R_candidates_2025_05_14.fits'
if field_name == 'COSMOS':
    cat_name = 'COSMOS_5sig_HSC_Z_nonDet_HSC_G_nonDet_HSC_R_candidates_2025_06_06.fits'
if field_name == 'CDFS':
    cat_name = 'CDFS_5sig_HSC_Z_nonDet_HSC_G_nonDet_r_candidates_2026_04_08_with_euclid.fits'


t = Table.read(cat_dir / cat_name)
print(len(t) , 'objects in catalogue')
# print(t.colnames)
# exit()

# Remove errant objects
t = t[t['Vmax'] > 0]
t = t[t['Muv'] < 0]

# In CDFS, remove VOICEless data
if field_name == 'CDFS':
    t = t[t['err_r'] > 0]
    print('In CDFS, ', len(t), ' objects after removing VOICEless')
# Restrict redshift range
# t = t[t['Zphot'] > 5.7]
# t = t[t['Zphot'] < 6.3]
# print('In CDFS, ', len(t), ' objects after restricting to 5.7<z<6.3')

# In CDFS, restrict to z=5.5-6.5
if field_name == 'CDFS':
    t = t[t['Zphot'] > 5.5]
    t = t[t['Zphot'] < 6.5]
    print('In CDFS, ', len(t), ' objects after restricting to 5.5<z<6.5')
    

# Restrict CDFS to > 8sigma in HSC-Z
# mask = (t['flux_HSC-Z'] / t['err_HSC-Z'] > 8)
# t = t[mask]
# print('In CDFS, ', len(t), ' objects after restricting to 8sigma')

# plt.hist(t['Muv'], bins=np.arange(-24, -19, 0.1))
# plt.show()
# exit()

# plt.scatter(t['Muv'], t['Vmax'], c=t['Zphot'], cmap='viridis', s=10)    
# plt.show()
# exit()

# Read in completeness matrix
if do_completeness:
    if field_name == 'COSMOS':
        completeness_dir = Path.cwd().parent / 'injection_recovery'
        completeness_name = f'completeness_matrix_z6_COSMOS.npy'
        completeness_matrix = np.load(completeness_dir / completeness_name)

        # Flip in y-axis to get correct Muv ordering
        completeness_matrix = np.flip(completeness_matrix, axis=0)

    # XMM compeletness matrices split into XMM1,2,3
    if field_name == 'XMM':
        completeness_dir = Path.cwd().parent / 'injection_recovery'
        completeness_name1 = f'completeness_matrix_z6_XMM1.npy'
        completeness_name2 = f'completeness_matrix_z6_XMM2.npy'
        completeness_name3 = f'completeness_matrix_z6_XMM3.npy'

        completeness_matrix1 = np.load(completeness_dir / completeness_name1)
        completeness_matrix2 = np.load(completeness_dir / completeness_name2)
        completeness_matrix3 = np.load(completeness_dir / completeness_name3)

        # Flip in y-axis to get correct Muv ordering
        completeness_matrix1 = np.flip(completeness_matrix1, axis=0)
        completeness_matrix2 = np.flip(completeness_matrix2, axis=0)
        completeness_matrix3 = np.flip(completeness_matrix3, axis=0)

        # Make a dict to store the matrices for each XMM subfield
        completeness_matrices = {
            'XMM1': completeness_matrix1,
            'XMM2': completeness_matrix2,
            'XMM3': completeness_matrix3
        }

    # Bins to snap galaxy zphot and Muv to
    Muv_completeness_bins = np.arange(-23, -20., 0.1)
    z_completeness_bins = np.arange(5.5, 6.5, 0.05)

#! ################## Muv BINNING ####################

Muv_bins = [-23., -22.5, -22.25, -22., -21.75, -21.5, -21.25, -21., -20.75, -20.5, -20.25, -20., -19.75, -19.5]
bin_widths = np.abs(np.diff(Muv_bins))
bin_centres = 0.5 * (np.array(Muv_bins[:-1]) + np.array(Muv_bins[1:]))

if field_name == 'CDFS':
    Muv_bins = [-24, -23.5, -23, -22.5, -22., -21.5, -21.]
    bin_widths = np.abs(np.diff(Muv_bins))
    bin_centres = 0.5 * (np.array(Muv_bins[:-1]) + np.array(Muv_bins[1:]))


print(f'Muv bins: {np.round(Muv_bins, 2)}')
print(f'Bin widths: {np.round(bin_widths, 2)}')
print(f'Bin centres: {np.round(bin_centres, 2)}')

# Split the table into these bins
binned_tables = []
for i in range(len(Muv_bins)-1):
    mask = (t['Muv'] >= Muv_bins[i]) & (t['Muv'] < Muv_bins[i+1])
    binned_tables.append(t[mask])

# Initialise LF sum
phi = np.zeros(len(Muv_bins)-1)
delta_phi = np.zeros(len(Muv_bins)-1)

# Get number of galaxies in each bin
n_gals = np.array([len(sub_table) for sub_table in binned_tables])

#! Loop through the bins
for i, sub_table in enumerate(binned_tables):
    
    # Print number of galaxies in each bin
    print(f'Bin {i} has {len(sub_table)} galaxies')

    Muv_range_here = [Muv_bins[i], Muv_bins[i+1]]
    print(f'Muv range: {Muv_range_here}')

    #! Go through all the objects in the sub-table
    for j, obj in enumerate(sub_table):

        ID = obj['ID']
        Vmax = obj['Vmax']
        z = obj['Zphot']
        Muv = obj['Muv']
        
        if do_completeness:
            # Find completeness for this Muv,z
            Muv_bin = np.digitize(['Muv'], Muv_completeness_bins)[0] - 1
            z_bin = np.digitize([z], z_completeness_bins)[0] - 1

            #? COSMOS case
            if field_name == 'COSMOS':
                # Ensure indices stay within bounds
                z_bin = max(0, min(z_bin, completeness_matrix.shape[1] - 1))
                Muv_bin = max(0, min(Muv_bin, completeness_matrix.shape[0] - 1))
                completeness = completeness_matrix[Muv_bin, z_bin]
                print('Completeness at Muv: ', Muv, ' z: ', z, ' is: ', completeness)

            #? XMM case
            if field_name == 'XMM':
                # Determine which XMM subfield the object is in based on its ID
                video_tile = obj['VISTA_tile_used']
                completeness_matrix = completeness_matrices[video_tile]
                # Ensure indices stay within bounds
                z_bin = max(0, min(z_bin, completeness_matrix.shape[1] - 1))
                Muv_bin = max(0, min(Muv_bin, completeness_matrix.shape[0] - 1))
                completeness = completeness_matrix[Muv_bin, z_bin]
                print('Completeness at Muv: ', Muv, ' z: ', z, ' is: ', completeness, 'in tile: ', video_tile)

        # Summand
        if do_completeness:
            summand = 1 / (Vmax * completeness)
        else:
            summand = 1 / Vmax
        phi[i] += summand

        # Compute error term
        if do_completeness:
            err_summand = 1 / (Vmax * completeness) ** 2
            delta_phi[i] += err_summand
        else:
            err_summand = 1 / Vmax ** 2
            delta_phi[i] += err_summand


##############################! ERROR, INCLUDING COSMIC VARIANCE ##############################
for i, lf in enumerate(phi):

    phi[i] /= bin_widths[i]  
    delta_phi[i] = np.sqrt(delta_phi[i]) / bin_widths[i]

print('LF values: ', phi)
print('LF errors: ', delta_phi)

#############################! LITERATURE COMPARISON ##############################
M = np.arange(-25, -18, 0.1)

bowler15_z6_sch = schechter(M, -20.77, 5.7e-4, -1.88)
bowler15_z6_dpl = dpl(M, -21.2, 1.9e-4, -2.1, -5.1)

bowler15_M_weighted = [-22.52, -22.08, -21.74, -21.49, -21.22]
bowler15_dM = [0.5, 0.5, 0.25, 0.25, 0.25]
bowler15_M = [-22.625, -22.125, -21.75, -21.5, -21.25]
bowler15_phi = [1.16e-6, 5.98e-6, 1.9e-5, 3.92e-5, 9.14e-5]
bowler15_dphi = [0.67e-6, 1.64e-6, 0.41e-5, 0.7e-5, 1.39e-5]

#! JT Schindler QSO LF z=6
schindler22_qso = dpl(M, -26.38, 1.778e-9, -1.70, -3.84)

#? MY VALUES IN XMM, WITH COMPLETENESS CORRECTION APPLIED.
xmm_LF = [1.52232581e-07, 2.27746489e-06, 4.13161771e-06, 1.83714466e-05,
 3.90277980e-05, 7.45208718e-05, 1.02450728e-04, 1.14245945e-04,
 1.52522595e-04, 6.64156917e-05, 3.47985444e-05, 8.52913347e-07,
 0.00000000e+00]
 
xmm_LF_err = [1.07651737e-07, 8.63741352e-07, 1.13299825e-06, 2.64674604e-06,
 4.11314503e-06, 6.27607487e-06, 7.49907283e-06, 9.27034989e-06,
 5.61527856e-05, 9.48876798e-06, 1.04010892e-05, 8.52913347e-07,
 0.00000000e+00]


if field_name == 'COSMOS':
    full_zspan_LF = [3.18107849e-06, 7.72650799e-06, 1.65541048e-05, 3.33368956e-05, 7.87426461e-05, 1.21480570e-04, 
                     1.93707969e-04, 2.78364345e-04, 3.12493734e-04, 4.73882743e-04, 4.13148837e-05]

    full_zspan_LF_err = [1.06290027e-06, 1.69434097e-06, 2.53083832e-06, 3.80276015e-06, 6.08074184e-06, 7.82026245e-06, 
                         1.17891825e-05, 1.83488356e-05, 2.59915168e-05, 2.75335959e-04, 8.89131186e-06]

#? Harikane+23, GOLDRUSH, HSC-only
Muv_h23 = np.array([
    -25.02, -24.52, -24.02, -23.52, -23.12, -22.82,
    -22.52, -22.22, -21.92, -21.62, -21.32, -21.02
])

LF_h23 = np.array([
    1.05e-8, 2.13e-8, 2.77e-8, 8.51e-8, 3.34e-7, 1.24e-6,
    2.67e-6, 4.48e-6, 1.10e-5, 3.69e-5, 7.35e-5, 1.77e-4
])

LFerrlo_h23 = np.array([
    1.05e-8, 2.13e-8, 2.23e-8, 2.85e-8, 0.72e-7, 0.14e-6,
    0.39e-6, 0.53e-6, 0.09e-5, 0.48e-5, 0.85e-5, 0.21e-4
])

LFerrup_h23 = np.array([
    4.11e-8, 4.21e-8, 4.19e-8, 5.38e-8, 0.72e-7, 0.15e-6,
    0.39e-6, 0.53e-6, 0.09e-5, 0.48e-5, 0.85e-5, 0.21e-4
])

#? Bouwens+21, CANDELS
Muv_b21 = np.array([
    -22.52, -22.02, -21.52, -21.02, -20.52,
    -20.02, -19.52, -18.77, -17.77, -16.77
])

# Shift Muv to right by 0.02 mag
Muv_b21 += 0.05

LF_b21 = np.array([
    2.0e-6, 1.4e-5, 5.1e-5, 1.69e-4, 3.17e-4,
    7.24e-4, 1.147e-3, 2.82e-3, 8.36e-3, 1.71e-2
])

LFerr_b21 = np.array([
    2.0e-6, 5.0e-6, 1.1e-5, 2.4e-5, 4.1e-5,
    8.7e-5, 1.57e-4, 4.40e-4, 1.66e-3, 5.26e-3
])

#############! PLOT MY DATA ###############
plt.figure(figsize=(10, 10))

label= field_name + r', ($\chi^2_{\mathrm{BD}} < \chi^2_{\mathrm{high-}z}$)'
label = f'{field_name}' + r', $6.5 < z < 7.5$'
label = f'{field_name}'

plt.errorbar(bin_centres, phi, yerr=delta_phi, fmt='o', color='red', 
             markersize=14, label=label, 
             elinewidth=3, markeredgecolor='black', zorder=10)

if field_name == 'COSMOS':
    # Also plot XMM
    plt.errorbar(bin_centres-0.05, xmm_LF, yerr=xmm_LF_err, fmt='s', color='orange',
                 markersize=14, label='XMM-LSS', elinewidth=3, markeredgecolor='black', zorder=10)

# plot my sample_chi2_10 data with offset in M
# sample_chi2_10_cut_LF = np.array(sample_chi2_10_cut_LF)
# sample_chi2_10_cut_LF_err = np.array(sample_chi2_10_cut_LF_err)

#? PLOTTING DIFFERENT BD CUTS IN XMM
# plt.errorbar(bin_centres + 0.1, sample_chi2_10_cut_LF, yerr=sample_chi2_10_cut_LF_err, fmt='o', color='orange', 
#              markersize=12, label=r'XMM-LSS ($\chi^2_\mathrm{BD} < 10$)', elinewidth=3, markeredgecolor='black',
#              marker='h')

#? PLOTTING DIFFERENT REDSHIFT RANGES IN COSMOS
# plt.errorbar(bin_centres[2:] + 0.07, full_zspan_LF, yerr=full_zspan_LF_err, fmt='o', color='orange',
#              markersize=12, label=r'COSMOS, $6.5 < z < 7.5$', elinewidth=3, markeredgecolor='black',
#              marker='h')


###############! PLOT LITERATURE DATA ###############
plt.plot(M, bowler15_z6_sch, color='blue', label='Bowler+15, Schechter')
plt.plot(M, bowler15_z6_dpl, color='blue', label='Bowler+15, DPL', linestyle='--')

plt.errorbar(bowler15_M_weighted, bowler15_phi, yerr=bowler15_dphi, fmt='o', color='blue',
             markersize=10, label=r'Bowler+15, $5.7<z<6.3$', elinewidth=2)

plt.errorbar(Muv_h23, LF_h23, yerr=[LFerrlo_h23, LFerrup_h23], fmt='h', color='green',
                markersize=10, label='Harikane+23', elinewidth=2, alpha=0.8)

plt.errorbar(Muv_b21, LF_b21, yerr=LFerr_b21, fmt='d', color='purple',
                markersize=10, label='Bouwens+21', elinewidth=2, alpha=0.8)

plt.plot(M, schindler22_qso, color='gray', label='Schindler+22, QSO', linestyle=':')

plt.tick_params(which='major', length=10, width=3)
plt.tick_params(axis='both', which='minor', length=5, width=2)

plt.xlabel(r'$M_{\mathrm{UV}}$', fontsize=25)
plt.ylabel(r'$\mathrm{Number \ of \ objects \  / \ mag \ / \ Mpc^{3}}$', fontsize=21.5)

# Increase size of tick labels
plt.xticks(fontsize=20)
plt.yticks(fontsize=20)


plt.ylim([1e-8, 1e-2])
if field_name == 'XMM':
    plt.ylim([1e-8, 1e-3])
plt.xlim(-24, -19)


plt.legend(loc='upper left', fontsize=17)


plt.tight_layout()  # Leaves extra space at the bottom
plt.yscale('log')
plt.savefig(plot_dir / f'z6_LF_{field_name}.pdf', dpi=100, bbox_inches='tight')
plt.show()
