"""
The code that computes the LF! Beautiful.

Created: Thursday 14th November 2024.
"""

import numpy as np
from astropy.table import Table, vstack, Column
import matplotlib.pyplot as plt
from pathlib import Path
from scipy.optimize import curve_fit
import emcee
import corner

plt.rcParams.update({'font.size': 15})
plt.rcParams['axes.linewidth'] = 4
plt.rcParams['figure.dpi'] = 100

plot_dir = Path.cwd().parents[1] / 'plots' / 'LF'

#! Switches
run_type = ''
#run_type = 'with_euclid'   

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
    Concatenates multiple sets of (x, y, dy) arrays.
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
cat_name = 'XMM_5sig_HSC_Z_nonDet_HSC_G_nonDet_HSC_R_candidates_2025_05_14.fits'



t = Table.read(cat_dir / cat_name)

# Remove the Lya emitters which have z>7.5 with no emission line.
t = t[t['Vmax'] > 0]
t = t[t['Muv'] < 0]

# plt.scatter(t['Muv'], t['Vmax'], c=t['Zphot'], cmap='viridis', s=10)    
# plt.show()
# exit()

# Restrict to Muv < -20.5
# print('Number of galaxies before Muv cut: ', len(t))
# #t = t[t['Muv'] < -20.95]
# print('Number of galaxies after Muv cut: ', len(t))

# Read in completeness matrix
completeness_dir = Path.cwd().parent / 'injection_recovery'
completeness_name = 'completeness_matrix_2.npy'
completeness_matrix = np.load(completeness_dir / completeness_name)

# Flip in y-axis to get correct Muv ordering
completeness_matrix = np.flip(completeness_matrix, axis=0)

# Bins to snap galaxy zphot and Muv to
Muv_completeness_bins = np.arange(-23, -20., 0.1)
z_completeness_bins = np.arange(6.5, 7.5, 0.05)

#! ################## Muv BINNING ####################

Muv_bins = [-23., -22.5, -22.25, -22., -21.75, -21.5, -21.25, -21., -20.75, -20.5, -20.25, -20.]
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
        
        # Find completeness for this Muv,z
        Muv_bin = np.digitize([Muv], Muv_completeness_bins)[0] - 1
        z_bin = np.digitize([z], z_completeness_bins)[0] - 1

        # Ensure indices stay within bounds
        z_bin = max(0, min(z_bin, completeness_matrix.shape[1] - 1))
        Muv_bin = max(0, min(Muv_bin, completeness_matrix.shape[0] - 1))
    

        completeness = completeness_matrix[Muv_bin, z_bin]

        # Summand
        summand = 1 / (Vmax) #* completeness)
        phi[i] += summand

        # Compute error term
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

#? MY VALUES, VARIOUS CUTS
sample_chi2_10_cut_LF = [1.33157075e-07, 1.18241778e-06, 2.52909446e-06, 9.08111135e-06, 
                         1.65676294e-05, 2.97647415e-05, 3.56246948e-05, 3.88829869e-05,
                         2.74801812e-05, 1.70064150e-05, 1.04995291e-05]
sample_chi2_10_cut_LF_err = [9.41586427e-08, 4.18552119e-07, 6.35377120e-07, 1.24680926e-06,
                             1.79523078e-06, 2.52202978e-06, 2.70418308e-06, 3.84325176e-06,
                             3.62746107e-06, 2.65738666e-06, 3.28051799e-06]


#############! PLOT MY DATA ###############
plt.figure(figsize=(10, 10))

plt.errorbar(bin_centres, phi, yerr=delta_phi, fmt='o', color='red', 
             markersize=14, label=r'XMM-LSS ($\chi^2_{\mathrm{BD}} < \chi^2_{\mathrm{high-}z}$)', 
             elinewidth=3, markeredgecolor='black')

# plot my sample_chi2_10 data with offset in M
sample_chi2_10_cut_LF = np.array(sample_chi2_10_cut_LF)
sample_chi2_10_cut_LF_err = np.array(sample_chi2_10_cut_LF_err)

plt.errorbar(bin_centres + 0.1, sample_chi2_10_cut_LF, yerr=sample_chi2_10_cut_LF_err, fmt='o', color='orange', 
             markersize=12, label=r'XMM-LSS ($\chi^2_\mathrm{BD} < 10$)', elinewidth=3, markeredgecolor='black',
             marker='h')



###############! PLOT LITERATURE DATA ###############
plt.plot(M, bowler15_z6_sch, color='blue', label='Bowler+15, Schechter')
plt.plot(M, bowler15_z6_dpl, color='blue', label='Bowler+15, DPL', linestyle='--')

plt.errorbar(bowler15_M_weighted, bowler15_phi, yerr=bowler15_dphi, fmt='o', color='blue',
             markersize=10, label=r'Bowler+15, $5.7<z<6.3$', elinewidth=2)

plt.tick_params(which='major', length=10, width=3)
plt.tick_params(axis='both', which='minor', length=5, width=2)

plt.xlabel(r'$M_{\mathrm{UV}}$', fontsize=25)
plt.ylabel(r'$\mathrm{Number \ of \ objects \  / \ mag \ / \ Mpc^{3}}$', fontsize=21.5)

# Increase size of tick labels
plt.xticks(fontsize=20)
plt.yticks(fontsize=20)


plt.ylim([1e-8, 4e-4])
plt.xlim(-23.5, -19.5)


plt.legend(loc='lower right', fontsize=17)


plt.tight_layout()  # Leaves extra space at the bottom
plt.yscale('log')
plt.savefig(plot_dir / f'z6_LF_XMM.pdf', dpi=100, bbox_inches='tight')
plt.show()
