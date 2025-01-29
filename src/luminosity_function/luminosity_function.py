"""
The code that computes the LF! Beautiful.

Created: Thursday 14th November 2024.
"""

import numpy as np
from astropy.table import Table, vstack
import matplotlib.pyplot as plt
from pathlib import Path
from scipy.optimize import curve_fit

plt.rcParams.update({'font.size': 15})
plt.rcParams['axes.linewidth'] = 4
plt.rcParams['figure.dpi'] = 100


def dpl(phiStar, alpha, beta, M, Mstar):

    ''' Produces DPL given parameters alpha, beta, normalisation phi*, char. mag. M* and mag array M. '''

    numerator = np.log(10) * phiStar / 2.5

    denomA = 10 ** (0.4 * (alpha + 1) * (M - Mstar))
    denomB = 10 ** (0.4 * (beta + 1)  * (M - Mstar))

    denominator = denomA + denomB

    phi = numerator/denominator

    return phi




def schechter(phiStar, alpha, M, Mstar):

    ''' Produces Schechter function given parameters: normalisation phi*, faint slope alpha, mag array M and char. mag M*. '''

    coeff = np.log(10) / 2.5

    faint = (10 ** (0.4 * (Mstar - M))) ** (alpha+1)


    bright_exponent = -10 ** (0.4 * (Mstar - M))
    bright = np.exp(bright_exponent)

    phi = coeff * phiStar * faint * bright

    return phi



# Read in the catalogue
cat_dir = Path.cwd().parents[1] / 'data' / 'catalogues' / 'candidates'
cat_name = 'COSMOS_5sig_Y_J_nonDet_HSC_G_nonDet_HSC_R_nonDet_HSC_I_candidates_2024_11_20.fits'
t = Table.read(cat_dir / cat_name)

# Restrict to 6.5 < z < 7.5
t = t[(t['Zphot'] > 6.49) & (t['Zphot'] < 7.51)]

# Restrict to Muv < -20.5
print('Number of galaxies before Muv cut: ', len(t))
#t = t[t['Muv'] < -20.6]
print('Number of galaxies after Muv cut: ', len(t))

# Read in completeness matrix
completeness_dir = Path.cwd().parent / 'injection_recovery'
completeness_name = 'completeness_matrix.npy'
completeness_matrix = np.load(completeness_dir / completeness_name)

# Flip in the y/Muv direction
#completeness_matrix = np.flip(completeness_matrix, axis=0)


Muv_bins = np.arange(-23, -20., 0.1)
z_bins = np.arange(6.5, 7.5, 0.05)

# Find the minimum and maximum Muv
Muv_min = np.min(t['Muv'])
Muv_max = np.max(t['Muv'])
print(f'Minimum Muv: {Muv_min}, Maximum Muv: {Muv_max}')

# Make bins of width 0.5 within these limits
bin_width = 0.5
Muv_bins = np.arange(Muv_min-0.25, Muv_max+0.25, bin_width)
print(Muv_bins)

# Split the table into these bins
binned_tables = []
for i in range(len(Muv_bins)-1):
    mask = (t['Muv'] >= Muv_bins[i]) & (t['Muv'] < Muv_bins[i+1])
    binned_tables.append(t[mask])

# Initialise LF sum
LF_sum = np.zeros(len(Muv_bins)-1)
LF_error = np.zeros(len(Muv_bins)-1)

#! Loop through the bins
for i, sub_table in enumerate(binned_tables):
    
    # Print number of galaxies in each bin
    print(f'Bin {i} has {len(sub_table)} galaxies')
    #plt.hist(sub_table['Muv'], bins=np.arange(Muv_min, Muv_max, 0.1), alpha=0.5)

    #! Go through all the objects in the sub-table
    for j, obj in enumerate(sub_table):

        Vmax = obj['Vmax']

        # Get redshift and Muv for this object
        z = obj['Zphot']
        Muv = obj['Muv']

        #print(f'z: {z}, Muv: {Muv}, Vmax: {Vmax}')
        
        # Find completeness for this Muv,z
        Muv_bin = np.digitize([Muv], Muv_bins)[0] - 1
        z_bin = np.digitize([z], z_bins)[0] - 1

        # if z_bin == 19:
        #     z_bin = 18

        #print(f'Muv bin: {Muv_bin}, z bin: {z_bin}')

        completeness = completeness_matrix[z_bin, Muv_bin]

        #print(f'Completeness: {completeness}')

        # Summand
        #summand = 1 / (Vmax * completeness)
        summand = 1 / (Vmax) # * completeness)
        LF_sum[i] += summand

        # Compute error term
        err_summand = 1 / Vmax ** 2
        LF_error[i] += err_summand

# Histogram of Muv split by bins
#plt.show()

#! Error
for i, lf in enumerate(LF_sum):

    # Account for bin width
	lf = lf / bin_width
	LF_error[i] =  np.sqrt(LF_error[i]) / bin_width

print(LF_sum)
print(LF_error)

#! Existing data

# Bouwens+21 LF points
b21x = [-22.19, -21.69, -21.19, -20.69, -20.19, -19.69, -19.19, -18.69, -17.94, -16.94]
b21y = [1e-6, 4.1e-5, 4.7e-5, 1.98e-4, 2.83e-4, 5.89e-4, 1.172e-3, 1.433e-3, 5.760e-3, 8.320e-3]
b21dy = [2e-6, 1.1e-5, 1.5e-5, 3.6e-5, 6.6e-5, 1.26e-4, 3.36e-4, 4.19e-4, 1.440e-3, 2.9e-3]

# McLure+13
m13x = [-21, -20.5, -20, -19.5, -19, -18.5, -18, -17.5, -17]
m13y = [0.00003, 0.00012, 0.00033, 0.00075, 0.0011, 0.0021, 0.0042, 0.0079, 0.011]
m13dy = [0.00001, 0.00002, 0.00005, 0.00009, 0.0002, 0.0006, 0.0009, 0.0019, 0.0025]

# Varadaraj+23 LF points
bins_mag = np.array([-22.175, -22.925, -23.675, -24.425])
phi_values = np.array([2.70e-6, 2.81e-7, 2.37e-8, np.nan])  
phi_errors = np.array([0.66e-6, 1.54e-7, 2.50e-8, np.nan])  

# Harikane+24 points
h24x = [-23.2, -22.7]
h24y = [1.6e-7, 4.8e-7]
h24dy_up = [3.7e-7, 4.6e-7]
h24dy_lo = [1.3e-7, 2.6e-7]

#! Fit the LF!

def dpl_fit(M, phiStar, alpha, beta, Mstar):
    """ For fitting the LFs"""
    return np.log(10) * phiStar / (2.5 * (10 ** (0.4 * (alpha + 1) * (M - Mstar)) + 10 ** (0.4 * (beta + 1) * (M - Mstar))))


# Combine my data with McLure+13
Muv_combined = np.concatenate((Muv_bins[:-1], m13x))  # Muv
LF_combined = np.array(LF_sum.tolist() + m13y)  # LF
LF_error_combined = np.array(LF_error.tolist() + m13dy) # LF error

# Combine wih Bouwens+21
# Muv_combined = np.concatenate((Muv_bins[:-1], b21x))
# LF_combined = np.concatenate((LF_sum, b21y))
# LF_error_combined = np.concatenate((LF_error, b21dy))

print('Combined Muv:', Muv_combined)
print('Combined LF:', LF_combined)

# Initial parameter guess (phi*, alpha, beta, M*)
p0 = [3.6e-4, -2.1, -4.2, -20.3]

# Fit the function
popt, pcov = curve_fit(dpl_fit, Muv_combined, LF_combined, sigma=LF_error_combined, p0=p0, maxfev=10000)

# Extract best-fit parameters
phiStar_fit, alpha_fit, beta_fit, Mstar_fit = popt

print(f"Best-fit parameters:\nPhi*: {phiStar_fit:.3e}, Alpha: {alpha_fit:.2f}, Beta: {beta_fit:.2f}, M*: {Mstar_fit:.2f}")


# mag range to plot over
M = np.arange(-25, -19, 0.1)

# Bowler+17
z7_gal = dpl(2.3*10**(-4.), -2.19, -4.60, M, -20.60)

# Harikane+24
z7_harikane = dpl(10**(-3.74), -2.08, -4.81, M, -21.01)

# Harikane GOLDRUSH
DPLy_h22 = dpl(10**(-3.05), -1.89, -3.81, M, -20.12) + dpl(10**(-8.49), -1.23, -2.73, M, -24.9)



# Varadaraj+23 LF points
bins_mag = np.array([-22.175, -22.925, -23.675, -24.425])
phi_values = np.array([2.70e-6, 2.81e-7, 2.37e-8, np.nan])  
phi_errors = np.array([0.66e-6, 1.54e-7, 2.50e-8, np.nan])  

# Bowler+17 LF points
b17x = [-22.86, -22.40, -21.85]
b17y = [3.59e-7, 1.16e-6, 2.75e-6]
b17dy = [2.54e-7, 0.58e-6, 1.04e-6]

b21dx = np.ones(len(b21x)) * 0.25

# UVISTA DR6 LF points
muv = [-22.79501799, -22.29501799, -21.79501799, -21.29501799, -20.79501799, -20.29501799]
phi = [4.64101649e-07, 2.33761186e-06, 4.46459736e-06, 1.09373428e-05, 6.40936287e-06]
phi_err = [4.15105134e-07, 9.35645234e-07, 1.28888147e-06, 2.14176670e-06, 1.75013185e-06]

# Print ratio between my LF points and the Harikane+24 LF function at that magnitude
ratios = []
for i, M_ in enumerate(Muv_bins[:-1]):

    ratio = LF_sum[i] / dpl(10**(-3.74), -2.08, -4.81, M_, -21.01)
    ratios.append(round(ratio, 2))
print('###### Ratios #######')
#print(ratios)

#! Plot!
plt.figure(figsize=(10, 10))
plt.plot(M, z7_gal, color='black', linewidth=3, label='Bowler+17')
plt.plot(M, z7_harikane, color='blue', linewidth=3, label='Harikane+24')
plt.plot(M, DPLy_h22, color='gray', linewidth=3, label='Harikane+22', alpha=0.5)

#! My best fit DPL function
# M_fit = np.linspace(-25, -18, 100)
# LF_fit = dpl_fit(M_fit, *popt)
# plt.plot(M_fit, LF_fit, color='red', linewidth=3, label="Best-fit DPL")

# McLure+13
plt.errorbar(m13x, m13y, yerr=m13dy, color='purple', label='McLure+13', marker='D', markersize=10, alpha=0.8, linestyle='none')

# Bowuens et al. 2021
plt.errorbar(b21x, b21y, color='orange', yerr=b21dy, label='Bouwens+21', marker='o', markerfacecolor='none', markersize=10, alpha=0.8, linestyle='none')

# Bowler+17
plt.errorbar(b17x, b17y, color='green', yerr=b17dy, label='Bowler+17', marker='s', markersize=10, alpha=0.8, linestyle='none', zorder=3)

# Varadaraj+23
plt.errorbar(bins_mag[:3], phi_values[:3], yerr=phi_errors[:3], fmt='o', color='black', 
             ecolor='black', elinewidth=2, label='Varadaraj+23', markersize=13, markeredgecolor='black', zorder=2)

# Harikane+24 LF points
plt.errorbar(h24x, h24y, yerr=[h24dy_lo, h24dy_up], fmt='D', color='blue', alpha=0.8, label='Harikane+24', markersize=10, zorder=2)


# Plot the LF new points
plt.errorbar(Muv_bins[:-1], LF_sum, yerr=LF_error, fmt='o', color='red', 
             ecolor='red', elinewidth=4, label=r'UltraVISTA$+Euclid$', markersize=17, markeredgecolor='black')

# Plot the UVISTA DR6 LF points
# plt.errorbar(muv[:-1], phi, yerr=phi_err, fmt='o', color='red', 
#              ecolor='red', elinewidth=4, label='UltraVISTA DR6', markersize=17, markeredgecolor='black')

plt.tick_params(which='major', length=10, width=3)
plt.tick_params(axis='both', which='minor', length=5, width=2)

plt.xlabel(r'$M_{\mathrm{UV}}$', fontsize=20)
plt.ylabel(r'$\mathrm{Number \ of \ objects \  / \ mag \ / \ Mpc^{3}}$', fontsize=20)

plt.xlim(-25, -19)
plt.ylim(1e-10, 1e-2)

plt.ylim([10**(-10), 0.03])
plt.legend(loc='upper left', fontsize=15)
plt.tight_layout()
plt.yscale('log')
plot_dir = Path.cwd().parents[1] / 'plots' / 'LF'
plt.savefig(plot_dir / 'LF_UVISTA.pdf')
#plt.savefig(plot_dir / 'LF_UVISTA_complete.pdf')
plt.show()


