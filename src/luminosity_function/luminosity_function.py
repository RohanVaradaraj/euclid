"""
The code that computes the LF! Beautiful.

Created: Thursday 14th November 2024.
"""

import numpy as np
from astropy.table import Table, vstack, Column
import matplotlib.pyplot as plt
from pathlib import Path
from scipy.optimize import curve_fit


plt.rcParams.update({'font.size': 15})
plt.rcParams['axes.linewidth'] = 4
plt.rcParams['figure.dpi'] = 100


#! ############### FUNCTIONS ####################
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



def dpl_fit(M, phiStar, alpha, beta, Mstar):
    """ For fitting the LFs"""
    return np.log(10) * phiStar / (2.5 * (10 ** (0.4 * (alpha + 1) * (M - Mstar)) + 10 ** (0.4 * (beta + 1) * (M - Mstar))))



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



def adaptive_muv_binning(muv_values, min_bin_size=10):
    """
    Bins Muv values starting from the brightest end (most negative) such that each bin contains at least `min_bin_size` galaxies.

    Parameters:
        muv_values (array-like): The Muv values of the galaxies.
        min_bin_size (int): Minimum number of galaxies per bin.

    Returns:
        bin_centers (numpy array): The central value of each bin.
        bin_widths (numpy array): The width of each bin.
        bin_edges (numpy array): The bin edges.
    """
    # Sort Muv values in descending order (brightest first)
    muv_sorted = np.sort(muv_values)
    
    bin_edges = [muv_sorted[0]]  # Start from the brightest galaxy
    i = 0

    while i < len(muv_sorted):
        # Ensure at least `min_bin_size` galaxies per bin
        next_i = min(i + min_bin_size, len(muv_sorted))
        
        # Define bin edge
        bin_edges.append(muv_sorted[next_i - 1])  # Include the faintest galaxy in this bin
        i = next_i  # Move to the next bin

    bin_edges = np.array(bin_edges)
    bin_centers = 0.5 * (bin_edges[:-1] + bin_edges[1:])  # Midpoints
    bin_widths = np.abs(np.diff(bin_edges))  # Bin width

    return bin_centers, bin_widths, bin_edges



#!################## READ IN DATA ####################

#? Read catalogue
cat_dir = Path.cwd().parents[1] / 'data' / 'catalogues' / 'candidates'
cat_name = 'COSMOS_5sig_Y_J_nonDet_HSC_G_nonDet_HSC_R_nonDet_HSC_I_candidates_2025_01_31.fits'
#cat_name = 'COSMOS_5sig_Y_J_nonDet_HSC_G_nonDet_HSC_R_nonDet_HSC_I_candidates_2024_11_20.fits' # inclusive
t = Table.read(cat_dir / cat_name)


# Remove the Lya emitters which have z>7.5 with no emission line.
t = t[t['Vmax'] > 0]
t = t[t['Muv'] < 0]

# Remove object with ID 1109577
# t = t[t['ID'] != 1109577]

# plt.scatter(t['Muv'], t['Vmax'], c=t['Zphot'], cmap='viridis', s=10)    
# plt.show()
# exit()

# Restrict to Muv < -20.5
print('Number of galaxies before Muv cut: ', len(t))
t = t[t['Muv'] < -20.95]
print('Number of galaxies after Muv cut: ', len(t))

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
# Find the minimum and maximum Muv of the sample
Muv_min = np.min(t['Muv'])
Muv_max = np.max(t['Muv'])
print(f'Sample Minimum Muv: {Muv_min}, Maximum Muv: {Muv_max}')

#? Original binning
bin_width = 0.4
Muv_bins = np.arange(Muv_min-0.25, Muv_max, bin_width)
#print(Muv_bins)

# bin_edges = np.array([-22.5, -22.1, -21.7, -21.3, -20.9, -20.8, -20.7, -20.6, -20.5])
# Muv_bins = 0.5 * (bin_edges[:-1] + bin_edges[1:])
# bin_widths = np.abs(np.diff(bin_edges))

#? Experimenting with different binning schemes
#Muv_bins = np.array([-20.8, -21.2, -21.6, -22.0, -22.4, -22.8])[::-1]
# Muv_bins = np.array([-22.55, -22.15, -21.75, -21.35, -20.95])

#? Final binning scheme for UltraVISTA
#Muv_bins = np.array([-22.6, -22.2, -21.8, -21.4, -20.8]) # This one seems to do well

bin_widths = np.abs(np.diff(Muv_bins))
bin_centres = 0.5 * (np.array(Muv_bins[:-1]) + np.array(Muv_bins[1:]))

# Muv_bins = np.arange(Muv_min-0.1, Muv_max, 0.5)
# bin_widths = np.abs(np.diff(Muv_bins))
# bin_centres = 0.5 * (Muv_bins[:-1] + Muv_bins[1:])

# bin_edges = [-22.65, -22.25, -21.85, -21.45, -21.25, -21.00]
# bin_widths = [0.4, 0.4, 0.4, 0.4, 0.2, 0.25]
# bin_centres = 0.5 * (np.array(bin_edges[:-1]) + np.array(bin_edges[1:]))
# Muv_bins = bin_edges

#Muv_bins, bin_widths, bin_edges = adaptive_muv_binning(t['Muv'], min_bin_size=30)
print(f'Muv bins: {np.round(Muv_bins, 2)}')
print(f'Bin widths: {np.round(bin_widths, 2)}')
print(f'Bin centres: {np.round(bin_centres, 2)}')
#exit()

# plt.hist(t['Muv'], bins=np.arange(Muv_min, Muv_max, 0.01), alpha=0.5)
# for Muv in Muv_bins:
#     plt.axvline(Muv, color='black', linestyle='--')

# plt.show()
# exit()


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
    # plt.hist(sub_table['Muv'], bins=np.arange(Muv_min, Muv_max, 0.1), alpha=0.5)
    # plt.axvline(Muv_bins[i], color='black', linestyle='--')
    # plt.axvline(Muv_bins[i+1], color='black', linestyle='--')
    # plt.show()

    Muv_range_here = [Muv_bins[i], Muv_bins[i+1]]
    completeness_here = completeness_matrix[:, np.digitize([Muv_range_here], Muv_completeness_bins)[0] - 1]
    print(completeness_here.shape)

    # Get mean, median completeness
    print(f'Mean completeness: {np.mean(completeness_here):.4f}, Median completeness: {np.median(completeness_here):.4f}')

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
        
        #print(f'Object {ID} in Muv bin {i+1} at z={z:.4f}, Muv={Muv:.4f} with Vmax {Vmax:.4f} has completeness {completeness_matrix[Muv_bin, z_bin]:.4f}')

        completeness = completeness_matrix[Muv_bin, z_bin]

        # Summand
        summand = 1 / (Vmax * completeness)
        phi[i] += summand

        # Compute error term
        err_summand = 1 / Vmax ** 2
        delta_phi[i] += err_summand

# Histogram of Muv split by bins
#plt.show()

print('LF values:', phi)

##############################! ERROR, INCLUDING COSMIC VARIANCE ##############################
for i, lf in enumerate(phi):

    # Account for bin width
	lf = lf / bin_widths[i]
	delta_phi[i] =  np.sqrt(delta_phi[i]) / bin_widths[i]


"""
 ------------------------------------------------- 
| Total fractional error on number counts = 0.168 | 
 ------------------------------------------------- 
 
Error Budget: 
Poisson uncertainty (relative): 0.141 
Cosmic variance (relative): 0.091 

Observed counts: 50 +/- 9 
"""

# For each bin at its median completeness
cv_per_bin = np.array([0.118, 0.110, 0.104, 0.096])

# Add fractional cosmic variance error in quadrature
#cosmic_variance_frac_err = 0.091

# Add this much percentage error to the LF in quadrature
delta_phi = np.sqrt(delta_phi**2 + (cv_per_bin * phi)**2)

print('LF errors:', delta_phi)


#! ############################## Existing data points ##############################

# Bouwens+21 LF points
b21x = [-22.19, -21.69, -21.19, -20.69, -20.19, -19.69, -19.19, -18.69, -17.94, -16.94]
b21y = [1e-6, 4.1e-5, 4.7e-5, 1.98e-4, 2.83e-4, 5.89e-4, 1.172e-3, 1.433e-3, 5.760e-3, 8.320e-3]
b21dy = [2e-6, 1.1e-5, 1.5e-5, 3.6e-5, 6.6e-5, 1.26e-4, 3.36e-4, 4.19e-4, 1.440e-3, 2.9e-3]
b21x, b21y, b21dy = np.array(b21x), np.array(b21y), np.array(b21dy)


# McLure+13
m13x = [-21, -20.5, -20, -19.5, -19, -18.5, -18, -17.5, -17]
m13y = [0.00003, 0.00012, 0.00033, 0.00075, 0.0011, 0.0021, 0.0042, 0.0079, 0.011]
m13dy = [0.00001, 0.00002, 0.00005, 0.00009, 0.0002, 0.0006, 0.0009, 0.0019, 0.0025]
m13x, m13y, m13dy = np.array(m13x), np.array(m13y), np.array(m13dy)

# Finkelstein+15
f15x = [-22.0, -21.5, -21.0, -20.5, -20.0, -19.5, -19.0, -18.5, -18.0]
f15y = [0.0046e-3, 0.0187e-3, 0.0690e-3, 0.1301e-3, 0.2742e-3, 0.3848e-3, 0.5699e-3, 2.5650e-3, 3.0780e-3]
f15y_up = [0.0049e-3, 0.0085e-3, 0.0156e-3, 0.0239e-3, 0.0379e-3, 0.0633e-3, 0.2229e-3, 0.8735e-3, 1.0837e-3]
f15y_lo = [0.0028e-3, 0.0067e-3, 0.0144e-3, 0.0200e-3, 0.0329e-3, 0.0586e-3, 0.1817e-3, 0.7161e-3, 0.8845e-3]
f15x, f15y, f15y_up, f15y_lo = np.array(f15x), np.array(f15y), np.array(f15y_up), np.array(f15y_lo)

# Harikane+22 LF points
yh22x = [-24.92, -24.42, -23.92, -23.42, -22.92, -22.42, -21.92]
#yh22y = [5e-10, 1.31e-8, 4.39e-8, 1.83e-7, 1.06e-6, 2.75e-6] # 'Galaxy LF'
yh22y = [8.89e-9, 2.41e-8, 9.02e-8, 1.62e-7, 4.63e-7, 1.95e-6, 3.47e-6] # Dropout LF
yh22dyU = [15.29e-9, 1.75e-8, 2.13e-8, 0.3e-7, 7.53e-7, 1.13e-6, 1.7e-6]
yh22dyL = [7.41e-9, 0.98e-8, 1.66e-8, 0.23e-7, 2.71e-7, 0.67e-6, 1.03e-6]
yh22x, yh22y, yh22dyU, yh22dyL = np.array(yh22x), np.array(yh22y), np.array(yh22dyU), np.array(yh22dyL)

# Varadaraj+23 LF points
v23x = np.array([-22.175, -22.925, -23.675])#, -24.425])
v23y = np.array([2.70e-6, 2.81e-7, 2.37e-8])#, np.nan])  
v23dy = np.array([0.66e-6, 1.54e-7, 2.50e-8])#, np.nan])  
v23x, v23y, v23dy = np.array(v23x), np.array(v23y), np.array(v23dy)

# Harikane+24 points
h24x = [-23.2, -22.7]
h24y = [1.6e-7, 4.8e-7]
h24dy_up = [3.7e-7, 4.6e-7]
h24dy_lo = [1.3e-7, 2.6e-7]
h24x, h24y, h24dy_up, h24dy_lo = np.array(h24x), np.array(h24y), np.array(h24dy_up), np.array(h24dy_lo)

# Bowler+17 LF points
b17x = [-22.86, -22.40, -21.85]
b17y = [3.59e-7, 1.16e-6, 2.75e-6]
b17dy = [2.54e-7, 0.58e-6, 1.04e-6]
b17x, b17y, b17dy = np.array(b17x), np.array(b17y), np.array(b17dy)

# # UVISTA DR6 LF points, without completeness, on full inclusive sample of 243 objects
# muv = [-22.79501799, -22.29501799, -21.79501799, -21.29501799, -20.79501799, -20.29501799]
# phi = [4.64101649e-07, 2.33761186e-06, 4.46459736e-06, 1.09373428e-05, 6.40936287e-06]
# phi_err = [4.15105134e-07, 9.35645234e-07, 1.28888147e-06, 2.14176670e-06, 1.75013185e-06]

####################! Fit the LF! ####################

#? Combine my points with other points:

Muv_combined, LF_combined, delta_phi_combined = concatenate_arrays(
    (bin_centres, phi, delta_phi), # My new UltraVISTA points
    #(Muv_bins[:-2], phi[:-1], delta_phi[:-1]), # My new UltraVISTA points, minus last bin affected by BD cut removing galaxies
    (v23x[:-1], v23y[:-1], v23dy[:-1]), # Varadaraj+23, minus CDFS BD
    #(v23x, v23y, v23dy), # Varadaraj+23
    (f15x, f15y, f15y_up), # Finkelstein+15
    #(m13x, m13y, m13dy), # McLure+13
    #(b17x, b17y, b17dy), # Bowler+17
)

#? Curve fitting
# Initial parameter guess (phi*, alpha, beta, M*)
p0 = [3.6e-4, -2.1, -4.8, -20.3]

# Fit the function
popt, pcov = curve_fit(dpl_fit, Muv_combined, LF_combined, sigma=delta_phi_combined, p0=p0, maxfev=10000)

# Extract best-fit parameters
phiStar_fit, alpha_fit, beta_fit, Mstar_fit = popt

print(f"Best-fit parameters:\nphi*: {phiStar_fit:.3e}, Alpha: {alpha_fit:.2f}, Beta: {beta_fit:.2f}, M*: {Mstar_fit:.2f}")

# Get errors on the parameters
perr = np.sqrt(np.diag(pcov))
print(f"Errors on parameters:\nphi*: {perr[0]:.3e}, Alpha: {perr[1]:.2f}, Beta: {perr[2]:.2f}, M*: {perr[3]:.2f}")
err_phiStar_fit, err_alpha_fit, err_beta_fit, err_Mstar_fit = perr


# mag range to plot over
M = np.arange(-25, -19, 0.1)

# Bowler+17
z7_gal = dpl(2.3*10**(-4.), -2.19, -4.60, M, -20.60)

# Harikane+24
z7_harikane = dpl(10**(-3.74), -2.08, -4.81, M, -21.01)

# Harikane GOLDRUSH
DPLy_h22 = dpl(10**(-3.05), -1.89, -3.81, M, -20.12) + dpl(10**(-8.49), -1.23, -2.73, M, -24.9)


# Print ratio between my LF points and the Harikane+24 LF function at that magnitude
ratios = []
for i, M_ in enumerate(Muv_bins[:-1]):

    ratio = phi[i] / dpl(10**(-3.74), -2.08, -4.81, M_, -21.01)
    ratios.append(round(ratio, 2))
print('###### Ratios #######')
#print(ratios)

#! Plot!
plt.figure(figsize=(10, 10))
plt.plot(M, z7_gal, color='green', linewidth=3, label='Bowler+17', alpha=0.7, linestyle='--')
plt.plot(M, z7_harikane, color='blue', linewidth=3, label='Harikane+24', alpha=0.7, linestyle=':')
# plt.plot(M, DPLy_h22, color='gray', linewidth=3, label='Harikane+22', alpha=0.5)

#! My best fit DPL function
M_fit = np.linspace(-25, -18, 100)
LF_fit = dpl_fit(M_fit, *popt)
plt.plot(M_fit, LF_fit, color='red', linewidth=5, label="Best-fit DPL", alpha=0.8)



# McLure+13
plt.errorbar(m13x, m13y, yerr=m13dy, color='magenta', label='McLure+13', marker='D', markersize=10, alpha=0.8, linestyle='none', markerfacecolor='none')

# Bowuens et al. 2021
plt.errorbar(b21x, b21y, color='orange', yerr=b21dy, label='Bouwens+21', marker='o', markerfacecolor='none', markersize=10, alpha=0.8, linestyle='none')

# Bowler+17
plt.errorbar(b17x-0.05, b17y, color='green', yerr=b17dy, label='Bowler+17', marker='s', markersize=10, alpha=0.8, linestyle='none', zorder=3)

# Varadaraj+23
plt.errorbar(v23x-0.05, v23y, yerr=v23dy, fmt='o', color='black', 
             ecolor='black', elinewidth=2, label='Varadaraj+23', markersize=10, markeredgecolor='black', zorder=2)

# Harikane+24 LF points
plt.errorbar(h24x, h24y, yerr=[h24dy_lo, h24dy_up], fmt='D', color='blue', alpha=0.8, label='Harikane+24', markersize=10, zorder=2)

# Plot Finkelstein+15
plt.errorbar(f15x, f15y, yerr=[f15y_lo, f15y_up], fmt='p', color='purple', alpha=0.8, label='Finkelstein+15', markersize=12, zorder=4)

# Harikane+22 LF points
plt.errorbar(yh22x-0.05, yh22y, yerr=[yh22dyL, yh22dyU], fmt='s', color='gray', alpha=0.8, label='Harikane+22', markersize=10, zorder=2, markerfacecolor='none')

# Plot the LF new points
plt.errorbar(bin_centres, phi, yerr=delta_phi, fmt='o', color='red', xerr=bin_widths/2,
             ecolor='red', elinewidth=4, label='This work', markersize=13, markeredgecolor='black', zorder=5)


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
#plt.savefig(plot_dir / 'LF_UVISTA.pdf')
plt.savefig(plot_dir / 'LF_UVISTA_complete.pdf')
plt.show()


