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
from matplotlib.lines import Line2D
import corner

plt.rcParams.update({'font.size': 15})
plt.rcParams['axes.linewidth'] = 4
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

plot_dir = Path.cwd().parents[1] / 'plots' / 'LF'

#! Switches
#run_type = ''
run_type = 'with_euclid'   
fit = True
kron = False # Use Muv updates based on Kron photometry

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



def dpl_fit(M, phiStar, Mstar, alpha, beta):
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


# Log-likelihood function
def log_likelihood(theta, M, phi, phi_err):
    #!DPL
    # M_star, phi_star, alpha, beta = theta
    # model = dpl_fit(M, phi_star, M_star, alpha, beta)
    # return -0.5 * np.sum(((phi - model) / phi_err) ** 2)
    #!Schechter
    M_star, phi_star, alpha = theta
    model = schechter_fit(M, phi_star, M_star, alpha)
    return -0.5 * np.sum(((phi - model) / phi_err) ** 2)

#* Flat priors
def log_prior(theta):
    M_star, phi_star, alpha, beta = theta
    if -22.5 < M_star < -19 and 1e-5 < phi_star < 9e-4 and -3 < alpha < -1 and -6 < beta < -3:
        return 0.0
    return -np.inf

#* Gaussiean priors on alpha, beta
# def log_prior(theta):
#     M_star, phi_star, alpha, beta = theta

#     # Flat priors for M_star and phi_star (keep reasonable limits)
#     if not (-22.5 < M_star < -19. and 1e-6 < phi_star < 1e-3):
#         return -np.inf  # Return -inf if out of bounds

#     # Gaussian priors for alpha and beta
#     logp_alpha = -0.5 * ((alpha + 2.) / 0.5) ** 2  # Mean = -2.0, Sigma = 0.5
#     logp_beta = -0.5 * ((beta + 4.5) / 0.7) ** 2    # Mean = -4.5, Sigma = 0.7

#     return logp_alpha + logp_beta  # Sum of log priors (log-likelihood adds them)


#* Gaussian priors on all params?
def log_prior(theta):
    #!DPL
    #M_star, phi_star, alpha, beta = theta
    # Gaussian priors (mean, sigma)
    # logp_M_star  = -0.5 * ((M_star  + 21) / 1.0) ** 2   # Mean = -20.5, Sigma = 1.0
    # logp_phi_star = -0.5 * ((phi_star - 4e-4) / 1e-3) ** 2  # Mean = 4e-4, Sigma = 1e-3
    # logp_alpha   = -0.5 * ((alpha + 2.0) / 0.5) ** 2   # Mean = -2.0, Sigma = 0.5
    # logp_beta    = -0.5 * ((beta + 5.0) / 0.5) ** 2   # Mean = -4.0, Sigma = 0.5
    #return logp_phi_star + logp_M_star + logp_alpha + logp_beta

    #! DPL
    M_star, phi_star, alpha = theta
    # Gaussian priors (mean, sigma)
    logp_M_star  = -0.5 * ((M_star  + 21.15) / 1.0) ** 2   # Mean = -20.5, Sigma = 1.0
    logp_phi_star = -0.5 * ((phi_star - 0.19e-3) / 1e-3) ** 2  # Mean = 4e-4, Sigma = 1e-3
    logp_alpha   = -0.5 * ((alpha + 2.06) / 0.5) ** 2   # Mean = -2.0, Sigma = 0.5

    return logp_phi_star + logp_M_star + logp_alpha


# Full probability function
def log_probability(theta, M, phi, phi_err):
    lp = log_prior(theta)
    if not np.isfinite(lp):
        return -np.inf
    return lp + log_likelihood(theta, M, phi, phi_err)

#!################## READ IN DATA ####################

#? Read catalogue
cat_dir = Path.cwd().parents[1] / 'data' / 'catalogues' / 'candidates'
#cat_name = 'COSMOS_5sig_Y_J_nonDet_HSC_G_nonDet_HSC_R_nonDet_HSC_I_candidates_2025_01_31.fits'
#cat_name = 'COSMOS_5sig_Y_J_nonDet_HSC_G_nonDet_HSC_R_nonDet_HSC_I_candidates_2024_11_20.fits' # inclusive

if run_type == '':
    cat_name = 'COSMOS_5sig_Y_J_nonDet_HSC_G_nonDet_HSC_R_nonDet_HSC_I_candidates_2025_02_14.fits' # just vista
if run_type == 'with_euclid':
    cat_name = 'COSMOS_5sig_Y_J_nonDet_HSC_G_nonDet_HSC_R_nonDet_HSC_I_candidates_2025_02_14_with_euclid.fits' # with euclid
    cat_name = 'Euclid_UltraVISTA_z7_sample.fits'
    #cat_name = 'COSMOS_5sig_Y_J_nonDet_HSC_G_nonDet_HSC_R_nonDet_HSC_I_candidates_2025_02_14_with_euclid_WITH_MAG_AUTO.fits' # with euclid and mag_auto

t = Table.read(cat_dir / cat_name)

# Remove the Lya emitters which have z>7.5 with no emission line.
t = t[t['Vmax'] > 0]
t = t[t['Muv'] < 0]

# Remove certain IDs
IDs_to_remove = [345668, 488168, 605134, 729259, 796957]
for remove in IDs_to_remove:
    t = t[t['ID'] != remove]

# print('Min/max Muv:')
# print(np.min(t['Muv']), np.max(t['Muv']))
# exit()

# Remove object with ID 1109577
# t = t[t['ID'] != 1109577]

# plt.scatter(t['Muv'], t['Vmax'], c=t['Zphot'], cmap='viridis', s=10)    
# plt.show()
# exit()

# Restrict to Muv < -20.5
print('Number of galaxies before Muv cut: ', len(t))
#t = t[t['Muv'] < -20.95]
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

#* VISTA
if run_type == '':
    Muv_bins = [-22.8, -22.4, -22., -21.8,  -21.6, -21.4, -21.2, -21.0, -20.8, -20.6, -20.4, -20.2]
    bin_widths = np.abs(np.diff(Muv_bins))
    bin_centres = 0.5 * (np.array(Muv_bins[:-1]) + np.array(Muv_bins[1:]))

#* VISTA + Euclid
if run_type == 'with_euclid':
    Muv_bins = [-22.4, -22., -21.8,  -21.6, -21.4, -21.2, -21.0, -20.8, -20.6, -20.4, -20.2]
    bin_widths = np.abs(np.diff(Muv_bins))
    bin_centres = 0.5 * (np.array(Muv_bins[:-1]) + np.array(Muv_bins[1:]))
    #? Adding brighter bin after adjusting Muv with kron photometry
    # Muv_bins = [-22.8, -22.4, -22., -21.8,  -21.6, -21.4, -21.2, -21.0, -20.8, -20.6, -20.4, -20.2]
    # bin_widths = np.abs(np.diff(Muv_bins))
    # bin_centres = 0.5 * (np.array(Muv_bins[:-1]) + np.array(Muv_bins[1:]))

print(f'Muv bins: {np.round(Muv_bins, 2)}')
print(f'Bin widths: {np.round(bin_widths, 2)}')
print(f'Bin centres: {np.round(bin_centres, 2)}')

# plt.hist(t['Muv'], bins=np.arange(Muv_min, Muv_max, 0.01), alpha=0.5)
# for Muv in Muv_bins:
#     plt.axvline(Muv, color='black', linestyle='--')

# plt.show()
# exit()

#! Depending on whether we want to use Kron-adjusted magnitudes, we use different columns for Muv
if kron == False:
    Muv_str = 'Muv'
if kron:
    Muv_str = 'Muv_linear_fit'

# Split the table into these bins
binned_tables = []
for i in range(len(Muv_bins)-1):
    mask = (t[Muv_str] >= Muv_bins[i]) & (t[Muv_str] < Muv_bins[i+1])
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
        Muv = obj[Muv_str]
        
        # Find completeness for this Muv,z
        Muv_bin = np.digitize([Muv], Muv_completeness_bins)[0] - 1
        z_bin = np.digitize([z], z_completeness_bins)[0] - 1

        # Ensure indices stay within bounds
        z_bin = max(0, min(z_bin, completeness_matrix.shape[1] - 1))
        Muv_bin = max(0, min(Muv_bin, completeness_matrix.shape[0] - 1))
    

        completeness = completeness_matrix[Muv_bin, z_bin]

        # Summand
        summand = 1 / (Vmax * completeness)
        phi[i] += summand

        # Compute error term
        err_summand = 1 / Vmax ** 2
        delta_phi[i] += err_summand



##############################! ERROR, INCLUDING COSMIC VARIANCE ##############################
for i, lf in enumerate(phi):

    phi[i] /= bin_widths[i]  
    delta_phi[i] = np.sqrt(delta_phi[i]) / bin_widths[i]


#? cosmic variance from https://www.ph.unimelb.edu.au/~mtrenti/cvc/CosmicVariance.html
if run_type == '':
    cv_per_bin = [0.118, 0.110, 0.108, 0.110, 0.102, 0.097, 0.095, 0.097, 0.096, 0.103, 0.123]
if run_type == 'with_euclid':
    cv_per_bin = [0.175, 0.171, 0.164, 0.157, 0.148, 0.144, 0.142, 0.144, 0.149, 0.171]
    #cv_per_bin = [0.118, 0.110, 0.108, 0.110, 0.102, 0.097, 0.095, 0.097, 0.096, 0.103, 0.123]

# Add this much percentage error to the LF in quadrature
delta_phi = np.sqrt(delta_phi**2 + (cv_per_bin * phi)**2)

print('LF values:', phi)
print('LF errors:', delta_phi)

# for i, LF in enumerate(phi):
#     print(LF, delta_phi[i])
#! ############################## Existing data points ##############################

# Bouwens+21 
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

# Harikane+22 
yh22x = [-24.92, -24.42, -23.92, -23.42, -22.92, -22.42, -21.92]
#yh22y = [5e-10, 1.31e-8, 4.39e-8, 1.83e-7, 1.06e-6, 2.75e-6] # 'Galaxy LF'
yh22y = [8.89e-9, 2.41e-8, 9.02e-8, 1.62e-7, 4.63e-7, 1.95e-6, 3.47e-6] # Dropout LF
yh22dyU = [15.29e-9, 1.75e-8, 2.13e-8, 0.3e-7, 7.53e-7, 1.13e-6, 1.7e-6]
yh22dyL = [7.41e-9, 0.98e-8, 1.66e-8, 0.23e-7, 2.71e-7, 0.67e-6, 1.03e-6]
yh22x, yh22y, yh22dyU, yh22dyL = np.array(yh22x), np.array(yh22y), np.array(yh22dyU), np.array(yh22dyL)

# Varadaraj+23 LF points
v23x = [-22.175, -22.925, -23.675]#, -24.425])
v23y = [2.70e-6, 2.81e-7, 2.37e-8]#, np.nan])  
v23dy = [0.66e-6, 1.54e-7, 2.50e-8]#, np.nan])  
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

#Rojas-ruiz+25 LF points
rr25x = [-21, -22]
rr25y = [107.31e-6, 4.6464e-6]
rr25dy_up = [54.76e-6, 5.231e-6]
rr25dy_lo = [40.77e-6, 2.918e-6]

# Franco+25 JWST COSMOS-Web points at z=6.4 (5.1 < z < 8.8)
fr25x = [-22.75, -22.25, -21.75, -21.25, -20.75, -20.25, -19.75]
fr25y = [3.72e-6, 7.1e-6, 26.58e-6, 69.76e-6, 151.07e-6, 232.66e-6, 260.12e-6]
fr25dy_up = [1.16e-6, 1.55e-6, 3.91e-6, 7.72e-6, 13.25e-6, 17.84e-6, 18.83e-6]
fr25dy_lo = [1.16e-6, 1.61e-6, 3.98e-6, 7.81e-6, 13.73e-6, 18.53e-6, 19.35e-6]

# # UVISTA DR6 LF points, without completeness, on full inclusive sample of 243 objects
# muv = [-22.79501799, -22.29501799, -21.79501799, -21.29501799, -20.79501799, -20.29501799]
# phi = [4.64101649e-07, 2.33761186e-06, 4.46459736e-06, 1.09373428e-05, 6.40936287e-06]
# phi_err = [4.15105134e-07, 9.35645234e-07, 1.28888147e-06, 2.14176670e-06, 1.75013185e-06]

####################! Fit the LF! ####################

#? Combine my points with other points:

if run_type == '':
    # Leave out last four bins
    bin_centres_to_fit = bin_centres[:-3]
    phi_to_fit = phi[:-3]
    delta_phi_to_fit = delta_phi[:-3]
if run_type == 'with_euclid':
    # Leave out last three bins
    bin_centres_to_fit = bin_centres[:-3]
    phi_to_fit = phi[:-3]
    delta_phi_to_fit = delta_phi[:-3]

Muv_combined, LF_combined, delta_phi_combined = concatenate_arrays(
    (bin_centres_to_fit, phi_to_fit, delta_phi_to_fit), # My new UltraVISTA points
    (v23x[:-1], v23y[:-1], v23dy[:-1]), # Varadaraj+23, minus CDFS BD
    #(v23x, v23y, v23dy), # Varadaraj+23
    (f15x, f15y, f15y_up), # Finkelstein+15
    #(m13x, m13y, m13dy), # McLure+13
    #(b17x, b17y, b17dy), # Bowler+17
)

#?####################### SCIPY Curve fitting #################################
# Initial parameter guess (phi*, alpha, beta, M*)
# p0 = [3.6e-4, -2.1, -4.8, -20.3]

# # Fit the function
# popt, pcov = curve_fit(dpl_fit, Muv_combined, LF_combined, sigma=delta_phi_combined, p0=p0, maxfev=10000)

# # Extract best-fit parameters
# phi_star_best, alpha_best, beta_best, M_star_best = popt

# # print(f"Best-fit parameters:\nphi*: {phiStar_fit:.3e}, Alpha: {alpha_fit:.2f}, Beta: {beta_fit:.2f}, M*: {Mstar_fit:.2f}")

# # Get errors on the parameters
# perr = np.sqrt(np.diag(pcov))
# print(f"Errors on parameters:\nphi*: {perr[0]:.3e}, Alpha: {perr[1]:.2f}, Beta: {perr[2]:.2f}, M*: {perr[3]:.2f}")
# err_phiStar_fit, err_alpha_fit, err_beta_fit, err_Mstar_fit = perr

#? ####################### EMCEE ############################
if fit:
    ##!! DPL
    # # Initial parameter guesses
    # initial = [-21.2, 2.5e-4, -2.0, -4.5]
    # ndim, nwalkers = len(initial), 50

    # #pos = initial + 1e-3 * np.random.randn(nwalkers, ndim)
    # pos = np.array([
    #     np.random.uniform(-22, -19, nwalkers),  # M_star
    #     np.random.uniform(1e-4, 1e-3, nwalkers),  # phi_star
    #     np.random.uniform(-2.5, -1.5, nwalkers),  # alpha
    #     np.random.uniform(-5, -3.5, nwalkers)  # beta
    # ]).T

    # # Run MCMC
    # sampler = emcee.EnsembleSampler(nwalkers, ndim, log_probability, args=(Muv_combined, LF_combined, delta_phi_combined))
    # sampler.run_mcmc(pos, 50000, progress=True)

    # # Extract best-fit parameters
    # samples = sampler.get_chain(discard=1000, thin=15, flat=True)
    # M_star_best, phi_star_best, alpha_best, beta_best = np.median(samples, axis=0)
    # print(f"Best-fit parameters:\nM*: {M_star_best:.4f}, phi*: {phi_star_best:.4e}, Alpha: {alpha_best:.4f}, Beta: {beta_best:.4f}")

    # # Corner plot of fitting
    # labels = [r"$M^*$", r"$\phi^*$", r"$\alpha$", r"$\beta$"]
    # fig = corner.corner(samples, labels=labels, quantiles=[0.16, 0.5, 0.84], show_titles=True)

    # # Convert phi* label values to scientific notation
    # for ax in fig.get_axes():
    #     if ax.get_xlabel() == r"$\phi^*$":
    #         ax.get_xaxis().get_offset_text().set_visible(False)

    # plt.savefig(plot_dir / 'corner_plot.pdf')
    # plt.close()

    # tau = sampler.get_autocorr_time()
    # #print("Autocorrelation Time Estimate:", np.mean(tau))

    # # Get 1sigma upper and lower errors
    # M_star_best_err_up, phi_star_best_err_up, alpha_best_err_up, beta_best_err_up = np.percentile(samples, 84, axis=0) - [M_star_best, phi_star_best, alpha_best, beta_best]
    # M_star_best_err_lo, phi_star_best_err_lo, alpha_best_err_lo, beta_best_err_lo = [M_star_best, phi_star_best, alpha_best, beta_best] - np.percentile(samples, 16, axis=0)
    # print('Upper errors on best-fit parameters:', M_star_best_err_up, phi_star_best_err_up, alpha_best_err_up, beta_best_err_up)
    # print('Lower errors on best-fit parameters:', M_star_best_err_lo, phi_star_best_err_lo, alpha_best_err_lo, beta_best_err_lo)    

    ##!! SCHECHTER
    # # Initial parameter guesses
    # initial = [-21.15, 0.19e-3, -2.06]
    # ndim, nwalkers = len(initial), 50
    # #pos = initial + 1e-3 * np.random.randn(nwalkers, ndim)
    # pos = np.array([
    #     np.random.uniform(-22, -19, nwalkers),  # M_star
    #     np.random.uniform(1e-4, 1e-3, nwalkers),  # phi_star
    #     np.random.uniform(-2.5, -1.5, nwalkers),  # alpha
    # ]).T

    # # Run MCMC
    # sampler = emcee.EnsembleSampler(nwalkers, ndim, log_probability, args=(Muv_combined, LF_combined, delta_phi_combined))
    # sampler.run_mcmc(pos, 50000, progress=True)

    # # Extract best-fit parameters
    # samples = sampler.get_chain(discard=1000, thin=15, flat=True)
    # M_star_best_sch, phi_star_best_sch, alpha_best_sch = np.median(samples, axis=0)
    # print(f"Best-fit parameters:\nM*: {M_star_best_sch:.4f}, phi*: {phi_star_best_sch:.4e}, Alpha: {alpha_best_sch:.4f}")

    # # Corner plot of fitting
    # labels = [r"$M^*$", r"$\phi^*$", r"$\alpha$", r"$\beta$"]
    # fig = corner.corner(samples, labels=labels, quantiles=[0.16, 0.5, 0.84], show_titles=True)

    # # Convert phi* label values to scientific notation
    # for ax in fig.get_axes():
    #     if ax.get_xlabel() == r"$\phi^*$":
    #         ax.get_xaxis().get_offset_text().set_visible(False)

    # plt.savefig(plot_dir / 'corner_plot_sch.pdf')
    # plt.close()

    # tau = sampler.get_autocorr_time()
    # #print("Autocorrelation Time Estimate:", np.mean(tau))

    # # Get 1sigma upper and lower errors
    # M_star_best_err_up, phi_star_best_err_up, alpha_best_err_up = np.percentile(samples, 84, axis=0) - [M_star_best_sch, phi_star_best_sch, alpha_best_sch]
    # M_star_best_err_lo, phi_star_best_err_lo, alpha_best_err_lo  = [M_star_best_sch, phi_star_best_sch, alpha_best_sch] - np.percentile(samples, 16, axis=0)
    # print('Upper errors on best-fit parameters:', M_star_best_err_up, phi_star_best_err_up, alpha_best_err_up)
    # print('Lower errors on best-fit parameters:', M_star_best_err_lo, phi_star_best_err_lo, alpha_best_err_lo)    


    if run_type == 'with_euclid':
        # Best-fit parameters from emcee]
        # M*: -21.1304, phi*: 9.1146e-05, Alpha: -2.1085, Beta: -4.5987
        # Upper errors on best-fit parameters: 0.27047277962511274 6.654507642353416e-05 0.21369413904090906 0.32097933349011587
        # Lower errors on best-fit parameters: 0.24948880288177477 3.8252481522955595e-05 0.17156830499011377 0.3703468864815038
        M_star_best = -21.1304
        phi_star_best = 9.1146e-05
        alpha_best = -2.1085
        beta_best = -4.5987

        # Schechter
        # M*: -20.9796, phi*: 1.6091e-04, Alpha: -1.9935
        # Upper errors on best-fit parameters: 0.1891394070272625 7.51883613921671e-05 0.16454132162330004
        # Lower errors on best-fit parameters: 0.20335279382762295 5.686421729145396e-05 0.15175426010573934
        M_star_best_sch = -20.9796
        phi_star_best_sch = 1.6091e-04
        alpha_best_sch = -1.9935

    if run_type == '':
        # M*: -20.8885, phi*: 1.3985e-04, Alpha: -2.0130, Beta: -4.2563
        # Upper errors on best-fit parameters: 0.27409088522794534 9.835983881875717e-05 0.25148557544353145 0.24448325643043045
        # Lower errors on best-fit parameters: 0.2850048058690291 6.376084300199492e-05 0.21107831155776502 0.2876519425081936
        M_star_best = -20.8885
        phi_star_best = 1.3985e-04
        alpha_best = -2.0130
        beta_best = -4.2563

#############################! Plotting ##########################§
# mag range to plot over
M = np.arange(-26, -18, 0.1)

# Bowler+17
z7_gal = dpl(2.3*10**(-4.), -2.19, -4.60, M, -20.60)

# Harikane+24
z7_harikane = dpl(10**(-3.74), -2.08, -4.81, M, -21.01)

# Harikane GOLDRUSH
DPLy_h22 = dpl(10**(-3.05), -1.89, -3.81, M, -20.12) + dpl(10**(-8.49), -1.23, -2.73, M, -24.9)

# Bouwens+21
b21_schechter = schechter(0.19e-3, -2.06, M, -21.15)

# Franco+25 DPL
fr25_dpl = dpl(23.39e-5, -2.10, -3.42, M, -20.88)

# Print ratio between my LF points and the Harikane+24 LF function at that magnitude
ratios = []
for i, M_ in enumerate(Muv_bins[:-1]):

    ratio = phi[i] / dpl(10**(-3.74), -2.08, -4.81, M_, -21.01)
    ratios.append(round(ratio, 2))

#print('###### Ratios #######')
#print(ratios)

#! Plot!
if fit == True:
    plt.figure(figsize=(10, 10))
if fit == False:
    plt.figure(figsize=(10, 7))

#! Literature functional fits
plt.plot(M, z7_gal, color='green', linewidth=3, alpha=0.6, linestyle=':') # label='Bowler+17', 
plt.plot(M, z7_harikane, color='blue', linewidth=3, alpha=0.9, linestyle='-.') #label='Harikane+25', 
plt.plot(M, DPLy_h22, color='gray', linewidth=3, alpha=0.8, linestyle=':') #label='Harikane+22', 
plt.plot(M, b21_schechter, color='orange', linewidth=3, alpha=0.9, linestyle='--') #label='Bouwens+21', 
if fit:
    plt.plot(M, fr25_dpl, color='tab:orange', linewidth=3, alpha=0.9, linestyle=':') # label=Franco+25

#! My best fit DPL function
if fit:
    M_fit = np.linspace(-25, -18, 100)
    # LF_fit = dpl_fit(M_fit, *popt)
    # plt.plot(M_fit, LF_fit, color='red', linewidth=5, label="Best-fit DPL", alpha=0.8)

    # Plot the best-fit DPL function from emcee
    LF_fit_sch = schechter_fit(M_fit, phi_star_best_sch, M_star_best_sch, alpha_best_sch)
    plt.plot(M_fit, LF_fit_sch, color='tab:red', linewidth=5, label="Best-fit Schechter", alpha=0.9, linestyle='--')

    LF_fit = dpl_fit(M_fit, phi_star_best, M_star_best, alpha_best, beta_best)
    plt.plot(M_fit, LF_fit, color='tab:red', linewidth=7, label="Best-fit DPL", alpha=0.9)



#! Literature LF points
# McLure+13
plt.errorbar(m13x, m13y, yerr=m13dy, color='yellowgreen', label='McLure+13', marker='D', markersize=10, alpha=0.8, linestyle='none', markerfacecolor='none')

# Bowuens et al. 2021
plt.errorbar(b21x, b21y, color='orange', yerr=b21dy, marker='o', markerfacecolor='none', markersize=10, alpha=0.8, linestyle='none') #label='Bouwens+21',

# Bowler+17
plt.errorbar(b17x-0.05, b17y, color='green', yerr=b17dy, marker='s', markersize=10, alpha=0.8, linestyle='none', zorder=3) #label='Bowler+17'

# Varadaraj+23
if fit == False:
    plt.errorbar(v23x+0.05, v23y, yerr=v23dy, fmt='o', color='black', 
                ecolor='black', elinewidth=2, label='Varadaraj+23', markersize=10, markeredgecolor='black', zorder=4)
if fit == True:
    plt.errorbar(v23x+0.05, v23y, yerr=v23dy, fmt='o', color='black', 
                ecolor='black', elinewidth=2, label='Varadaraj+23', markersize=15, markeredgecolor='black', zorder=4)    
# plt.errorbar(v23x+0.05, v23y, yerr=v23dy, fmt='s', color='black', 
#              ecolor='red', elinewidth=4, label='Varadaraj+23', markersize=13, markeredgecolor='black', zorder=4)

# Harikane+24 LF points
plt.errorbar(h24x, h24y, yerr=[h24dy_lo, h24dy_up], fmt='D', color='blue', alpha=0.8, markersize=10, zorder=2) #label='Harikane+25',

# Plot Finkelstein+15
plt.errorbar(f15x-0.05, f15y, yerr=[f15y_lo, f15y_up], fmt='p', color='tab:purple', alpha=0.9, label='Finkelstein+15', markersize=12, zorder=4)

# Harikane+22 LF points
plt.errorbar(yh22x-0.05, yh22y, yerr=[yh22dyL, yh22dyU], fmt='s', color='gray', alpha=0.8, markersize=10, zorder=2, markerfacecolor='none') #label='Harikane+22',

# Rojas-ruiz+25 points
plt.errorbar(rr25x, rr25y, yerr=[rr25dy_lo, rr25dy_up], fmt='P', color='tab:pink', alpha=0.8, label=r'Rojas-Ruiz+25', markersize=10, zorder=2) #, markerfacecolor='none')

# Franco+25 CWEB points
plt.errorbar(fr25x, fr25y, yerr=[fr25dy_lo, fr25dy_up], fmt='H', color='tab:orange', alpha=0.8, markersize=13, zorder=2)

# Plot the LF new points
if run_type == 'with_euclid':
    if fit == False:
        plt.errorbar(bin_centres[:-3], phi[:-3], yerr=delta_phi[:-3], fmt='o', color='tab:red', xerr=bin_widths[:-3]/2,
                    ecolor='tab:red', elinewidth=4, label=r'UltraVISTA + $Euclid$', markersize=15, markeredgecolor='black', zorder=5)
        
        # Add an upper limit at a -22.6 bin
        up_lim = 1.841 / 4090450.2883619294 / 0.4
        print('2sigma upper limit:' , up_lim)
        x = -22.6
        y = up_lim           # your y value
        yerr = up_lim * 0.4  # your y error
        xerr = 0.2           # your x error

        # 1. Plot y error bar with upper limit arrow and capsize
        plt.errorbar(
            x, y,
            yerr=yerr,
            fmt='o',
            color='tab:red',
            ecolor='tab:red',
            elinewidth=4,
            markersize=15,
            markeredgecolor='black',
            zorder=5,
            uplims=True,
            capsize=9  # Controls the size of the upper limit arrow
        )

        # 2. Add x error bar separately with no capsize (capsize=0 removes caps)
        plt.errorbar(
            x, y,
            xerr=xerr,
            fmt='none',
            ecolor='tab:red',
            elinewidth=4,
            capsize=0  # No caps for x-error
        )

        #? Put the VISTA LF on the same plot to compare the scatter directly
        LF_vista = [1.30515748e-06, 3.32642337e-06, 7.34422140e-06, 6.27708092e-06,
        1.57390423e-05, 2.74080836e-05, 4.12950686e-05, 3.83475493e-05,
        5.04532349e-05, 3.29772144e-05, 4.49456387e-06]
        LF_vista_err = [4.42753787e-07, 7.49204721e-07, 1.55601274e-06, 1.41177516e-06,
        2.48501893e-06, 3.60833675e-06, 5.01404632e-06, 4.72413805e-06,
        5.90207334e-06, 4.52903972e-06, 9.73488143e-07]
        Muv_bins_vista = [-22.8, -22.4, -22., -21.8,  -21.6, -21.4, -21.2, -21.0, -20.8, -20.6, -20.4, -20.2]
        bin_widths_vista = np.abs(np.diff(Muv_bins_vista))
        bin_centres_vista = 0.5 * (np.array(Muv_bins_vista[:-1]) + np.array(Muv_bins_vista[1:]))

        plt.errorbar(bin_centres_vista[:-3], LF_vista[:-3], yerr=LF_vista_err[:-3], fmt='o', color='dodgerblue',
            ecolor='dodgerblue', elinewidth=2, label='UltraVISTA-only', markersize=11, markeredgecolor='black', zorder=5, markerfacecolor='dodgerblue', alpha=1)

    if fit == True:
        #We are leaving out the last three points in the fitting, so make them faint
        plt.errorbar(bin_centres[:-3], phi[:-3], yerr=delta_phi[:-3], fmt='o', color='tab:red', xerr=bin_widths[:-3]/2,
                    ecolor='tab:red', elinewidth=4, label='This work', markersize=20, markeredgecolor='black', zorder=5)
if run_type == '':
    # We are leaving out the last four points in the fitting, so make them faint
    plt.errorbar(bin_centres[:-3], phi[:-3], yerr=delta_phi[:-3], fmt='o', color='red', xerr=bin_widths[:-3]/2,
                ecolor='red', elinewidth=4, label='This work', markersize=15, markeredgecolor='black', zorder=5)
    
    # plt.errorbar(bin_centres[-4:], phi[-4:], yerr=delta_phi[-4:], fmt='o', color='red', xerr=bin_widths[-4:]/2,
    #             ecolor='red', elinewidth=3, markersize=8, markeredgecolor='red', markeredgewidth=1, markerfacecolor='white', zorder=5, alpha=1)

# Plotting all points without making incomplete ones fainter
# plt.errorbar(bin_centres, phi, yerr=delta_phi, fmt='o', color='red', xerr=bin_widths/2,
#             ecolor='red', elinewidth=4, label='This work', markersize=15, markeredgecolor='black', zorder=5)

#Draw an up arrow at the value of M*
if fit:
    plt.annotate('', xy=(M_star_best, 2e-5), xytext=(M_star_best, 4e-7),
                    arrowprops=dict(facecolor='black',  lw=1), fontsize=40)
    plt.text(M_star_best-1, 2e-7, r'$M^*=-21.13^{+0.27}_{-0.25}$', fontsize=20)


plt.tick_params(which='major', length=10, width=3)
plt.tick_params(axis='both', which='minor', length=5, width=2)

if fit == False:
    plt.xlabel(r'$M_{\mathrm{UV}}$', fontsize=25)
    # plt.ylabel(r'$\mathrm{Number \ of \ objects \  / \ mag \ / \ Mpc^{3}}$', fontsize=21.5)
    plt.ylabel(r'$\log_{10}(\phi\,/\,\rm{mag}^{-1}\,\rm{Mpc}^{-3})$', fontsize=25)
if fit:
    plt.xlabel(r'$M_{\mathrm{UV}}$', fontsize=28)
    # plt.ylabel(r'$\mathrm{Number \ of \ objects \  / \ mag \ / \ Mpc^{3}}$', fontsize=21.5)
    plt.ylabel(r'$\log_{10}(\phi\,/\,\rm{mag}^{-1}\,\rm{Mpc}^{-3})$', fontsize=28)
# 

# Increase size of tick labels
plt.xticks(fontsize=20)
plt.yticks(fontsize=20)

# Add minor ticks in between x axis ticks
plt.minorticks_on()

#plt.xlim(-25, -19)
# plt.ylim(1e-10, 1e-2)

if fit == True:
    if run_type == '':
        plt.text(-21, 1.5e-10, 'UltraVISTA only', size=25)
    if run_type == 'with_euclid':
        #plt.text(-21.93, 1.3e-9, r'UltraVISTA + $Euclid$', size=31)
        plt.text(-24.4, 2e-3, r'UltraVISTA + $Euclid$', size=31)
    plt.text(-19.9, 2e-3,  r'$z=7$', size=31)
if fit == False:
    if run_type == '':
        plt.text(-21.6, 2.4e-7, 'UltraVISTA only', size=25)
    if run_type == 'with_euclid':
        #plt.text(-21.85, 2.4e-7, r'UltraVISTA + $Euclid$', size=25)
        plt.text(-21.85, 2.4e-7, r'', size=25)
    plt.text(-22.9, 9e-5,  r'$z=7$', size=25)

if fit == True:
    plt.ylim([1e-9, 5e-3])
    plt.xlim(-24.5, -19)
if fit == False:
    plt.ylim([3e-7, 2e-4])
    plt.xlim([-23, -20.7])

if fit == True:

    # custom legend handles
    mclure13_handle = Line2D([], [], color='yellowgreen', linestyle='none', marker='D', markersize=10, label='McLure+13', markerfacecolor='none')
    finkelstein15_handle = Line2D([], [], color='tab:purple', linestyle='none', marker='p', markersize=12, label='Finkelstein+15')
    bowler17_handle = Line2D([], [], color='green', linestyle=':', marker='s', markersize=10, label='Bowler+17')
    bouwens21_handle = Line2D([], [], color='orange', linestyle='--', marker='o', markersize=10, label='Bouwens+21', markerfacecolor='none')
    harikane22_handle = Line2D([], [], color='gray', linestyle=':', marker='s', markersize=10, label='Harikane+22', markerfacecolor='none')
    varadaraj23_handle = Line2D([], [], color='black', linestyle='none', marker='o', markersize=15, label='Varadaraj+23', markeredgecolor='black')
    harikane25_handle = Line2D([], [], color='blue', linestyle='-.', marker='D', markersize=10, label='Harikane+25')
    this_work_handle = Line2D([], [], color='tab:red', linestyle='none', marker='o', markersize=15, label='This work', markeredgecolor='black')
    this_work_dpl_handle = Line2D([], [], color='tab:red', linestyle='-', linewidth=7, label='DPL fit', alpha=0.9)
    this_work_sch_handle = Line2D([], [], color='tab:red', linestyle='--', linewidth=5, label='Schechter fit', alpha=0.9)
    rojasruiz25_handle = Line2D([], [], marker='P', color='tab:pink', alpha=0.8, label=r'Rojas-Ruiz+25', markersize=10, zorder=2, linestyle='none')
    franco25_handle = Line2D([], [], marker='h', color='tab:orange', alpha=0.8, label='Franco+25', markersize=13, zorder=2, linestyle=':', linewidth=3)

    # Now manually create the legend using only the desired handles
    handles = [
        mclure13_handle,
        finkelstein15_handle,
        bowler17_handle,
        bouwens21_handle,
        harikane22_handle,
        varadaraj23_handle,
        harikane25_handle,
        rojasruiz25_handle,
        franco25_handle,
        this_work_handle,
        this_work_dpl_handle,
        this_work_sch_handle,
]

    #plt.legend(loc='upper left', fontsize=17)
    plt.legend(handles=handles, loc='lower right', fontsize=17.5 if fit else 12, ncol=2 if not fit else 2, frameon=True)
if fit == False:
    # Dummy plots for studies with LF points and a best fit function, to combine the two
    plt.plot([], [], color='green', linestyle=':', marker='s', markersize=10, label='Bowler+17')
    plt.plot([], [], color='gray', linestyle=':', marker='s', markersize=10, label='Harikane+22', markerfacecolor='none')
    plt.plot([], [], color='blue', linestyle='-.', marker='D', markersize=10, label='Harikane+25')
    plt.plot([], [], color='orange', linestyle='--', marker='o', markersize=10, label='Bouwens+21', markerfacecolor='none')
    plt.plot([], [], marker='h', color='tab:orange', alpha=0.8, label='Franco+25', markersize=13, zorder=2, linestyle='none')
    plt.legend(fontsize=17.5, ncol=3, bbox_to_anchor=(0.45, -0.15), loc='upper center')



plt.yscale('log')

if fit == False:
    yticks = [1e-6, 1e-5, 1e-4]
    plt.yticks(yticks)
    plt.gca().set_yticklabels([r"$-6\,$", r"$-5\,$", r"$-4\,$"])
if fit == True:
    yticks = [1e-9, 1e-8, 1e-7, 1e-6, 1e-5, 1e-4, 1e-3]
    plt.yticks(yticks)
    plt.gca().set_yticklabels([r"$-9\,$", r"$-8\,$", r"$-7\,$", r"$-6\,$", r"$-5\,$", r"$-4\,$", r"$-3\,$"])
    # Make axis labels larger
    plt.gca().tick_params(axis='both', which='major', labelsize=24)


plt.tight_layout()  # Leaves extra space at the bottom
#plt.savefig(plot_dir / 'LF_UVISTA.pdf')
#plt.savefig(plot_dir / 'LF_UVISTA_complete_emcee_finer_bins.pdf')
if run_type == '':
    plt.savefig(plot_dir / 'LF_UVISTA_newbdcut_points.pdf', bbox_inches='tight')
    #plt.savefig(plot_dir / 'LF_UVISTA_newbdcut_fit.pdf')
if run_type == 'with_euclid':
    if fit == False:
        plt.savefig(plot_dir / 'LF_UVISTA_newbdcut_with_euclid_points_with_CWEB_kron.pdf', bbox_inches='tight')
    #plt.savefig(plot_dir / 'LF_UVISTA_bologna_interview.pdf')
    if fit == True:
        plt.savefig(plot_dir / 'LF_UVISTA_newbdcut_with_euclid_fit_with_CWEB_kron.pdf', bbox_inches='tight')
        #plt.savefig(plot_dir / 'LF_UVISTA_newbdcut_with_euclid_UMASS_interview.pdf', bbox_inches='tight')
plt.show()


