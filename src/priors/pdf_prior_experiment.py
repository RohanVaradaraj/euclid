"""
Experimenting with placing a luminosity function prior on the P(z) of an SED.

Just trying on one object or a handful of objects for now.

Created: Wednesday 14th May 2025.
"""
from pathlib import Path
import sys
sed_fitting_path = Path.cwd().parents[0] / 'sed_fitting'
sys.path.append(str(sed_fitting_path))

from sed_fitting_codes import parse_spec_file
import matplotlib.pyplot as plt
import glob
import numpy as np
from astropy.table import Table
from astropy.cosmology import FlatLambdaCDM
import matplotlib.cm as cm
import matplotlib.colors as mcolors

plt.rcParams['axes.linewidth'] = 2.5
plt.rcParams.update({'font.size': 15})
plt.rcParams['figure.dpi'] = 100

#! SWITCH TO PLOT P(z) FOR ALL OBJECTS, TO TAKE A LOOK
plot_all_pdfs = False

#! INTERESTING IDs
IDs = [376795, 905849, 1053436]

sed_dir = Path.cwd().parents[1] / 'data' / 'sed_fitting' / 'zphot' / 'XMM' / 'best_fits' / 'det_HSC-Z_DR3_z7'
spec_files = glob.glob(str(sed_dir / '*.spec'))

#! Define a cosmological model
H = 70
omegaM = 0.3
omegaV = 0.7
cosmo = FlatLambdaCDM(H0=H, Om0=omegaM)

#! Bouwens+21 LF parameter evolution
def alpha_evolution(z):
    """
    Evolution of the faint end slope of the LF.
    """
    return -1.94 - 0.11 * (z - 6)

def phi_star_evolution(z):
    """
    Evolution of the characteristic value of the LF.
    """
    exponent = -0.33 * (z - 6) - 0.024 * (z - 6)**2
    return 0.4e-3 * 10 ** exponent

def Mstar_evolution(z):
    """
    Evolution of the characteristic magnitude of the LF.
    """
    z_t = 2.46
    return np.where(
        z < z_t,
        -20.89 - 1.09 * (z - z_t),
        -21.03 - 0.04 * (z - z_t)
    )


#! Schechter function from the parameters
def Schecter_function(M, alpha, phi_star, M_star):
    """
    Produces Schechter function given parameters: normalisation phi*, faint slope alpha, mag array M and char. mag M*. 
    """

    coeff = np.log(10) / 2.5

    faint = (10 ** (0.4 * (M_star - M))) ** (alpha+1)


    bright_exponent = -10 ** (0.4 * (M_star - M))
    bright = np.exp(bright_exponent)

    phi = coeff * phi_star * faint * bright

    return phi


#! Double schechter function GSMF
def Double_Schechter(M, M_star, alpha1, alpha2, phi_star1, phi_star2):
    """
    Produces double Schechter function given parameters: normalisation phi*, faint slope alpha, mag array M and char. mag M*. 
    """

    exp_term = np.exp( -10 ** (M - M_star) ) * 10 ** (M - M_star)
    
    high_mass = phi_star1 * 10 ** ( (M - M_star) * alpha1 )

    low_mass = phi_star2 * 10 ** ( (M - M_star) * alpha2 )

    phi = np.log(10) * exp_term * ( high_mass + low_mass )

    return phi


#! McLeod+21 GSMF parameter evolution
z_mcleod = [0.5, 1, 1.5, 2, 2.5, 3.25]

M_star_mcleod = [10.8, 10.72, 10.72, 10.77, 10.77, 10.84]
dM_star_mcleod = [0.06, 0.07, 0.05, 0.06, 0.10, 0.18]

log_phi1_mcleod = [-2.77, -2.80, -2.94, -3.18, -3.39, -4.30]
dlog_phi1_mcleod = [0.07, 0.09, 0.05, 0.08, 0.11, 0.23]

alpha1_mcleod = [-0.61, -0.46, -0.55, -0.68, -0.62, -0.00]
dalpha1_mcleod = [0.23, 0.34, 0.22, 0.29, 0.50, 1.03]

log_phi2_mcleod = [-3.26, -3.26, -3.54, -3.84, -3.78, -3.94]
dlog_phi2_mcleod = [0.17, 0.23, 0.22, 0.46, 0.50, 0.37]

alpha2_mcleod = [-1.52, -1.53, -1.65, -1.73, -1.74, -1.79]
dalpha2_mcleod = [0.05, 0.07, 0.07, 0.12, 0.13, 0.09]

# Fit quadratic polynomials to the data
coeffs_M_star = np.polyfit(z_mcleod, M_star_mcleod, 2)
coeffs_log_phi1 = np.polyfit(z_mcleod, log_phi1_mcleod, 2) 
coeffs_alpha1 = np.polyfit(z_mcleod, alpha1_mcleod, 2)
coeffs_log_phi2 = np.polyfit(z_mcleod, log_phi2_mcleod, 2)
coeffs_alpha2 = np.polyfit(z_mcleod, alpha2_mcleod, 2)

# Define a function that evaluates the GSMF at z
def McLeod_GSMF(z):
    alpha1 = np.polyval(coeffs_alpha1, z)
    phi_star1 = 10 ** np.polyval(coeffs_log_phi1, z)
    M_star = np.polyval(coeffs_M_star, z)
    alpha2 = np.polyval(coeffs_alpha2, z)
    phi_star2 = 10 ** np.polyval(coeffs_log_phi2, z)

    M = np.arange(7.5, 12, 0.1)
    gsmf = Double_Schechter(M, M_star, alpha1, alpha2, phi_star1, phi_star2)
    return M, gsmf

#! Convert GSMF to rest-UV LF using M/L ratio of 0.5
def GSMF_to_LF(M, gsmf):
    """
    Convert GSMF to LF using M/L ratio of 0.5.
    """
    # Convert GSMF to LF
    # M/L = 0.5
    M_L = 0.5
    L = gsmf * M_L

    # Convert to absolute magnitude
    M_abs = -19 - 2.5 * (M - 9)

    return M_abs, L

# Plot the resulting LF
z_grid = np.linspace(0.25, 3.75, 100)
cmap = cm.coolwarm
norm = mcolors.Normalize(vmin=min(z_grid), vmax=max(z_grid))

masses = np.arange(7.5, 12, 0.1)
for z in z_grid:
    M, gsmf = McLeod_GSMF(z)
    M_abs, LF = GSMF_to_LF(masses, gsmf)
    plt.plot(M_abs, LF, color=cmap(norm(z)), linewidth=2)

plt.yscale('log')
plt.show()





# Create a grid of z values for plotting
# z_grid = np.linspace(0.25, 3.75, 100)

# cmap = cm.coolwarm
# norm = mcolors.Normalize(vmin=min(z_grid), vmax=max(z_grid))

# # Plot the GSMF at each z of the z_grid
# for z in z_grid:
#     alpha1 = np.polyval(coeffs_alpha1, z)
#     phi_star1 = 10 ** np.polyval(coeffs_log_phi1, z)
#     M_star = np.polyval(coeffs_M_star, z)
#     alpha2 = np.polyval(coeffs_alpha2, z)
#     phi_star2 = 10 ** np.polyval(coeffs_log_phi2, z)

#     M = np.arange(7.5, 12, 0.1)
#     gsmf = Double_Schechter(M, M_star, alpha1, alpha2, phi_star1, phi_star2)
#     plt.plot(M, gsmf, color=cmap(norm(z)), linewidth=2)

# plt.yscale('log')

# # Colorbar legend

# plt.xlabel(r'$\mathrm{log}(M^{*}/M_{\odot})$')
# plt.ylabel(r'$\Phi$')


# plt.ylim(1e-7, 0)
# plt.xlim(7., 12.5)
# plt.show()

# # Plot all of these in subplots
# fig, axs = plt.subplots(3, 2, figsize=(10, 10))
# axs[0, 0].errorbar(z_mcleod, M_star_mcleod, yerr=dM_star_mcleod, fmt='o', color='purple', label='M*')
# axs[0, 0].set_xlabel('z')
# axs[0, 0].set_ylabel('M*')

# axs[0, 1].errorbar(z_mcleod, log_phi1_mcleod, yerr=dlog_phi1_mcleod, fmt='o', color='purple', label='log phi1')
# axs[0, 1].set_xlabel('z')
# axs[0, 1].set_ylabel('log phi1')

# axs[1, 0].errorbar(z_mcleod, alpha1_mcleod, yerr=dalpha1_mcleod, fmt='o', color='purple', label='alpha1')
# axs[1, 0].set_xlabel('z')
# axs[1, 0].set_ylabel('alpha1')

# axs[1, 1].errorbar(z_mcleod, log_phi2_mcleod, yerr=dlog_phi2_mcleod, fmt='o', color='purple', label='log phi2')
# axs[1, 1].set_xlabel('z')
# axs[1, 1].set_ylabel('log phi2')

# axs[2, 0].errorbar(z_mcleod, alpha2_mcleod, yerr=dalpha2_mcleod, fmt='o', color='purple', label='alpha2')
# axs[2, 0].set_xlabel('z')
# axs[2, 0].set_ylabel('alpha2')
# axs[2, 1].axis('off')  # Hide the last subplot

# # Plot the fitted curves
# axs[0, 0].plot(z_grid, np.polyval(coeffs_M_star, z_grid), color='purple', linewidth=2)
# axs[0, 1].plot(z_grid, np.polyval(coeffs_log_phi1, z_grid), color='purple', linewidth=2)
# axs[1, 0].plot(z_grid, np.polyval(coeffs_alpha1, z_grid), color='purple', linewidth=2)
# axs[1, 1].plot(z_grid, np.polyval(coeffs_log_phi2, z_grid), color='purple', linewidth=2)
# axs[2, 0].plot(z_grid, np.polyval(coeffs_alpha2, z_grid), color='purple', linewidth=2)

exit()

#! Experimenting on an ID
ID = IDs[0]
ID = 376795

#! Open the parent catalogue ID
cat_dir = Path.cwd().parents[1] / 'data' / 'catalogues' / 'candidates'
parent_cat_name = 'XMM_5sig_HSC_Z_nonDet_HSC_G_nonDet_HSC_R_candidates_2025_05_14.fits'
parent_cat = Table.read(cat_dir / parent_cat_name)

# Get Muv of this ID
ID_row = parent_cat['ID'] == ID
Muv = parent_cat['Muv'][ID_row][0]


# pAD ID with zeros to left so that there are nine digits
ID_str = str(ID).zfill(9)
spec_file = f'{sed_dir}/Id{ID_str}.spec'
file = parse_spec_file(spec_file)

phot = file.get('phot')
zpdf = file.get('zpdf')
sed = file.get('sed')
params = file.get('model')

plt.figure(figsize=(10, 6))
plt.plot(zpdf['z'], zpdf['P(z)'], color='purple', linewidth=3, label='Flat prior')


# Open the .chi file as ascii commented header
chi_file = spec_file.replace('.spec', '.chi').split('/')[-1]
chi = Table.read(sed_dir.parents[1] / 'det_HSC-Z_DR3' / chi_file, format='ascii.commented_header')

#! Compute the comoving volume at each redshift in z_grid increment
comoving_volume = cosmo.comoving_volume(zpdf['z']).value

#! Compute LF * comoving volume as function of redshift
M_values = np.arange(-23, -19, 0.1)
cmap = cm.coolwarm
norm = mcolors.Normalize(vmin=min(M_values), vmax=max(M_values))

# # Create figure and axis
# fig, ax = plt.subplots(figsize=(8, 5))

# Compute the comoving volume in the 0.01 interval of each z_grid value
del_comoving_volume = np.diff(cosmo.comoving_volume(zpdf['z']).value)
del_comoving_volume = np.append(del_comoving_volume, del_comoving_volume[-1])  # Append the last value to match the length of z_grid

# for M in M_values:
#     alpha = alpha_evolution(z_grid)
#     phi_star = phi_star_evolution(z_grid)
#     M_star = Mstar_evolution(z_grid)
    
#     schechter = Schecter_function(M, alpha, phi_star, M_star)
#     #schechter_volume = schechter * del_comoving_volume
#     schechter_volume = schechter * comoving_volume
#     schechter_volume /= np.max(schechter_volume) #, z_grid)

#     color = cmap(norm(M))
#     ax.plot(z_grid, schechter_volume, linewidth=2.5, color=color)

# # Add colorbar
# sm = cm.ScalarMappable(cmap=cmap, norm=norm)
# sm.set_array([])
# cbar = fig.colorbar(sm, ax=ax)
# cbar.set_label(r'$M_{\mathrm{UV}}$')

# ax.set_xlabel(r'$z$')
# #ax.set_ylabel(r'$\Phi \times \frac{dV_{c}}{dz}$')
# ax.set_ylabel(r'$\Phi \times V_{c}$')
# fig.tight_layout()
# plt.show()

#plt.plot(zpdf['z'], del_comoving_volume, color='blue', linewidth=3, label='Delta Comoving volume')
#plt.plot(zpdf['z'], comoving_volume, color='green', linewidth=3, label='Comoving volume')
# plt.show()
# plt.clf()
# plt.close()

alpha = alpha_evolution(zpdf['z'])
phi_star = phi_star_evolution(zpdf['z'])
M_star = Mstar_evolution(zpdf['z'])

schechter = Schecter_function(Muv, alpha, phi_star, M_star)
schechter_volume = schechter * del_comoving_volume
#schechter_volume = zpdf['z'] * np.exp(-zpdf['z'])

# Use gsmf for z<4
#schechter_volume[zpdf['z'] < 4] = 


# Multiply P(z) by the Schechter volume prior
pdf_with_prior = zpdf['P(z)'] * schechter_volume

# Normalize the PDF with prior
pdf_with_prior /= np.max(pdf_with_prior)

# Plot the PDF with prior
plt.plot(zpdf['z'], pdf_with_prior, color='red', linewidth=3, label='LF prior')
plt.xlabel('z')
plt.ylabel('P(z)')
plt.legend()
plt.show()


#? Independently computing P(z) from chi2 values
# # Extract redshifts and chi2 values
# z_vals = np.unique(chi['Z'])  # shape (N_z,), sorted
# n_z = len(z_vals)

# # Check total number of rows to get number of models
# n_models = len(chi) // n_z

# # Reshape Chi2 values into (n_models, n_z)
# chi2_grid = chi['Chi2'].data.reshape((n_models, n_z))

# # Get the minimum Chi2 at each redshift over all models
# chi2_min_z = np.min(chi2_grid, axis=0)

# # Compute unnormalized P(z)
# ln_likelihood = -0.5 * chi2_min_z
# p_z_unnorm = np.exp(ln_likelihood - np.max(ln_likelihood))  # subtract max for numerical stability

# # Normalize P(z)
# #p_z = p_z_unnorm / np.trapz(p_z_unnorm, z_vals)  # Trapezoidal integration over redshift
# p_z = p_z_unnorm / np.max(p_z_unnorm)  # Normalize to max value

# # Plot P(z)
# plt.plot(z_vals, p_z, label='P(z)', color='red', linewidth=2)
# plt.xlabel('Redshift z')
# plt.ylabel('Probability')
# plt.title('Photometric Redshift Probability Distribution P(z)')
# plt.grid(True)
# plt.legend()
# plt.show()

#? Plot the LF parameter evolution
# Plot the LF evolution, each on own panel
# if not plot_all_pdfs:
#     fig, axs = plt.subplots(1, 3, figsize=(15, 5))

#     # Faint end slope
#     axs[0].plot(z_grid, alpha_evolution(z_grid), color='purple', linewidth=3)
#     axs[0].set_xlabel('z')
#     axs[0].set_ylabel(r'$\alpha$')

#     # Characteristic value
#     axs[1].plot(z_grid, phi_star_evolution(z_grid), color='purple', linewidth=3)
#     axs[1].set_xlabel('z')
#     axs[1].set_ylabel(r'$\Phi^*$')
#     # log y axis
#     axs[1].set_yscale('log')

#     # Characteristic magnitude
#     axs[2].plot(z_grid, Mstar_evolution(z_grid), color='purple', linewidth=3)
#     axs[2].set_xlabel('z')
#     axs[2].set_ylabel(r'$M^*$')
#     # Flip y axis
#     axs[2].invert_yaxis()

#     # Increase tick size
#     for ax in axs:
#         ax.tick_params(axis='both', which='major', labelsize=15)
#         ax.tick_params(axis='both', which='minor', labelsize=15)

#     plt.tight_layout()
#     plt.show()

#? Plot the PDFs
# if plot_all_pdfs:
#     for spec_file in spec_files:

#         file = parse_spec_file(spec_file)
#         ID = spec_file.split('/')[-1].split('Id')[-1].lstrip('0').split('.spec')[0]

#         zpdf = file.get('zpdf')

#         # Plot P(z) for the object
#         plt.figure(figsize=(10, 6))
#         plt.plot(zpdf['z'], zpdf['P(z)'], color='purple', linewidth=3)
#         plt.xlabel('z')
#         plt.ylabel('P(z)')
#         plt.title(ID)
#         plt.show()
#         plt.close()
