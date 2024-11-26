"""
Draw different euclid fields as boxes on the z=3-9 LFs

Created: Wednesday 30th October 2024.
"""

import numpy as np
from scipy.integrate import quad
from astropy.cosmology import WMAP9 as cosmo
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.colors import Normalize
from pathlib import Path

plt.rcParams['axes.linewidth'] = 2.5
plt.rcParams.update({'font.size': 15})
plt.rcParams['figure.dpi'] = 100

# Define the DPL LF in terms of magnitude
def dpl_function(M, Phi_star, M_star, alpha, beta):
    term_faint = 10 ** (0.4 * (alpha + 1) * (M - M_star))
    term_bright = 10 ** (0.4 * (beta + 1) * (M - M_star))
    return Phi_star / (term_faint + term_bright)

# Comoving volume element in Mpc^3 per steradian per dz
def dV_dz(z):
    return cosmo.differential_comoving_volume(z).value

# Calculate absolute magnitude limit at each redshift
def absolute_magnitude_limit(z, m_lim):
    distance_modulus = cosmo.distmod(z).value
    return m_lim - distance_modulus + 2.5 * np.log10(1 + z)

# Function to find limits of a survey box
def survey_box_limits(m_lim, area, z_min, z_max):
    M_lim = absolute_magnitude_limit(np.mean([z_min, z_max]), m_lim)

    # Convert area to steradians
    area_sr = area * (np.pi / 180) ** 2
    Omega = area_sr

    # Comoving volume in Mpc^3 between z_min and z_max per steradian
    volume_per_sr = cosmo.comoving_volume(z_max).value - cosmo.comoving_volume(z_min).value

    # Total comoving volume for the survey area
    V_survey = area_sr * volume_per_sr

    # Minimum detectable number density (1 galaxy per survey volume)
    Phi_min = 1 / V_survey  # in galaxies per Mpc^3

    dmin = cosmo.comoving_distance(z_min).value
    dmax = cosmo.comoving_distance(z_max).value

    V = Omega/3 * (dmax**3 - dmin**3)

    Phi_min = 1 / V

    return M_lim, Phi_min



M = np.arange(-35, -10, 0.1)

# Adams+23
z3_gal = dpl_function(M, 10**(-3.20), -21.18, -1.85, -4.95)
z3_agn = dpl_function(M, 10**(-6.57), -25.57, -1.37, -4.91)

z4_gal = dpl_function(M, 10**(-3.62), -21.68, -2.10, -5.29)
z4_agn = dpl_function(M, 10**(-7.77), -27.18, -2.02, -4.34)

z5_gal = dpl_function(M, 10**(-3.57), -21.60, -1.94, -5.73)
z5_agn = dpl_function(M, 10**(-8.65), -27.58, -2.12, -6.31)

z3_lf = z3_gal + z3_agn
z4_lf = z4_gal + z4_agn
z5_lf = z5_gal + z5_agn

# Bowler+15
z6_gal = dpl_function(M, 1.9*10**(-4.), -21.20, -2.10, -5.10)

# Schindler+23
z6_agn = dpl_function(M, 10**(-8.75), -26.38, -1.70, -3.84)

z6_lf = z6_gal + z6_agn

# Bowler+17
z7_gal = dpl_function(M, 2.3*10**(-4.), -20.60, -2.19, -4.60)

# Matsuoka+23
z7_agn = dpl_function(M, 1.3*10**(-9.), -25.6, -1.20, -3.34)

z7_lf = z7_gal + z7_agn

# Donnan+22
z8_lf = dpl_function(M, 3.3*10**(-4.), -20.02, -2.04, -4.26)
z9_lf = dpl_function(M, 2.1*10**(-4.), -19.93, -2.10, -4.29)

# List of redshift values and LFs
redshifts = [3, 4, 5, 6, 7, 8, 9]
lfs = [z3_lf, z4_lf, z5_lf, z6_lf, z7_lf, z8_lf, z9_lf]

# Normalize the colormap for the redshift range
norm = Normalize(vmin=min(redshifts), vmax=max(redshifts))
cmap = cm.viridis

# Plot the LFs with the spectral color map
fig, ax = plt.subplots(figsize=(8, 6))
for z, lf in zip(redshifts, lfs):
    plt.plot(M, lf, color=cmap(norm(z)), lw=4)

# Varadaraj+23 LF points
bins_mag = np.array([-22.175, -22.925, -23.675, -24.425])
phi_values = np.array([2.70e-6, 2.81e-7, 2.37e-8, np.nan])  
phi_errors = np.array([0.66e-6, 1.54e-7, 2.50e-8, np.nan])  

plt.errorbar(bins_mag[:3], phi_values[:3], yerr=phi_errors[:3], fmt='o', color='red', 
             ecolor='red', elinewidth=5, label=r'Varadaraj+23, $z=7$', markersize=17, markeredgecolor='black', zorder=10)

# plt.errorbar(bins_mag[3], 2.07e-8, yerr=0, marker=r'$\downarrow$', color='red', ecolor='red', 
#              elinewidth=5, markersize=17)

#! Survey boxes

# # Eventual full Euclid auxiliary fields
# M_lim, Phi_min = survey_box_limits(26.3, 4.5, 6.5, 7.5)
# plt.plot([M_lim, M_lim], [Phi_min, 1], color='gray', lw=2.5, linestyle='dashed') # vertical line
# plt.plot([M_lim, -30], [Phi_min, Phi_min], color='gray', lw=2.5, linestyle='dashed') # horizontal line

# EDFs, eventually
M_lim, Phi_min = survey_box_limits(26, 59, 6.5, 7.5)
plt.plot([M_lim, M_lim], [Phi_min, 1], color='black', lw=3.5, alpha=0.9)
plt.plot([M_lim, -30], [Phi_min, Phi_min], color='black', lw=3.5, alpha=0.9)
#plt.text(-20.76, 4e-9, 'Final EDFs', color='gray')
plt.text(-20.76, 1e-9, 'EDFs', color='black')


# COSMOS
# M_lim, Phi_min = survey_box_limits(26.4, 4.5, 6.5, 7.5)
# plt.plot([M_lim, M_lim], [Phi_min, 1], 'k', lw=3.5, alpha=0.9) # vertical line
# plt.plot([M_lim, -30], [Phi_min, Phi_min], 'k', lw=3.5, alpha=0.9) # horizontal line
# plt.text(-20.4, 3.1e-8, 'EAFs')

# EDFs, Q1
# M_lim, Phi_min = survey_box_limits(24, 43, 6.5, 7.5)
# plt.plot([M_lim, M_lim], [Phi_min, 1], 'k', lw=3.5, alpha=0.9) # vertical line
# plt.plot([M_lim, -30], [Phi_min, Phi_min], 'k', lw=3.5, alpha=0.9) # horizontal line
# plt.text(-23.5, 5.4e-10, 'EDFs, Q1')

# CWEB
# M_lim, Phi_min = survey_box_limits(28.2, 0.54, 6.5, 7.5)
# plt.plot([M_lim, M_lim], [Phi_min, 1], color='orange', lw=3.5, alpha=0.8, zorder=3) #, linestyle='dotted') # lw = 2 for Euclid version
# plt.plot([M_lim, -30], [Phi_min, Phi_min], color='orange', lw=3.5, alpha=0.8, zorder=3) #, label = 'COSMOS-Web, 0.54 deg' + r'$^{2}$', linestyle='dotted')
# plt.text(-18.7, 1e-5, 'COSMOS-Web \n0.54 deg' + r'$^{2}$', color='orange', fontsize=15)

# CWEB + MIRI
M_lim, Phi_min = survey_box_limits(28.2, 0.19, 6.5, 7.5)
plt.plot([M_lim, M_lim], [Phi_min, 1], color='orange', lw=3.5, alpha=0.8, zorder=3) #, linestyle='dotted') # lw = 2 for Euclid version
plt.plot([M_lim, -30], [Phi_min, Phi_min], color='orange', lw=3.5, alpha=0.8, zorder=3) #, label = 'COSMOS-Web, 0.54 deg' + r'$^{2}$', linestyle='dotted')
plt.text(-18.7, 2e-8, 'COSMOS-Web \n +MIRI \n0.19 deg' + r'$^{2}$', color='orange', fontsize=15)

# UltraVISTA + VIDEO
M_lim, Phi_min = survey_box_limits(26.2, 10., 6.5, 7.5)
plt.plot([M_lim, M_lim], [Phi_min, 1], color='deepskyblue', lw=3.5, alpha=0.8, zorder=3) #, linestyle='dotted')
plt.plot([M_lim, -30], [Phi_min, Phi_min], color='deepskyblue', lw=3.5, alpha=0.8, zorder=3) #, linestyle='dotted', label='VISTA, 10 deg' + r'$^{2}$')
plt.text(-20.5, 1e-8, 'VISTA \n10 deg' + r'$^{2}$', color='deepskyblue', fontsize=15)

# MIDIS
M_lim, Phi_min = survey_box_limits(29.49, 24.9 / 3600, 6.5, 7.5)
plt.plot([M_lim, M_lim], [Phi_min, 1], color='maroon', lw=3.5, alpha=0.8, zorder=3) #, linestyle='dotted') # lw = 2 for Euclid version
plt.plot([M_lim, -30], [Phi_min, Phi_min], color='maroon', lw=3.5, alpha=0.8, zorder=3) #, label = 'COSMOS-Web, 0.54 deg' + r'$^{2}$', linestyle='dotted')
plt.text(-18.7, 2e-6, 'JADES \n24.9 arcmin' + r'$^{2}$', color='maroon', fontsize=15)

plt.xlim(-25.5, -16.5)
plt.yscale('log')
plt.ylim(1e-11, 5e-2)

plt.legend()
plt.xlabel(r"$M_{\mathrm{UV}}$")
plt.ylabel(r'$\Phi \ [N \  / \ \mathrm{mag} \ / \ \mathrm{Mpc}^{3}]$')
plt.tight_layout()

plot_dir = Path.cwd().parents[1] / 'plots' / 'LF'
plt.savefig(plot_dir / 'LF_STScI.pdf')
#plt.show()
