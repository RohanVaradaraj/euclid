"""
estimate_sample_from_lf.py

Use Rebecca's LF determination to estimate the number of galaxies we should find in COSMOS.

Created: Tuesday 29th October 2024.
"""

import numpy as np
from scipy.integrate import quad
from astropy.cosmology import WMAP9 as cosmo

# DPL function parameters in units of mag^-1 Mpc^-3
Phi_star = 2.3e-4  # per mag per Mpc^3
M_star = -20.6     # Characteristic absolute magnitude
alpha = -2.19      # Faint-end slope
beta = -4.6        # Bright-end slope

# Define the DPL LF in terms of magnitude
def dpl_function(M, Phi_star, M_star, alpha, beta):
    term_faint = 10 ** (0.4 * (alpha + 1) * (M - M_star))
    term_bright = 10 ** (0.4 * (beta + 1) * (M - M_star))
    return Phi_star / (term_faint + term_bright)

# Field area conversion
A_sq_deg = 1.51  # Field area in square degrees
A_sr = A_sq_deg * (np.pi / 180) ** 2

# Redshift range and apparent magnitude limit
z_min, z_max = 6.5, 7.5
m_lim = 26.0  # Apparent magnitude limit

# Comoving volume element in Mpc^3 per steradian per dz
def dV_dz(z):
    return cosmo.differential_comoving_volume(z).value

# Calculate absolute magnitude limit at each redshift
def absolute_magnitude_limit(z, m_lim):
    distance_modulus = cosmo.distmod(z).value
    return m_lim - distance_modulus + + 2.5*np.log10(1+z)

# Main function to estimate expected number of galaxies
def expected_number_galaxies_dpl():
    # Discretize redshift range
    z_values = np.linspace(z_min, z_max, 50)
    total_galaxies = 0

    for z in z_values:
        M_lim = absolute_magnitude_limit(z, m_lim)
        # Integrate DPL up to M_lim, from faint limit -17
        phi_integral, _ = quad(dpl_function, -25, M_lim, args=(Phi_star, M_star, alpha, beta))
        # Sum over volume in each redshift slice
        total_galaxies += phi_integral * dV_dz(z) * A_sr * (z_values[1] - z_values[0])

    return total_galaxies

# Calculate expected number of galaxies
N_galaxies = expected_number_galaxies_dpl()
print(f"Expected number of galaxies: {N_galaxies:.2f}")
