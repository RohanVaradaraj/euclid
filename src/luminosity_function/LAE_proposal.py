"""
Computing number counts for LAEs for Bologna job.

Created: Tuesday 4th Febraury 2025.
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
def dpl_function(L, Phi_star, L_star, alpha, beta):
    term_faint = 10 ** (0.4 * (alpha + 1) * (L - L_star))
    term_bright = 10 ** (0.4 * (beta + 1) * (L - L_star))
    return Phi_star / (term_faint + term_bright)

L = np.arange(41, 45, 0.1)


LAE_LF_z7p3 = dpl_function(L, 10**(-5.03), 42.98, -2.51, -5.02)
LAE_LF_z7p0 = dpl_function(L, 10**(-4.49), 43.23, -2.68, -4.51)

plt.plot(L, LAE_LF_z7p0, color='green', lw=4)
plt.plot(L, LAE_LF_z7p3, color='deepskyblue', lw=4)

plt.legend()
plt.xlabel(r"$\rm{log}L_{\alpha} \ \rm{[erg s]}^{-1}$")
plt.ylabel(r'$\Phi(L_{\alpha}) \ [N \  / \ \Delta\rm{log}L_{\alpha} \ / \ \mathrm{cMpc}^{3}]$')
plt.tight_layout()
plt.show()