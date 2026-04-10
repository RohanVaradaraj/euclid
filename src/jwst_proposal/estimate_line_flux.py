"""
Read in SEDs and estimate line fluxes given an equivalent width.

Created: Wednesday 8th October 2025.
"""

from astropy.table import Table
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
from astropy.constants import c
from scipy.optimize import curve_fit

# Paths & constants
sed_dir = Path.cwd().parent / 'sed_fitting' / 'pipes'

calzetti_windows = [
    (1268, 1284),
    (1309, 1316),
    (1342, 1371),
    (1435, 1496),
    (1562, 1583),
    (1677, 1740),
    (1760, 1830),
]

IDs = [178396, 387777, 615679, 878786, 879369]
redshifts = [7.1105, 6.7614, 6.8484, 7.0975, 7.0975]

# --- Line parameters ---
lam_rest = 1908.73  # Å (CIII] rest wavelength; average of doublet)
EW = 10             # Å (rest-frame equivalent width)
# Note: For a doublet, this can be refined per component

# NVλ1243
lam_rest = 1243
EW = 10

# lam_rest = 1216
# EW = 25

# Power-law continuum function
def powerlaw(lam, logA, beta):
    return 10**logA * lam**beta

for ID, z in zip(IDs, redshifts):

    # --- Load data ---
    t = Table.read(sed_dir / f"{ID}_microJy.fits")
    wlen_obs = t['wave']       # observed wavelength [Å]
    flux = t['flux']           # fν [μJy]

    # Convert to rest-frame
    lam_rest_arr = wlen_obs / (1 + z)

    # Mask base continuum region: 1250–2500 Å rest-frame
    mask_base = (lam_rest_arr > 1250) & (lam_rest_arr < 2500)

    # Mask Calzetti windows
    mask_cal = np.zeros_like(lam_rest_arr, dtype=bool)
    for wmin, wmax in calzetti_windows:
        mask_cal |= (lam_rest_arr >= wmin) & (lam_rest_arr <= wmax)

    mask = mask_base & ~mask_cal

    # Fit continuum in fν space
    popt, _ = curve_fit(powerlaw, lam_rest_arr[mask], flux[mask])
    logA, beta = popt

    # Observed wavelength of CIII]
    lam_obs = lam_rest * (1 + z)

    # fν at line position
    fnu_line = powerlaw(lam_rest, logA, beta)  # still in μJy

    print(f"ID {ID}: z={z:.4f}")
    # print(f"  Observed λ(CIII]) = {lam_obs:.1f} Å")
    # #print(f"  Continuum fν({lam_rest:.1f} Å rest) = {fnu_line:.4f} μJy")
    # # Print in erg/s/cm²/Hz
    # print(f"  Continuum fν({lam_rest:.1f} Å rest) = {fnu_line*1e-29:.4e} erg/s/cm²/Hz")
    print(f" Observed frame EW = {EW * (1 + z):.1f} Å")
    print(f" Observe lambda = {lam_obs:.1f} Å")

    f_line = fnu_line * 1e-29 * EW * (1 + z) * 2.998e18 / lam_obs**2  # erg/s/cm²
    print(f"  Estimated line flux = {f_line:.4e} erg/s/cm²")

    # Find width of line in Å assuming rest-frame FWHM = 200 km/s
    FWHM_rest = 200  # km/s
    print(f"  FWHM = {FWHM_rest*(1+z)} km/s")
    FWHM_obs = FWHM_rest * (1 + z) * lam_rest / 3e5  # Å
    sigma_obs = FWHM_obs / 2.355  # Å
    #print(f"  Observed FWHM = {FWHM_obs:.1f} Å (σ = {sigma_obs:.1f} Å)")