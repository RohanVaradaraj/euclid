"""
Set Spitzer errors to minimum of 5% to see if it makes a difference.

Created: Wednesday 3rd August 2025.
"""

from pathlib import Path
from astropy.io import fits
from astropy.table import Table
from astropy.nddata.utils import Cutout2D
from astropy.wcs import WCS
import sep
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Ellipse, Circle

t = Table.read(Path.cwd().parents[1] / 'data' / 'catalogues' / 'COSMOS_5sig_Y_J_nonDet_HSC_G_nonDet_HSC_R_nonDet_HSC_I.fits')

err1 = t['err_ch1cds']
err2 = t['err_ch2cds']
flux1 = t['flux_ch1cds']
flux2 = t['flux_ch2cds']

plt.hist(err1/flux1, bins=np.arange(0, 1, 0.01), label='Old ratio', histtype='step', color='blue', lw=2, alpha=0.5)

# Spitzer fluxes are under colnames flux_ch1cds, flux_ch2cds with errors err_ch1cds, err_ch2cds
# These errors are set to minimum 20% of the flux. Change to 5% of the flux.

# If the ratio of err/flux < 0.2, set err to 0.05 * flux
def set_spitzer_err_floor(table, flux_col, err_col, floor=0.05):
    """
    Set the error floor for Spitzer fluxes to a minimum percentage of the flux.
    
    Parameters:
    - table: Astropy Table containing the flux and error columns.
    - flux_col: Name of the column containing the flux values.
    - err_col: Name of the column containing the error values.
    - floor: Minimum percentage of the flux to set as the error floor (default is 0.05).
    
    Returns:
    - Updated table with modified error values.
    """
    mask = (table[err_col] / table[flux_col] < 0.21) & (table[err_col] / table[flux_col] > 0.)
    table[err_col][mask] = floor * table[flux_col][mask]
    return table

# Apply the error floor to both Spitzer channels
t = set_spitzer_err_floor(t, 'flux_ch1cds', 'err_ch1cds', floor=0.05)
t = set_spitzer_err_floor(t, 'flux_ch2cds', 'err_ch2cds', floor=0.05)

# Save table
t.write(Path.cwd().parents[1] / 'data' / 'catalogues' / 'COSMOS_5sig_Y_J_nonDet_HSC_G_nonDet_HSC_R_nonDet_HSC_I_5percent_IRACfloor.fits', overwrite=True)

# plt.hist(t['err_ch1cds']/t['flux_ch1cds'], bins=np.arange(0, 1, 0.01), label='New ratio', histtype='step', color='red', lw=2, alpha=0.5)
# plt.xlabel('Error / Flux')

# plt.legend()
# plt.show()
