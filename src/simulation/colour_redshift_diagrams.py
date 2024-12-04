"""
colour_redshift_diagrams.py

plot how colour changes with redshift
"""

import matplotlib.pyplot as plt
from astropy.table import Table
from pathlib import Path
import numpy as np
from brown_dwarf_colours import *

# Set up the grid
redshifts = np.arange(5.5, 10.01, 0.01)
EWs = np.arange(0, 250, 10)

Avs = np.arange(0, 0.7, 0.1)
ages = np.array([50, 100, 150, 200, 300, 400, 500])

# Broad redshift range of interest
redshift = '7'

z_lower = 4.9
z_upper = 10.0

table_dir = Path.cwd().parent.parent / 'data' / 'simulations' / 'LAEs' / 'tables'

# Subplots to do colour1 and colour2 under the same figure
fig, ax = plt.subplots(2, 1, figsize=(10, 8))

# Use common shared x axis
fig.subplots_adjust(hspace=0)
plt.setp(ax[0].get_xticklabels(), visible=False)

# Define color maps
cmap = plt.get_cmap('viridis')

for Av in Avs:
    for age in ages:

        # Read in the tables
        fluxes = Table.read(table_dir / f'LAE_fluxes_Av_{round(Av, 1)}_age_{age}Myr.fits')

        # Get magnitudes
        mag_i = flux_to_mag(fluxes['i'])
        mag_Ye = flux_to_mag(fluxes['Ye'])
        mag_Je = flux_to_mag(fluxes['Je'])
        mag_y = flux_to_mag(fluxes['y'])
        mag_Y = flux_to_mag(fluxes['Y'])
        mag_J = flux_to_mag(fluxes['J'])
        mag_Je = flux_to_mag(fluxes['Je'])
        mag_He = flux_to_mag(fluxes['He'])
        mag_VIS = flux_to_mag(fluxes['VIS'])
        mag_z = flux_to_mag(fluxes['z'])
        mag_nb921 = flux_to_mag(fluxes['nb921'])    

        # Calculate colors

        if redshift == '6':
            y_colour_1 = mag_i - mag_VIS
            y_colour_2 = mag_VIS - mag_z

        if redshift == '7':
            y_colour_1 = mag_y - mag_Ye
            y_colour_2 = mag_Y - mag_Ye

        if redshift == '8':
            y_colour_1 = mag_Ye - mag_Je
            y_colour_2 = mag_J - mag_Je

        ###########! HIGH EW LAEs ###########
        # Limit to desired redshift range and EW
        redshift_mask = (fluxes['redshift'] > z_lower) & (fluxes['redshift'] < z_upper)
        high_EW = redshift_mask & (fluxes['EW'] > 80)

        # Find the most extreme values at each redshift
        unique_redshifts = np.unique(fluxes['redshift'][high_EW])
        max_y_colour_1 = np.zeros_like(unique_redshifts)
        min_y_colour_1 = np.zeros_like(unique_redshifts)
        max_y_colour_2 = np.zeros_like(unique_redshifts)
        min_y_colour_2 = np.zeros_like(unique_redshifts)

        for i, z in enumerate(unique_redshifts):
            idx = fluxes['redshift'][high_EW] == z
            max_y_colour_1[i] = np.max(y_colour_1[high_EW][idx])
            min_y_colour_1[i] = np.min(y_colour_1[high_EW][idx])
            max_y_colour_2[i] = np.max(y_colour_2[high_EW][idx])
            min_y_colour_2[i] = np.min(y_colour_2[high_EW][idx])

        # Plot filled areas
        ax[0].fill_between(unique_redshifts, min_y_colour_1, max_y_colour_1, color='deepskyblue', alpha=0.7)
        ax[1].fill_between(unique_redshifts, min_y_colour_2, max_y_colour_2, color='deepskyblue', alpha=0.7)

        ###########! NON-EMITTING LBGs ###########
        non_emit = redshift_mask & (fluxes['EW'] < 5)

        # Find the most extreme values at each redshift
        unique_redshifts = np.unique(fluxes['redshift'][non_emit])
        max_y_colour_1 = np.zeros_like(unique_redshifts)
        min_y_colour_1 = np.zeros_like(unique_redshifts)
        max_y_colour_2 = np.zeros_like(unique_redshifts)
        min_y_colour_2 = np.zeros_like(unique_redshifts)

        for i, z in enumerate(unique_redshifts):
            idx = fluxes['redshift'][non_emit] == z
            max_y_colour_1[i] = np.max(y_colour_1[non_emit][idx])
            min_y_colour_1[i] = np.min(y_colour_1[non_emit][idx])
            max_y_colour_2[i] = np.max(y_colour_2[non_emit][idx])
            min_y_colour_2[i] = np.min(y_colour_2[non_emit][idx])

        # Plot filled areas
        ax[0].fill_between(unique_redshifts, min_y_colour_1, max_y_colour_1, color='gray', alpha=0.7)
        ax[1].fill_between(unique_redshifts, min_y_colour_2, max_y_colour_2, color='gray', alpha=0.7)

# Set plot limits and labels
for i in range(2):

    if redshift == '8':
        ax[i].set_xlim(8.2, 9.)
        ax[i].set_ylim(-0.5, 2)
    if redshift == '7':
        ax[i].set_xlim(6.5, 7.6)
        ax[i].set_ylim(-2, 2)
    if redshift == '6':
        ax[i].set_xlim(5.5, 6.5)
        ax[i].set_ylim(-1, 3)



    ax[i].set_xlabel(r'$z$')

    if redshift == '8':
        ax[i].set_ylabel(['$Y_{E} - J_{E}$', '$J_{\mathrm{VISTA}} - J_{E}$'][i])
    if redshift == '7':
        ax[i].set_ylabel(['$y_{\mathrm{HSC}} - Y_{E}$', '$Y_{\mathrm{VISTA}} - Y_{E}$'][i])
    if redshift == '6':
        ax[i].set_ylabel(['$i_{\mathrm{HSC}} - VIS$', '$VIS - z_{\mathrm{HSC}}$'][i])

  
ax[0].plot([], [], color='deepskyblue', label=r'EW > 80 $\AA$')
ax[0].plot([], [], color='gray', label=r'No Lyman-$\alpha$ emission')
ax[0].legend(loc='upper left')
plot_dir = Path.cwd().parent.parent / 'plots' / 'LAEs'
#plt.savefig(plot_dir / f'colour_redshift_evolution_z{redshift}.pdf')

# Add the LAE object
lae_z = 7.19309
lae_dz_sup = 7.30626
lae_dz_inf = 7.03676
dz_sup = lae_dz_sup - lae_z
dz_inf = lae_z - lae_dz_inf

# get colours
cat_dir = Path.cwd().parent.parent / 'data' / 'catalogues'
t = Table.read(cat_dir / 'COSMOS_5sig_Ye_2sig_Y_nonDet_HSC_G_nonDet_HSC_R_nonDet_HSC_I.fits')

t = t[t['ID'] == 178396]

# Fluxes
y_flux = t['flux_HSC-Y_DR3']
Ye_flux = t['flux_Ye']
Y_flux = t['flux_Y']

y_err = t['err_HSC-Y_DR3']
Ye_err = t['err_Ye']
Y_err = t['err_Y']

# Flux to magnitude conversion
def flux_to_mag(flux):
    return -2.5 * np.log10(flux)

# Magnitude errors from flux and flux errors
def mag_error(flux, flux_err):
    return (2.5 / np.log(10)) * (flux_err / flux)

# Calculate magnitudes
y_mag = flux_to_mag(y_flux)
Ye_mag = flux_to_mag(Ye_flux)
Y_mag = flux_to_mag(Y_flux)

# Calculate magnitude errors
y_mag_err = mag_error(y_flux, y_err)
Ye_mag_err = mag_error(Ye_flux, Ye_err)

# Magnitude difference and its error
y_minus_Ye = y_mag - Ye_mag
Y_minus_Ye = Y_mag - Ye_mag

# Error on the difference
Y_minus_Ye_err = np.sqrt(Ye_mag_err**2 + Ye_mag_err**2)
y_minus_Ye_err = np.sqrt(y_mag_err**2 + Ye_mag_err**2)

# Plot the LAE object
ax[0].errorbar(lae_z, y_minus_Ye, yerr=y_minus_Ye_err, xerr=[[dz_inf], [dz_sup]], fmt='o', color='red', markersize=12, markeredgecolor='black', elinewidth=4)
ax[1].errorbar(lae_z, Y_minus_Ye, yerr=Y_minus_Ye_err, xerr=[[dz_inf], [dz_sup]], fmt='o', color='red', markersize=12, markeredgecolor='black', elinewidth=4)

plt.savefig(plot_dir / f'colour_redshift_evolution_z{redshift}_WITH_LAE.pdf')

plt.show()
