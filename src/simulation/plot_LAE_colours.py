"""
plot_LAE_colours.py

This script reads in the output of the Lyman-alpha emitter simulation and plots various filter colours as a function of redshift and EW.

Created: Wednesday 15th May 2024.
"""

from astropy.table import Table
from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np
from brown_dwarf_colours import *
from lyman_alpha_utils import getSignalToNoise

plt.rcParams['axes.linewidth'] = 2.5
plt.rcParams.update({'font.size': 15})
plt.rcParams['figure.dpi'] = 100

# Function to scale plotting values from one range to another
def scale_values(values, min_old, max_old, min_new, max_new):
    return ((values - min_old) / (max_old - min_old)) * (min_new - max_new) + max_new

# Plot z vs EW with filter colours?
plot_z_vs_EW = False

# Plot colour colour diagrams?
plot_colour_colour = True

# If colour colour, what redshift range? 7 or 8?
z_LAE = 8

# Set up the grid
redshifts = np.arange(5.5, 10.01, 0.01)
EWs = np.arange(0, 250, 10)

table_dir = Path.cwd().parent.parent / 'data' / 'simulations' / 'LAEs' / 'tables'
dwarf_dir = Path.cwd().parent.parent / 'data' / 'simulations' / 'tables'

# Read in the tables
fluxes = Table.read(table_dir / 'LAE_fluxes.fits')
errors = Table.read(table_dir / 'LAE_errors.fits')

# get Ye and Je
Ye = fluxes['Ye']
Je = fluxes['Je']
y = fluxes['y']
Y = fluxes['Y']
J = fluxes['J']
Je = fluxes['Je']
He = fluxes['He']
VIS = fluxes['VIS']
z = fluxes['z']

# Get magnitudes
mag_Ye = flux_to_mag(Ye)
mag_Je = flux_to_mag(Je)
mag_y = flux_to_mag(y)
mag_Y = flux_to_mag(Y)
mag_J = flux_to_mag(J)
mag_Je = flux_to_mag(Je)
mag_He = flux_to_mag(He)
mag_VIS = flux_to_mag(VIS)
mag_z = flux_to_mag(z)

# Plot Ye-Je on imshow as a function of redshift and EW
zmin = np.min(fluxes['redshift'])
zmax = np.max(fluxes['redshift'])
EWmin = np.min(fluxes['EW'])
EWmax = np.max(fluxes['EW'])

colour = mag_y - mag_Ye


###########################! Plot filter-filter colours as a function of redshift and EW ###########################
if plot_z_vs_EW:
        
    # Reshape colour to match the grid shape
    colour = colour.reshape(len(EWs), len(redshifts))

    plt.figure(figsize=(10, 10))
    plt.imshow(colour, origin='lower', aspect='auto', extent=[zmin, zmax, EWmin, EWmax], cmap='coolwarm')

    # Limit the colourbar 
    plt.clim(-1.0, 5)

    plt.colorbar(label=r'$y - Y_{E}$')
    plt.xlabel(r'$z$')
    plt.ylabel(r'$\mathrm{EW}_{o} (\AA)$')
    plt.show()

###########################! Plot colour-colour diagrams ###########################
if plot_colour_colour:

    if z_LAE == 6:

        c1 = mag_z-mag_y
        c2 = mag_VIS-mag_z

        # Limit the redshift
        z_lower = 6.2
        z_upper = 6.6
        mask = (fluxes['redshift'] > z_lower) & (fluxes['redshift'] < z_upper)

        # Apply mask
        c1 = c1[mask]
        c2 = c2[mask]

        # Scale redshifts to map to size
        scaled_sizes = scale_values(fluxes['redshift'][mask], z_lower, z_upper, 1, 50)

        plt.figure(figsize=(10, 8))

        plt.scatter(c1, c2, c=fluxes['EW'][mask], cmap='viridis', s=scaled_sizes, alpha=0.7)
        plt.colorbar(label=r'$\mathrm{EW}_{0} \ (\AA)$')
        #plt.clim(0, 70)

        #! Plotting commands
        plt.xlabel(r'$z_{\mathrm{HSC}} - y_{\mathrm{HSC}}$')
        plt.ylabel(r'$I_{E} - z_{\mathrm{HSC}}$')


        #plt.xlim(-0.5, 0.7)
        #plt.ylim(-0.9, 6)

        plt.title(f'{z_lower}' + r'$< z < $' + f'{z_upper}')

        plt.show()

    if z_LAE == 7:

        c1 = mag_Ye-mag_Y
        c2 = mag_y-mag_Ye

        # Limit to where redshift is less than 7.3
        z_lower = 7.1
        z_upper = 7.5
        mask = (fluxes['redshift'] > z_lower) & (fluxes['redshift'] < z_upper)

        # Apply mask
        c1 = c1[mask]
        c2 = c2[mask]

        # Scale redshifts to map to size
        scaled_sizes = scale_values(fluxes['redshift'][mask], z_lower, z_upper, 1, 50)

        plt.figure(figsize=(10, 8))

        plt.scatter(c1, c2, c=fluxes['EW'][mask], cmap='viridis', s=scaled_sizes, alpha=0.7)
        plt.colorbar(label=r'$\mathrm{EW}_{0} \ (\AA)$')
        plt.clim(0, 70)

        #! Draw the selection function
        x = [0.3, 0.3]
        y = [-2, 2]
        plt.plot(x, y, marker='none', linestyle='--', color='black', lw=2.5, alpha=0.8)

        x = [0.3, 1.]
        y = [2, 2]
        plt.plot(x, y, marker='none', linestyle='--', color='black', lw=2.5, alpha=0.8)  


        #? ----- Load brown dwarf spectra ------
        bd_table = Table.read(dwarf_dir / 'dwarf_spectra_mags.fits') #! Spectra
        spex_table = Table.read(dwarf_dir / 'spex_template_mags.fits') #! Templates

        # Get the colours from spectra
        mag_y_dwarf = bd_table['y']
        mag_Ye_dwarf = bd_table['Ye']
        mag_Y_dwarf = bd_table['Y']

        # Get colours from templates
        mag_Ye_dwarf_template = spex_table['Ye']
        mag_Y_dwarf_template = spex_table['Y']
        mag_y_dwarf_template = spex_table['y']
        spectral_type = spex_table['Spectral Type']

        # Calculate the colours
        c1_dwarf = mag_Ye_dwarf - mag_Y_dwarf
        c2_dwarf = mag_y_dwarf - mag_Ye_dwarf
        c1_dwarf_template = mag_Ye_dwarf_template - mag_Y_dwarf_template
        c2_dwarf_template = mag_y_dwarf_template - mag_Ye_dwarf_template

        # Plot the brown dwarf colours
        plt.scatter(c1_dwarf, c2_dwarf, c='red', s=10, marker='*', label='Brown dwarfs', alpha=0.7)

        plt.scatter(c1_dwarf_template, c2_dwarf_template, c='black', s=13, marker='*', label='Spex templates', alpha=0.7)

        # Add text of spectral type next to the templates
        for i in range(len(spectral_type)):
            plt.text(c1_dwarf_template[i]+0.01, c2_dwarf_template[i]+0.01, f'{spectral_type[i]}', fontsize=8, color='black')

        #? ------ Draw typical LBG evolution tracks ------
        # lbg_table = Table.read(dwarf_dir / 'lbg_spectra_mags.fits')

        # # Apply redshift mask
        # lbg_table = lbg_table[(lbg_table['Redshift'] > z_lower) & (lbg_table['Redshift'] < z_upper)]

        # # Get the colours
        # mag_y_lbg = lbg_table['y']
        # mag_Ye_lbg = lbg_table['Ye']
        # mag_Y_lbg = lbg_table['Y']

        # # Calculate the colours
        # c1_lbg = mag_Ye_lbg - mag_Y_lbg
        # c2_lbg = mag_y_lbg - mag_Ye_lbg

        # # Plot the LBG colours
        # plt.scatter(c1_lbg, c2_lbg, c='deepskyblue', s=14, marker='o', label='LBGs', alpha=0.7)

        #! Plotting commands
        plt.xlabel(r'$y_{\mathrm{HSC}} - Y_{E}$')
        plt.ylabel(r'$Y_{E} - Y_{\mathrm{VISTA}}$')


        plt.xlim(-0.5, 0.7)
        plt.ylim(-0.9, 6)

        plt.title(f'{z_lower}' + r'$< z < $' + f'{z_upper}')

        plt.show()

    if z_LAE == 8:
        c1 = mag_Je-mag_J
        c2 = mag_Ye-mag_Je

        # Limit to where redshift is less than 7.3
        z_lower = 8.3
        z_upper = 8.5
        mask = (fluxes['redshift'] > z_lower) & (fluxes['redshift'] < z_upper)

        # Apply mask
        c1 = c1[mask]
        c2 = c2[mask]

        # Scale redshifts to map to size
        scaled_sizes = scale_values(fluxes['redshift'][mask], z_lower, z_upper, 1, 50)

        plt.figure(figsize=(10, 8))
        plt.scatter(c1, c2, c=fluxes['EW'][mask], cmap='viridis', s=scaled_sizes, alpha=0.7)

        #! Draw the selection function
        # x = [0.3, 0.3]
        # y = [-2, 2]
        # plt.plot(x, y, marker='none', linestyle='--', color='black', lw=2.5, alpha=0.8)

        # x = [0.3, 1.]
        # y = [2, 2]
        # plt.plot(x, y, marker='none', linestyle='--', color='black', lw=2.5, alpha=0.8)  


        #? ----- Load brown dwarf spectra ------
        bd_table = Table.read(dwarf_dir / 'dwarf_spectra_mags.fits') #! Spectra
        spex_table = Table.read(dwarf_dir / 'spex_template_mags.fits') #! Templates

        # Get the colours from spectra
        mag_Ye_dwarf = bd_table['Ye']
        mag_J_dwarf = bd_table['J']
        mag_Je_dwarf = bd_table['Je']

        # Get colours from templates
        mag_Ye_dwarf_template = spex_table['Ye']
        mag_J_dwarf_template = spex_table['J']
        mag_Je_dwarf_template = spex_table['Je']
        spectral_type = spex_table['Spectral Type']

        # Calculate the colours
        c1_dwarf = mag_Je_dwarf - mag_J_dwarf
        c2_dwarf = mag_Ye_dwarf - mag_Je_dwarf
        c1_dwarf_template = mag_Je_dwarf_template - mag_J_dwarf_template
        c2_dwarf_template = mag_Ye_dwarf_template - mag_Je_dwarf_template

        # Plot the brown dwarf colours
        plt.scatter(c1_dwarf, c2_dwarf, c='red', s=10, marker='*', label='Brown dwarfs', alpha=0.7)

        plt.scatter(c1_dwarf_template, c2_dwarf_template, c='black', s=13, marker='*', label='Spex templates', alpha=0.7)

        # Add text of spectral type next to the templates
        for i in range(len(spectral_type)):
            plt.text(c1_dwarf_template[i]+0.01, c2_dwarf_template[i]+0.01, f'{spectral_type[i]}', fontsize=8, color='black')

        #! Plotting commands
        plt.xlabel(r'$J_{E} - J_{\mathrm{VISTA}}$')
        plt.ylabel(r'$Y_{E} - J_{E}$')

        plt.colorbar(label=r'$\mathrm{EW}_{0} \ (\AA)$')
        plt.clim(0, 240)

        # plt.xlim(-0.5, 0.7)
        # plt.ylim(-0.9, 6)

        plt.title(f'{z_lower}' + r'$< z < $' + f'{z_upper}')

        plt.show()   