#!/usr/bin/env python3

"""
lyman_alpha_emitters.py

Simulate Lyman-alpha emitters (LAEs) passing through the Euclid+VISTA filter set.

Created: Friday 10th May 2024.
"""

from lyman_alpha_utils import *
from pathlib import Path
from matplotlib.animation import FuncAnimation
import matplotlib.animation as animation
import matplotlib.pyplot as plt
from brown_dwarf_colours import getFilters

table_dir = Path.cwd().parent.parent / 'data' / 'simulations' / 'LAEs' / 'tables'
plot_dir = Path.cwd().parent.parent / 'plots' / 'LAEs'

#! ####################### Making use of numpy vectorisation ##############################################
# # Choose an Muv
# Muv = -21.5

# # Set up the grid
# redshifts = np.arange(5.5, 10.01, 0.01)
# EWs = np.arange(0, 250, 10)
# # redshifts = np.arange(5.5, 5.7, 0.1) #* test values
# # EWs = np.arange(0, 40, 10)

# # Also run through an age-attenuation grid
# Avs = np.arange(0, 0.7, 0.1)
# ages = np.array([50, 100, 150, 200, 300, 400, 500])

# for i, Av in enumerate(Avs):
#     for j, age in enumerate(ages):

#         print('----------------------------')
#         print(f'Av: {Av}, age: {age} Myr')
#         print(f'Iteration {i+1} of {len(Avs)}, {j+1} of {len(ages)}')
#         print(f'Total iteration {i*len(ages) + j + 1} of {len(Avs)*len(ages)}')
#         print('----------------------------')


        
#         # Initialize a 2D grid to store filter_colour values
#         filter_colours = np.zeros((len(EWs), len(redshifts)))

#         # Generate redshift and EW grids
#         redshift_grid, EW_grid = np.meshgrid(redshifts, EWs)

#         # print number of combinations of redshift and EW
#         print(f'Number of EW iterations: {len(EWs)}')

#         print(redshift_grid.ravel(), EW_grid.ravel())

#         # Run Lyman-alpha emitter simulation for all combinations, at the given dust and age.
#         mags, errors = runLymanAlphaEmitter_vectorised(redshift_grid.ravel(), Muv, EW_grid.ravel(), Av, age)

#         # mags is now a list of dictionaries for each combination of redshift and EW. 
#         # Get the values for each filter and reshape them to match the grid shape
#         mags = {k: np.array([v[k] for v in mags]).reshape(len(EWs), len(redshifts)) for k in mags[0].keys()}

#         # For each key, create a variable with the same name and assign the corresponding values
#         for k, v in mags.items():
#             exec(f'{k} = v')

#             # Reshape each variable to match the grid shape
#             exec(f'{k} = {k}.reshape(len(EWs), len(redshifts))')


#         # Save an astropy table with redshift, EW and filter colours for all the above variables
#         t = Table([redshift_grid.ravel(), EW_grid.ravel(), y.ravel(), Ye.ravel(), Je.ravel(), He.ravel(), VIS.ravel(), Y.ravel(), J.ravel(), H.ravel(), Ks.ravel(), g.ravel(), r.ravel(), i.ravel(), nb816.ravel(), z.ravel(), nb921.ravel()], 
#         names=['redshift', 'EW', 'y', 'Ye', 'Je', 'He', 'VIS', 'Y', 'J', 'H', 'Ks', 'g', 'r', 'i', 'nb816', 'z', 'nb921'])

#         print(t)

#         # Convert errors dictionary to astropy table.
#         t_err = Table([[errors[k]] for k in errors.keys()], names=errors.keys())
#         print(t_err)

#         # Save the tables
#         t.write(table_dir / f'LAE_fluxes_Av_{Av}_age_{age}Myr.fits', overwrite=True)
#         t_err.write(table_dir / f'LAE_error_Av_{Av}_age_{age}Myr.fits', overwrite=True)

# exit()

# #Plot the imshow
# plt.figure(figsize=(10, 8))
# im = plt.imshow(y-Ye, extent=[redshifts[0], redshifts[-1], EWs[0], EWs[-1]], cmap='coolwarm', aspect='auto', origin='lower')
# plt.colorbar(im, label='mag')

# # Add labels
# plt.xlabel(r'$z$')
# plt.ylabel('EW ' + r'$(\AA)$')
# plt.title(r'$y - Y_{e}$')
# plt.savefig(plot_dir / f'z_vs_EW_Muv{Muv}.png')
# plt.show()

#!################################## MAKE AN ANIMATION OF THE REDSHIFTING SED #################

# Get filter sets
euclid_filters = getFilters('Euclid')
vista_filters = getFilters('VISTA')
hsc_filters = getFilters('HSC')

redshifts = np.arange(6.6, 7.5, 0.01)
plt.figure(figsize=(10, 6))

Muv=-23
EW=200

# Create the animation
ani = FuncAnimation(plt.gcf(), update_plot, fargs=(Muv, EW), frames=redshifts, interval=200)

# Save the animation as a GIF
#ani.save(plot_dir / f'redshifting_lyman_alpha_Muv{Muv}_EW_{EW}A_yOverlap.gif', writer='imagemagick')
ani.save(plot_dir / f'LAE_z7_Muv{Muv}_EW_{EW}A_yOverlap.gif', writer='imagemagick')

# Show the animation
plt.show()

#####################! INDIVIDUAL SEDs ############################

# redshift = 7.0
# EW = 100
# Muv = -21.5

# euclid_filters = getFilters('Euclid')
# vista_filters = getFilters('VISTA')
# hsc_filters = getFilters('HSC')

# for Av in np.arange(0, 0.7, 0.1):

#     wlen, flux_sed = makeLBG(redshift=redshift, SFH_component='constant', age=(0, 13.8), massformed=10., metallicity=0.2, 
#                                 dust_type='Calzetti', Av=Av, nebular=True, logU=-2.)

#     wlen, flux_sed = set_Muv(z=redshift, Muv_target=Muv, wlen=wlen, flux=flux_sed)

#     wlen, flux_sed = add_emission_line(EW=EW, z=redshift, wlen=wlen, flux=flux_sed)

#     fluxes = convolveFilters([euclid_filters, vista_filters, hsc_filters], (wlen, flux_sed), magnitudes=False)

#     depths, _ = simulateDepths([euclid_filters, vista_filters, hsc_filters])

#     errors = getErrors(depths)

#     signal_to_noise = getSignalToNoise(fluxes, errors)

#     scattered_fluxes = scatterFluxes(fluxes, depths)

#     # Plotting
#     #plt.figure(figsize=(10, 6))
#     plt.plot(wlen, flux_sed, color='deepskyblue', lw=2.5, alpha=0.8)

# plt.xlabel(r'Wavelength ($\AA$)')
# plt.ylabel('Flux')
# plt.legend()
# plt.xlim(5000, 30000)
# plt.yscale('log')
# plt.ylim(3e-32, 1e-27)
# plt.tight_layout()

# plt.show()
# exit()


    # # Get the centre of the Euclid, VISTA and HSC filters and then put into one big dictionary
    # filter_centres = {}
    # filter_widths = {}
    # for filter_name in euclid_filters.keys():
    #     filter_key = filter_name.split('e')[0]
    #     centre, width = filterCentreAndWidth(filter_key, 'Euclid')
    #     filter_centres[filter_name] = centre
    #     filter_widths[filter_name] = width

    # for filter_name in vista_filters.keys():
    #     centre, width = filterCentreAndWidth(filter_name, 'VISTA')
    #     filter_centres[filter_name] = centre
    #     filter_widths[filter_name] = width

    # for filter_name in hsc_filters.keys():
    #     filter_keys = {'g': 'HSC-G_DR3', 'r': 'HSC-R_DR3', 'i': 'HSC-I_DR3', 'nb816': 'HSC-NB0816_DR3', 'z': 'HSC-Z_DR3', 'nb921': 'HSC-NB0921_DR3', 'y': 'HSC-Y_DR3'}
    #     centre, width = filterCentreAndWidth(filter_keys[filter_name], 'HSC')
    #     filter_centres[filter_name] = centre
    #     filter_widths[filter_name] = width

    # # convert filter centres to angstroms from microns
    # filter_centres = {k: v*1e4 for k, v in filter_centres.items()}
    # filter_widths = {k: v*1e4 for k, v in filter_widths.items()}

    # # If signal to noise is larger than 2, plot the fluxes as errorbars
    # for filter_name, flux in scattered_fluxes.items():

    #     # Set colour based on filter name
    #     if filter_name in ['VIS', 'Ye', 'Je', 'He']:
    #         colour = 'red'
    #         fmt = 'o'
    #     elif filter_name in ['Y', 'J', 'H', 'Ks']:
    #         colour = 'blue'
    #         fmt = 's'
    #     elif filter_name in ['g', 'r', 'i', 'nb816', 'z', 'nb921', 'y']:
    #         colour = 'black'
    #         fmt='D'

    #     if signal_to_noise[filter_name] > 2:
    #         plt.errorbar(filter_centres[filter_name], flux, yerr=errors[filter_name], xerr=filter_widths[filter_name]/2, fmt=fmt, color=colour)

    #     # Else, plot a 2sigma upper limit
    #     else:
    #         plt.errorbar(filter_centres[filter_name], flux+2*errors[filter_name], fmt='v', uplims=True, color=colour)

    # # Dummy plots for labels for each instrument
    # plt.errorbar(0, 1, yerr=0, xerr=0, fmt='o', label='Euclid', color='red')
    # plt.errorbar(0, 1, yerr=0, xerr=0, fmt='s', label='VISTA', color='blue')
    # plt.errorbar(0, 1, yerr=0, xerr=0, fmt='D', label='HSC', color='black')

    # plt.xlabel(r'Wavelength ($\AA$)')
    # plt.ylabel('Flux')
    # plt.legend()
    # plt.xlim(5000, 30000)
    # plt.yscale('log')
    # plt.ylim(3e-32, 1e-27)
    # plt.tight_layout()

    # plt.title(r'$z = $' + f'{round(redshift, 2)}, ' + r'$EW_{0} = $' + f'{EW} $\AA$' + r', $M_{\mathrm{UV}} = $' + f'{Muv}')

    # plt.show()












