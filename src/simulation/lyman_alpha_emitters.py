#!/usr/bin/env python3

"""
lyman_alpha_emitters.py

Simulate Lyman-alpha emitters (LAEs) passing through the Euclid+VISTA filter set.

Created: Friday 10th May 2024.
"""

from lyman_alpha_utils import *
from pathlib import Path

table_dir = Path.cwd().parent.parent / 'data' / 'simulations' / 'LAEs' / 'tables'



#! ####################### Making use of numpy vectorisation ##############################################
# Choose an Muv
Muv = -21.5

# Set up the grid
redshifts = np.arange(5.5, 10.01, 0.01)
EWs = np.arange(0, 250, 10)


# redshifts = np.arange(5.5, 5.7, 0.1) #* test values
# EWs = np.arange(0, 40, 10)

# Initialize a 2D grid to store filter_colour values
filter_colours = np.zeros((len(EWs), len(redshifts)))

# Generate redshift and EW grids
redshift_grid, EW_grid = np.meshgrid(redshifts, EWs)

# print number of combinations of redshift and EW
print(f'Number of EW iterations: {len(EWs)}')

print(redshift_grid.ravel(), EW_grid.ravel())

# Run Lyman-alpha emitter simulation for all combinations
mags, errors = runLymanAlphaEmitter_vectorised(redshift_grid.ravel(), Muv, EW_grid.ravel())

# mags is now a list of dictionaries for each combination of redshift and EW. 
# Get the values for each filter and reshape them to match the grid shape
mags = {k: np.array([v[k] for v in mags]).reshape(len(EWs), len(redshifts)) for k in mags[0].keys()}

# For each key, create a variable with the same name and assign the corresponding values
for k, v in mags.items():
    exec(f'{k} = v')

    # Reshape each variable to match the grid shape
    exec(f'{k} = {k}.reshape(len(EWs), len(redshifts))')


# Save an astropy table with redshift, EW and filter colours for all the above variables
t = Table([redshift_grid.ravel(), EW_grid.ravel(), y.ravel(), Ye.ravel(), Je.ravel(), He.ravel(), VIS.ravel(), Y.ravel(), J.ravel(), H.ravel(), Ks.ravel(), g.ravel(), r.ravel(), i.ravel(), nb816.ravel(), z.ravel(), nb921.ravel()], 
names=['redshift', 'EW', 'y', 'Ye', 'Je', 'He', 'VIS', 'Y', 'J', 'H', 'Ks', 'g', 'r', 'i', 'nb816', 'z', 'nb921'])

print(t)

# Convert errors dictionary to astropy table.
t_err = Table([[errors[k]] for k in errors.keys()], names=errors.keys())
print(t_err)

# Save the tables
t.write(table_dir / 'LAE_fluxes.fits', overwrite=True)
t_err.write(table_dir / 'LAE_errors.fits', overwrite=True)

# Plot the imshow
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

# # Get filter sets
# euclid_filters = getFilters('Euclid')
# vista_filters = getFilters('VISTA')
# hsc_filters = getFilters('HSC')

# redshifts = np.arange(5.6, 10., 0.02)
# plt.figure(figsize=(10, 6))

# Muv=-21.5
# EW=100

# # Create the animation
# ani = FuncAnimation(plt.gcf(), update_plot, fargs=(Muv, EW), frames=redshifts, interval=50)

# # Save the animation as a GIF
# ani.save(plot_dir / f'redshifting_lyman_alpha_Muv{Muv}_EW_{EW}A.gif', writer='imagemagick')

# # Show the animation
# plt.show()

#!##################### PLOT 2D GRID OF Z vs EW COLOURED BY DIFFERENT COLOURS, LOOPING SLOWLY #################

# # # Choose an Muv
# Muv = -21.5

# # Set up the grid
# redshifts = np.arange(5.5, 5.7, 0.1)
# EWs = np.arange(50, 100, 10)

# # Initialize a 2D grid to store filter_colour values
# filter_colours = np.zeros((len(EWs), len(redshifts)))

# # Loop through the redshifts
# for i, z in enumerate(redshifts):
#     print('z = ', round(z, 2))

#     # Make the Lyman-break galaxy and set the absolute magnitud

#     # Loop through the EWs
#     for j, EW in enumerate(EWs):
#         print('EW = ', EW)
#         mags = runLymanAlphaEmitter(z, Muv, EW)
#         # Calculate the filter_colour
#         filter_colour = mags['y'] - mags['Ye']
#         # Store the filter_colour in the grid
#         filter_colours[j, i] = filter_colour

# # Plot the imshow
# plt.figure(figsize=(10, 8))
# im = plt.imshow(filter_colours, extent=[redshifts[0], redshifts[-1], EWs[0], EWs[-1]], cmap='coolwarm', aspect='auto', origin='lower')
# plt.colorbar(im, label='mag')

# # Add labels
# plt.xlabel(r'$z$')
# plt.ylabel('EW ' + r'$(\AA)$')
# plt.title(r'$y - Y_{e}$')
# #plt.savefig(plot_dir / f'z_vs_EW_Muv{Muv}.png')
# plt.show()
# exit()













