#!/usr/bin/env python3

"""
Adds lyman alpha emission lines to SED models.

Modified from some code used from my z=7 paper. We want to know how these sources appear in our overlapping filter set with Euclid.

Created: Thursday 19th July 2022. Hottest day ever!
Modified: Thursday 18th April 2024. 

"""

import numpy as np
import matplotlib.pyplot as plt
from astropy.io import ascii
import os
from pathlib import Path

plt.rcParams['axes.linewidth'] = 2.5
plt.rcParams.update({'font.size': 15})
plt.rcParams['figure.dpi'] = 100

# Set up directories
lya_dir = Path.cwd().parent.parent / 'data' / 'simulations' / 'LAEs'
sed_dir = Path.home() / 'lephare' / 'lephare_dev' / 'sed' / 'GAL' / 'BC03' / 'ASCII'

#sed_dir = '/mnt/zfsusers/varadaraj/lephare/lephare_dev/sed/GAL/LYMAN_ALPHA/ASCII/'
os.chdir(sed_dir)

for f in os.listdir(sed_dir)[::10]:

    # Make sure we only take the .ascii files
    if f[-5:] == 'ascii':

        # Print file name
        #print(f)

        # Open the file!
        sed = ascii.read(f)

        # plt.plot(sed['col1'], sed['col2'])
        # plt.yscale('log')
        # plt.show()
        # continue

        # Get lower and upper limits in wlen for measuring continuum

        # Flux value
        low = next(i for i in sed['col1'] if i > 1249.0)
        # Index of the above flux value
        idx_low = np.where(sed['col1'] == 1255.0)[0][0]

        upp = next(i for i in sed['col1'] if i > 1299.0)
        idx_upp = np.where(sed['col1'] == 1305.0)[0][0]

        # Print this information
        #print('Lower continuum bound: ', low, idx_low)
        #print('Upper continuum bound: ', upp, idx_upp)
        #print('Unnormalised continuum flux: ', np.sum(sed['col2'][idx_low:idx_upp]))

        # Lines to see where we are measuring the continuum
        # plt.vlines(low, 0, 0.002, color='red')
        # plt.vlines(upp, 0, 0.002, color='red')
        # plt.show()
        # continue

        # Find SED resolution
        h = sed['col1'][300] - sed['col1'][299]
        #print('h: ', h)

        # Compute continuum flux
        f_cont = np.sum(sed['col2'][idx_low:idx_upp]) / (idx_upp - idx_low)
        #print('Continuum flux: ', f_cont)

        # Create EW grid from 10 to 240 Angstroms, in steps of 10.
        EWs = np.arange(0, 250, 10)

        tmp = []

        # Iterate through the EWs
        for i, ew in enumerate(EWs):

            # Find intensity of line
            f_lya = f_cont * ( (EWs[i]/h) - 1)

            #print('EW and flux: ', ew, f_lya)

            # Find where Lyman alpha line needs to go
            idx_lya = np.where( (sed['col1'] > 1214.0) & (sed['col1'] < 1225.0) )[0][0]

            # Make sed copy
            tmp = sed

            # Change to 1216A and add Lya flux.
            tmp['col1'][idx_lya] = 1216.0
            tmp['col2'][idx_lya] = f_lya

            # Save the SED to the Lyman alpha folder
            lya_str = f[:-6] + '_EW{0}A'.format(EWs[i]) + '.ascii'
            ascii.write(sed, lya_dir / lya_str, overwrite=True, format='no_header')

            # Plot
            # plt.plot(sed['col1'], sed['col2'], label='Continuum')
            # plt.plot(tmp['col1'], tmp['col2'], label='Lyman alpha')
            # plt.yscale('log')
            # plt.show()

#            f.append()

