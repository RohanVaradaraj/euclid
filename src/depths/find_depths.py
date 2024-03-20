#!/usr/bin/env python3

""" find_depths.py

Aim: Cleaner python 3 version of my depth codes, make this nice!

Created: Fri 29th Nov 2019

"""

###################### Import useful modules #########################
import numpy as np
from new_depth_codes import get_depths
#from jwst_depth_codes import get_depths

import os

######################### Set-up ####################################
#fields = ['COSMOS']
fields = ['CDFS1', 'CDFS2', 'CDFS3']
#fields = ['XMM1'] #, 'XMM2', 'XMM3']

# required aperture diameter
#reqAper = 'all' # or specific one
#reqFilters = ['Y', 'J', 'H', 'Ks']
#reqFilters = ['HSC-R_DR3', 'HSC-Z_DR3', 'HSC-NB0816_DR3', 'HSC-NB0921_DR3']
#reqFilters = ['HSC-G_DR3'] #, 'HSC-R_DR3', 'HSC-I_DR3', 'HSC-Z_DR3', 'HSC-Y_DR3', 'HSC-NB0816_DR3', 'HSC-NB0921_DR3']
reqFilters = ['ch1cds', 'ch2cds']

#reqFilters = ['f090w', 'f115w', 'f150w', 'f200w', 'f277w', 'f356w', 'f410m', 'f444w']
#reqFilters = ['f444w']


# enter the prefered queue.
queue = 'berg'
overwrite = True # False is the default anyway

####################################################################
# required aperture diameters to run through
# change at your peril!
#apDiametersAS = [1.8, 2.0, 3.0, 4.0, 5.0]

# for IRAC [2.8, 3.8, 5.8, 9.8, 11.6]
#apDiametersAS = np.array([0.1, 0.2, 0.3, 0.4, 0.5])
#apDiametersAS = np.array([1.0, 1.8, 2.0, 3.0, 4.0])
apDiametersAS = [2.8, 3.8, 5.8, 9.8, 11.6]

# For JWST
#apDiametersAS = [0.2, 0.3, 0.5, 0.9, 1.0]

############################### Loop ################################
## Loop through the different fields
for ff, fieldName in enumerate(fields):
    
    print('#############################################')
    print("Analysing field ", fieldName)
    outputDir = '/mnt/vardy/vardygroupshare/data/depths/{0}/'.format(fieldName)
    #outputDir = '../../../../hoy/temporaryFilesROHAN/depths/{0}/'.format(fieldName) # vardy is full

    if os.path.isdir(outputDir) == False:
        os.system('mkdir ' + outputDir)
    
    # get the depths, nice clean code
    get_depths(fieldName, queue = queue, reqFilters = reqFilters, \
               overwrite = overwrite, outputDir = outputDir, apDiametersAS = apDiametersAS)
