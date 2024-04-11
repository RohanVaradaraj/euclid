#!/usr/bin/env python3

""" 
find_depths.py

Calculates 5sigma depths. Specify the fields, filters and aperture sizes you want to run through the depth code.

Modified from original to work for Euclid data and to be be more general.

Originally created: Fri 29th Nov 2019
Modified: Wednesday 20th March 2024.

"""

###################### Import useful modules #########################
import numpy as np
from new_depth_codes import get_depths
from pathlib import Path
import os

os.system('module load sextractor')
######################### Set-up ####################################
fields = ['COSMOS']

# required aperture diameter
reqFilters = ['Y', 'J', 'H', 'VIS']

# enter the prefered queue.
queue = 'normal'
overwrite = True # False is the default anyway

####################################################################
# required aperture diameters to run through
# change at your peril!
#apDiametersAS = [1.8, 2.0, 3.0, 4.0, 5.0]
apDiametersAS = [0.3, 0.5, 0.6, 1.0, 1.2, 1.8] # For euclid.

############################### Loop ################################
## Loop through the different fields
for ff, fieldName in enumerate(fields):
    
    print('#############################################')
    print("Analysing field ", fieldName)
    outputDir = Path.home().parent.parent / 'vardy' / 'vardygroupshare'/ 'rohan' / 'euclid' / 'data' / 'depths' / fieldName

    if os.path.isdir(outputDir) == False:
        os.system('mkdir ' + outputDir)
    
    # get the depths, nice clean code
    get_depths(fieldName, queue = queue, req_filters = reqFilters, \
               overwrite = overwrite, output_dir = outputDir, ap_diametersAS = apDiametersAS)
