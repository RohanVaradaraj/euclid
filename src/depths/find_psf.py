#!/usr/bin/env python3

""" find_psf.py

Aim: cleaner version of code run_psf.py from before.  Batch run this like the depths.

Created: Wed 22nd Jan 2020

"""

###################### Import useful modules #########################
import numpy as np
from new_psf_codes import get_psf
import os
from pathlib import Path

######################### Set-up ####################################
fields = ['COSMOS']

reqFilters = ['Y']

# enter the prefered queue.
queue = 'berg'
overwrite = True # this is the default anyway

# to run the first stage for each filter, set to true
# then check the pdf, update the star_param file
# then run again set to false
stars = True

############################### Loop ################################
## Loop through the different fields
for ff, fieldName in enumerate(fields):
    
    print('#############################################')
    print("Analysing field ", fieldName)

    outputDir = Path.home().parent.parent / 'vardy' / 'vardygroupshare' / 'rohan' / 'euclid' / 'data' / 'psf' / fieldName

    if os.path.isdir(outputDir) == False:
        os.system('mkdir ' + outputDir)
    
    # get the depths, nice clean code
    get_psf(fieldName, queue = queue, reqFilters = reqFilters, \
               overwrite = overwrite, outputDir = outputDir, starsOnly = stars)
