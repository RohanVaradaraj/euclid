#!/usr/bin/env python3

"""
pipeline_lephare.py

Takes the final catalogue and puts it through all the lephare preparation and running.

First calls convert_fits_txt.py to convert the catalogue into the desired LePhare input format.

Then runs the LePhare SED fitting.

Then calls good_seds.py to take all objects that are suitable high-z candidates.

Then runs plot_spec.py to produce a pdf of the best fitting SEDs and their P(z)

The code is run in stages.

Created: Tuesday 10th May 2022
Modified: Wednesday 17th April 2022.

"""

'''SETUP'''

import os
from astropy.table import Table
from pathlib import Path

'''SWITCHES'''

# Memory and node to run on in the queue.
mem = 3
node = 'cmb'

# If true, this deletes all previous .spec files in the good folder.
overwrite = False

# May be best to run convert by itself first in case it takes a while on the full catalogues.
convert = False

# All of these can be run at the same time.

# Build library from templates
build_lib =   True

# Build filters
build_filts = True

# Build theoretical magnitudes
mags =        True

# Run photometric redshifts!
z_phot =      True

# Extract SEDs with desired properties (e.g. z>6)
good_seds =   False

# Plot SEDs
plots =       False

# Parameters for the process are in lephare_params.txt
# Lephare config is in ~/lephare/lephare_dev/config/rohan_test.para

# Read catalogue name and other things in from lephare input txt file
lephareFile = 'lephare_params.txt'

lephareParams = Table.read(lephareFile, format='ascii.commented_header')

# Catalogue
cat = str(lephareParams['catalogue'][0])

# Directory to run lephare zphot into
zphotDir = str(lephareParams['dir_source'][0])
#zphotDir = str(lephareParams['dir_target'][0])

# Field
field = str(lephareParams['field'][0])
print(field)

# Parameter file
para = 'euclid.para'

'''RUN convert_fits_txt.py'''

if convert:

    # Read catlaogue name in from lephare input txt file

    lephareParams = Table.read('lephare_params.txt', format='ascii.commented_header')
    cat = str(lephareParams['catalogue'][0])

    print('#######################################################################')
    print('Running conversion to txt file on {0}'.format(cat))

    # Create .sh file
    shell = 'tmp_run_txt_conv.sh'
    f = open(shell, 'w')
    f.write('#!/bin/bash\n')
    f.write('./convert_fits_txt.py')
    f.close()

    # make executable
    os.system('chmod u+x ' + shell)

    # Add to queue
    os.system("addqueue -c 'convert_fits_txt.py' -q {0} -m {1} -d ./{2}".format(node, mem, shell))


    print('CONVERSION DONE ', '#######################################################################')

'''BUILD TEMPLATES FROM SED LIBRARIES'''

if build_lib:

    print('####### BUILDING LIBRARIES######')

    # Move into lephare directories.
    os.chdir('/mnt/zfsusers/varadaraj/lephare/lephare_dev/config')

    # Build stars
    os.system('$LEPHAREDIR/source/sedtolib -t S -c $LEPHAREDIR/config/{0}'.format(para))

    # Build quasars
    os.system('$LEPHAREDIR/source/sedtolib -t Q -c $LEPHAREDIR/config/{0}'.format(para))

    # Build galaxies
    os.system('$LEPHAREDIR/source/sedtolib -t G -c $LEPHAREDIR/config/{0}'.format(para))


    print('#######################################################################')


'''BUILD FILTERS'''

if build_filts:

    # Move into lephare directories.
    os.chdir('/mnt/zfsusers/varadaraj/lephare/lephare_dev/config')

    os.system('$LEPHAREDIR/source/filter  -c $LEPHAREDIR/config/{0}'.format(para))



'''COMPUTE THEORETICAL MAGNITUDES'''

if mags:

    # Move into lephare directories.
    os.chdir('/mnt/zfsusers/varadaraj/lephare/lephare_dev/config')

    # Stars
    os.system('$LEPHAREDIR/source/mag_star -c  $LEPHAREDIR/config/{0}'.format(para))

    # Quasars
    os.system('$LEPHAREDIR/source/mag_gal  -t Q -c $LEPHAREDIR/config/{0}'.format(para))

    # Galaxies
    os.system('$LEPHAREDIR/source/mag_gal  -t G -c $LEPHAREDIR/config/{0}'.format(para))



'''RUN PHOTO-Zs'''

if z_phot:

    # Change into the directory where we will output the .spec files
    os.chdir(zphotDir)

    # Remove previous files, if desired.
    if overwrite:
        os.system('rm ./*')

    # Run the command
    os.system('$LEPHAREDIR/source/zphota -c $LEPHAREDIR/config/{0}'.format(para))



'''EXTRACT INITIAL HIGH REDSHIFT SEDS'''

if good_seds:

    # Move into code directories
    os.chdir('/mnt/vardy/vardygroupshare/HSC_SSP_DR3/codes/')

    # There are parameters within lephare_params.txt that we can change here.
    os.system('python3 good_seds.py')


'''SAVE SED PLOTS'''

if plots:

#    print('SENDING PLOTS TO QUEUE')

    # Move into code directories
    os.chdir('/mnt/vardy/vardygroupshare/HSC_SSP_DR3/codes/')

    os.system('python3 plot_spec.py')
    print('Done!')

    # Create .sh file
#    shell = 'tmp_plot_seds.sh'
#    f = open(shell, 'w')
#    f.write('#!/bin/bash\n')
#    f.write('./plot_spec.py')
#    f.close()

    # make executable
#    os.system('chmod u+x ' + shell)

    # Add to queue
#    os.system("addqueue -c 'plot_seds.py' -q {0} -m {1} -d ./{2}".format(node, mem, shell))

#    print('Done! Check progress of plots in the queue.')
