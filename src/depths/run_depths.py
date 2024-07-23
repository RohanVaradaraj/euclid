#!/usr/bin/env python4

"""
run_depths.py

Run the depth code (noise_calc.py) in all filters and fields.

Created: Wed 5th January 2022

"""

import os

#########################################################################################

# Switches for which field to run on.
# Change same switches in noise_calc.py

dr3 =    True
dr2 =    False
cdfs =   False
cosmos = False
uvista = False

# Run irac on dr3 or cosmos
irac =   False
servs =  False
cds =    False

video =  True
vstack = True
###################################### XMM DR3 ###########################################

if dr3:

    # Base directory
    baseDir = '/mnt/hoy/temporaryFilesROHAN/HSCpatchesDR3/'

    # Fields
    fields = ['XMM1', 'XMM2', 'XMM3']

    # Filters
    #filts =  ['G'] #, 'R', 'I', 'Z', 'Y', 'NB0816', 'NB0921']
    #filts =  ['R', 'Z', 'NB0816', 'NB0921']
    #filts = ['Y', 'J', 'Ks', 'H']
    filts = ['YJ']

    # Node
    node = 'normal'

    for i in range(0, len(fields)):
        for j in range(0, len(filts)):
            # Enter relevant directory
            if video == False:
                os.chdir(baseDir+'HSC-{0}_{1}_DEEP'.format(filts[j], fields[i]))
            if video == True and vstack == False:
                os.chdir(baseDir+'{0}_{1}_DEEP'.format(filts[j], fields[i]))
            if video and vstack:
                os.chdir('/mnt/hoy/temporaryFilesROHAN/VIDEO/stacking/{0}'.format(fields[i]))

            # Add to queue
            if video == False:
                os.system("addqueue -c 'HSC-{0}_{1} DR3 depth' -q '{2}' -m 10 -d ./depthfile.sh".format(filts[j], fields[i], node))
            if video == True and vstack == False:
                os.system("addqueue -c '{0}_{1} depth' -q '{2}' -m 10 -d ./depthfile.sh".format(filts[j], fields[i], node))
            if video and vstack:
                os.system("addqueue -c '{0}_{1} depth' -q '{2}' -m 10 -d ./depthfile.sh".format(filts[j], fields[i], node))

###################################### CDFS #########################################

# Run on CDFS:

if cdfs:
    baseDir = '/mnt/hoy/temporaryFilesROHAN/CDFS/'

    # Fields
    fields = ['CDFS1', 'CDFS2', 'CDFS3']

    # Filters
#    filts = ['VOICE-u', 'VOICE-g', 'VOICE-r', 'VOICE-i'] #, 'HSC-G', 'HSC-R', 'HSC-I', 'HSC-Z']
#    filts = ['Y', 'J'] #, 'H', 'Ks']
    filts = ['YJ']

    # Node
    node = 'normal'

    for i in range(0, len(fields)):
        for j in range(0, len(filts)):
            # Enter relevant directory
            if vstack == False:
                os.chdir(baseDir+'{0}_{1}'.format(fields[i], filts[j]))
            if video and vstack:
                os.chdir('/mnt/hoy/temporaryFilesROHAN/VIDEO/stacking/{0}'.format(fields[i]))


            # Add to queue
            if vstack == False:
                os.system("addqueue -c '{0}_{1} depth' -q '{2}' -m 10 -d ./depthfile.sh".format(fields[i], filts[j], node))
            if video and vstack:
                os.system("addqueue -c '{0}_{1} depth' -q '{2}' -m 10 -d ./depthfile.sh".format(filts[j], fields[i], node))

#################################### XMM DR2 #############################################

# Run on XMM DR2:
if dr2:

    # Base directory
    baseDir = '/mnt/hoy/temporaryFilesROHAN/HSCpatchesDR2/'

    # Fields
    fields = ['XMM1', 'XMM2', 'XMM3']

    # Filters
    #filts =  ['G', 'R', 'I', 'Z', 'Y', 'NB0816', 'NB0921']
    filts = ['G', 'I', 'Y']

    # Node
    node = 'normal' # 'cmb'

    for i in range(0, len(fields)):
        for j in range(0, len(filts)):
            # Enter relevant directory
            os.chdir(baseDir+'HSC-{0}_{1}_DEEP'.format(filts[j], fields[i]))

            # Add to queue
            os.system("addqueue -c 'HSC-{0}_{1} DR2 depth' -q '{2}' -m 10 -d ./depthfile.sh".format(filts[j], fields[i], node))

###################################### COSMOS ###########################################

# Run on COSMOS
if cosmos == True and cds == False:

    # Base directory
    baseDir = '/mnt/hoy/temporaryFilesROHAN/HSCpatchesDR3/'

    # Filters
    filts =  ['G', 'R', 'I', 'Z', 'Y', 'NB0816', 'NB0921']

    # Node
    node = 'berg'

    for i in range(0, len(filts)):
        # Enter relevant directory
        os.chdir(baseDir+'HSC-{0}_COSMOS'.format(filts[i]))

        # Add to queue
        os.system("addqueue -c 'HSC-{0}_COSMOS DR2 depth' -q '{1}' -m 10 -d ./depthfile.sh".format(filts[i], node))

###################################### Spitzer/IRAC ###########################################

# Run on Spitzer/IRAC Cosmic Dawn Survey in CDFS
if irac and cds:

    # Base directory
    baseDir = '/mnt/hoy/temporaryFilesROHAN/Spitzer/CDFS/'

    # Fields
    #fields = ['CDFS1'] #, 'CDFS2', 'CDFS3']
    fields = ['COSMOS']

    # Filters
    filts =  ['ch1'] #, 'ch2']

    # Node
    node = 'berg'

    for i in range(0, len(filts)):
        for j in range(0, len(fields)):
            # Enter relevant directory
            os.chdir(baseDir+'{0}_{1}'.format(fields[j], filts[i]))

            # Add to queue
            os.system("addqueue -c '{0} CDS {1} depth' -q '{2}' -m 10 -d ./depthfile.sh".format(fields[j], filts[i], node))

# Run on Spitzer/IRAC Cosmic Dawn Survey in COSMOS
if cds and cosmos:

    # Base directory
    baseDir = '/mnt/hoy/temporaryFilesROHAN/Spitzer/COSMOS/'

    # Fields
    fields = ['COSMOS']

    # Filters
    filts =  ['ch1', 'ch2']

    # Node
    node = 'normal'

    for i in range(0, len(filts)):
        for j in range(0, len(fields)):
            # Enter relevant directory
            os.chdir(baseDir+'{0}_{1}'.format(fields[j], filts[i]))

            # Add to queue
            os.system("addqueue -c '{0} CDS {1} depth' -q '{2}' -m 10 -d ./depthfile.sh".format(fields[j], filts[i], node))


# Run on Spitzer/IRAC SERVS in XMM
if irac and servs:

    # Base directory
    baseDir = '/mnt/hoy/temporaryFilesROHAN/Spitzer/XMM/'

    # Fields
    fields = ['XMM1', 'XMM2', 'XMM3']

    # Filters
    filts =  ['ch1', 'ch2']

    # Node
    node = 'planet'

    for i in range(0, len(filts)):
        for j in range(0, len(fields)):
            # Enter relevant directory
            os.chdir(baseDir+'{0}_{1}'.format(fields[j], filts[i]))

            # Add to queue
            os.system("addqueue -c '{0} SERVS {1} depth' -q '{2}' -m 10 -d ./depthfile.sh".format(fields[j], filts[i], node))

###################################### UltraVISTA ###########################################

# Run on UltraVISTA
if uvista:

    # Base directory
    baseDir = '/mnt/hoy/temporaryFilesROHAN/COSMOS/'

    # Filters
    filts =  ['YJ_dr5'] #, 'Ks_dr5', 'H_dr5']
    #filts =  ['Ks_dr4'] #, 'J_dr4', 'H_dr4', 'Ks_dr4']

    # Node
    node = 'cmb'

    for i in range(0, len(filts)):
        # Enter relevant directory
        os.chdir(baseDir+'{0}'.format(filts[i]))

        # Add to queue
        os.system("addqueue -c 'UltraVISTA {0} depth' -q '{1}' -m 10 -d ./depthfile.sh".format(filts[i], node))
