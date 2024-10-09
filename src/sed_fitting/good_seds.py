#!/usr/bin/env python3

"""
good_seds.py

Extracts good SED fits from LePhare .out file and the corresponding .spec files.

This will later be incorporated into a larger python wrapper for lePhare.

Created: Friday 22nd April 2022

"""

import os
from astropy.io import fits, ascii
from astropy.table import Table
import matplotlib.pyplot as plt
import numpy as np

'''SETUP'''

# Switch for clearing all previous .spec files in the folder
overwrite = False

#lephareFile = 'lephare_params_DR3.txt'
#lephareFile = 'lephare_params.txt'

lephareFile = 'euclid_para.txt'

# Read parameters in from lephare input txt file
lephareParams = Table.read(lephareFile, format='ascii.commented_header')

z_cut = lephareParams['z_cut'][0]
chi_cut = lephareParams['chi_cut'][0]
field = lephareParams['field'][0]
zphotDir = str(lephareParams['dir_source'][0])
goodDir = str(lephareParams['dir_target'][0])
inFile = str(lephareParams['inputName'][0])

print(field)
# .out file directory
outDir = '/users/varadaraj/lephare/lephare_dev/test/'

# .spec file directory
specDir = zphotDir
#specDir = '/mnt/hoy/temporaryFilesROHAN/lephare/secondFitting/{0}/'.format(field) # Brown dwarf fitting
#specDir = '/mnt/hoy/temporaryFilesROHAN/lephare/secondFitting/5{0}_emlines/'.format(field) # irac fitting


# Output file name
#outFile = '{0}FULL.out'.format(field)
#outFile = '2{0}vis.out'.format(field)
outFile = inFile[:-3] + '.out'
print(outFile)

# Read in the ascii table
data = ascii.read(outDir + outFile, format='no_header')

# PLOT THE PHOTOMETRIC REDSHIFTS

#plt.hist(data['col2'], bins=np.arange(0, 10, 0.05), color='orange')
plt.hist(data['col2'], color='orange', bins=np.arange(0, 10, 0.1))
plt.xlabel('z')
plt.ylabel('Count')
plt.title('LePhare photo-z output ({0})'.format(field))
plt.xticks(np.arange(0,10, step=1), fontsize=10)
plt.savefig('/mnt/vardy/vardygroupshare/HSC_SSP_DR3/simulations/plots/initial_ZPhot_{0}.pdf'.format(field))
#plt.savefig(goodDir + '../../plots/z_distn_{0}.pdf'.format(field))
plt.show()

# PLOT THE HIGH-Z MODEL CHI-SQUARE

print(data['col6'])
plt.hist(data['col6'], color='orange', bins=np.arange(0, 50, 0.1))
plt.xlabel('chi2')
plt.ylabel('Count')
plt.show()

#plt.savefig(zphotDir + '../../plots/bdchi_distn_{0}.pdf'.format(field))

# Following Bowler et al. 2012 for selection.

# Cut at z > 6 or 6.5
#print('##################')
#print('CUTTING AT z > ', z_cut)
#data = data[:][data['col2'] > z_cut]
#data = data[:][data['col2'] > 6]

#print('{0} objects left'.format(len(data)))

# Cut at at chi_squared < 11.3, equivalent to 2 sigma.
#print('###################')
#print('CUTTING AT Chi2 < ', chi_cut)
#data = data[:][data['col6'] < chi_cut]
#exit()
#print('{0} objects left'.format(len(data)))

#delta = data['col15'] - data['col6']
#plt.hist(delta, bins=np.arange(0, 10, 0.2))
#plt.show()


# Cut at chi_squared_1 - chi_squared_2 > 4. Also ensure secondary model is at z<5.

#cond1 = (data['col15'] - data['col6'] > 4) # & (data['col14'] < 3) # FOR THE SIMULATION, CHECK LOW-Z SECOND SOLUTION
#cond2 = (data['col14'] > 3) # If secondary solution at high-z, keep it.

#data = data[:][cond1] # | cond2] 

#data = data[:][(data['col15'] - data['col6'] > 4)]  # This step actually from rebecca's z=8-10 paper




''' QUASARS '''
# Cut so quasar is delta chi2 =4 better than galaxy.
#data = data[:][ ( data['col6'] - data['col18'] > 4)]

#print(len(data))

# Cut at at chi_squared < 11.3, same for galaxies but probably arbitrary.
#print('###################')
#print('CUTTING AT Chi2_QSO < ', 200)
#data = data[:][data['col18'] < 200]


# Require stellar fit bad.

#plt.hist(data['col21'], bins=np.arange(0, 110, 2))
#plt.show()


data = data[:][(data['col21'] > 10) | ( (data['col21'] < 10) & (data['col21'] - data['col6'] > 4) )] # Allowing good stars
#data = data[:][(data['col21'] > 10)]
#data = data[:][(data['col21'] < 10)] # taking a closer look at good stellar fits...

print(len(data))

# Copy corresponding surviving candidates into a separate folder

# Change to .spec dir
os.chdir(specDir)

if overwrite:
    os.system('rm {0}'.format(goodDir) + '*')

for i, id in enumerate(data['col1']):

    os.system('cp Id{0}.spec {1}'.format(str(id).zfill(9), goodDir))
    print('Moved {0}'.format(str(id).zfill(9)))
