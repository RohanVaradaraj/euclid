#!/usr/bin/env python3

"""
UV_magnitude.py

Compute UV absolute magnitude! Step 1 of the luminosity function.

Created: Monday 15th August 2022

"""

import os
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import ascii, fits
from scipy.integrate import simps
from cosmology_calculator import cosmology_calculator
from astropy.cosmology import FlatLambdaCDM
from scipy import stats

def mag_to_flux(m):
    '''Convert mags to flux count'''
    flux = 10**(-0.4*(m+48.6))
    return flux

def flux_to_mag(f):
    mag = -2.5 * np.log10(f) - 48.6
    return mag

# Filter directory
filtDir = '/mnt/zfsusers/varadaraj/lephare/lephare_dev/filt/myfilters/artificial/'

# Filter to use
filtName = '1500_100.txt'
#filtName = '1600_100.txt'

# Define the subdirectory and field for the following SED directory.
#subDir = 'vis_check1'
#subDir = 'orig'
subDir = 'inclusive'
#subDir = 'secure'
#subDir = 'strict'
#subDir = 'vis_check2'
#subDir = 'DR3'

# Error subdirs
#subDir = 'primary'
#subDir = 'lya'
#subDir = 'goldrush'
#
field = 'XMM'
#field = 'CDFS'

# SED directory
#sedDir = '/mnt/zfsusers/varadaraj/sed_fitting/{0}/fitting/final/{1}/'.format(field, subDir)
if field == 'XMM':
#    sedDir = '/mnt/zfsusers/varadaraj/sed_fitting/XMM/fitting/fixed_DUD/6primary/'
#    sedDir = '/mnt/zfsusers/varadaraj/sed_fitting/XMM/fitting/fixed_DUD/6inclusive/'
#    sedDir = '/mnt/zfsusers/varadaraj/sed_fitting/XMM/fitting/fixed_DUD/5.3lyaONLY/'
    sedDir =  '/mnt/zfsusers/varadaraj/sed_fitting/XMM/fitting/fixed_DUD/insecure/primary/'
#    sedDir = '/mnt/zfsusers/varadaraj/sed_fitting/XMM/fitting/fixed_DUD/TEST/'

#    sedDir = '/mnt/hoy/temporaryFilesROHAN/lephare/goldrush/'
#    sedDir = '/mnt/zfsusers/varadaraj/sed_fitting/XMM/fitting/fixed_DUD/4irac/'

if field == 'CDFS':
    #sedDir =  '/mnt/zfsusers/varadaraj/sed_fitting/{0}/fitting/4irac/'.format(field)
#    sedDir =  '/mnt/zfsusers/varadaraj/sed_fitting/{0}/fitting/final/vis_check2/'.format(field)
    #sedDir =  '/mnt/zfsusers/varadaraj/sed_fitting/{0}/fitting/5.2lyaONLY/'.format(field)
    sedDir =  '/mnt/zfsusers/varadaraj/sed_fitting/CDFS/fitting/insecure/primary/'

# Other directories
#sedDir = '/mnt/hoy/temporaryFilesROHAN/lephare/goldrush/'

# Read in model stuff, so that we can access the rest frame model.
modDir = '/mnt/zfsusers/varadaraj/lephare/lephare_work/lib_bin/'
modDoc = 'LIB_GAL.doc'

# Read in catalogue for plotting against fwhm
catDir = '/mnt/vardy/vardygroupshare/data/catalogues/{0}FULL/cutting/'.format(field)
catDir = '/mnt/vardy/vardygroupshare/data/catalogues/{0}FULL/'.format(field)

if field == 'XMM':
#    cat = 'XMMFULL_DR2_MASKVISTADET_YJ_1.8as_IRAC2.8as_2022_03_24_inclusive_FWHM.fits'
    cat = 'XMMFULL_DR2_MASKVISTADET_YJ_1.8as_IRAC2.8as_2022_10_26.fits'

if field == 'CDFS':
#    cat = 'CDFSFULL_DR2_MASKVISTADET_YJ_1.8as_IRAC2.8as_2022_08_18_inclusive_FWHM.fits'
    cat = 'CDFSFULL_DR2_MASKVISTADET_YJ_1.8as_IRAC2.8as_2022_08_18.fits'

hdu = fits.open(catDir + cat)
header = hdu[0].header
data = hdu[1].data

# Write to file switch
write = False

# Plot switch
plot = False

# Are we saving errors?
errors = False

# Define plotting colours based on field
if field == 'XMM':
        colour = 'C0'
if field == 'CDFS':
        colour = 'C1'

# Need to prepare all the other stupid LePhare output stuff
# Define start and end points of each section in .spec file. Depends on field because of no. of filters.
if field == 'XMM':
    ds_phot = 9
    de_phot = 20

    ds_mod = 201
    de_mod = -1

    ds_pz = 20
    de_pz = 201

    ds_param = 3
    de_param = 9

if field == 'CDFS':
    ds_phot = 9
    de_phot = 21

    ds_mod = 202
    de_mod = -1

    ds_pz = 21
    de_pz = 202

    ds_param = 3
    de_param = 9


# Column names
namesPhot=['phot', 'yerr', 'wlen', 'xerr', 'modelPhot', 'col6', 'col7']
namesModel = ['wlen', 'flux']
namesPDF = ['z', 'P']
namesParam = ['Type', 'Nline', 'Model', 'Library', 'Nband', 'Zphot', 'Zinf', 'Zsup', 'Chi2', 'PDF', 'Extlaw', 'EB-V', 'Lir', 'Age', 'Mass', 'SFR', 'SSFR']


# Load tophat filter
filtFile = np.genfromtxt(filtDir+filtName)
tophat = filtFile[:,1]
tophat_wlen = filtFile[:,0]

# Empty mag list to append to
M = []
zphots = []
m = []
IDs = []
fwhms = []

dust = []

Ymag = []
Jmag = []
YJ = []
YJsynth = []

# Error empty arrays
zsup = []
zinf = []
Msup = []
Minf = []

# Define the cosmology
H = 70
omegaM = 0.3
omegaV = 0.7

# Alternate cosmology (Planck)
#H = 69.6
#omegaM = 0.286

# Set up the astropy cosmology
cosmo = FlatLambdaCDM(H0=H, Om0=omegaM) #, Ode0=omegaV)

# Open file to write the mags to a table.
file = open('errorsMin_{0}_{1}.txt'.format(subDir, field), 'w')
#file = open('errors_{0}_{1}_inc.txt'.format(subDir, field), 'w')
#file = open('Muv_{0}_{1}_inc.txt'.format(subDir, field), 'w')
#file = open('Muv_altcos_{0}_lya.txt'.format(field), 'w')

#file = open('Muv_GOLDRUSH_{0}_{1}.txt'.format(subDir, field), 'w')
#file = open('MuvErrors_GOLDRUSH_{0}_{1}.txt'.format(subDir, field), 'w')

if errors == False:
    file.write('#ID	Muv Av')

# File for errors
if errors == True:
    file.write('#ID zsup zinf Muv_sup Muv_inf')

sedDir = '/mnt/zfsusers/varadaraj/sed_fitting/XMM/fitting/fixed_DUD/1initial/endsley/'

# Loop through files
for f in sorted(os.listdir(sedDir)):

    # Get the ID
    ID = f[2:-5].lstrip('0')

    # Print this
    print('--------------')
    print(ID)
    print('--------------')

    IDs = np.append(IDs, ID)

    # Find this object in the catalogue and its fwhm.
    print(data['ID' == int(ID)])

    # Get FWHM for plotting
    #fwhm = data['FWHM_IMAGE'][data['ID'] == int(ID)]
    #print('FWHM: ', fwhm)
    #fwhms = np.append(fwhms, fwhm)

    # And get the flux in Y to convert to an apparent magnitude
    phot = ascii.read(sedDir+f, format='basic', data_start=ds_phot, data_end=de_phot, delimiter=' ', names=namesPhot)

    yflux = data['flux_Y'][data['ID'] == int(ID)]
    jflux = data['flux_J'][data['ID'] == int(ID)]

    yjflux = (yflux + jflux) / 2

    yj = flux_to_mag(yjflux)
    YJ = np.append(YJ, yj)

    yjSynth = (phot['phot'][-4] + phot['phot'][-3]) / 2
    YJsynth = np.append(YJsynth, yjSynth)

    print(yflux)
    yflux = flux_to_mag(yflux)
    jflux = flux_to_mag(jflux)
    Ymag = np.append(Ymag, yflux)
    Jmag = np.append(Jmag, jflux)

    # Load model and parameters
    model = ascii.read(sedDir+f, format='basic', data_start=ds_mod, data_end=de_mod, delimiter=' ', names=namesModel)

    params = ascii.read(sedDir+f, format='basic', data_start=ds_param, data_end=de_param, delimiter=' ', names=namesParam)


    

    # Get Y and J mags
#    modY = phot['modelPhot'][-4]
#    modJ = phot['modelPhot'][-3]

    if ID == '595119':
        # Get redshifts and errors
        z = params[1][5]
        zInf = params[1][6]
        zSup = params[1][7]

        zlo = z - zInf
        zup = zSup - z

    else:

        # Get redshifts and errors
        z = params[0][5]
        zInf = params[0][6]
        zSup = params[0][7]

        zlo = z - zInf
        zup = zSup - z

    print('Redshift {0} + {1} - {2}'.format(z, zlo, zup))

    # Append these errors to array
    zinf = np.append(zinf, zlo)
    zsup = np.append(zsup, zup)

    # Get dust extinction
    Av = params['Type' == 'GAL-1']['EB-V'].round(1)
    dust = np.append(dust, Av)

#    if f == 'Id000138952.spec':
#        z = 7.177
#    if f == 'Id001548625.spec':
#        z = 7.084

    zphots = np.append(zphots, z)
    print(round(z, 2))

    # Extract each model by finding where wavelength jumps.
    model_wlens = np.split(model['wlen'], np.where(np.diff(model['wlen']) < 0)[0] + 1)
    model_fluxes = np.split(model['flux'], np.where(np.diff(model['wlen']) < 0)[0] + 1)


    # Convert model to flux
    SED = mag_to_flux(model_fluxes[0])
    wlen = model_wlens[0]

    if ID == '595119':
        SED = mag_to_flux(model_fluxes[1])
        wlen = model_wlens[1]

    # Get rest frame wavelength grid
    wlen_rest = wlen / (1+z)

    #---------------------------------------------------------------------------

    #                Compute M1500 with distance modulus

    #---------------------------------------------------------------------------

    # Luminosity distance
    DL = cosmo.luminosity_distance(z).value * 10 ** 6 # put into pc

    # And the uncertainty on DL
    upDL = cosmo.luminosity_distance(zSup).value * 10 ** 6
    loDL = cosmo.luminosity_distance(zInf).value * 10 ** 6

    # Filter in observed frame
    filterO = np.zeros(len(wlen))
    for i, lam in enumerate(wlen):
        if (lam > 1450.0*(1+z)) and (lam < 1550.0*(1+z)):
            filterO[i] = 1.0

    #plt.plot(wlen, filterO)
    #plt.show()
    #exit()

    # Filter in rest frame
    filterE = np.zeros(len(wlen))
    for i, lam in enumerate(wlen):
        if (lam > 1450.0) and (lam < 1550.0):
            filterE[i] = 1.0

    L = SED*(1+z)*4*np.pi*DL**2

    # Integrals in K correction
    I1 = simps(SED * filterO, wlen)
    I2 = simps(L * filterE, wlen_rest)

    # Each k correction component
    K1 = 2.5*np.log10(1+z)
    K2 = 2.5*np.log10(I2/I1)

    # Convolve SED with tophat filter at rest frame 1500. Gets flux density in T1500(1+z).
    conv = np.sum(SED * filterO) / np.sum(filterO)

    # Impose minimum error of 5%
    convUp = conv + (conv * 0.05)
    convLo = conv - (conv * 0.05)

    # Compute apparent magnitude
    m1500Up = -2.5*np.log10(convUp)-48.6
    m1500Lo = -2.5*np.log10(convLo)-48.6


    m1500 = -2.5*np.log10(conv)-48.6


    # Save apparent mag
    m = np.append(m, m1500)

    # print this as a check
    print(m1500.round(1))

    # Compute abs mag using distance modulus and the redshift correction
    M1500 = m1500 - 5*np.log10(DL/10) + K1 # - K2

    # K1 at each distance
    K1up = 2.5*np.log10(1+zSup)
    K1lo = 2.5*np.log10(1+zInf)

    # uncertainties
    Mup = m1500Up - 5*np.log10(upDL/10) + K1up
    Mlo = m1500Lo - 5*np.log10(loDL/10) + K1lo

    plus = Mlo - Mup
    minus = Mup - Mlo

    print('+/- {0} {1}'.format(round(plus,2), round(minus,2)))

    Msup = np.append(Msup, plus)
    Minf = np.append(Minf, minus)

    # print this!
    print(M1500.round(1))

    print('Plus/minus ', max([plus, minus]).round(3))

    # Add to M array for plotting below
    M = np.append(M, M1500)

    # Write ID and Muv to file
    if write:
        if errors == False:
            file.write('\n')
            file.write('{0}	{1} {2}'.format(ID, M1500.round(2), Av))

        if errors:
            file.write('\n')
            file.write('{0} {1} {2} {3} {4}'.format(ID, zup.round(2), zlo.round(2), plus.round(2), minus.round(2)))

# If I don't want to plot
if plot == False:
    exit()

###################################

# Print the range and mean, Muv
print('------------------ Muv information ---------------------------')
print(str(min(M).round(1)) + ' < Muv < ' + str(max(M).round(1)))
print('Mean value: ', np.mean(M).round(1))

# Print the range and mean, muv
print('------------------ muv information ---------------------------')
print(str(min(m).round(1)) + ' < muv < ' + str(max(m).round(1)))
print('Mean value: ', np.mean(m).round(1))

####################################

# Plot Muv distributions

diff = m - yj
diffSynth = m - YJsynth

n, bins, patches = plt.hist(diff, bins=12)
#plt.xlabel(r'$m_{1500}$ - YJ')
#plt.show()
plt.clf()
print(np.mean(m-yj))


idx = n.argmax()

mode = bins[idx]
print(mode)

plt.scatter(zphots, diff, label='true')
plt.scatter(zphots, diffSynth, label='synthetic')
plt.ylabel(r'$m_{1500}$ - $m_{YJ}$')
plt.xlabel('z')
plt.legend()
plt.show()

exit()
# Muv - z
#plt.scatter(zphots, M)
#plt.ylabel(r'$M_{1500}$')
#plt.xlabel(r'$z_{\mathrm{phot}}$')
#plt.show()


plt.scatter(M, fwhms, color=colour)
plt.xlabel(r'$M_{1500}$')
plt.ylabel('FWHM')
plt.show()

plt.scatter(Ymag, fwhms, color=colour)
plt.xlabel(r'$m_{Y}$')
plt.ylabel('FWHM')
plt.show()
#exit()

###############################
