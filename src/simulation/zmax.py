#!/usr/bin/env python3

"""
zmax.py

Shift the SED of each of my final sample of galaxies to find its zmax, to then get its Vmax.

Step 2 of the luminosity function!

This is the first piece of code to date where I will start using tabs instead of spaces.

Created: Wednesday 24th August 2022

"""
# ----------------------
#   IMPORT LIBRARIES
# ----------------------

import numpy as np
import os
from astropy.io import fits, ascii
from astropy.cosmology import WMAP9
from astropy.cosmology import FlatLambdaCDM
from astropy.table import Table
import matplotlib.pyplot as plt
from new_depth_codes import grid_depths
from astropy.constants import c
from scipy.interpolate import splev, splrep, interp1d
from scipy.integrate import simps
# ------------------
#  FUNCTION SETUP
# ------------------

def mag_to_flux(m):
	'''Convert mags and their errors to flux'''
	flux = 10**(-0.4*(m+48.6))
	return flux

def flux_to_mag(flux):
	'''Convert flux to mag'''
	mag = -2.5*np.log10(flux)-48.6
	return mag



# ------------------
#  DIRECTORY SETUP
# ------------------

# ------ Switches ------

# Check plots
plot_check = False

# Print checks
print_check = False

# Write to file switch
write = True

# Secure objects or all?
secure = True

# Test code on a flat spectrum
test = False
testFlux = 1e-28

# Which sample to run
pri = False
inc = True
lya = False

# ------ FIELD ------

field = 'XMM'
#field = 'CDFS'

#subDir = 'inclusive'
#subDir = 'orig'

if field == 'XMM':
    subDir = 'DR3' # xmm objects
if field == 'CDFS':
    subDir = 'vis_check2' # Where the final cdfs objects are

# --- COVERING FRACTION ---
covFrac = 0.83


# ------ ROOT DIRECTORIES ------

# Root directory for my home folder
root = '/mnt/zfsusers/varadaraj/'

# Directory for the SEDs, dependent on field
sedBase = 'sed_fitting/{0}/fitting/'.format(field)

if field == 'XMM':

	# Primary candidates
	if pri:
		sedDir = root + sedBase + 'fixed_DUD/' + '6primary/'

	# Inclusive objects
	if inc:
#		sedDir = root + sedBase + 'fixed_DUD/' + '6inclusive/'
		sedDir = root + sedBase + 'fixed_DUD/' + 'insecure/primary/'

	# Lyman alpha objects
	if lya:
		sedDir = root + sedBase + 'fixed_DUD/' + '5.3lyaONLY/'

if field == 'CDFS':

	# Primary candidates
	if pri:
		sedDir = root + sedBase + 'final/' + subDir + '/'

	# Inclusive objects
	if inc:
#		sedDir = root + sedBase + 'final/' + '6inclusive/'
		sedDir = root + sedBase + 'insecure/primary/'

        # Lyman alpha objects
	if lya:
		sedDir = root + sedBase + '5.2lyaONLY/'


# ------ .SPEC FILE SETUP ------

# Define start and end points of each section in .spec file. Depends on field because of no. of filters. TAKEN FROM plot_spec.py
if field == 'XMM':
	ds_mod = 201
	de_mod = -1

	ds_phot = 9
	de_phot = 20

if field == 'CDFS':
	ds_mod = 202
	de_mod = -1

	ds_phot = 9
	de_phot = 21

# Define names of columns in the .spec file
namesModel = ['wlen', 'flux']

# The start and end of parameter part is the same for either field. We need this to read in the redshift. Also from plot_spec.py
ds_param = 3
de_param = 9
namesParam = ['Type', 'Nline', 'Model', 'Library', 'Nband', 'Zphot', 'Zinf', 'Zsup', 'Chi2', 'PDF', 'Extlaw', 'EB-V', 'Lir', 'Age', 'Mass', 'SFR', 'SSFR']
namesPhot=['phot', 'yerr', 'wlen', 'xerr', 'modelPhot', 'col6', 'col7']

# ------ FILTERS ------

# Directory for filters, base
filtDir = '/lephare/lephare_dev/filt/myfilters/'

# Load the filters

#---VISTA---
Y = np.genfromtxt(root + filtDir + 'VISTA/' + 'VISTA_Y.txt')
J = np.genfromtxt(root + filtDir + 'VISTA/' + 'VISTA_J.txt')
H = np.genfromtxt(root + filtDir + 'VISTA/' + 'VISTA_H.txt')
Ks = np.genfromtxt(root + filtDir + 'VISTA/' + 'VISTA_Ks.txt')

#---HSC---
G = np.genfromtxt(root + filtDir + 'HSC/' + 'g_HSC.txt')
R = np.genfromtxt(root + filtDir + 'HSC/' + 'r_HSC.txt')
I = np.genfromtxt(root + filtDir + 'HSC/' + 'i_HSC.txt')
Z = np.genfromtxt(root + filtDir + 'HSC/' + 'z_HSC.txt')
y = np.genfromtxt(root + filtDir + 'HSC/' + 'y_HSC.txt')
nb816 = np.genfromtxt(root + filtDir + 'HSC/' + 'hsc_nb816.txt')
nb921 = np.genfromtxt(root + filtDir + 'HSC/' + 'hsc_nb921.txt')

#---VST---
u = np.genfromtxt(root + filtDir + 'VST/' + 'vst_u.txt')
g = np.genfromtxt(root + filtDir + 'VST/' + 'vst_g.txt')
r = np.genfromtxt(root + filtDir + 'VST/' + 'vst_r.txt')
i = np.genfromtxt(root + filtDir + 'VST/' + 'vst_i.txt')


# Filter arrays
if field == 'CDFS':
	filters = [u, g, r, i, G, R, I, Z, Y, J, H, Ks]
	filterNames = ['u', 'g', 'r', 'i', 'HSC-G', 'HSC-R', 'HSC-I', 'HSC-Z', 'Y', 'J', 'H', 'Ks']
#			0    1	  2    3     4        5        6        7       8    9    10   11
if field == 'XMM':
	filters = [G, R, I, nb816, Z, nb921, y, Y, J, H, Ks]
	filterNames = ['HSC-G', 'HSC-R_DR3', 'HSC-I', 'HSC-NB0816_DR3', 'HSC-Z_DR3', 'HSC-NB0921_DR3', 'HSC-Y', 'Y', 'J', 'H', 'Ks']
#			 0	   1	       2	      3              4	           5             6       7    8    9    10

Xfilt = [] # Filter wavelengths
Yfilt = [] # Filter transmission/response



# ------ FILTER PREPARATION ------

# Unpack the np genfromtxt filters, then pack them into a convenient array. E.g. (Xfilt[0], Yfilt[0]) in CDFS is VST-u.
for filter in filters:
	x = []
	y = []

	# Open the filter and append the values to these x and y arrays
	for i in range(0, len(filter)):
		x.append(filter[i][0])
		y.append(filter[i][1])

	# Append to the data structure
	Xfilt.append(x)
	Yfilt.append(y)


# The VISTA filters need to be normalised so that the maximum value is 1.
for i, y in enumerate(Yfilt):
	# Check to see if the greatest value goes above 1.
	if np.max(Yfilt[i]) > 1:
		Yfilt[i] = [value/100.0 for value in Yfilt[i]]

	# Normalise all filters to 1???
#	Yfilt[i] = Yfilt[i]/np.max(Yfilt[i])



# ------ CATALOGUES ------

# Load in catalogues
#catDir = '/mnt/vardy/vardygroupshare/data/catalogues/{0}FULL/cutting/'.format(field)
catDir = '/mnt/vardy/vardygroupshare/data/catalogues/{0}FULL/'.format(field)

if field == 'XMM':
#	cat = catDir + 'XMMFULL_DR2_MASKVISTADET_YJ_1.8as_IRAC2.8as_2022_10_26_6primary.fits' #.format(subDir)
	cat = catDir + 'XMMFULL_DR2_MASKVISTADET_YJ_1.8as_IRAC2.8as_2022_10_26.fits'

if field == 'CDFS':
#	cat = catDir + 'CDFSFULL_DR2_MASKVISTADET_YJ_1.8as_IRAC2.8as_2022_08_18_{0}.fits'.format(subDir)
	cat = catDir + 'CDFSFULL_DR2_MASKVISTADET_YJ_1.8as_IRAC2.8as_2022_08_18.fits'.format(subDir)

hdu = fits.open(cat)

header = hdu[0].header
cat    = hdu[1].data



# ------ DEPTHS ------

# Load in depth grids for Y+J

# Make list of tiles based on field name
tiles = ['{0}1'.format(field), '{0}2'.format(field), '{0}3'.format(field)]

# Depth directory
depthDir = '/mnt/vardy/vardygroupshare/HSC_SSP_DR3/depths/'

# Loop through tiles and load each grid.
for i, tileName in enumerate(tiles):

	# Grid directory
	gridDir = depthDir + tileName + '/grids/VIDEO/'

	# Grid file name
	grid = '{0}_YJ_2.0as_grid_depths.fits'.format(tileName)

        # Assign depth table to variable based on where we are in the loop.
	if i == 0:
		depth1 = Table.read(gridDir + grid)
	if i == 1:
		depth2 = Table.read(gridDir + grid)
	if i == 2:
		depth3 = Table.read(gridDir + grid)



# ------ PSF ------

psfDir = '/mnt/vardy/vardygroupshare/data/psf/'




# ------ COSMOLOGY ------

# Define the cosmology
H = 70
#H = 69.6

omegaM = 0.3
#omegaM = 0.286

# From omegaM we get omegaVacuum = 0.7

# Set up the astropy cosmology
cosmo = FlatLambdaCDM(H0=H, Om0=omegaM)



# ------ Vmax VARIABLES ------

# Survey area in square degrees
#if field == 'XMM':
#	fieldArea = 4.26
#	fieldArea = 3.39

#if field == 'CDFS':        	# ------ NOW DONE WITHIN THE LOOP
#	fieldArea = 4.08
#	fieldArea = 3.89

# Covering fraction of foreground objects?
#covFrac = 0.8
#fieldArea = fieldArea * covFrac

# Area of whole sky
skyArea = 41252.96


# ------ READ IN Muv ------

# Read
#if secure == False:
#uvMags = np.genfromtxt('Muv_{0}.txt'.format(field))
#if secure == True:
#uvMags = np.genfromtxt('Muv_secure_{0}.txt'.format(field))

if pri:
	uvMags = np.genfromtxt('Muv/Muv_{0}_{1}.txt'.format(subDir, field)) # primary

if inc:
	uvMags = np.genfromtxt('Muv/Muv_orig_{1}_inc.txt'.format(subDir, field)) # inclusive

if lya:
	uvMags = np.genfromtxt('Muv/Muv_orig_{1}_lya.txt'.format(subDir, field)) # lyman alpha

#uvMags = np.genfromtxt('Muv_altcos_{0}.txt'.format(field)) # Planck cosmology
#uvMags = np.genfromtxt('Muv_altcos_{0}_lya.txt'.format(field))


# Convert to useable format
uvIDs = uvMags[:,0]
uvIDs = [int(j) for j in uvIDs]
uvVals = uvMags[:,1]

# Get dust extinction as well
Av = uvMags[:,2]

# ------ OTHER STUFF ------

# Some arrays for saving things.

# Checks and plots
delta_z = []
sig2noise = []
magComp = []
magCat = []
magDiff = []

# Final values we want in the LF
IDarray = []
zArray = []
zmaxArray = []
MuvArray = []
VmaxArray = []

# Checking my magnitude issues
dY = []
dJ = []
dYJ = []

# Filters
Y = []
J = []

modelmagsY = []
modelmagsJ = []

# ----------------------------------
#   START ZMAX CALCULATION LOOPING
# ----------------------------------


# Loop through SED directory and begin the measuring and shifting
for f in sorted(os.listdir(sedDir)):

	# ------ GET IDs ------

	# Print the file name
	print('------------')

	# Get the ID from the filename
	ID = f[2:-5].lstrip('0')
	print(ID)


	# ------ LOAD SEDs ------
	model = ascii.read(sedDir+f, format='basic', data_start=ds_mod, data_end=de_mod, delimiter=' ', names=namesModel)

	# This file contains the stellar model and each galaxy model. Find each by where wavelength jumps
	model_wlens = np.split(model['wlen'], np.where(np.diff(model['wlen']) < 0)[0] + 1)
	model_fluxes = np.split(model['flux'], np.where(np.diff(model['wlen']) < 0)[0] + 1)

	# LePhare magnitudes to compare my convolved mags against
	phot = ascii.read(sedDir+f, format='basic', data_start=ds_phot, data_end=de_phot, delimiter=' ', names=namesPhot)

	# Model photometry magnitudes to compare to.
	modelMagY = phot['modelPhot'][-4]
	modelMagJ = phot['modelPhot'][-3]

	print('Model mag in Y: ', modelMagY)
	print('Model mag in J: ', modelMagJ)

	modelmagsY.append(modelMagY)
	modelmagsJ.append(modelMagJ)

	# The high-z SED is the first model
	sed = model_fluxes[0]
	wlen = model_wlens[0]

	# ------ TESTING AREA ------
	if test:

		# Convert SEDs into flat spectra to check magnitudes.

		sed = np.ones(len(sed))
		sed = sed * testFlux
		testMag = flux_to_mag(testFlux)
		print('TESTING MAGNITUDE: {0} in cgs, {1} mag'.format(testFlux, testMag))

	# Convert the SED into flux
	sed = mag_to_flux(sed)

	# Plot the SED and filters to check normalizations.
	if plot_check == True:
		plt.plot(wlen, sed/np.max(sed)) # Normalised for plotting.
		plt.plot(Xfilt[0], Yfilt[0])
		plt.plot(Xfilt[1], Yfilt[1])
		plt.plot(Xfilt[2], Yfilt[2])
		plt.plot(Xfilt[3], Yfilt[3])
		plt.plot(Xfilt[4], Yfilt[4])
		plt.plot(Xfilt[5], Yfilt[5])
		plt.plot(Xfilt[6], Yfilt[6])
		plt.plot(Xfilt[7], Yfilt[7])
		plt.plot(Xfilt[8], Yfilt[8])
		plt.plot(Xfilt[9], Yfilt[9])
		plt.plot(Xfilt[10], Yfilt[10])
		if field == 'CDFS':
			plt.plot(Xfilt[11], Yfilt[11])
		plt.show()



	# ------ LOAD SED PARAMETERS -------

	params = ascii.read(sedDir+f, format='basic', data_start=ds_param, data_end=de_param, delimiter=' ', names=namesParam)

	# Get the redshift
	z_obj = params[0][5]


	# ------ LOAD CATALOGUE DATA ------

	# Find the index of this current ID (the catalogue isn't necessarily in ID order)
	idxID = np.where(cat['ID'] == int(ID))[0][0]

	# Get the x and y of this object
	x = cat[idxID]['X_IMAGE']
	y = cat[idxID]['Y_IMAGE']

#	if print_check:
#		print(x, y)

	# Get the VISTA tile the object lies in
	tile = cat[idxID]['VISTA_tile_used']

#	if print_check:
#		print(tile)


	# Get the HSC pointing the circle lies in for the area! Only for CMM.
	if field == 'XMM':
		circle = cat[idxID]['HSC_circle']

		# Ultradeep pointing
		if (circle == '1') | (circle == '2U'):
			fieldArea = 4.33 #1.79
		if (circle == '2L') | (circle == '3') | (circle == '2L3') | (circle == '2U3'):
			fieldArea = 4.33 # Objects in DEEP will be seen in UDEEP

	# CDFS area
	if field == 'CDFS':
		#fieldArea = (3.29 + 4.33) # Objects in CDFS will be seen in XMM, they are bright enough!
		fieldArea = (3.89 + 4.33) # Objects in CDFS will be seen in XMM, they are bright enough!
		#fieldArea = 3.89 # Only cdfs

	fieldArea = fieldArea * covFrac


	# ------ DEPTHS ------

	# Now find the correct depth table we loaded earlier.
	if tile[-1] == '1':
		depthTable = depth1

	if tile[-1] == '2':
		depthTable = depth2

	if tile[-1] == '3':
		depthTable = depth3

	# Get the depth!
	depth = grid_depths(depthTable, x, y)[0]


	# ------ INTERPLOTING SED TO FILTER WAVELENGTH GRID ------

	# We need to interpolate the filter to the SED wavelength grid.
	# The filters have a much higher resolution than the SED. Can see this in the plot above.

	# ^^^^NOOOOOO!!!!! IT'S THE OTHER WAY AROUND!!!!!

	# Y and J are the 4th and 3rd from last filters in both fields.

	interpY = np.interp(Xfilt[-4], wlen, sed)
	interpJ = np.interp(Xfilt[-3], wlen, sed)

	# Plot these lower-res filters.
	if plot_check == True:
		plt.plot(Xfilt[-4], interpY/np.max(interpY), label='SED, interpY')
		plt.plot(Xfilt[-3], interpJ/np.max(interpJ), label='SED, interpJ')

		plt.plot(Xfilt[-4], Yfilt[-4], label='Y')
		plt.plot(Xfilt[-3], Yfilt[-3], label='J')

		plt.plot(wlen, 75*sed/np.max(sed), label='SED, arbitrarily scaled') # Normalised for plotting
		plt.legend()
		plt.show()


	# ------ SETTING UP REDSHIFT SHIFTING ------

	# Initial redshift
	z0 = z_obj

	# Redshift step
	dz = 0.01

	# Maximum redshift we will go up to
	zmax = 7.50
#	if field == 'CDFS':
#		zmax = 7.30

	# Minimum redshift we are interested in
	zmin = 6.50

	# Get luminosity distance for z0
	DL0 = cosmo.luminosity_distance(z0)

	# ----- MEASURING INITIAL FLUX IN FILTER ------


	# Initial convolving of filter and SED, normalised by filter area.

	# Convert to frequency space. This is also used in the loop.
	YfreqX = [ (c.value / (p * 1e-10)) for p in Xfilt[-4]]
	JfreqX = [ (c.value / (p * 1e-10)) for p in Xfilt[-3]]

	# Frequency space
#	Yconv0 = np.trapz((interpY * Yfilt[-4]), YfreqX) / np.trapz(Yfilt[-4], YfreqX) # Trapezium rule
#	Jconv0 = np.trapz((interpJ * Yfilt[-3]), JfreqX) / np.trapz(Yfilt[-3], JfreqX)

	Yconv0 = simps((interpY * Yfilt[-4]), YfreqX) / simps(Yfilt[-4], YfreqX) # Simpson's rule
	Jconv0 = simps((interpJ * Yfilt[-3]), JfreqX) / simps(Yfilt[-3], JfreqX)

	# Print initial magnitudes to compare against LePhare model.
	mY0 = flux_to_mag(Yconv0)
	mJ0 = flux_to_mag(Jconv0)

	print('mY0: ', mY0.round(2))
	print('mJ0: ', mJ0.round(2))

	# --- Getting mYJ: add in flux space ---
	flux = (Yconv0 + Jconv0) / 2

	mYJ0 = flux_to_mag(flux)

	if print_check:
		print('mYJ0: ', mYJ0.round(2))
#		print('Depth here: ', depth.round(2))



	# Print the starting redshift and mYJ.
	print(str(z_obj.round(2)) + ': ' + str(round(mYJ0, 2)))


	############## HERE IS WHERE WE CHECK IF XMM OBJECT CAN BE FOUND IN CDFS #############
	if field == 'XMM' and mYJ0 <=24.71:
		fieldArea = fieldArea + 3.98
		print('Including CDFS area')

	# ------ CHECKING OUTPUTTED MAGNITUDES ------

	# ------ CHECK: mYJ from the catalogue flux.

	# Get flux and convert to mag
	catfluxYJ = cat[idxID]['flux_YJ']
	catmagYJ = flux_to_mag(catfluxYJ)


	# APPEND FOR CHECKING
	deltaMag = catmagYJ - mYJ0
	magDiff.append(deltaMag)

	# Get error
	caterrYJ = cat[idxID]['err_YJ']

	# Get signal-to-noise
	SN = catfluxYJ / caterrYJ

	sig2noise.append(SN)

	# My convolution minus LePhare
	dY.append(mY0 - modelMagY)
	dJ.append(mJ0 - modelMagJ)

	# ------------- REDSHIFT LOOPING --------------

	# First, define initial mYJ as mYJ0
	mYJ = mYJ0

	# Change accuracy of z to 2 d.p.
	precision = 0.01
	z = z0.round(2)

	# Store some values
	seds = []
	wlens = []
	redshifts = []
	appMags = []

	# Begin the loop!

	while (mYJ < depth) and (z <= zmax):

		# Increment z (we have already computed mYJ for z0, and don't need to shift anything for it).
		z = z + dz

		# First get the luminosity distance
		DL = cosmo.luminosity_distance(z)

		# Define a redshift ratio
		#zratio = (1+z0) / (1+z)

		# Shift the wavelength grid
		wlen = wlen * ((1+z) / (1+z0))

		# Shift the SED! Need to convert back into magnitude here. Same equation from Muv, but ratios of z_(i) and z_(i+1)
		sed = flux_to_mag(sed)
		sed = sed - (5 * np.log10(DL0/DL)) + (2.5 * np.log10( (1+z) / (1+z0) ))

		seds.append(sed)
		wlens.append(wlen)

		# Convert SED back into flux to convolve
		sed = mag_to_flux(sed)


		# Interpolate SED to filter wavelength grid
		interpY = np.interp(Xfilt[-4], wlen, sed)
		interpJ = np.interp(Xfilt[-3], wlen, sed)


		# ----- MEASURING FLUX IN FILTER ------

		# Initial convolving of filter and SED, normalised by filter area.

#		Yconv = np.trapz((interpY * Yfilt[-4]), YfreqX) / np.trapz(Yfilt[-4], YfreqX) # Trapezium rule in frequency space
#		Jconv = np.trapz((interpJ * Yfilt[-3]), JfreqX) / np.trapz(Yfilt[-3], JfreqX)

		Yconv = simps((interpY * Yfilt[-4]), YfreqX) / simps(Yfilt[-4], YfreqX) # Simpson's rule 
		Jconv = simps((interpJ * Yfilt[-3]), JfreqX) / simps(Yfilt[-3], JfreqX)
		# ------ COMPUTE MAGNITUDE ------


		# --- Method 1: add in mag space ---

		# And compute the observed magnitudes at this redshift

		'''
		mY = flux_to_mag(Yconv)
		mJ = flux_to_mag(Jconv)


		# We need to get mYJ. See notes in my notebook for derivation.
		mYJ = -2.5 * np.log10( 10 ** (-0.4 * mY) + 10 ** (-0.4 * mJ) )
		'''

		# --- Method 2: add in flux space ---
		flux = (Yconv + Jconv) / 2

		mYJ = flux_to_mag(flux)


		# Print redshift and magnitude
		#print(str(round(z, 2)) + ': ' + str(mYJ))

		# Save redshift and magnitude
		redshifts.append(z)
		appMags.append(mYJ)


	# Print information after exiting the while loop
#	print('------')
	print('Final mYJ: ', mYJ.round(2))
	print('Depth: ', round(depth, 2))
	print('Zfinish: ',round(z/precision)*precision)
	print('delta-z: ', round(z-z0, 2))
	print('SO ZMAX is ', redshifts[-2].round(2), ' with mYJ ', appMags[-2].round(2), ' against depth ', depth.round(2))
	print('------')

	# Save redshift delta
	delta_z.append(z-z0)



	# ------ COMPUTE Vmax -------

	# V = (4pi/3) * (r1 ** 3 - r2 ** 3) # I THINK THIS IS WRONG
	# IT SHOULD BE:
	# V = (solid angle of field) / 3 * ( (comoving distance at zmax)**3 - (comoving distance at z=6.5)**3 )

	# Convert field area to steradians
	fieldArea = fieldArea * (np.pi/180)**2
	skyArea = skyArea * (np.pi/180)**2

	# Object's Vmax
	#Vmax =      (4 * np.pi /3) * (fieldArea / skyArea) * ( cosmo.comoving_distance(z) ** 3 - cosmo.comoving_distance(6.5) ** 3 )
	Vmax =      fieldArea/3 * ( cosmo.comoving_distance(z) ** 3 - cosmo.comoving_distance(6.5) ** 3 )

	# Maximum Vmax an object can have
	#Vuniverse = (4 * np.pi /3) * (fieldArea / skyArea) * ( cosmo.comoving_distance(7.5) ** 3 - cosmo.comoving_distance(6.5) ** 3 )

	Vuniverse = fieldArea/3 * ( cosmo.comoving_distance(7.5) ** 3 - cosmo.comoving_distance(6.5) ** 3 )

	V = 8.22 * (np.pi/180)**2 / 3 * ( cosmo.comoving_distance(7.5) ** 3 - cosmo.comoving_distance(6.5) ** 3 )

	print('VUNIVERSE', V)

	#Vcheck = 4*np.pi/3 * ( cosmo.comoving_distance(7.5) ** 3 - cosmo.comoving_distance(6.5) ** 3 )

	# Make sure objec Vmax doesn't exceed the maximum for our redshift range.
	Vmax = min(Vmax, Vuniverse)

	print('Vmax: ', Vmax)

	#print('Vuniverse: ', Vcheck)

	# ------ COLLECT Muv ------

	idxUV = uvIDs.index(int(ID))

	Muv = uvVals[idxUV]


	# ------ SAVE VALUES ------
	IDarray.append(ID)
	zArray.append(z0)
	zmaxArray.append(round(z, 2))
	MuvArray.append(Muv)
	VmaxArray.append(Vmax.value) # In Mpc3


# Save data to file
data = {'ID':IDarray, 'z':zArray, 'zmax':zmaxArray, 'Muv':MuvArray, 'Vmax':VmaxArray, 'Av':Av}

if write:

#	ascii.write(data, 'vmax_{0}_lya.txt'.format(field), format='commented_header', overwrite=False) # For each field, separate calc

	if pri:
		ascii.write(data, 'vmax_{0}_{1}.txt'.format(subDir, field), format='commented_header', overwrite=True)

	if inc:
		ascii.write(data, 'vmax_inclusive_{0}.txt'.format(field), format='commented_header', overwrite=True) # For inclusive objects

	if lya:
		ascii.write(data, 'vmax_{0}_lya.txt'.format(field), format='commented_header', overwrite=True) # For lya objects

#	ascii.write(data, 'vmax_altcos_{0}_lya.txt'.format(field), format='commented_header', overwrite=True)


exit()

print(zArray)
# Do some plotting
#plt.hist(delta_z, bins=5)
#plt.show()

plt.scatter(sig2noise, delta_z)
plt.xlabel('Signal to noise ratio')
plt.ylabel('zmax - z0')
plt.show()
plt.clf()

#plt.hist(magDiff)
#plt.xlabel('mag computed from CGS - mag computed from SED')
#plt.show()

#plt.clf()

plt.scatter(MuvArray, VmaxArray)
plt.xlabel('Muv')
plt.ylabel('Vmax [Mpc3]')
plt.title('{0}'.format(field))
plt.show()
plt.clf()


plt.hist(dY)
plt.xlabel('convolved mag - LePhare mag')
plt.title('Y {0}'.format(field))
plt.show()
plt.clf()

plt.hist(dJ)
plt.xlabel('convolved mag - LePhare mag')
plt.title('J {0}'.format(field))
plt.show()
plt.clf()

