#!/usr/bin/env python3

"""
luminosity_function.py

Takes the output of zmax.py and uv_magnitude.py to compute the luminosity function!!

Created: Friday 26th August 2022

"""

# ------ IMPORT LIBRARIES ------

import numpy as np
from astropy.table import Table, vstack
import matplotlib.pyplot as plt

import matplotlib.ticker as tck
from matplotlib.ticker import MultipleLocator, AutoMinorLocator, LogLocator, FormatStrFormatter, LogFormatter

plt.rcParams.update({'font.size': 15})
plt.rcParams.update({'font.size': 40}) # FOR z Muv plot
plt.rcParams['axes.linewidth'] = 4
plt.rcParams['figure.dpi'] = 100

# ------ SETUP AND DIRECTORIES ------

field = 'XMM'

subDir = 'vis_check2'
#subDir = 'inclusive'

# plot upper limits? Set to false when doing fields separately or for the inclusive sample
uplim = True

# ---- DOUBLE POWER LAW ----
def dpl(phiStar, alpha, beta, M, Mstar):

    ''' Produces DPL given parameters alpha, beta, normalisation phi*, char. mag. M* and mag array M. '''

    numerator = phiStar

    denomA = 10 ** (0.4 * (alpha + 1) * (M - Mstar))
    denomB = 10 ** (0.4 * (beta + 1)  * (M - Mstar))

    denominator = denomA + denomB

    phi = numerator/denominator

    return phi

def schechter(phiStar, alpha, M, Mstar):

    ''' Produces Schechter function given parameters: normalisation phi*, faint slope alpha, mag array M and char. mag M*. '''

    coeff = np.log(10) / 2.5

    faint = (10 ** (0.4 * (Mstar - M))) ** (alpha+1)


    bright_exponent = -10 ** (0.4 * (Mstar - M))
    bright = np.exp(bright_exponent)

    phi = coeff * phiStar * faint * bright

    return phi

# ------ READ IN VMAX stuff ------

# Primary candidates
#dataXMM  = Table.read('vmax/vmax_DR3_XMM.txt'.format(subDir), format='ascii.commented_header')
#dataCDFS = Table.read('vmax/vmax_{0}_CDFS.txt'.format(subDir), format='ascii.commented_header')

# Each field
#dataXMM  = Table.read('vmax/vmax_XMM_ONLY.txt'.format(subDir), format='ascii.commented_header')
#dataCDFS = Table.read('vmax/vmax_CDFS_ONLY.txt'.format(subDir), format='ascii.commented_header')

# Inclusive candidates (3.98) NOW WITH COMPLETENESS!!!
#dataXMM  = Table.read('vmax_inclusive_XMM_comp.txt', format='ascii.commented_header')
#dataCDFS = Table.read('vmax_inclusive_CDFS_comp.txt', format='ascii.commented_header')

# Primary + lyman alpha (3.98)!!!!!!!!!!!!!! # NOW WITH COMPLETENESS! yay
dataXMM  = Table.read('vmax_XMM_comp.txt', format='ascii.commented_header')
dataCDFS = Table.read('vmax_CDFS_comp.txt', format='ascii.commented_header')

# Alternate cosmology
#dataXMM  = Table.read('vmax/vmax_altcos_XMM.txt'.format(subDir), format='ascii.commented_header')
#dataCDFS = Table.read('vmax/vmax_altcos_CDFS.txt'.format(subDir), format='ascii.commented_header')

# ------ Reading in some other data --------

# Bouwens 2021
dataBouwens = Table.read('../ref_catalogues/bouwens21_z7.dat', format='ascii.commented_header')

# Errors on my candidate Muv-z
errorXMM = Table.read('errorsMin_primary_XMM.txt', format='ascii.commented_header')
errorCDFS = Table.read('errorsMin_primary_CDFS.txt', format='ascii.commented_header')

#data = dataCDFS
#data = dataXMM

data = vstack([dataXMM, dataCDFS])

print(data)

print(data['C'])
#print(len(data[0][:]))

# Take the Muv data to make bins
Muv = data['Muv']
#print(np.min(Muv), np.max(Muv))

# MEAN REDSHIFT
print(np.mean(data['z']))
#exit()
'''
# ------ Some first year report plotting ------

# All candidates, rather than just the secure ones
moreDir = 'inclusive'

# All candidates
moreXMM  = Table.read('vmax/vmax_{0}_XMM.txt'.format(moreDir), format='ascii.commented_header')
moreCDFS = Table.read('vmax/vmax_{0}_CDFS.txt'.format(moreDir), format='ascii.commented_header')

more = vstack([moreXMM, moreCDFS])
moreMuv = more['Muv']
print(np.min(moreMuv), np.max(moreMuv))

# Secure histogram
plt.figure(figsize=(23,20))
n, bins, patches = plt.hist(Muv, bins=9)
binW = abs(bins[1]-bins[0]).round(2)
#print(binW)
plt.xlabel(r'$M_{UV}$')
plt.ylabel('Count')
plt.show()
plt.clf()

# Inclusive histogram
plt.figure(figsize=(23,20))
n, bins, patches = plt.hist(moreMuv, bins=9, color='red')
binW = abs(bins[1]-bins[0]).round(2)
print(binW)
plt.xlabel(r'$M_{UV}$')
plt.ylabel('Count')
#plt.show()
plt.clf()


# Secure scatter
plt.figure(figsize=(18,12))

# Get errors
XMM_zerr = np.array(list(zip(errorXMM['zinf'], errorXMM['zsup']))).T
CDFS_zerr = np.array(list(zip(errorCDFS['zinf'], errorCDFS['zsup']))).T

XMM_Merr = np.array(list(zip(errorXMM['Muv_inf'], errorXMM['Muv_sup']))).T
CDFS_Merr = np.array(list(zip(errorCDFS['Muv_inf'], errorCDFS['Muv_sup']))).T

# Candidates
#plt.errorbar(dataXMM['z'], dataXMM['Muv'], xerr=XMM_zerr, yerr=XMM_Merr, label='XMM-LSS', color='purple', alpha=0.8, linestyle='none', marker='o', markersize=10)
#plt.errorbar(dataCDFS['z'], dataCDFS['Muv'], xerr=CDFS_zerr, yerr=CDFS_Merr, label='ECDF-S', color='orange', marker='s', alpha=0.8, linestyle='none', markersize=10)

#print(errorXMM['zsup'][:-4])

plt.errorbar(dataXMM['z'][:-4], dataXMM['Muv'][:-4], xerr=(errorXMM['zinf'][:-4], errorXMM['zsup'][:-4]), yerr=(np.abs(errorXMM['Muv_inf'][:-4]), errorXMM['Muv_sup'][:-4]), label='XMM-LSS', color='red', alpha=0.8, linestyle='none', marker='o', markersize=20, elinewidth=3.5)
plt.errorbar(dataXMM['z'][-4:], dataXMM['Muv'][-4:], xerr=(errorXMM['zinf'][-4:], errorXMM['zsup'][-4:]), yerr=(np.abs(errorXMM['Muv_inf'][-4:]), errorXMM['Muv_sup'][-4:]), label='XMM-LSS', color='red', alpha=0.8, linestyle='none', marker='d', markersize=20, elinewidth=3.5) #, mfc='white')

#plt.scatter(dataXMM['z'], dataXMM['Muv'], s=250, alpha=1, color='red')

print(errorCDFS['Muv_inf'], errorCDFS['Muv_sup'])

plt.errorbar(dataCDFS['z'], dataCDFS['Muv'], xerr=(errorCDFS['zinf'], errorCDFS['zsup']), yerr=(np.abs(errorCDFS['Muv_inf']), errorCDFS['Muv_sup']), label='ECDF-S', color='blue', marker='s', alpha=0.8, linestyle='none', markersize=20, elinewidth=3.5)

# Adding Bouwens
plt.scatter(dataBouwens['col11'], dataBouwens['Muv'], label='Bouwens+21', s=150, color='gray', alpha=0.8)

# Adding Harikane
goldrushGals = {6.47657:-22.1, 6.33213:-21.7, 7.34258:-24.1, 7.42788:-24.0, 6.67466:-21.6, 6.55:-22.6, 6.67154:-21.7, 6.66697:-22.2}
lists = sorted(goldrushGals.items()) # sorted by key, return a list of tuples
goldrushz, goldrushMuv = zip(*lists) # unpack a list of pairs into two tuples

plt.scatter(goldrushz, goldrushMuv, label='Harikane+22', s=500, color='gold', alpha=0.9, marker='p')

# Inclusive scatter
#plt.scatter(moreXMM['z'], moreXMM['Muv'], facecolors='none', edgecolors='green', s=100)
#plt.scatter(moreCDFS['z'], moreCDFS['Muv'], facecolors='none', edgecolors='blue', s=100)
plt.ylim([-24.5, -20])
plt.xlim([6.1, 7.5])
plt.ylabel(r'$M_{\mathrm{UV}}$')
plt.xlabel(r'$z$')
plt.gca().invert_yaxis()
#plt.legend(fontsize=)
#plt.show()
plt.savefig('/mnt/vardy/vardygroupshare/HSC_SSP_DR3/plots/Muv_v_zphot.pdf', bbox_inches='tight')
exit()


# Secure z-hist
plt.figure(figsize=(15,15))
plt.hist(data['z'])
plt.xlabel(r'$z$')
plt.ylabel('Count')
plt.show()
plt.clf()

# Inclusive z-hist
plt.figure(figsize=(15,15))
plt.hist(more['z'], color='red')
plt.xlabel(r'$z$')
plt.ylabel('Count')
plt.show()
plt.clf()


# Inclusive dust scatter
Av = more['Av']
plt.scatter(moreMuv, Av)
plt.show()
exit()
'''
#----------------------------------------------

# Number of bins
nbins = 4

# Uniform 0.75 mag width, 4 bins. THIS IS THE ONE
binEnds = [-24.8, -24.05, -23.3, -22.55, -21.8]
binErr = [0.375, 0.375, 0.375, 0.375]


# Inclusive candidate bins. THIS IS THE ONE
#binEnds = [-24.3, -23.55, -22.8, -22.3, -21.8]
#binErr = [0.375, 0.375, 0.25, 0.25]

#binEnds = [-24.5, -23.75, -23.0, -22.5, -22.0]
#binErr = [0.375, 0.375, 0.25, 0.25]

# Inclusive candidate bins, making sure the -24.1 object is in middle of its bin
#binEnds = [-24.475, -23.725, -22.975, -22.225, -21.725]
#binErr = [0.375, 0.375, 0.375, 0.25]

# Different sizes
#binEnds = [-24.65, -23.65, -22.9, -22.15, -21.4]


# CDFS bins
#binEnds = [-24.0, -23.0, -22.0]
#binErr = [0.5, 0.5]

# XMM bins
#binEnds = [-24.6, -23.6, -22.6, -22.1, -21.6]
#binErr = [0.5, 0.5, 0.25, 0.25]

#n, bins, patches = plt.hist(Muv, bins=nbins)
n, bins, patches = plt.hist(Muv, bins=binEnds)
#plt.show()

print(bins)
#exit()

# Split data into these bins
bin1 = np.where((bins[0] <= data['Muv']) & (data['Muv'] < bins[1]))[0]
bin2 = np.where((bins[1] <= data['Muv']) & (data['Muv'] < bins[2]))[0]
bin3 = np.where((bins[2] <= data['Muv']) & (data['Muv'] < bins[3]))[0]
bin4 = np.where((bins[3] <= data['Muv']) & (data['Muv'] < bins[4]))[0]
if nbins == 5:
	bin5 = np.where((bins[4] <= data['Muv']) & (data['Muv'] < bins[5]))[0]

#print(bin1, bin2, bin3)

# Split
brightest = data[:][bin1]
bright = data[:][bin2]
middle = data[:][bin3]
faint = data[:][bin4]
if nbins == 5:
	faintest = data[:][bin5]

databins = [brightest, bright, middle, faint] # Same order as on plot
#databins = [brightest, bright] # CDFS

if nbins == 5:
	databins = [brightest, bright, middle, faint, faintest] # Same order as on plot

bin_width = abs(bins[1] - bins[0]).round(2)
#print(bin_width)

# Initialise LF sum
LFsum = np.zeros(len(databins))
LFerror = np.zeros(len(databins))

# And a Harikane contaminated version
LFcontam = np.zeros(len(databins))

# Dummy completeness value
#completeness = [0.5, 0.5, 0.5, 0.5]
completeness = np.ones(len(databins))
#completeness = [1, 1, 1, 1, 0.5]



# Approximate Harikane contamination factor
contamination = [0.5, 0.5, 0.5, 0.5] # 4/8
#contamination = [0.5, 0.5, 0.5, 0.5, 0.5] # 4/8
#contamination = [0.375, 0.375, 0.375, 0.375] # 5/8

for i, split in enumerate(databins):

	for j, obj in enumerate(split):

#		print('{0}, {1}'.format(j, obj))

		# Compute LF term for this object
		#val = 1 / (completeness[i] * obj['Vmax'])
		val = 1 / (obj['C'] * obj['Vmax'])

		# Contamination factor
		contamVal = 1 / (obj['C'] * contamination[i] * obj['Vmax'])

		# Add to sum
		LFsum[i] = LFsum[i] + val
		LFcontam[i] = LFcontam[i] + contamVal

		# Compute error term too
		errval = 1 / (obj['Vmax']) ** 2

		# Add to sum
		LFerror[i] = LFerror[i] + errval

	# Upper lims
	if uplim:
		if len(split) == 0:
			LFsum[i] = 1/64304431.93531033/(binErr[i]*2) # Standard cosmology
			#LFsum[i] = 1/69236712.74110535/(binErr[i]*2) # Planck
			LFcontam[i] = 1/64304431.93531033/(binErr[i]*2)

# Multiply LF sum by 1/bin width
#LFsum = [(1/bin_width) * val for val in LFsum]

# Square root error term and divide by bin width
#LFerror = [ (1/bin_width) * np.sqrt(val) for val in LFerror]

#print(LFsum)
#print(LFerror)

for i, lf in enumerate(LFsum):

	lf = 1/(binErr[i]*2) * lf
#	lf = 1/(bin_width) * lf
	LFcontam[i] = 1/(binErr[i]*2) * LFcontam[i]

	LFerror[i] = 1/(binErr[i]*2) * np.sqrt(LFerror[i])
#	LFerror[i] = 1/(bin_width) * np.sqrt(LFerror[i])

# UPPER LIM
if uplim:
	LFerror[0] = 1e-8


# Use crude contamination rate from Harikane
#contam = LFsum * 8/3
#contam = LFsum * 2
contam = LFcontam

# PRINTING VALUES
print('LF VALUES')
print('POINTS: ', LFsum)
print('ERRORS: ', LFerror)
print('CONTAMINATED LF VALUES')
print(contam)


# Muv plotting values
xbin_plot = [ (bins[1] + bins[0])/2, (bins[2] + bins[1])/2, (bins[3] + bins[2])/2, (bins[4] + bins[3])/2 ]
#xbin_plot = [ (bins[1] + bins[0])/2, (bins[2] + bins[1])/2] # CDFS

if nbins == 5:
	xbin_plot = [ (bins[1] + bins[0])/2, (bins[2] + bins[1])/2, (bins[3] + bins[2])/2, (bins[4] + bins[3])/2, (bins[5] + bins[4])/2 ]
print(xbin_plot)

# ------ EXISTING DATA ------

# --- Bowler+ 2014 ---

rb14x = [-22.66, -22.17, -21.75]
rb14y = [3.59e-7, 1.69e-6, 4.10e-6]

rb14dx = [0.25, 0.25, 0.25]
rb14dyU = [6.13e-7, 2.38e-6, 5.47e-6]
rb14dyL = [1.05e-7, 1e-6, 2.73e-6]

#rb14err = np.array(list(zip(rb14dyL, rb14dyU)))
rb14err = [rb14dyL, rb14dyU]

# --- Bowler+ 2017 ---
rb17x = [-22.86, -22.40, -21.85]
rb17y = [3.59e-7, 1.16e-6, 4.25e-6]

rb17dx = [0.25, 0.25, 0.25]
rb17dyU = [2.54e-7, 0.58e-6, 1.04e-6]
rb17dyL = rb17dyU

rb17err = [rb17dyL, rb17dyU]

# --- Hairkane+ 2022 ---

yh22x = [-24.92, -24.42, -23.92, -23.42, -22.92, -22.42, -21.92]
#yh22y = [5e-10, 1.31e-8, 4.39e-8, 1.83e-7, 1.06e-6, 2.75e-6] # 'Galaxy LF'
yh22y = [8.89e-9, 2.41e-8, 9.02e-8, 1.62e-7, 4.63e-7, 1.95e-6, 3.47e-6] # Dropout LF

yh22dx = [0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25]
yh22dyU = [15.29e-9, 1.75e-8, 2.13e-8, 0.3e-7, 7.53e-7, 1.13e-6, 1.7e-6]
yh22dyL = [7.41e-9, 0.98e-8, 1.66e-8, 0.23e-7, 2.71e-7, 0.67e-6, 1.03e-6]

yh22y_contam = [val * 0.5 for val in yh22y]

#yh22err = np.array(list(zip(yh22dyL, yh22dyU)))
yh22err = [yh22dyL, yh22dyU]


# --- Bouwens+ 2021 ---

b21x = [-22.19, -21.69, -21.19, -20.69, -20.19, -19.69, -19.19, -18.69, -17.94, -16.94]

b21y = [1e-6, 4.1e-5, 4.7e-5, 1.98e-4, 2.83e-4, 5.89e-4, 1.172e-3, 1.433e-3, 5.760e-3, 8.320e-3]

b21dy = [2e-6, 1.1e-5, 1.5e-5, 3.6e-5, 6.6e-5, 1.26e-4, 3.36e-4, 4.19e-4, 1.440e-3, 2.9e-3]

b21dx = np.ones(len(b21x)) * 0.25


# --- Finkelstein+ 2015 ---
f15x = [-22.0, -21.5, -21.0, -20.5, -20.0, -19.5, -19.0, -18.5, -18.0]

f15y = [0.0046, 0.0187, 0.0690, 0.1301, 0.2742, 0.3848, 0.5699, 2.5650, 3.0780]

f15dyU = [0.0049, 0.0085, 0.0156, 0.0239, 0.0379, 0.0633, 0.2229, 0.8735, 1.0837]

f15dyL = [0.0028, 0.0067, 0.0144, 0.02, 0.0329, 0.0586, 0.1817, 0.7161, 0.8845]

f15y = np.asarray(f15y) * 10**(-3)
f15dyU = np.asarray(f15dyU) * 10**(-3)
f15dyL = np.asarray(f15dyL) * 10 ** (-3)

f15dx = np.ones(len(f15x)) * 0.25

#f15uplims = [1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0]


# --- McLure+ 2013 ---
mcl13x = [-21, -20.5, -20, -19.5, -19, -18.5, -18, -17.5, -17]

mcl13y = [3e-5, 1.2e-4, 3.3e-4, 7.5e-4, 1.1e-3, 2.1e-3, 4.2e-3, 7.9e-3, 1.1e-2]

mcl13dy = [1e-5, 2e-5, 5e-5, 9e-5, 2e-4, 6e-4, 9e-4, 1.9e-3, 2.4e-3]


# --- Ono et al. 2018
o18x = [-24.165, -23.665, -23.165]

o18y = [0.00001e-4, 0.00010e-4, 0.00091e-4]

o18dx = [0.25, 0.25, 0.25]

o18dyU = [0.00019e-4, 0.00039e-4, 0.0008e-4]

o18dyL = [0.00001e-4, 0.00009e-4, 0.00044e-4]

# ------ FITTING FUNCTIONS ------

# --- Double power law ---
DPLx = np.arange(-26, -16, 0.1)

# Bowler+17
DPLy_b17 = dpl(2.3e-4, -2.19, -4.6, DPLx, -20.60)

sch_b17 = schechter(4.2e-4, -2.07, DPLx, -20.49)

# Harikane+22
#				phi			alpha	beta  M      M*				AGN component
DPLy_h22 = dpl(10**(-3.05), -1.89, -3.81, DPLx, -20.12) + dpl(10**(-8.49), -1.23, -2.73, DPLx, -24.9)

######
DPLval1 = dpl(10**(-3.05), -1.89, -3.81, -22.925, -20.12)
DPLval2 = dpl(10**(-3.05), -1.89, -3.81, -23.675, -20.12)
print(DPLval1)
print(DPLval2)
#exit()
######
sch_h22 = schechter(10**(-3.14), -1.88, DPLx, -20.49)

# Adams+22 AGN LF extrapolated to z=6.


DPLy_adams22 = dpl(10**(-8.65), -2.12, -6.31, DPLx, -27.58) * 10**(-2.3 / (7.0-6.0))

# 							phi*						  alpha   beta   M       M*
DPLy_adams22 = dpl(10**(-8.65) * 10**(0.82 * (4.8-6.74)), -2.12, -6.31, DPLx, -27.58)



# Matsuoka z=6 AGN
def matsuoka_LF(M, M_star, alpha, beta, phi_star):

    a = alpha + 1
    b = beta + 1

    denom1 = 10 ** (0.4 * a * (M-M_star))
    denom2 = 10 ** (0.4 * b * (M-M_star))

    denom = denom1 + denom2

    phi = phi_star / denom

    return phi

matsuoka18 = matsuoka_LF(M=DPLx, M_star=-24.9, alpha=-1.23, beta=-2.73, phi_star = 10.9 * 1e-9)

# Leethochawalit z=8
DPLy_leet22 = dpl(1.00e-5, -2.23, -4.16, DPLx, -22.33)

# Donnan z=8
DPLy_donnan22 = dpl(3.30e-4, -2.04, -4.26, DPLx, -20.02)

# Bouwens Schechter function
sch_b21 = schechter(0.19e-3, -2.06, DPLx, -21.15)

# ------------- REDSHIFT 8 DATA --------------

# Bowler+20
rb20x = [-22.90, -22.15, -21.65]
rb20y = [0.14e-6, 0.58e-6, 2.95e-6]

rb20dx = [0.5, 0.25, 0.25]
rb20dy = [0.06e-6, 0.33e-6, 0.98e-6]

# Bouwens+21
b21_x = [-21.85, -21.35, -20.85, -20.10, -19.35, -18.60, -17.60]
b21_y = [0.000003, 0.000012, 0.000041, 0.00012, 0.000657, 0.0011, 0.00302]

b21_dx = np.diff(b21_x)
b21_dx = np.insert(b21_dx, 0, 0.5)
b21_dy = [0.000002, 0.000004, 0.000011, 0.000040, 0.000233, 0.00034, 0.00114]

# Donnan+ 22
cd22x = [-22.17, -21.42]
cd22y = [0.63e-6, 3.92e-6]

cd22dx = [0.5, 0.25]
cd22dyL = [0.30e-6, 1.56e-6]
cd22dyU = [0.50e-6, 2.34e-6]

# Rojas-Ruiz+20
rr20x = [-23, -22, -21]
rr20y = [5.6904e-6, 7.0378e-6, 40.326e-6]

rr20dx = [0.5, 0.5, 0.5]

rr20dyL = [2.413e-6, 3.047e-6, 18.802e-6]
rr20dyU = [11.184e-6, 13.413e-6, 81.276e-6]

# Bagley+22

bag22x = np.arange(-23.15, -20.65, 0.1)
bag22x = np.flip(bag22x)

print(len(bag22x))

bag22yU = [76.55, 56.79, 43.52, 38.65, 34.63, 30.36, 24.73, 19.12, 14.83, 11.7, 9.46, 7.71, 6.32, 5.16, 4.37, 3.76, 3.27, 2.85, 2.49, 2.27, 2.16, 2.1, 1.95, 1.72, 1.49]
bag22yU = [val * 1e-6 for val in bag22yU]

bag22yL = [31.65, 23.7, 18.64, 16.75, 15.17, 13.16, 11., 8.37, 6.48, 5.07, 4.09, 3.33, 2.74, 2.26, 1.91, 1.63, 1.41, 1.25, 1.08, 1., 0.97, 0.96, 0.88, 0.77, 0.68]
bag22yL = [val * 1e-6 for val in bag22yL]

# ---------- REDSHIFT z>8 JWST DATA ----------------

# Harikane+23 z=9-17
yh23x = [-18.03, -19.03, -20.03, -21.03]
yh23y = [1.12e-3, 2.24e-4, 4.08e-5, 4e-5]

yh23dx = [0.5, 0.5, 0.5, 0.5]
yh23dyU = [1.03e-3, 1.87e-4, 9.6e-5, 9.42e-5]
yh23dyL = [0.9e-3, 1.46e-4, 3.92e-5, 3.85e-5]

yh23dy = [yh23dyL, yh23dyU]

# Perez-Gonzalez z=9-12
pg23x = [-16.5, -17.5, -18.5, -19.5]
pg23y = [10**(-2.3), 10**(-2.77), 10**(-3.53), 10**(-4.17)]

pg23dx = [0.5, 0.5, 0.5, 0.5]
pg23dyL = [10**(-2.3-0.27), 10**(-2.77-0.15), 10**(-3.53-0.20), 10**(-4.17+0.25)]
pg23dyU = [10**(-2.3+0.16), 10**(-2.77+0.11), 10**(-3.53+0.14), 10**(-4.17-0.61)]

pg23dy = [pg23dyL, pg23dyU]

# ------ PLOT ------
plt.clf()
plt.close()
#plt.figure(figsize=(20, 10)) # full
plt.figure(figsize=(10, 10)) # full


#plt.figure(figsize=(12, 10)) # indiv fields

# MY OLD DATA BEFORE COMPLETENESS CORRECTION
#plt.errorbar(xbin_plot, [2.07e-8, 1.87e-8, 2.12e-7, 1.98e-6], yerr=LFerror, xerr=binErr, marker='o', color='orange', label='This work', markersize=5, elinewidth=2, uplims=[1,0,0,0], zorder=13, linestyle='none')

# My data
#plt.errorbar(xbin_plot, LFsum, yerr=LFerror/2, xerr=bin_width/2, fmt='ro', label='This work', markersize=12) # Equal bin widths
#plt.errorbar(xbin_plot, LFsum, yerr=LFerror, xerr=binErr, fmt='ro', label='This work', markersize=15, elinewidth=3.5, uplims=[1,0,0,0,0], zorder=10) # Custom bins, 5
plt.errorbar(xbin_plot, LFsum, yerr=LFerror, xerr=binErr, fmt='ro', label='This work', markersize=13, elinewidth=3.5, uplims=[1,0,0,0], zorder=10) # Custom bins, 4
#plt.errorbar(xbin_plot, LFsum, yerr=LFerror, xerr=binErr, fmt='ro', label='This work', markersize=15, elinewidth=3.5, zorder=10) # Individual fields / inclusive cands

#plt.errorbar(xbin_plot, LFsum, yerr=LFerror/2, xerr=binErr, fmt='ro', label='This work', markersize=15, elinewidth=3.5, uplims=[1,0,0,0,0]) # Halving the y error?

# Contam factor
#plt.plot(xbin_plot[1:], contam[1:], marker='o', ms=13, color='red', alpha=1.0, ls='none', mew=3, fillstyle='none', zorder=10, label='This work, contamination factor') # Normal
#plt.plot(xbin_plot, contam, marker='o', ms=13, color='red', alpha=1.0, ls='none', mew=3, fillstyle='none', zorder=10, label='This work, contamination factor') # Inclusive

# Bowler et al. 2017
#plt.errorbar(rb17x, rb17y, xerr=rb17dx, yerr=rb17err, color='green', label='Bowler+17', marker='s', linestyle='none', alpha=0.8, markersize=10)

# Harikane et al. 2022
plt.errorbar(yh22x, yh22y, xerr=yh22dx, yerr=yh22err, label='Harikane+22', marker='v', linestyle='none', markersize=10, alpha=0.8, color='blue')

# CONTAMINATED H22
#plt.errorbar(yh22x, yh22y_contam, xerr=yh22dx, yerr=yh22err, label='Harikane+22, z=7 contam', marker='v', linestyle='none', markerfacecolor='none', markersize=12, alpha=0.8, color='blue')

# Bowuens et al. 2021
#plt.errorbar(b21x, b21y, xerr=b21dx, color='orange', yerr=b21dy, label='Bouwens+21', marker='o', markerfacecolor='none', markersize=10, alpha=0.8, linestyle='none')

# Finkelstein et al. 2015
#plt.errorbar(f15x, f15y, xerr=f15dx, color='magenta', yerr=[f15dyL, f15dyU], label='Finkelstein+15', markersize=10, alpha=0.8, marker='s', markerfacecolor='none', linestyle='none')

# McLure et al. 2013
#plt.errorbar(mcl13x, mcl13y, xerr=0.25, yerr=mcl13dy, color='brown', markersize=10, marker='p', markerfacecolor='none', linestyle='none', alpha=0.8, label='McLure+13')

# Ono et al. 2018
#plt.errorbar(o18x, o18y, xerr=0.25, yerr=[o18dyL, o18dyU], marker='D', color='darkturquoise', markersize=6, alpha=0.6, markerfacecolor='none', linestyle='none', label='Ono+18')

# Double power law from Bowler+17
#plt.plot(DPLx, DPLy_b17, color='green', label='DPL fit, Bowler+17', alpha=0.4, linewidth=5)
plt.plot(DPLx, DPLy_b17, color='green', label='Bowler+17, z=7 DPL fit', alpha=0.4, linewidth=5)
#plt.plot(DPLx, DPLy_b17, color='green', label='Bowler+17', alpha=0.4, linewidth=5)

# Schechter function from Bowler+17
#plt.plot(DPLx, sch_b17, color='green', label='Bowler+17, Schechter fit', alpha=0.4, linestyle='dotted', linewidth=5)

# Schechter function from Bouwens+21
#plt.plot(DPLx, sch_b21, color='orange', label='Schechter fit, Bouwens+21', alpha=0.4, linestyle='dashdot', linewidth=5)

# DPL + DPL from Harikane+22
#plt.plot(DPLx, DPLy_h22, color='blue', label='Harikane+22, z=7 DPL+DPL fit', alpha=0.4, linewidth=5)

# Galaxy component Schechter function from Harikane+22
#plt.plot(DPLx, sch_h22, color='blue', label='Schechter fit, Harikane+22', alpha=0.4, linestyle='dotted', linewidth=3)

# ---------- REDSHIFT 8 ------------

# Adams+ 22 AGN LF extrapolated to z=7
#plt.plot(DPLx, DPLy_adams22, color='purple', label='Adams+22, z=7 AGN LF', alpha=0.4, linestyle='dotted', linewidth=5)

# Matusoka+18 AGN LF at z=6
#plt.plot(DPLx, matsuoka18, color='red', label='Matsuoka+18', alpha=0.4, linestyle='dashed', linewidth=5)

# leethochawalit+22 z=8
#plt.plot(DPLx, DPLy_leet22, color='gray', label='Leethochawalit+22', alpha=0.4, linestyle='dashed', linewidth=5)

# Adams extrapolated point AGN
#plt.scatter(-24, 5e-10, color='purple', label='Adams+22 z=6 AGN estimation', alpha=0.8, marker='o', facecolor='purple', s=90)

# Bouwens+21 z=8
#plt.errorbar(b21_x, b21_y, xerr=b21_dx/2, color='orange', yerr=b21_dy, label='Bouwens+21, z=8', marker='o', markersize=10, alpha=0.8, linestyle='none')

# Bowler+20 z=8
#plt.errorbar(rb20x, rb20y, xerr=rb20dx, yerr=rb20dy, color='deepskyblue', label='Bowler+20, z=8', marker='s', linestyle='none', alpha=0.8, markersize=10)

# Donnan+22 z=8
#plt.errorbar(cd22x, cd22y, xerr=cd22dx, color='magenta', yerr=[cd22dyL, cd22dyU], label='Donnan+22, z=8', markersize=10, alpha=0.8, marker='H', linestyle='none')

# Rojas-Ruiz+20
#plt.errorbar(rr20x, rr20y, xerr=rr20dx, color='deepskyblue', yerr=[rr20dyL, rr20dyU], label='Rojas-Ruiz+20, $7.0<z<8.3$', markersize=10, alpha=0.8, marker='D', linestyle='none')

# Donnan+22 DPL
#plt.plot(DPLx, DPLy_donnan22, color='magenta', label='Donnan+22, z=8 DPL fit', alpha=0.4, linewidth=5, linestyle='dashed')

#plt.text(-20.43, 5.e-8, "z=7 AGN LF extrapolation", rotation=20, color='purple', fontsize=15)
#plt.text(-20.33, 5.306e-5, "z=8 UV LF, SuperBoRG", rotation=21, color='gray', fontsize=15)
#plt.text(-19.535, 0.00103, "z=7 UV LF", rotation=21, color='green', fontsize=15)

# Bagley+22
#plt.fill_between(bag22x, bag22yU, bag22yL, label='Bagley+22, $8.5 \leq z \leq 11$', alpha=0.5, color='orange')

# ----------- JWST RESULTS ----------
#plt.errorbar(yh23x, yh23y, xerr=yh23dx, yerr=yh23dy, label='Harikane+23, z=9', marker='s', linestyle='none', markersize=10, alpha=0.8, color='black')

#plt.errorbar(pg23x, pg23y, xerr=pg23dx, yerr=pg23dy, label='Pérez-González+23, z=9', marker='o', linestyle='none', markersize=10, alpha=0.4, color='black')



# ------ style ------
#ax = plt.gca()

plt.yscale('log')

locmin = tck.LogLocator(base=10.0, subs=np.arange(2, 10) * .1, numticks=100)
plt.axes().yaxis.set_minor_locator(locmin)
#plt.axes().yaxis.set_minor_formatter(tck.NullFormatter())

plt.axes().xaxis.set_minor_locator(tck.AutoMinorLocator(10))
#plt.axes().yaxis.set_minor_locator(tck.MultipleLocator(10))

plt.tick_params(which='major', length=10, width=3)
plt.tick_params(axis='both', which='minor', length=5, width=2)

#ax.yaxis.set_minor_formatter(FormatStrFormatter("%.1f"))

#plt.xlim([-25, -16.6])
plt.xlim([-25, -18])

plt.xlabel(r'$M_{\mathrm{UV}}$', fontsize=20)
plt.ylabel(r'$\mathrm{Number \ of \ objects \  / \ mag \ / \ Mpc^{3}}$', fontsize=20)

#plt.text(-24.7, 0.008, 'ECDF-S', fontsize=25)
#plt.text(-24.7, 0.008, 'XMM-LSS', fontsize=25)

plt.xticks(fontsize=20)
plt.yticks(fontsize=20)

#plt.vlines(-24, 1e-10, 1e-6)
#plt.vlines(-25, 1e-10, 1e-6)

#plt.ylim([10**(-16), 10**(-10)])
plt.ylim([10**(-10), 0.03])
#plt.legend(fontsize=16, loc='lower right')
plt.legend(fontsize=16, loc='upper left')
#plt.title(' Luminosity Function')

#plt.savefig('../plots/LF_{0}_only.png'.format(field))

#plt.savefig('../plots/LF_data.pdf', bbox_inches='tight')
#plt.savefig('../plots/LF_H22.pdf', bbox_inches='tight')
#plt.savefig('../plots/LF_z8.pdf', bbox_inches='tight')
#plt.savefig('../plots/LF_leet.pdf', bbox_inches='tight')

#plt.savefig('../plots/LF_AGN.pdf', bbox_inches='tight')

#plt.savefig('../plots/LF_inclusive.pdf', bbox_inches='tight')

#plt.savefig('../plots/LF_talk.pdf', bbox_inches='tight')

#plt.savefig('../plots/LF_ESA_proposal.pdf', bbox_inches='tight')


plt.show()



