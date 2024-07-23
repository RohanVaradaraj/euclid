#!/usr/bin/env python3

"""
plot_spec.py

Plot the .spec files produced by LePhare
/mnt/vardy/vardygroupshare/HSC_SSP_DR3/codes/./tmp_plot_seds.sh-513909.out
Created: Monday 25th April 2022

"""

'''SETUP'''

# Import libraries

# Use backend that won't display when running in the queue
import matplotlib as mpl
##mpl.use('Agg')

import os
from astropy.table import Table
from astropy.io import ascii, fits
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.ticker as mticker

mpl.rcParams['figure.dpi'] = 100
plt.rcParams['axes.linewidth'] = 2
plt.rcParams.update({'font.size': 16})

def mag_to_Jy(m):
    '''Convert mags and their errors to Jy'''
    flux = 10**(-0.4*(m+48.6))
    return flux

def flux_to_mag(flux):
    '''Convert flux to mag, for the secondary y axis'''
    mag = -2.5*np.log10(flux*1e-29)-48.6
    return mag

def magJy_error(flux, error):
    error_upper = flux * (10.0 ** (0.4 * error) - 48.6)
    error_lower = flux * (48.6 - 10.0 ** (-0.4 * error))
    error_flux = (error_lower + error_upper) / 2.0
    return error_flux

# Choose whether second model plotted is stars or not
irac =       True

# Different params required for stellar plots (when avoiding bluest bands)
star =       False

# Plot bagpipes on top?
pipes =      True

# Save individual plots as png
save_indiv = False

# Show plots?
display =    True

# Need slightly different plotting for lya, XMM
lya =        False

# Add specz for rebels
REBELS =     False

# Read parameters in from lephare input txt file
lephareParams = Table.read('lephare_params_DR3.txt', format='ascii.commented_header')
#lephareParams = Table.read('lephare_params.txt', format='ascii.commented_header')

#lephareParams = Table.read('lephare_params_sim.txt', format='ascii.commented_header') # SIMULATION

# Get field
field = lephareParams['field'][0]
catName = lephareParams['inputName'][0]

# Get directory
# ------------------- TAKE CARE HERE -------------------------------
#goodDir = str(lephareParams['dir_target'][0]) # Use this for Lyman alpha, when irac=True, self reference.
#goodDir = str(lephareParams['dir_source'][0])
#listDir = str(lephareParams['dir_target'][0])
#listDir = str(lephareParams['dir_source'][0]) # For plotting LePhare withour IRAC



if irac:
    goodDir = str(lephareParams['dir_source'][0])
    listDir = str(lephareParams['dir_source'][0])

if star and pipes:
    goodDir = str(lephareParams['dir_source'][0])
    listDir = str(lephareParams['dir_target'][0])

if lya and star == False:
    goodDir = str(lephareParams['dir_target'][0])
    listDir = str(lephareParams['dir_target'][0])
if lya and star == True:
    goodDir = str(lephareParams['dir_source'][0])
    listDir = str(lephareParams['dir_target'][0])

goodDir = str(lephareParams['dir_target'][0])
listDir = str(lephareParams['dir_target'][0])

# Title from this name
title = goodDir.split('/')[-2]
print(title)
# Read in the catalogue to plot upper limits
catDir = '/mnt/hoy/temporaryFilesROHAN/lephare/inputs/'
cat = 	Table.read(catDir+catName, format='ascii.commented_header')

# Directory to save individual plots to.

plotDir = '/mnt/zfsusers/varadaraj/sed_fitting/{0}/plots/finalPlots/'.format(field) # PRIMARY
#plotDir = '/mnt/zfsusers/varadaraj/sed_fitting/{0}/plots/insecure/'.format(field) # INSECURE

if lya:
    plotDir = '/mnt/zfsusers/varadaraj/sed_fitting/{0}/plots/lya/sed/'.format(field) # Lyman alpha

baseDir = '/mnxt/zfsusers/varadaraj/sed_fitting/{0}/'.format(field)

###### BAGPIPES ########

bagDir = '/mnt/vardy/vardygroupshare/HSC_SSP_DR3/codes/pipes/sed/{0}/'.format(field) # PRIMARY
#bagDir = '/mnt/vardy/vardygroupshare/HSC_SSP_DR3/codes/pipes/sed/{0}/insecure/'.format(field) # INSECURE
#bagDir = '/mnt/vardy/vardygroupshare/HSC_SSP_DR3/codes/pipes/sed/{0}/low_Av/'.format(field) # Av < 0.5

if lya:
    bagDir = '/mnt/vardy/vardygroupshare/HSC_SSP_DR3/codes/pipes/sed/{0}/lya/'.format(field) # LYMAN ALPHA
    bagDir = '/mnt/vardy/vardygroupshare/HSC_SSP_DR3/codes/pipes/sed/{0}/low_Av_lya/'.format(field) # Av < 0.5


# Read in the parent catalogue for a few things

if field == 'XMM':
    parDir = '/mnt/vardy/vardygroupshare/data/catalogues/{0}FULL/cutting/'.format(field)
    #parDir = '/mnt/vardy/vardygroupshare/HSC_SSP_DR3/simulations/mockCats/{0}/'.format(field) # SIMULATION

if field == 'CDFS':
    parDir = '/mnt/vardy/vardygroupshare/data/catalogues/{0}FULL/cutting/'.format(field)
    #parDir = '/mnt/vardy/vardygroupshare/HSC_SSP_DR3/simulations/mockCats/{0}/'.format(field) # SIMULATION
par = str(lephareParams['catalogue'][0])

parentCat = fits.open(parDir + par)
parentCat = parentCat[1].data

# Define start and end points of each section in .spec file. Depends on field because of no. of filters.
if field == 'XMM' or field == 'REBELS':
    ds_phot = 9
    de_phot = 20 #20 #18 if no GR
    if star:
        de_phot = 18

    ds_mod = 201 #201 # 199 if no GR
    if star:
        ds_mod = 199
    de_mod = -1

    ds_pz = 20 #20 #18 if no GR
    de_pz = 201 #201 # 199 if no GR
    if star:
        ds_pz = 18
        de_pz = 199

    ds_param = 3
    de_param = 9

if field == 'CDFS': # No ugG(rR): 21 -> 18(16), 202 -> 199(197)
    ds_phot = 9
    if star == False:
        de_phot = 21
    if star == True:
        de_phot = 16

    if star == False:
        ds_mod = 202
    if star == True:
        ds_mod = 197
    de_mod = -1

    if star == False:
        ds_pz = 21
    if star == True:
        ds_pz = 16
    if star == False:
        de_pz = 202
    if star == True:
        de_pz = 197

    ds_param = 3
    de_param = 9

'''
# The above need to change if we include IRAC
if irac:

    if field == 'XMM' or field == 'REBELS':
        ds_phot = 9
        de_phot = 22

        ds_mod = 203
        de_mod = -1

        ds_pz = 22
        de_pz = 203

        ds_param = 3
        de_param = 9

    if field == 'CDFS':
        ds_phot = 9
        de_phot = 23

        ds_mod = 204
        de_mod = -1

        ds_pz = 23
        de_pz = 204

        ds_param = 3
        de_param = 9
'''
# From these values, we can get the number of filters which will be useful for extracting fluxes from the .in file.
n_bands = de_phot - ds_phot

# Double the value to correspond to each flux and error
n_in = 2 * n_bands

# Restrict input file to these columns
#cat = cat.columns[0:n_in]

####### CONFUSED OBJECTS ########
conf = ['1510053', '219612', '1548625', '165441', '661623', '1482093', '1017553', '125919', '1090797', '1025296']
conf = ['1510053', '219612', '1548625', '165441', '661623', '1017553', '125919', '1090797', '1025296'] # For main paper plot.

# Begin pdf file:
#with PdfPages(baseDir + 'plots/' + title + '.pdf') as pdf:
#with PdfPages('/mnt/hoy/temporaryFilesROHAN/lephare/results/tmp.pdf') as pdf: # For testing
#with PdfPages('/mnt/hoy/temporaryFilesROHAN/lephare/results/8XMM_smallPhot_0.0.pdf') as pdf:
#with PdfPages('/mnt/hoy/temporaryFilesROHAN/lephare/results/8XMM_dwarf_smallPhotSE.pdf') as pdf:
with PdfPages('tmp.pdf') as pdf: # For testing

    print('Saving to ', baseDir + 'plots/' + title + '.pdf')

    # Spectrum directory
    specDir = goodDir

    # Initialise an array for the photo-z, chi2 and other things.
    z_phots = []
    chis = []
    mags = []

    for f in sorted(os.listdir(listDir)):

        print('#######################################')
        print(f)
        print('#######################################')

        # Column names
        namesPhot=['phot', 'yerr', 'wlen', 'xerr', 'modelPhot', 'col6', 'col7']
        namesModel = ['wlen', 'flux']
        namesPDF = ['z', 'P']
        namesParam = ['Type', 'Nline', 'Model', 'Library', 'Nband', 'Zphot', 'Zinf', 'Zsup', 'Chi2', 'PDF', 'Extlaw', 'EB-V', 'Lir', 'Age', 'Mass', 'SFR', 'SSFR']

        # Load the .spec files
        # These data_start/end values were found just by looking at the .spec files.
        phot = ascii.read(specDir+f, format='basic', data_start=ds_phot, data_end=de_phot, delimiter=' ', names=namesPhot)

        model = ascii.read(specDir+f, format='basic', data_start=ds_mod, data_end=de_mod, delimiter=' ', names=namesModel)

        prob_z = ascii.read(specDir+f, format='basic', data_start=ds_pz, data_end=de_pz, delimiter=' ', names=namesPDF)

        params = ascii.read(specDir+f, format='basic', data_start=ds_param, data_end=de_param, delimiter=' ', names=namesParam)

        # If plotting REBELS catalogue, find index of each object to print z_phot and z_spec.
        if REBELS:
            rebels = Table.read('/mnt/vardy/vardygroupshare/data/catalogues/XMMFULL/REBELS.fits')
            obj =  f[2:-5].lstrip('0')
            idx = np.where(rebels['ID_2'] == float(obj))

        # For everything else, find the ID from the input file.
        if field == 'XMM' or field == 'CDFS':
            obj = f[2:-5].lstrip('0')
            idx = np.where(cat['ID'] == int(obj))

        print(obj)

        # Initialize flux and error arrays, to add .in data.
        flux = []
        error = []

        # Match ID with the input file to find the object.
        cat_obj = cat[idx]

        # Save flux and error
        for i in range(0, n_bands):
            flux = np.append(flux, cat_obj.columns[2*i+1])
            error = np.append(error, cat_obj.columns[2*i+2])

        print(cat_obj)

        print('INPUT DATA')
        print(input)

        print('OUTPUT')
        print(phot)


        ############ BAGPIPES ############

        if pipes:
            # Load in the model SED
            bagSED = Table.read(bagDir + '{0}_microJy.fits'.format(obj))

            # Convert from uJy to cgs
            bagSED['flux'] = bagSED['flux'] * 1e-29 #(1/3.631) * 10 ** -28.44
#            if obj == '308281':
#                bagSED['flux2'] = bagSED['flux2'] * 1e-29
#                bagSED['flux3'] = bagSED['flux3'] * 1e-29

            # Load model photometry
            bagPhot = Table.read(bagDir + '{0}_microJy_synthpoints.fits'.format(obj))
            bagPhot['flux'] = bagPhot['flux'] * 1e-29

            # Load P(z)
            bagPz = Table.read(bagDir + '{0}_Pz.fits'.format(obj))
            redshifts = np.percentile(bagPz['z_phot'], (16, 50, 84), axis=0).T

            zBest = round(redshifts[1], 2)

        # Also grab extra bands not used in BD fitting, since we need them for bagpipes
        par_idx = np.where(parentCat['ID'] == int(obj))[0]
        print('##############', par_idx)

        par_obj = parentCat[par_idx]

        print(par_obj)

        # Get name of object
        #name = par_obj['name'][0]
        name = par_obj['ID'][0]
        print(name)

        # Define the extra bands
        if field == 'XMM':
            extraFilt = ['HSC-G', 'HSC-R_DR3', 'HSC-I', 'HSC-NB0816_DR3', 'HSC-Z_DR3', 'HSC-NB0921_DR3', 'HSC-Y', 'Y', 'J', 'H', 'Ks', 'ch1servs', 'ch2servs']
            #extraFilt = [string.replace('-', '_').split('_DR3')[0] for string in extraFilt] # SIMULATION
            extraX = [ 0.481612e4, 0.623411e4,  0.774058e4, 0.817688e4,  0.891151e4, 0.921420e4, 0.977993e4, 0.102141e5,  0.125441e5, 0.164650e5, 0.214842e5,  0.355726e5,  0.450487e5]
            extraXerr = np.array([0.138600e4,  0.150410e4, 0.155180e4,  0.113000e3, 0.773000e3, 0.134000e3, 0.783000e3, 0.926000e3, 0.172500e4, 0.291600e4, 0.309200e4, 0.743200e4, 0.100970e5])

        if field == 'CDFS':
            extraFilt = ['u', 'g', 'r', 'i', 'HSC-G', 'HSC-R', 'HSC-I', 'HSC-Z', 'Y', 'J', 'H', 'Ks', 'ch1cds', 'ch2cds']
            #extraFilt = [string.replace('-', '_') for string in extraFilt] # SIMULATION
            extraXerr = np.array([ 0.558720e3, 0.132981e4, 0.136217e4, 0.156046e4, 0.138600e4,  0.150410e4, 0.155180e4, 0.773000e3,  0.926000e3, 0.172500e4, 0.291600e4, 0.309200e4, 0.743200e4, 0.100970e5])
            extraX = [0.357369e4, 0.472813e4, 0.629929e4, 0.756951e4, 0.481612e4, 0.623411e4, 0.774058e4, 0.891151e4, 0.102141e5, 0.125441e5, 0.164650e5, 0.214842e5, 0.355726e5, 0.450487e5]

            # Define the extra bands
#            if field == 'XMM':
#                extraFilt = ['HSC-G_DR3', 'HSC-R_DR3', 'ch1servs', 'ch2servs']
#                extraXerr = np.array([0.138600e4, 0.150410e4, 0.743200e4,  0.100970e5])
#                extraX = [0.481612e4, 0.623411e4, 0.355726e5, 0.450487e5]

#            if field == 'CDFS':
#                extraFilt = ['u', 'g', 'r', 'HSC-G', 'HSC-R', 'ch1cds', 'ch2cds']
#                extraXerr = np.array([0.558720e3, 0.132981e4, 0.136217e4,  0.138600e4, 0.150410e4, 0.743200e4, 0.100970e05])
#                extraX = [0.357369e4, 0.472813e4, 0.629929e4, 0.481612e4, 0.623411e4, 0.355726e5, 0.450487e5]


            # Empty flux array to add fluxes to
        extraFlux = []
        extraYerr = []

        # Extract flux in each filter.
        for i, filt in enumerate(extraFilt):
            extraFlux = np.append(extraFlux, par_obj['flux_{0}'.format(filt)])
            extraYerr = np.append(extraYerr, par_obj['err_{0}'.format(filt)])

        if field == 'XMM' and (par_obj['flux_ch1'] != -99. or par_obj['flux_ch2'] != -99):
            extraFlux[-2] = par_obj['flux_ch1']
            extraFlux[-1] = par_obj['flux_ch2']

            extraYerr[-2] = par_obj['err_ch1']
            extraYerr[-1] = par_obj['err_ch2']

            print('Using UDS IRAC')
            # Get significance
#            extraSigma = extraFlux / extraYerr
#            extraLims = [int(i < 2.) for i in extraSigma]

            # Replace <2sigma values with 2sigma upper limit
#            for i, val in enumerate(extraFlux):
#                if np.abs(extraFlux[i] / extraYerr[i]) < 2. and extraFilt[i][0:2] != 'ch':
#                    extraFlux[i] = np.abs(extraFlux[i]) + np.abs(2*extraYerr[i])
#                    extraYerr[i] = extraFlux[i]*0.1
#                if np.abs(extraFlux[i] / extraYerr[i]) < 2. and extraFilt[i][0:2] == 'ch':
#                    extraFlux[i] = np.abs(2*extraYerr[i])
#                    extraYerr[i] = extraFlux[i]*0.1



        #################################

        # Extract each model by finding where wavelength jumps.
        model_wlens = np.split(model['wlen'], np.where(np.diff(model['wlen']) < 0)[0] + 1)
        model_fluxes = np.split(model['flux'], np.where(np.diff(model['wlen']) < 0)[0] + 1)

        # Convert magnitudes to fluxes. Also deal with the '******' for smaller photometry.
        if phot['yerr'][0] == '**********':
            phot['yerr'][0] = 4.0
        print(phot)

        phot_Jy = mag_to_Jy(phot['phot'])
        model_Jy = mag_to_Jy(phot['modelPhot'])

        best_Jy = mag_to_Jy(model_fluxes[0])
        if len(model_fluxes) > 1:
            second_Jy = mag_to_Jy(model_fluxes[1])
        # Check if the third model exists. Do the same for plotting later on.
        if len(model_fluxes) > 2:
            stellar_Jy = mag_to_Jy(model_fluxes[2])

        xerr = phot['xerr']
        yerr = magJy_error(phot_Jy, phot['yerr'])

        print(yerr)


        ######## DOING IT PROPERLY: ALL CATALOGUE FLUXES ################
        phot_Jy = np.array(extraFlux)
        yerr = np.array(extraYerr)
        x = np.array(extraX)
        xerr = np.array(extraXerr)

        # Create separate arrays for measurements with sigma < 2 to plot as an upper limit.

        sigma = phot_Jy/yerr

        print('SIGMA ARRAY')
        print(sigma)

        uplims = [int(i < 2.) for i in sigma]

        print('UPPER LIMIT BOOLEAN ARRAY')
        print(uplims)

        # Replace each <2sigma value with two sigma upper limit.
        for i, val in enumerate(phot_Jy):
            if np.abs(phot_Jy[i] / yerr[i]) < 2.:
                phot_Jy[i] = np.abs(phot_Jy[i]) + 2*np.abs(yerr[i])
                yerr[i] = phot_Jy[i]*0.25

        if lya and obj == '156132':
            phot_Jy[-3] = np.abs(phot_Jy[-3]) + 2*np.abs(yerr[-3])
            yerr[i] = phot_Jy[-3]*0.1
            uplims[-3] = 1


        #### Make sure IRAC 20% error in Lyman alpha ####
        if np.abs(yerr[-2]/phot_Jy[-2]) < 0.2:
            yerr[-2] = 0.2 * phot_Jy[-2]
        if np.abs(yerr[-1]/phot_Jy[-1]) < 0.2:
            yerr[-1] = 0.2 * phot_Jy[-1]


        ######### CONFUSED OBJECTS ########
        if obj in conf:
            irac_y = phot_Jy[-2:]
            irac_dy = yerr[-2:]
            irac_x = x[-2:]
            irac_dx = xerr[-2:]
            irac_lims = uplims[-2:]

        #### Plotting in Jy and microns #### There are factors of 10000 everywhere for wlen. Defining Jansky factor here.
        uJy = 1e-29


        ### CREATE FIGURE ###
        fig, ax1 = plt.subplots(figsize=(9,6))

        # Define size and location of inset P(z)
        left, bottom, width, height = [0.65, 0.23, 0.2, 0.2]
        if (field == 'CDFS' and obj == '1482093') or (field == 'CDFS' and obj == '243379') or (lya == True and obj == '156132' and pipes == True) or (field == 'CDFS' and obj == '705214'): #  or (pipes == False and obj == '16448')
            left, bottom, width, height = [0.65, 0.65, 0.2, 0.2]
        ax2 = fig.add_axes([left, bottom, width, height])

        # Plot model photometry
        #if irac:
            #ax1.scatter(phot['wlen']/10000, model_Jy/uJy, facecolors='none', edgecolors='gray', alpha=1.0, s=80, zorder=5) #, label='Model photometry')

        # Plot best model
        if irac:
            ax1.plot(model_wlens[0]/10000, best_Jy/uJy, linestyle='dashed', color='blue', alpha=1.0, label='z = {0}, {1} = {2}'.format(round(params['Type' == 'GAL-1']['Zphot'], 2), r'${\chi}^{2}$', round(params['Type' == 'GAL-1']['Chi2'], 1)), linewidth=3)

        # Save best model redshift, chi2 etc
        z_phots = np.append(z_phots, params['Type' == 'GAL-1']['Zphot'])
        chis = np.append(chis, params['Type' == 'GAL-1']['Chi2'])

        # Plot second model
        if irac:
            if len(model_fluxes)>2:
                ax1.plot(model_wlens[1]/10000, second_Jy/uJy, color='red', alpha=1.0, label='z = {0}, {1} = {2}'.format(round(params[1]['Zphot'], 2), r'${\chi}^{2}_{\mathrm{gal2}}$', round(params[1]['Chi2'], 1)), linewidth=3)

        if lya==True and star==False:
            ax1.plot(model_wlens[1]/10000, second_Jy/uJy, color='red', alpha=1.0, label='z = {0}, {1} = {2}'.format(round(params[1]['Zphot'], 2), r'${\chi}^{2}_{\mathrm{gal2}}$', round(params[1]['Chi2'], 1)), linewidth=3)


        #if star:
        #    if len(model_fluxes)==2:
        ax1.plot(model_wlens[1]/10000, second_Jy/uJy, color='red', alpha=1.0, label='{0} = {1}'.format( r'${\chi}^{2}_{\mathrm{BD}}$', round(params[5]['Chi2'], 1)), linewidth=2)

        # Plot stellar model
        if star:
            if len(model_fluxes)>2:
                ax1.plot(model_wlens[2]/10000, stellar_Jy/uJy, color='red', alpha=1.0, label='{0} = {1}'.format( r'${\chi}^{2}_{\mathrm{BD}}$', round(params[5]['Chi2'], 1)), linewidth=2)

        # ... or the quasar model if we're running that
#        if len(model_fluxes) > 2:
#            ax1.plot(model_wlens[2], stellar_Jy, color='purple', alpha=0.6, label='z = {0}, {1} = {2}'.format(params[4]['Zphot'], r'${\chi}^{2}$', params[4]['Chi2']))

        # This is the stellar model for goldrush, but there is no GAL-2 so need to spoof it.
        #ax1.plot(model_wlens[1]/10000, second_Jy/uJy, color='red', alpha=1.0, label='{0} = {1}'.format( r'${\chi}^{2}_{\mathrm{BD}}$', params[5]['Chi2'].round(1)), linewidth=2)

        ### PHOTOMETRY ###
#        ax1.errorbar(phot['wlen'], phot_Jy, xerr=xerr/2, yerr=yerr, ls='None', color='red', capsize=1, \

        if not obj in conf:
            ax1.errorbar(x/10000, phot_Jy/uJy, xerr=xerr/20000, yerr=yerr/uJy, ls='None', color='black', capsize=2.5, \
            elinewidth=2, markeredgewidth=1, alpha=1, uplims=uplims, zorder=10)

        if obj in conf:
            ax1.errorbar(x[:-2]/10000, phot_Jy[:-2]/uJy, xerr=xerr[:-2]/20000, yerr=yerr[:-2]/uJy, ls='None', color='black', capsize=2.5, \
            elinewidth=2, markeredgewidth=1, alpha=1, uplims=uplims[:-2], zorder=10)
            ax1.errorbar(irac_x/10000, irac_y/uJy, xerr=irac_dx/20000, yerr=irac_dy/uJy, ls='None', color='g', capsize=2.5, \
            elinewidth=2, markeredgewidth=1, alpha=1, uplims=irac_lims, zorder=10)[-1][0].set_linestyle('--')

        ### WIRDS ##
        ax1.errorbar(1.253, 1.83574e-30/uJy, xerr=0.158, yerr=1.83574e-30/uJy*0.1, ls='None', color='green', capsize=2.5, \
        elinewidth=2, markeredgewidth=1, alpha=1, zorder=10, marker='s', markersize=8)

        ########## BAGPIPES ############
        if pipes:
            # SED
            ax1.plot(bagSED['wave']/10000, bagSED['flux']/uJy, color='deepskyblue', alpha=1.0, label='{0} = {1}'.format(r'$z_{\mathrm{BP}}$', zBest), linewidth=2)
#            if obj == '308281':
#                ax1.plot(bagSED['wave'], bagSED['flux2'], color='darkred', alpha=0.7, label='z = {0}, 2'.format(zBest2))
#                ax1.plot(bagSED['wave'], bagSED['flux3'], color='brown', alpha=0.7, label='z = {0}, 3'.format(zBest3))

            # Model photometry
            ax1.scatter(bagPhot['wave']/10000, bagPhot['flux']/uJy, facecolors='none', edgecolors='gray', alpha=1.0, s=60, zorder=5) #, label='Model photometry')

            # Add photometry missing from BD fitting
#            ax1.errorbar(extraX, extraFlux, xerr=extraXerr/2, yerr=extraYerr, ls='None', color='red', capsize=1, \
#            elinewidth=1.3, markeredgewidth=1, alpha=1, uplims=extraLims)

        ###############################


        ##### Adding REBELS information #######
        if REBELS:
            # Fudge information into legend
            ax1.plot([], [], ' ', label='{0}'.format(rebels['REBELS_id'][idx][0]))
            ax1.plot([], [], ' ', label='$z_{phot}$ = ' + str(rebels['z_phot'][idx][0]))
            if str(rebels['z_spec'][idx][0]) == 'nan':
                ax1.plot([], [], ' ', label='No $z_{spec}$')
            else:
                ax1.plot([], [], ' ', label='$z_{spec}$ = ' + str(round(rebels['z_spec'][idx][0], 2)))

        ###################################




        ###### Plotting limits ########

#        ticks = plt.xticks() /10000
#        plt.set_xticklabels(ticks)

        if pipes == False:
#            ax1.set_xlim([3000, 26000])
            ax1.set_xlim([0.3, 2.6])
        if pipes:
            ax1.set_xlim([0.3, 5.5])
#            ax1.set_xlim([3000, 55000])

        # P(z) plotting limits
        #ax2.set_xlim([0, 8.1])

        if field == 'XMM':
            if irac == True:
#                ax1.set_ylim([1e-31/uJy, 2.577])
                ax1.set_ylim([1e-31/uJy, 10]) # GOLDRUSH
            if star == True:
                ax1.set_ylim([1e-31/uJy, 1e-28/uJy])
        if field == 'CDFS':
            if irac == True:
                ax1.set_ylim([4e-31/uJy, 9.4])
            if star == True:
                ax1.set_ylim([4e-31/uJy, 9.4])                
                #ax1.set_ylim([4e-31/uJy, 17.6])
#        ax1.set_ylim([10**(-31), np.max(bagSED['flux']) + 0.01*np.max(bagSED['flux'])])

        ##############################




        # ############## Other stuff (labels) ##################
        ax1.set_xlabel('{0} [{1}]'.format(r'$\lambda$', r'$\mu$m'), labelpad=0)
#        ax1.set_ylabel('Flux [erg $\mathrm{s}^{-1} \mathrm{cm}^{-2} \mathrm{Hz}^{-1}$]') # ergs
        ax1.set_ylabel('Flux [$\mu$Jy]', labelpad=0) # microjansky
        ax1.legend(loc='upper left', fontsize=14)

        # Move legend for goldrush plots
        #ax1.legend(loc='lower right', fontsize=15)

        # Set title as the object ID
#        ax1.set_title('ID ' + f[2:-5].lstrip('0'))
        ax1.set_title(name)

        # Log the y axis of course!
        ax1.set_yscale('log')

        # Add magnitude axis
        secax = ax1.secondary_yaxis('right', functions=(flux_to_mag, mag_to_Jy))
        secax.set_ylabel(r'$\mathrm{m_{AB}}$')
        secax.set_yticks([28, 27, 26, 25, 24, 23, 22])

        if obj == '1610530':
            ax1.set_ylim([1e-31/uJy, 1e-27/uJy])
            secax.set_yticks([28, 27, 26, 25, 24, 23, 22, 21, 20])


        if field == 'CDFS':
            secax.set_yticks([27, 26, 25, 24, 23, 22, 21])
        secax.yaxis.set_major_formatter(mticker.ScalarFormatter())

        ax1.tick_params(which='major', length=5, width=3)
        ax1.tick_params(axis='both', which='minor', length=5, width=2)

        secax.tick_params(which='major', length=5, width=3)
        secax.tick_params(axis='both', which='minor', length=5, width=2)
        ######################################################


        ############### Plot redshift PDF, inset onto SED plot ##################

        # LePhare P(z)
        if pipes == False:
            ax2.plot(prob_z['z'], prob_z['P'], color='green', linewidth=2, linestyle='dashed')


        # BAGPIPES P(z)

        zmin = np.min(prob_z['z'])
        zmax = np.max(prob_z['z'])

        if pipes:
#            ax2.tick_params(left = False, labelleft = False)
            bagbins = np.arange(zmin, zmax, 0.1)
            n,bins = np.histogram(bagPz['z_phot'], bins=bagbins, density=True)
            ax2.plot(bins[:-1], n/np.max(n), color='black', linewidth=2)
            ax2.set_xlim([zmin, zmax])
            ax2.set_ylabel('P(z)')
#            ax2.vlines(x=redshifts, ymin=0, ymax=np.max(bagPz['z_phot']+1), linestyle='dashed', color='red')


        # Labels
        ax2.set_xlabel('z', fontsize=14, labelpad=0)
        if pipes == False:
            ax2.set_ylabel('P(z)', fontsize=14, labelpad=0)


        # Keep z ticks fine
        #if irac == False:
        #ax2.set_xticks([0, 1, 2, 3, 4, 5, 6, 7, 8])

        # other stuff
        #ax2.set_yscale('log')
        #ax2.grid()

        #######################################################################

        if display:
            plt.show()
        pdf.savefig(fig)

        ######## SAVE FIGURE ########
        fig.savefig('../plots/1610530_SED.pdf', bbox_inches='tight')
        #exit()

        if save_indiv:
            if star:
#                plt.savefig(plotDir+'DWARF_{0}.pdf'.format(obj), bbox_inches='tight')
                plt.savefig(plotDir+'DWARF_{0}.pdf'.format(name), bbox_inches='tight')
                print('Saved to ' + plotDir+'DWARF_{0}.pdf'.format(name))
            if star == False:
#                plt.savefig(plotDir+'{0}.pdf'.format(obj), bbox_inches='tight')
                plt.savefig(plotDir+'{0}.pdf'.format(name), bbox_inches='tight')
                print('Saved to ' + plotDir+'{0}.pdf'.format(name))
#                plt.savefig(plotDir+'GOLDRUSH_{0}.pdf'.format(obj), bbox_inches='tight')
        plt.close()

        # Get the magnitude in YJ
#        mag = -2.5*np.log10(cat_obj['flux_YJ'])-48.6
#        mags = np.append(mags, mag)

exit()


print(z_phots)
plt.hist(z_phots, bins=np.arange(6, 8, 0.1))
#plt.savefig('z_phot_4.pdf')
plt.xlabel(r'$z_{phot}$')
plt.ylabel('Count')
plt.show()
plt.clf()

print(chis)

plt.hist(chis)
plt.xlabel(r'${\chi}^{2}$')
plt.ylabel('Count')
plt.show()
plt.clf()

print(mags)
plt.hist(mags)
plt.xlabel(r'$m_{YJ}$')
plt.ylabel('Count')
plt.show()

# Extra code not used atm.

# 1) READING IN ORIGINAL LEPHARE INPUT FILE
# 2) PLOTTING IN TERMS OF MAGNITUDE, WITH INTERACTIVE SED WINDOW

'''
# Input photometry directory
inputDir = '/mnt/vardy/vardygroupshare/HSC_SSP_DR3/lephare/inputs/'
input = 'REBELS.in'

# Read input filters
imDir = '/mnt/vardy/vardygroupshare/data/XMM1/'

inputs = Table.read(imDir + 'images.lis', format='ascii.commented_header')

availFilters = np.array(inputs['Name'])

# Load input files with their peculiar formats.
filters = ['HSC-G', 'HSC-R_DR3', 'HSC-I', 'HSC-Z_DR3', 'HSC-NB0816_DR3', 'HSC-Y', 'HSC-NB0921_DR3', 'Y', 'J', 'H', 'Ks']

namesFilt = ['ID']
plotFilt = []
errFilt = []

for i, filterName in enumerate(filters):

    namesFilt = namesFilt + ['flux_{0}'.format(filterName)]

    namesFilt = namesFilt + ['err_{0}'.format(filterName)]

    plotFilt = plotFilt + ['flux_{0}'.format(filterName)]

    errFilt = errFilt + ['err_{0}'.format(filterName)]

# Find unused filters
remainder = list(set(availFilters) - set(filters))

# Append unused filters to the end
for i, filterName in enumerate(remainder):

    namesFilt = namesFilt + ['flux_{0}'.format(filterName)]

    namesFilt = namesFilt + ['err_{0}'.format(filterName)]

namesFilt = namesFilt + ['context', 'spec_z']
'''


'''
# Plotting in Mag

# Plot photometry
plt.errorbar(phot['wlen'], phot['phot'], xerr=phot['xerr']/2, yerr=phot['yerr'], ls='None', color='red', capsize=2, \
elinewidth=1.3, markeredgewidth=1, alpha=1)

# Plot model photometry
plt.scatter(phot['wlen'], phot['modelPhot'], facecolors='none', edgecolors='blue', label='Model photmetry')

# Plot best model
plt.plot(model_wlens[0], model_fluxes[0], color='black', alpha=0.6, label='Best fit')

# Plot second model
plt.plot(model_wlens[1], model_fluxes[1], color='gray', alpha=0.6, label='Second fit', lw=0.8)

# Remove zeros from the stellar model
idx = np.where(model_fluxes[2] < 99.)

star_plot = model_fluxes[2][idx]
star_wlen = model_wlens[2][idx]

# Plot stellar model
plt.plot(star_wlen, star_plot, color='green', alpha=0.3, label='Stellar fit', lw=0.5)

# Remove 99s for plotting limits
idx = np.where(phot['phot'] != 99.)
xlims = phot['wlen'][idx]
ylims = phot['phot'][idx]

# Plotting limits
plt.xlim([min(xlims)-2000, max(xlims)+3000])
plt.ylim([min(ylims)-1.5, max(ylims)+4])

# Other stuff
plt.xlabel('Wavelength')
plt.ylabel('Magnitude')
plt.legend()

plt.gca().invert_yaxis()
plt.xscale('log')
#plt.show()
plt.clf()
'''
