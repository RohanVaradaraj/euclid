#!/usr/bin/env python3

"""
bagpipes_fit.py

Fit SEDs using BAGPIPES so that we can include nebular emission.

This will also give us the properties of galaxies.

Created: Tuesday 15th November 2022, modified for Euclid Tuesday 17th December

"""

import numpy as np
import bagpipes as pipes
import os
from astropy.io import fits, ascii
from astropy.table import Table
import matplotlib.pyplot as plt
import matplotlib as mpl
from pathlib import Path

# Use backend that doesn't display plots to user, to prevent issues when running on the queue.
#mpl.use('Agg')

# Load multinest
#os.system('module load multinest')

indiv = False
table = False


catDir = Path.cwd().parents[1] / 'data' / 'catalogues' #/ 'candidates'
outdir = Path.cwd() / 'pipes'
tableDir = Path.cwd() / 'pipes' / 'tables'
if not os.path.exists(outdir):
    os.makedirs(outdir)
if not os.path.exists(tableDir):
    os.makedirs(tableDir)

#catName = 'LAE.fits'
catName = 'Euclid_UltraVISTA_z7_sample_kron_piecewise.fits'
catName = 'COSMOS_5sig_Y_J_nonDet_HSC_G_nonDet_HSC_R_nonDet_HSC_I_with_irac_handbook.fits'

filtFile = 'filt_list.txt'

filters = ['HSC-G_DR3', 'HSC-R_DR3', 'HSC-I_DR3', 'HSC-NB0816_DR3', 'HSC-Z_DR3', 'HSC-NB0921_DR3', 'HSC-Y_DR3', 'Y', 'J', 'H', 'Ks', 'VIS', 'Ye', 'Je', 'He', 'f115w', 'f150w', 'f277w', 'f444w', 'ch1cds', 'ch2cds']

def load_obj(ID):

    """ Load objects from the xmm or cdfs catalogue in the format bagpipes needs."""

    # Open catalogue
    cat = Table.read(catDir / catName)

    # Convert ID column values to strings
    cat['ID'] = cat['ID'].astype(int)

    # Extract object with ID
    #cat['ID'].format

    # BAGPIPES ID is just the row number
    #ID = 178396
    #cat_here = cat[cat['ID'] == ID]
    cat = cat[cat['ID'] == int(ID)]
    print(cat)

    # # But also print our ID for this object.
    # print('########## OBJECT ID: ' + str(cat['ID']) + '##############')

    # Empty arrays for getting names of columns
    fluxes = []
    fluxerr = []

    # Loop through filters
    for i, filterName in enumerate(filters):

        # Append flux name
        fluxes = fluxes + ['flux_{0}'.format(filterName)]

        # Append flux error name
        fluxerr = fluxerr + ['err_{0}'.format(filterName)]

    # Extract fluxes and errors
    # Now get the flux
    flux = cat[fluxes][0][:]
    flux = np.array(flux)
    error = cat[fluxerr][0][:]
    error = np.array(error)



    # print the ch1 and ch2 magnitudes
    # mag_ch1 = -2.5*np.log10(cat['flux_ch1cds'][0]) - 48.6
    # mag_ch2 = -2.5*np.log10(cat['flux_ch2cds'][0]) - 48.6
    # print('ch1 mag: ', mag_ch1)
    # print('ch2 mag: ', mag_ch2)

    # Convert to microJanskys
    flux = flux * 1e29 #* 3.631 * 10**28.44
    error = error * 1e29 #* 3.631 * 10**28.44

    # If flux or error is nan, blow up the error
    # for i, flux in enumerate(fluxes):
    #     if np.isnan(flux) or np.isnan(errors[i]):
    #         fluxes[i] = 0
    #         errors[i] = 9.9e99

    # blow up the errors associated with any missing fluxes.
    for i, f in enumerate(flux):
        if (f == -99.) or (error[i] == -99.) or (error[i] < 0):
            flux[i] = 0
            error[i] = 9.9e99

    # Replace fluxes for ch1 and ch2 with tractor values from COSMOS2020
    # flux[-1] = 0.39101
    # error[-1] = 0.00924
    # flux[-2] = 0.2084
    # error[-2] = 0.00958


    # Turn these into a 2D array
    photometry = np.c_[flux, error]

    return photometry


############# DEFINE FITTING INSTRUCTIONS ##############

# Load filter list
filt_lis = np.loadtxt(filtFile, dtype='str')

exp = {}                                  # Tau-model star-formation history component
exp["age"] = (0.01, 15.)                   # Vary age between 10 Myr and 15 Gyr. In practice
                                          # the code automatically limits this to the age of
                                          # the Universe at the observed redshift.

exp["tau"] = (0.05, 10.)                   # Vary tau between 300 Myr and 10 Gyr
exp["massformed"] = (5., 15.)             # vary log_10(M*/M_solar) between 1 and 15
exp["metallicity"] = (0.1, 3)                  # vary Z between 0 and 1.0 Z_oldsolar

dust = {}                                 # Dust component
dust["type"] = "Calzetti"                 # Define the shape of the attenuation curve
#dust["Av"] = (0., 4.)                     # Vary Av between 0 and 4 magnitudes to allow for very dusty low-z interlopers.
dust["Av"] = (0., 7)                     # Vary Av between 0 and 0.5 magnitudes to avoid dust-age degeneracy.

fit_instructions = {}                     # The fit instructions dictionary
#fit_instructions["redshift"] = lae_z
fit_instructions["redshift"] = (0, 10)
fit_instructions["dust"] = dust

# Add a constant SFH option.
constant = {}
constant["massformed"] = (1., 15.)
constant["metallicity"] = 0.2
constant["age_max"] = (0., 13.8)
constant["age_min"] = 0.

delayed = {}                         # Delayed Tau model t*e^-(t/tau)
delayed["age"] = (0.01, 13.8)           # Time since SF began: Gyr
delayed["tau"] = (0.3, 10)           # Timescale of decrease: Gyr
delayed["metallicity"] = (0.1, 3)
delayed["massformed"] = (0., 15.)

# And bursty option
burst = {}                                   # A burst component
burst["age"] = 0.1                           # Fix age to 0.1 Gyr
burst["metallicity"] = 0.2          # Vary metallicity from 0 to 2.5 Solar
burst["massformed"] = (0., 13.)  


# Nebular emission component
nebular = {}
nebular["logU"] = (-4., -2.)
#nebular["logU"] = (-2., -1.) # For extreme Lya fitting.

# Add to fit instructions

# EXPONENTIALLY DECLINING SFH
#fit_instructions["exponential"] = exp

# CONSTANT SFH
# fit_instructions["constant"] = constant

# # BURSTY SFH
#fit_instructions["burst"] = burst

# # DELAYED TAU SFH
fit_instructions["delayed"] = delayed

# NEBULAR EMISSION LINES
fit_instructions["nebular"] = nebular




############ FITTING OBJECTS BY LOOPING ###########

# Load catalogue
t = Table.read(catDir / catName)

# Restrict to sample for LBG proposal
# IDs = [387777, 615679, 878786]
# z_grism = [6.7614, 6.8484, 7.0975]  

#? The extra component to the south of the very clumpy one...
IDs = [879369]
z_grism = [7.0975]

t = t[[id in IDs for id in t['ID']]]



# print(load_obj(IDs[0]))
# exit()

for i, ID in enumerate(IDs):

    lephareID = IDs[i]

    galaxy = pipes.galaxy(str(ID), load_obj, spectrum_exists=False, filt_list=filt_lis)

    # Update the redshift prior to the lephare zphot
    fit_instructions["redshift"] = (z_grism[i]-0.2, z_grism[i]+0.2)

    # fig = galaxy.plot()
    # exit()

    fit = pipes.fit(galaxy, fit_instructions, run='LBG_proposal')
    fit.fit(verbose=True) 
    
    fig = fit.plot_spectrum_posterior(save=True, show=False)
    fig = fit.plot_sfh_posterior(save=True, show=False)
    fig = fit.plot_corner(save=True, show=False)

    list(fit.posterior.samples)

    # Get best fit redshift
    # redshift = lae_z
    redshift = np.percentile(fit.posterior.samples["redshift"], (50))
    print(np.mean(redshift))

    # Get spectrum
    spec_post = fit.posterior.samples["spectrum_full"]
    post = np.percentile(spec_post, (16, 50, 84), axis=0).T
    wavs = fit.posterior.model_galaxy.wavelengths*(1.+redshift)

    # convert to microJy using the conversion from Adam
    # this is just Fnu = Flamda*lambda**2/c all in units of A
    conversion = (wavs**2)/(10**-29*2.9979*10**18)
    microJy = post[:,1]*conversion

#    microJyA = post[:,1]*conversion
#    microJyB = post[:,2]*conversion

    print(post.shape)

    # Then you can just save to a fits file
#    newtb1 = Table([wavs, wavs/(1.+redshift), microJy, microJyA, microJyB], names = ['wave', 'restwave', 'flux', 'flux2', 'flux3'])
    newtb1 = Table([wavs, wavs/(1.+redshift), microJy], names = ['wave', 'restwave', 'flux'])
    outName = outdir / '{0}_microJy.fits'.format(lephareID)
    newtb1.write(outName, overwrite = True)
    print('SED saved to ', outName)

    # also extract the integrated photometry points
    phot_post = np.percentile(fit.posterior.samples["photometry"],
                                  (16, 50, 84), axis=0).T

    filterwavs = fit.galaxy.filter_set.eff_wavs

    # convert the units
    conversion = (filterwavs**2)/(10**-29*2.9979*10**18)
    microJy = phot_post[:,1]*conversion

    newtb2 = Table([filterwavs, filterwavs/(1.+redshift), microJy], names = ['wave', 'restwave', 'flux'])
    outName = outdir / '{0}_microJy_synthpoints.fits'.format(lephareID)
    newtb2.write(outName, overwrite = True)
    print('Synthetic photometry saved to ', outName)

    # Also extract P(z)
    # z_phot = fit.posterior.samples["redshift"]
    # #print(z_phot)
    # newtb3 = Table([z_phot], names=['z_phot'])
    # outName = outdir / '{0}_Pz.fits'.format(lephareID)
    # newtb3.write(outName, overwrite = True)
    # print('P(z) saved to ', outName)

    #for key, value in fit.posterior.samples.items():
    #    print(key)

    # Also extract dust extinction
    Av = fit.posterior.samples["dust:Av"]
    newtb4 = Table([Av], names=['Av'])
    outName = outdir / '{0}_Av.fits'.format(lephareID)
    newtb4.write(outName, overwrite = True)
    print('Av saved to ', outName)

    # Also extract stellar mass
    mass = fit.posterior.samples["stellar_mass"]
    newtb5 = Table([mass], names=['M'])
    outName = outdir / '{0}_stellarMass.fits'.format(lephareID)
    newtb5.write(outName, overwrite = True)
    print('Stellar mass saved to ', outName)

    # Print available components of fit.posterior
    print(fit.posterior.samples.keys())




    # Print all the available fit.posterior.samples keys
    print(fit.posterior.samples.keys())

    # Print best value of each parameter
    for key in fit.posterior.samples.keys():
        print(key, np.percentile(fit.posterior.samples[key], (50)))

        # Get 1sigma errors
        low = np.percentile(fit.posterior.samples[key], (16))
        high = np.percentile(fit.posterior.samples[key], (84))
        print('1sigma errors:', low, high)



    ############! Measuring Oiii + Hbeta equivalent width ##############
    # Get the spectrum
    spec_post = fit.posterior.samples["spectrum_full"]
    post = np.percentile(spec_post, (16, 50, 84), axis=0).T
    wavs = fit.posterior.model_galaxy.wavelengths*(1.+redshift)

    # convert to microJy using the conversion from Adam
    # this is just Fnu = Flamda*lambda**2/c all in units of A
    conversion = (wavs**2)/(10**-29*2.9979*10**18)
    spec = post[:,1]*conversion

    # plot the spectrum
    #plt.plot(wavs, spec[0])
    #plt.yscale('log')
#    plt.ylim(1e-32, 1e-28)

    #plt.show()



    cont_mask = ((wavs > 4500) & (wavs < 4700)) | ((wavs > 5100) & (wavs < 5500))
    cont_level = np.median(spec[cont_mask])  # Use median flux as continuum level


    # Define the observed-frame wavelength range for [OIII] + Hbeta
    oiii_hb_mask = (wavs > 4825) & (wavs < 5030)

    # Extract flux values in the line region and sum them
    d_wav = np.mean(np.diff(wavs))  # Wavelength bin width in Å
    Oiii_Hb_flux = np.sum((spec[oiii_hb_mask] - cont_level) * d_wav)  # Integrated flux above continuum


    plt.plot(wavs, spec, label="Median Spectrum")
    plt.axhline(y=cont_level, color="r", linestyle="--", label="Continuum Level")
    plt.axvspan(4500, 4700, color="gray", alpha=0.3, label="Continuum Regions")
    plt.axvspan(5100, 5500, color="gray", alpha=0.3)
    plt.axvspan(4825, 5030, color="blue", alpha=0.3, label="[OIII] + Hbeta Region")
    plt.xlim(500, 6000)
    plt.legend()
    plt.xlabel("Wavelength (Å)")
    plt.ylabel("Flux")
    plt.yscale("log")
    plt.show()


    # Compute equivalent width in the observed frame
    EW_obs = Oiii_Hb_flux / cont_level  # Correct formula

    # Convert to rest-frame equivalent width
    #EW_rest = EW_obs / (1 + lae_z)
    EW_rest = EW_obs / (1 + z_grism[i])

    print('Equivalent width of Oiii+Hbeta (rest-frame):', EW_rest)


############ TABLE OF POSTERIORS ###############

# Save outputs to a LaTeX table. Exits before we get to the actual BAGPIPES fitting below.
# if table:

#     bagDir = '/mnt/vardy/vardygroupshare/HSC_SSP_DR3/codes/pipes/sed/{0}/low_Av/'.format(field)

#     if field == 'XMM':
#         IDs = np.arange(0, 21)
# #        IDs = np.arange(0, 5) # lya

#     if field == 'CDFS':
#         IDs = np.arange(0, 10)

#     ids = []
#     extinc = []
#     z = []
#     M = []
#     z_phots = []
#     z_photsBP = []


#     for i, id in enumerate(IDs):

#         lephareID = cat_ID(id)
#         print(lephareID)

#         z_phot = cat_z(id, field)
#         z_phots = np.append(z_phots, z_phot)

#         ids = np.append(ids, lephareID)

#         ### Redshift ###
#         bagPz = Table.read(bagDir + '{0}_Pz.fits'.format(lephareID))
#         redshifts = np.percentile(bagPz['z_phot'], (16, 50, 84), axis=0).T

#         # Get errors
#         z_ = round(redshifts[1], 2)
#         zU = str(round(redshifts[2]-redshifts[1], 2))
#         zL = str(round(redshifts[1]-redshifts[0], 2))

#         # Make table entry
#         z = np.append(z, '$' + str(z_) + '^{+' + zU + '}' + '_{-' + zL + '}$')
#         z_photsBP = np.append(z_photsBP, z_)

#         ### DUST ###
#         dust = Table.read(bagDir + '{0}_Av.fits'.format(lephareID))
#         Av = np.percentile(dust['Av'], (16, 50, 84), axis=0).T

#          # Erors
#         Av_ =  round(Av[1], 2)
#         AvU = str(round(Av[2]-Av[1], 2))
#         AvL = str(round(Av[1]-Av[0], 2))

#         # Table entry
#         extinc = np.append(extinc, '$' + str(Av_) + '^{+' + AvU + '}' + '_{-' + AvL + '}$')

#         ### Stellar Mass ###
#         stellarMass = Table.read(bagDir + '{0}_stellarMass.fits'.format(lephareID))
#         mass = np.percentile(stellarMass['M'], (16, 50, 84), axis=0).T

#         # Errors
#         mass_ = round(mass[1], 2)
#         massU = str(round(mass[2]-mass[1], 2))
#         massL = str(round(mass[1]-mass[0], 2))

#         # Mass
#         M = np.append(M, '$' + str(mass_) + '^{+' + massU + '}' + '_{-' + massL + '}$')

#     data = {'ID':ids, 'z':z, 'Av':extinc, 'stellar mass':M}

#     os.chdir(tableDir)
#     ascii.write(data, 'bagpipes_{0}_low_Av.txt'.format(field), format='latex', overwrite=True)

# #    plt.scatter(z_phots, z_photsBP)
# #    plt.xlabel('LePhare')
# #    plt.ylabel('BAGPIPES')
# #    plt.show()

#     print(np.mean(Av_))

#     exit()

# def cat_ID(ID):

#     """ get lephare ID. """

#     # Open catalogue
#     cat = Table.read(catDir / catName, hdu=1)

#     # Extract object with ID
#     cat['ID'].format

#     # BAGPIPES ID is just the row number
#     ID = int(ID)
#     print(ID)
#     cat = cat[:][ID]

#     return str(cat['ID'])

# def cat_z(ID, field):

#     """ get lephare redshift. """

#     # Open catalogue
#     if field =='XMM':
#         cat = Table.read('vmax/vmax_DR3_XMM.txt', format='ascii.commented_header')
#     if field =='CDFS':
#         cat = Table.read('vmax/vmax_vis_check2_CDFS.txt', format='ascii.commented_header')
#     # Extract object with ID

#     cat['ID'].format

#     # BAGPIPES ID is just the row number
#     ID = int(ID)
#     cat = cat[:][ID]

#     return cat['z']