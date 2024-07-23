#!/usr/bin/env python3

"""
Plot the outputs of the LePhare fitting.

Created: Friday 12th July 2024.
"""

from astropy.io import ascii
from astropy.table import Table
import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path
import glob
from sed_fitting_codes import *
import sys
from matplotlib.backends.backend_pdf import PdfPages

# Backend for running on queue
import matplotlib as mpl
mpl.use('Agg')

cutout_path = Path.cwd().parents[0] / 'cutouts'
sys.path.append(str(cutout_path))
from cutout_codes import *

plt.rcParams['axes.linewidth'] = 2.5
plt.rcParams.update({'font.size': 15})
plt.rcParams['figure.dpi'] = 100

# Function to convert magnitudes to fluxes
def mag_to_flux(mag):
    return 10**(-0.4 * (mag + 48.6))

#! Output PDF file
output_dir = Path.cwd().parents[1] / 'plots' / 'seds'
output_pdf = 'det_Je_J_LBG.pdf'

# Set up directories
#zphot_dir = Path.cwd().parents[1] / 'data' / 'sed_fitting' / 'zphot' / 'test_euclid'
zphot_dir = Path.cwd().parents[1] / 'data' / 'sed_fitting' / 'zphot' / 'det_Je_J_LBG'

# Read in the input .in file
input_dir = Path.home().parents[1] / 'hoy' / 'temporaryFilesROHAN' / 'lephare' / 'inputs' / 'euclid'
#input_name = 'euclid_test.in'
#input_name = 'euclid_lya.in'
input_name = 'det_Je_J.in'
flux_table = Table.read(input_dir / input_name, format='ascii.commented_header')

# Define start and end points of each section in the LePhare .spec file
ds_phot = 9
de_phot = 28

ds_mod = 209
de_mod = -1

ds_pz = 28
de_pz = 209

ds_param = 3
de_param = 9

# From these values, we can get the number of filters which will be useful for extracting fluxes from the .in file.
n_bands = de_phot - ds_phot

# Double the value to correspond to each flux and error
n_in = 2 * n_bands

# Column names
names_phot=['phot', 'yerr', 'wlen', 'xerr', 'modelPhot', 'col6', 'col7']
names_sed = ['wlen', 'flux']
names_zpdf = ['z', 'P']
names_param = ['Type', 'Nline', 'Model', 'Library', 'Nband', 'Zphot', 'Zinf', 'Zsup', 'Chi2', 'PDF', 'Extlaw', 'EB-V', 'Lir', 'Age', 'Mass', 'SFR', 'SSFR']

# Get filters
filter_dict = filter_widths()

# Load the crossmatched catalogue of known objects
crossmatch_dir = Path.cwd().parents[1] / 'data' / 'catalogues'
crossmatch_name = 'XMATCH_COSMOS_5sig_Je_3sig_J_nonDet_HSC_G_nonDet_HSC_R_nonDet_HSC_I_nonDet_HSC_Z_nonDet_HSC_Y_nonDet_Y.fits'
crossmatch = Table.read(crossmatch_dir / crossmatch_name, format='fits')
crossmatch['Redshift'] = crossmatch['Redshift'].astype(float)

# Load the parent catalogue to get RA,DEC
parent_cat = Table.read(crossmatch_dir / 'COSMOS_5sig_Je_3sig_J_nonDet_HSC_G_nonDet_HSC_R_nonDet_HSC_I_nonDet_HSC_Z_nonDet_HSC_Y_nonDet_Y.fits', format='fits')

# Collect all .spec files
spec_files = glob.glob(str(zphot_dir / '*.spec'))

with PdfPages(str(output_dir/output_pdf)) as pdf:
    #! Loop through files
    for i, spec_file in enumerate(spec_files):
        print(f'Object {i+1} of {len(spec_files)}')

        phot = ascii.read(spec_file, format='basic', data_start=ds_phot, data_end=de_phot, delimiter=' ', names=names_phot)
        zpdf = ascii.read(spec_file, format='basic', data_start=ds_pz, data_end=de_pz, delimiter=' ', names=names_zpdf)
        sed = ascii.read(spec_file, format='basic', data_start=ds_mod, data_end=de_mod, delimiter=' ', names=names_sed)
        params = ascii.read(spec_file, format='basic', data_start=ds_param, data_end=de_param, delimiter=' ', names=names_param)

        # Get the ID of the object from the file name
        ID = spec_file.split('/')[-1].split('Id')[-1].lstrip('0').split('.spec')[0]
        print(ID)

        # Index in input file of this ID
        if len(np.where(flux_table['ID'] == int(ID))[0]) == 0:
            raise ValueError('ID not found: check the correct directory is in use.')

        idx = np.where(flux_table['ID'] == int(ID))[0][0]

        # Match ID with the input file to find the object.
        object_flux = flux_table[idx]

        # Initialize flux and error arrays, to add .in data.
        flux = []
        error = []

        # In object_FLUX, ignore first and last columns. Get flux and errors, which are in format flux, err, flux, err,...
        for i in range(1, n_in, 2):
            flux.append(object_flux[i])
            error.append(object_flux[i+1])

        flux = np.array(flux)
        error = np.array([np.abs(e) for e in error])

        # Extract each model by finding where wavelength jumps.
        model_wlens = np.split(sed['wlen'], np.where(np.diff(sed['wlen']) < 0)[0] + 1)
        model_fluxes = np.split(sed['flux'], np.where(np.diff(sed['wlen']) < 0)[0] + 1)

        # Filter centres and widths
        central_wavelengths = [value[0] for value in filter_dict.values()]
        central_wavelengths = [c*10000 for c in central_wavelengths]

        fwhm_values = [value[1] for value in filter_dict.values()]
        fwhm_values = [f*10000 for f in fwhm_values]

        # Get each of the SEDs
        primary_sed = mag_to_flux(model_fluxes[0])
        primary_wlen = model_wlens[0]

        secondary_sed = mag_to_flux(model_fluxes[1])
        secondary_wlen = model_wlens[1]

        if len(model_fluxes) > 2:
            stellar_sed = mag_to_flux(model_fluxes[2])
            stellar_wlen = model_wlens[2]
        else:
            stellar_sed = None
            stellar_wlen = None

        # Get parameters
        secondary_is_stellar = False

        zphot_1 = round(params['Zphot'][0], 2)
        chi2_1 = round(params['Chi2'][0], 1)

        print('GAL 1 SOLUTION: ', zphot_1, chi2_1)

        zphot_2 = round(params['Zphot'][1], 2)
        chi2_2 = round(params['Chi2'][1], 1)

        print('GAL 2 SOLUTION: ', zphot_2, chi2_2)

        if int(zphot_2) == -1:
            secondary_is_stellar = True

        chi2_star = round(params['Chi2'][-1], 1)
        print('STELLAR SOLUTION: ', chi2_star)

        #! Sigma array
        sigma = flux / error

        # Create figure with multiple axes
        fig = plt.figure(figsize=(10, 6))
        ax1 = fig.add_subplot(111)
        left, bottom, width, height = [0.70, 0.23, 0.18, 0.18]
        ax2 = fig.add_axes([left, bottom, width, height])

        ############! PLOT RESULTS OF FITTING #############
        # Plot best model
        ax1.plot(primary_wlen, primary_sed, color='deepskyblue', alpha=1.0, label=f'z = {zphot_1}, ' + r'$\chi^{2} = $' + f'{chi2_1}', linewidth=3, zorder=4)

        # plot secondary galaxy solution and stellar model
        if not secondary_is_stellar:
            ax1.plot(secondary_wlen, secondary_sed, color='darkorange', label=f'z = {zphot_2}, ' + r'$\chi^{2} = $' + f'{chi2_2}', linewidth=2, alpha=0.7)
            ax1.plot(stellar_wlen, stellar_sed, color='red', label=r'$\chi^{2} = $' + f'{chi2_star}', linewidth=2, alpha=0.7)

        if secondary_is_stellar:
            ax1.plot(secondary_wlen, secondary_sed, color='red', label=r'$\chi^{2} = $' + f'{chi2_star}', linewidth=2, alpha=0.7)

        ax2.plot(zpdf['z'], zpdf['P'], color='black', linewidth=2)
        ax2.set_xlim(0, 10)
        ax2.set_xlabel(r'$z_{\mathrm{phot}}$')
        ax2.set_ylabel(r'$P(z)$')

        ###########! PLOT MODEL AND REAL PHOTOMETRY ##############

        # Model photometry
        model_photometry = mag_to_flux(phot['modelPhot'])
        ax1.scatter(central_wavelengths, model_photometry, marker='o', s=100, alpha=0.6, zorder=5, edgecolor='black', facecolor='none', linewidth=2)

        # Real photometry with upper limits for sigma < 2
        for i in range(len(flux)):
            if sigma[i] < 2:
                ax1.scatter(central_wavelengths[i], flux[i] + 2 * error[i], marker='v', color='black', s=100, zorder=6)
            else:
                ax1.errorbar(central_wavelengths[i], flux[i], yerr=error[i], fmt='o', color='black', markersize=8, zorder=6, elinewidth=2)

        ############! CHECK IF THIS IS AN EXISTING OBJECT ############
        # Check if this object is in the crossmatched catalogue
        if int(ID) in crossmatch['ID']:
            # Get the object
            obj = crossmatch[np.where(crossmatch['ID'] == int(ID))]

            # Get the redshift
            ref_id = obj['Object Name'][0]
            z_spec = round(obj['Redshift'][0], 2)

            ax1.set_title(f'ID {ID}, {ref_id}, z = {z_spec}')
        else:
            ax1.set_title(f'ID {ID}')

        ax1.set_yscale('log')
        ax1.legend(loc='upper right')
        ax1.set_ylim(3e-32, 1e-29)
        ax1.set_xlim(3000, 35000)

        ax1.set_ylabel(r'$f_{\nu}$ (erg s$^{-1}$ cm$^{-2}$ Hz$^{-1}$)')
        ax1.set_xlabel(r'$\lambda (\AA)$')

        pdf.savefig(fig)
        plt.close(fig)

        ############! CUTOUT ##############
        # Get ra, dec of object from parent catalogue
        obj = parent_cat[np.where(parent_cat['ID'] == int(ID))]
        ra = obj['RA'][0]
        dec = obj['DEC'][0]

        cutout_fig, cutout_axs = Cutout(ra, dec, size=10., save_cutout=False)

        pdf.savefig(cutout_fig)
        plt.close(cutout_fig)


