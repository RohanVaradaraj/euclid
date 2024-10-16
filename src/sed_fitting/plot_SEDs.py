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
from selection import generate_selection_name, generate_input_name
from astropy.coordinates import SkyCoord
import json

# Backend for running on queue
import matplotlib as mpl
mpl.use('Agg')

cutout_path = Path.cwd().parents[0] / 'cutouts'
sys.path.append(str(cutout_path))
from cutout_codes import *
from vista_cutouts import *

crosstalk_path = Path.cwd().parents[3] / 'HSC_SSP_DR3' / 'codes'
sys.path.append(str(crosstalk_path))
from full_catalogue_codes import label_ct

plt.rcParams['axes.linewidth'] = 2.5
plt.rcParams.update({'font.size': 15})
plt.rcParams['figure.dpi'] = 100

if len(sys.argv) > 1:
    filters_json = sys.argv[1]
    filters = json.loads(filters_json)
    bools_json = sys.argv[2]
    bools = json.loads(bools_json)
    all_filters_json = sys.argv[3]
    all_filters = json.loads(all_filters_json)

def stellar_type(model):
    stellar_dict = {
        1: 'M4',
        2: 'M5',
        3: 'M6',
        4: 'M7',
        5: 'M8',
        6: 'M9',
        7: 'L0',
        8: 'L1',
        9: 'L2',
        10: 'L3',
        11: 'L4',
        12: 'L5',
        13: 'L6',
        14: 'L7',
        15: 'L8',
        16: 'L9',
        17: 'T0',
        18: 'T1',
        19: 'T2',
        20: 'T3',
        21: 'T4',
        22: 'T5',
        23: 'T6',
        24: 'T7',
        25: 'T8',
    }
    return stellar_dict[model]

# Function to convert magnitudes to fluxes
def mag_to_flux(mag):
    return 10**(-0.4 * (mag + 48.6))

#! What is the main object we want to plot? highz, bd, dusty or lya?
object_type = 'highz'

#! Label crosstalk?
label_crosstalk = True

#! Output PDF file
output_dir = Path.cwd().parents[1] / 'plots' / 'seds'

#! Set up the output PDF name from the det and nondet filters
det_list = [f for f, t in filters.items() if t['type'] == 'detection']
if len(det_list) != 0:
    base_det = 'det_' + '_'.join(det_list)
if len(det_list) == 0:
    stack_list = [f for f, t in filters.items() if t['type'] == 'stacked-detection']
    stack_filters = stack_list[0].split('+')
    base_det = 'det_' + '_'.join(stack_filters)
    det_list = stack_filters

base_nondet = 'nonDet_' + '_'.join([f for f, t in filters.items() if t['type'] == 'non-detection'])
if object_type == 'highz':
    output_pdf = f'{base_det}_{base_nondet}.pdf'
else:
    output_pdf = f'{base_det}_{base_nondet}_{object_type}.pdf'
print('Saving to: ', output_pdf)

# Set up directories
# zphot_dir = Path.cwd().parents[1] / 'data' / 'sed_fitting' / 'zphot' / 'best_fits' / base_det
zphot_folder = base_det + f'_best_{object_type}'
zphot_dir = Path.cwd().parents[1] / 'data' / 'sed_fitting' / 'zphot' / 'best_fits' / zphot_folder
print(zphot_dir)

# Crossmatched catalogue name to get existing sources
crossmatch_name = 'all_COSMOS_highz.fits'

# Load the crossmatched catalogue of known objects
crossmatch_dir = Path.cwd().parents[1] / 'data' / 'ref_catalogues'
crossmatch = Table.read(crossmatch_dir / crossmatch_name, format='fits')
crossmatch['Redshift'] = crossmatch['Redshift'].astype(float)

# Read in the input .in file
input_dir = Path.home().parents[1] / 'hoy' / 'temporaryFilesROHAN' / 'lephare' / 'inputs' / 'euclid'
input_name = base_det + '.in'
flux_table = Table.read(input_dir / input_name, format='ascii.commented_header')

# Directory to read in the .out file
lephare_out_dir = Path.home() / 'lephare' / 'lephare_dev' / 'test'

# Determine the number of filters

# Define start and end points of each section in the LePhare .spec file
ds_phot = 9

print(det_list)

if 'e' in det_list[0]:
    de_phot = 27
if det_list == ['Y', 'J']:
    de_phot = 20

if 'e' in det_list[0]:
    ds_mod = 228
if det_list == ['Y', 'J']:
    ds_mod = 220

de_mod = -1

if 'e' in det_list[0]:
    ds_pz = 27
if det_list == ['Y', 'J']:
    ds_pz = 20

if 'e' in det_list[0]:
    de_pz = 228
if det_list == ['Y', 'J']:
    de_pz = 220

ds_param = 3
de_param = 9

#! Without euclid
# Check if 'no_euclid' is in input name
if 'no_euclid' in input_name:
    ds_phot = 9
    de_phot = 24

    ds_mod = 205
    de_mod = -1

    ds_pz = 24
    de_pz = 205

    ds_param = 3
    de_param = 9

# #! Brown dwarf fitting
# if 'bd' or 'lya' in output_pdf:
#     ds_phot = 9
#     de_phot = 24

#     ds_mod = 224
#     de_mod = -1

#     ds_pz = 24
#     de_pz = 224

#     ds_param = 3
#     de_param = 9

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

# Remove f444w
if 'dusty' not in output_pdf:
    filter_dict.pop('f444w')
    filter_dict.pop('ch1cds')
    filter_dict.pop('ch2cds')
    print('Removed f444w, ch1cds, ch2cds')

if det_list == ['Y', 'J']:
    filter_dict.pop('VIS')
    filter_dict.pop('Ye')
    filter_dict.pop('Je')
    filter_dict.pop('He')
    print('Removed VIS, Ye, Je, He')

if 'bd' or 'lya' in output_pdf:
    filter_dict.pop('f277w')
    filter_dict.pop('HSC-G_DR3')
    filter_dict.pop('HSC-R_DR3')
    print('Removed f277w, HSC-G_DR3, HSC-R_DR3')

# If running only VIISTA: Remove items with keys VIS, Ye, Je, He
if 'no_euclid' in input_name:
    filter_dict.pop('VIS')
    filter_dict.pop('Ye')
    filter_dict.pop('Je')
    filter_dict.pop('He')


# Load the parent catalogue to get RA,DEC
# Generate the name of the parent catalogue
parent_cat_dir = Path.cwd().parents[1] / 'data' / 'catalogues'
parent_cat_name = generate_selection_name('COSMOS', filters)

# Label crosstalk
if label_crosstalk:
    label_ct(str(parent_cat_dir / parent_cat_name), fieldName='COSMOS')

parent_cat = Table.read(parent_cat_dir / parent_cat_name, format='fits')

# Collect all .spec files
spec_files = glob.glob(str(zphot_dir / '*.spec'))

# Sort the files in increasing numerical order
spec_files = sorted(spec_files, key=lambda x: int(x.split('/')[-1].split('Id')[-1].lstrip('0').split('.spec')[0]))

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

        if len(model_fluxes) < 2:
            print('Only one solution')
            continue

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
        dz_sup = round(params['Zsup'][0], 2)
        dz_inf = round(params['Zinf'][0], 2)
        chi2_1 = round(params['Chi2'][0], 1)

        # Selection criterion for HSC+VISTA+JWST+Euclid
        # if chi2_1 > 20.4:
        #     print('Would not be selected as a high-z galaxy')
        #     continue

        print('GAL 1 SOLUTION: ', zphot_1, chi2_1)

        #! Read in the secondary solutions from the lephare_out tables
        bd_out_name = base_det + '_bd.out'
        bd_out = ascii.read(lephare_out_dir / bd_out_name, format='no_header')
        chi2_star = round(bd_out[bd_out['col1'] == int(ID)]['col21'][0], 1)

        dusty_out_name = base_det + '_dusty.out'
        dusty_out = ascii.read(lephare_out_dir / dusty_out_name, format='no_header')
        chi2_2 = round(dusty_out[dusty_out['col1'] == int(ID)]['col15'][0], 1)

        zphot_2 = round(params['Zphot'][1], 2)
        #chi2_2 = round(params['Chi2'][1], 1)

        print('GAL 2 SOLUTION: ', zphot_2, chi2_2)

        if int(zphot_2) == -1:
            secondary_is_stellar = True

        #chi2_star = round(params['Chi2'][-1], 1)
        print('STELLAR SOLUTION: ', chi2_star)

        stellar_model = params['Model'][-1]
        mlt_type = stellar_type(stellar_model)

        #! Sigma array
        sigma = flux / error

        # Create figure with multiple axes
        fig = plt.figure(figsize=(10, 6))
        ax1 = fig.add_subplot(111)
        left, bottom, width, height = [0.70, 0.23, 0.18, 0.18]
        ax2 = fig.add_axes([left, bottom, width, height])

        ############! PLOT RESULTS OF FITTING #############
        # Plot best model
        ax1.plot(primary_wlen, primary_sed, color='deepskyblue', alpha=1.0, label= rf'$z = {zphot_1}^{{+{round(dz_sup-zphot_1, 2)}}}_{{-{round(zphot_1-dz_inf, 2)}}}$, ' + r'$\chi^{2} = $' + f'{chi2_1}', linewidth=3, zorder=4)

        # plot secondary galaxy solution and stellar model
        if not secondary_is_stellar:
            ax1.plot(secondary_wlen, secondary_sed, color='darkorange', label=f'z = {zphot_2}, ' + r'$\chi^{2} = $' + f'{chi2_2}', linewidth=2, alpha=0.7)
            ax1.plot(stellar_wlen, stellar_sed, color='red', label=r'$\chi^{2} = $' + f'{chi2_star}, type = {mlt_type}', linewidth=2, alpha=0.7)

        if secondary_is_stellar:
            ax1.plot(secondary_wlen, secondary_sed, color='red', label=r'$\chi^{2} = $' + f'{chi2_star}, type = {mlt_type}', linewidth=2, alpha=0.7)

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

        ############! CROSSMATCH WITH EXISTING SOURCES ############
        # Get the RA and DEC of the object from the parent catalog
        obj_ra = parent_cat[np.where(parent_cat['ID'] == int(ID))]['RA'][0]
        obj_dec = parent_cat[np.where(parent_cat['ID'] == int(ID))]['DEC'][0]

        # Get the RA and DEC of all objects in the crossmatch catalog
        xmatch_ra = crossmatch['RA']
        xmatch_dec = crossmatch['DEC']

        # Create SkyCoord objects for both catalogs
        catalog1 = SkyCoord(ra=obj_ra * u.deg, dec=obj_dec * u.deg)  # Single object
        catalog2 = SkyCoord(ra=xmatch_ra * u.deg, dec=xmatch_dec * u.deg)  # Crossmatch catalog

        # Perform crossmatch using match_to_catalog_sky
        idx, d2d, d3d = catalog1.match_to_catalog_sky(catalog2)

        # Set the tolerance to 1 arcsecond
        tolerance = 1 * u.arcsec

        # side mission: check if there is a 'CLASS' column in the parent catalogue
        if 'CLASS' in parent_cat.colnames:
            ct = parent_cat[np.where(parent_cat['ID'] == int(ID))]['CLASS'][0]
            if ct > 0:
                print('Possible crosstalk')
        else:
            ct = 0.

        # Check if the match is within the tolerance
        if d2d < tolerance:
            # Get the index of the matched object
            matched_idx = idx  # Since idx is scalar when matching a single object

            # Retrieve the corresponding object from the crossmatch catalog
            obj = crossmatch[matched_idx]

            # Extract the object name and redshift
            name = obj['Object Name']
            redshift = round(obj['Redshift'], 2)

            # Update the title of the plot with the matched information. Add crosstalk if it is non-zero.
            title_string = f'ID {ID}, z = {redshift}'
            if ct > 0:
                title_string += f', POSSIBLE CROSSTALK'
            ax1.set_title(title_string)
        else:
            # If no match is found, just set the title with the ID
            title_string = f'ID {ID}'
            if ct > 0:
                title_string += f', POSSIBLE CROSSTALK'
            ax1.set_title(title_string)



        ax1.set_yscale('log')
        ax1.legend(loc='upper right')
        ax1.set_ylim(3e-32, 1e-29)
        ax1.set_xlim(3000, 40000)

        ax1.set_ylabel(r'$f_{\nu}$ (erg s$^{-1}$ cm$^{-2}$ Hz$^{-1}$)')
        ax1.set_xlabel(r'$\lambda (\AA)$')

        pdf.savefig(fig)
        plt.close(fig)

        ############! CUTOUT ##############
        # Get ra, dec of object from parent catalogue
        obj = parent_cat[np.where(parent_cat['ID'] == int(ID))]
        ra = obj['RA'][0]
        dec = obj['DEC'][0]

        if det_list == ['Y', 'J']:
            cutout_fig, cutout_axs = VistaCutout(ra, dec, size=10., save_cutout=False)
        else:
            cutout_fig, cutout_axs = Cutout(ra, dec, size=10., save_cutout=False)

        pdf.savefig(cutout_fig)
        plt.close(cutout_fig)


