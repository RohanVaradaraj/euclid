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
#from sed_fitting_codes import *
import sys
from matplotlib.backends.backend_pdf import PdfPages
#from make_comparison_tables import stellar_type
import matplotlib.ticker as mticker

cutout_path = Path.cwd().parents[0] / 'cutouts'
sys.path.append(str(cutout_path))
from cutout_codes import *

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

def mag_to_Jy(m):
    '''Convert mags and their errors to Jy'''
    flux = 10**(-0.4*(m+48.6))
    return flux

def flux_to_mag(flux):
    '''Convert flux to mag, for the secondary y axis'''
    mag = -2.5*np.log10(flux)-48.6
    return mag

# survey_style = {
#     'CFHT/HSC': {'marker': 'o', 'color': 'dimgray'},   
#     'VISTA' : {'marker': 'o', 'color': 'black'},       # VISTA photometry: blue squares
#     'Euclid': {'marker': 'D', 'color': 'magenta'}     # Euclid photometry: red triangles
# }

survey_style = {
    'y (HSC)': {'marker': 'D', 'color': 'deepskyblue'},   
    'Y (VISTA)' : {'marker': 'D', 'color': 'orange'},     
    r'$Y_{E} \ (\mathrm{Euclid})$': {'marker': 'D', 'color': 'deeppink'},
    'Other' : {'marker': 'o', 'color': 'black'}
}


# def get_survey_name(filter_name):
#     if ('HSC' in filter_name) or ('CFHT' in filter_name):
#         return 'CFHT/HSC'
#     if filter_name in ['VIS', 'Ye', 'Je', 'He']:
#         return 'Euclid'
#     if (filter_name in ['Y', 'J', 'H', 'Ks']):
#         return 'VISTA'

def get_survey_name(filter_name):
    if filter_name == 'HSC-Y_DR3':
        return 'y (HSC)'
    if filter_name == 'Y':
        return 'Y (VISTA)'
    if filter_name == 'Ye':
        return r'$Y_{E} \ (\mathrm{Euclid})$'
    else:
        return 'Other'

#! Output PDF file
output_dir = Path.cwd().parents[1] / 'plots' / 'seds'
output_pdf = 'poster_det_Ye_y_LAE.pdf'
#output_pdf = 'det_NB0921_Ye.pdf'

# Set up directories
#zphot_dir = Path.cwd().parents[1] / 'data' / 'sed_fitting' / 'zphot' / 'test_euclid'
zphot_dir = Path.cwd().parents[1] / 'data' / 'sed_fitting' / 'zphot' / 'det_Ye_y_LAE'

# Read in the input .in file
input_dir = Path.home().parents[1] / 'hoy' / 'temporaryFilesROHAN' / 'lephare' / 'inputs' / 'euclid'
#input_name = 'euclid_test.in'
#input_name = 'euclid_lya.in'
input_name = 'det_Ye_Y.in'
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

def filter_widths():
    """
    Returns dictionary of filter central wavelengths and FWHMs.
    """

    filt_dict = {
        'CFHT-u': (0.3783, 0.0704),
        'CFHT-g': (0.4858, 0.1440),
        'CFHT-r': (0.6253, 0.1219),
        #'CFHT-iy': (0.7696, 0.1368),
        'CFHT-z': (0.8904, 0.0907),
        'HSC-G_DR3': (0.4816, 0.1386),
        'HSC-R_DR3': (0.6234, 0.1504),
        'HSC-I_DR3': (0.7741, 0.1552),
        'HSC-NB0816_DR3': (0.8177, 0.0113),
        'HSC-Z_DR3': (0.8912, 0.0773),
        'HSC-NB0921_DR3': (0.9214, 0.0134),
        'HSC-Y_DR3': (0.9780, 0.0783),
        'Y': (1.0214, 0.0926),
        'J': (1.2544, 0.1725),
        'H': (1.6465, 0.2916),
        'Ks': (2.1484, 0.3092),
        'VIS': (0.7180, 0.3699),
        'Ye': (1.0812, 0.2626),
        'Je': (1.3670, 0.3991),
        'He': (1.7708, 0.4994)
    }
    return filt_dict

filter_dict = filter_widths()

# If running only VIISTA: Remove items with keys VIS, Ye, Je, He
if 'no_euclid' in input_name:
    filter_dict.pop('VIS')
    filter_dict.pop('Ye')
    filter_dict.pop('Je')
    filter_dict.pop('He')

# Load the crossmatched catalogue of known objects
crossmatch_dir = Path.cwd().parents[1] / 'data' / 'catalogues'
crossmatch_name = 'XMATCH_COSMOS_5sig_Ye_2sig_VISTA_Y_nonDet_HSC_G_nonDet_HSC_R_nonDet_HSC_I.fits'
crossmatch = Table.read(crossmatch_dir / crossmatch_name, format='fits')
crossmatch['Redshift'] = crossmatch['Redshift'].astype(float)

# Load the parent catalogue to get RA,DEC
parent_cat = Table.read(crossmatch_dir / 'COSMOS_5sig_Ye_2sig_HSC_Y_nonDet_HSC_G_nonDet_HSC_R_nonDet_HSC_I.fits', format='fits')
#parent_cat = Table.read(crossmatch_dir / 'XMATCH_HSC_NB.fits', format='fits')

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

        if str(ID) != '178396':
            continue

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
        dz_sup = round(params['Zsup'][0], 2)
        dz_inf = round(params['Zinf'][0], 2)
        chi2_1 = round(params['Chi2'][0], 1)

        print('GAL 1 SOLUTION: ', zphot_1, chi2_1)

        zphot_2 = round(params['Zphot'][1], 2)
        chi2_2 = round(params['Chi2'][1], 1)

        print('GAL 2 SOLUTION: ', zphot_2, chi2_2)

        if int(zphot_2) == -1:
            secondary_is_stellar = True

        chi2_star = round(params['Chi2'][-1], 1)
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
        ax1.plot(primary_wlen, primary_sed, color='blue', alpha=1.0, linewidth=3.5, zorder=4)
        print(rf'$z = {zphot_1}^{{+{round(dz_sup-zphot_1, 2)}}}_{{-{round(zphot_1-dz_inf, 2)}}}$, ' + r'$\chi^{2} = $' + f'{chi2_1}')

        # plot secondary galaxy solution and stellar model
        if not secondary_is_stellar:
            ax1.plot(secondary_wlen, secondary_sed, color='darkorange', linewidth=2.5, alpha=0.7, zorder=3)
            #ax1.plot(stellar_wlen, stellar_sed, color='red', linewidth=2.5, alpha=0.6, zorder=3)

        # if secondary_is_stellar:
        #     ax1.plot(secondary_wlen, secondary_sed, color='red', linewidth=3.5, alpha=0.6, zorder=3)

        print(r'$\chi^{2} = $' + f'{chi2_star}, type = {mlt_type}')
        print(f'z = {zphot_2}, ' + r'$\chi^{2} = $' + f'{chi2_2}')

        ax2.plot(zpdf['z'], zpdf['P'], color='black', linewidth=2)
        ax2.set_xlim(0, 10)
        ax2.set_xlabel(r'$z_{\mathrm{phot}}$')
        ax2.set_ylabel(r'$P(z)$')

        ###########! PLOT MODEL AND REAL PHOTOMETRY ##############

        # Model photometry
        model_photometry = mag_to_flux(phot['modelPhot'])
        ax1.scatter(central_wavelengths, model_photometry, marker='o', s=300, alpha=0.6, zorder=5, edgecolor='black', facecolor='none', linewidth=2, label='LAE model photometry')
        
        # Real photometry with upper limits for sigma < 2
        for i in range(len(flux)):
            survey = get_survey_name(list(filter_dict.keys())[i])
            print(survey)
            marker = survey_style[survey]['marker']
            color = survey_style[survey]['color']
            
            if sigma[i] < 2:
                ax1.scatter(central_wavelengths[i], flux[i] + 2 * error[i], marker='v', color=color, s=100, zorder=6)
            else:
                ax1.errorbar(central_wavelengths[i], flux[i], yerr=error[i], fmt=marker, color=color, markersize=12, zorder=6, elinewidth=2)

        # Add legend for photometry points
        for survey, style in survey_style.items():
            if survey == 'Other':
                ax1.scatter([], [], marker=style['marker'], color=style['color'])
            else:
                ax1.scatter([], [], marker=style['marker'], color=style['color'], label=survey)


        ############! CHECK IF THIS IS AN EXISTING OBJECT ############
        # Check if this object is in the crossmatched catalogue
        if int(ID) in crossmatch['ID']:
            # Get the object
            obj = crossmatch[np.where(crossmatch['ID'] == int(ID))]

            # Get the redshift
            ref_id = obj['Object Name'][0]
            z_spec = round(obj['Redshift'][0], 2)

            #ax1.set_title(f'ID {ID}, {ref_id}, z = {z_spec}')
        # else:
        #     ax1.set_title(f'ID {ID}')

        ax1.set_yscale('log')
        ax1.legend(loc='center right')
        ax1.set_ylim(1e-32, 1e-29)
        ax1.set_xlim(3000, 35000)

        # Convert the x axis into microns by dividing by 10000
        ax1.set_xticks([5000, 10000, 15000, 20000, 25000, 30000, 35000])
        ax1.set_xticklabels([0.5, 1, 1.5, 2, 2.5, 3, 3.5])

        ax1.set_ylabel(r'$f_{\nu}$ (erg s$^{-1}$ cm$^{-2}$ Hz$^{-1}$)')
        #ax1.set_xlabel(r'$\lambda (\AA)$')
        ax1.set_xlabel(r'$\lambda (\mu m)$')

        # Add magnitude axis
        secax = ax1.secondary_yaxis('right', functions=(flux_to_mag, mag_to_flux))

        mags = [31, 30, 29, 28, 27, 26, 25, 24, 23, 22]
        secax.set_yticks(mags)

        secax.yaxis.set_major_formatter(mticker.ScalarFormatter())

        ax1.tick_params(which='major', length=5, width=3)
        ax1.tick_params(axis='both', which='minor', length=3, width=2)

        secax.tick_params(which='major', length=5, width=3)
        secax.tick_params(axis='both', which='minor', length=5, width=2)

        secax.set_ylabel(r'$\mathrm{m_{AB}}$')

        plt.tight_layout()

        pdf.savefig(fig)
        plt.close(fig)

        ############! CUTOUT ##############
        #Get ra, dec of object from parent catalogue
        obj = parent_cat[np.where(parent_cat['ID'] == int(ID))]
        print(obj)
        ra = obj['RA'][0]
        dec = obj['DEC'][0]

        cutout_fig, cutout_axs = Cutout(ra, dec, size=10., save_cutout=False)

        pdf.savefig(cutout_fig)
        plt.close(cutout_fig)

