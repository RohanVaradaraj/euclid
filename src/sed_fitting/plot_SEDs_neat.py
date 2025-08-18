#!/usr/bin/env python3

"""
Plot the outputs of the LePhare fitting.
Created: Friday 12th July 2024.
"""

import sys
import glob
import json
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import matplotlib as mpl
from matplotlib.backends.backend_pdf import PdfPages

from pathlib import Path
from astropy.io import ascii
from astropy.table import Table
from astropy.coordinates import SkyCoord
import astropy.units as u
#plt.rcParams['text.usetex'] = True

# Configure matplotlib for non-interactive backend (for queues)
mpl.use('Agg')

# Custom imports
cutout_path = Path.cwd().parents[0] / 'cutouts'
sys.path.append(str(cutout_path))

crosstalk_path = Path.cwd().parents[3] / 'HSC_SSP_DR3' / 'codes'
sys.path.append(str(crosstalk_path))
from full_catalogue_codes import label_ct

sed_path = Path.cwd().parents[0] / 'sed_fitting'
sys.path.append(str(sed_path))

from sed_fitting_codes import *
from selection import generate_selection_name
from cutout_codes import *
from vista_cutouts import *
from all_cutouts import AllCutout
from xmm_cutouts import XMMCutout
from full_catalogue_codes import label_ct
from sed_fitting_codes import parse_spec_file

# Plotting configuration
plt.rcParams.update({
    'axes.linewidth': 2.5,
    'font.size': 15,
    'figure.dpi': 100
})

plt.rcParams.update({
    # Ticks on all sides, pointing inwards
    'xtick.top': True, 'xtick.bottom': True,
    # 'ytick.left': True, 'ytick.right': True,
    'xtick.direction': 'in', 'ytick.direction': 'in',

    # Major tick size and width
    'xtick.major.size': 6.5, 'ytick.major.size': 6.5,
    'xtick.major.width': 2, 'ytick.major.width': 2,

    # Minor tick size and width
    'xtick.minor.size': 3, 'ytick.minor.size': 3,
    'xtick.minor.width': 1.5, 'ytick.minor.width': 1.5,
})

# === Configuration ===
field_name = 'COSMOS'
individual_sed = False
save_indiv = True
indiv_ID = '381772'
indiv_pdf_name = '381772_no_euclid.pdf'  # Name of the individual SED PDF file
test = False # Run a limited number
N = 30
label_crosstalk = False
sort_by_Muv = False
Muv_avail = False
fontsize=22
remove_title = False

# === Handle command-line arguments ===
if len(sys.argv) > 1:
    filters = json.loads(sys.argv[1])
    bools = json.loads(sys.argv[2])
    all_filters = json.loads(sys.argv[3])
    run_type = json.loads(sys.argv[4])
    object_type = json.loads(sys.argv[5])
    field_name = json.loads(sys.argv[6])
else:
    filters, bools, all_filters, run_type, object_type = {}, [False]*3, {}, '', ''

euclid_blind = run_type == ''

# === Utilities ===
def stellar_type(model):
    return {
        1: 'M4', 2: 'M5', 3: 'M6', 4: 'M7', 5: 'M8', 6: 'M9',
        7: 'L0', 8: 'L1', 9: 'L2', 10: 'L3', 11: 'L4', 12: 'L5',
        13: 'L6', 14: 'L7', 15: 'L8', 16: 'L9', 17: 'T0', 18: 'T1',
        19: 'T2', 20: 'T3', 21: 'T4', 22: 'T5', 23: 'T6', 24: 'T7', 25: 'T8',
        -1: 'No fit'
    }.get(model, 'Unknown')

def mag_to_flux(mag):
    return 10**(-0.4 * (mag + 48.6))

def mag_to_Jy(m):
    return 10**(-0.4*(m + 48.6))

def flux_to_mag(flux):
    return -2.5 * np.log10(flux) - 48.6

# === Function to create IAU-compliant Euclid name ===
def truncate(value, decimals):
    factor = 10 ** decimals
    return int(value * factor) / factor

def make_euclid_name(ra_deg, dec_deg):
    c = SkyCoord(ra=ra_deg * u.deg, dec=dec_deg * u.deg, frame='icrs')
    ra_h = c.ra.hms
    dec_d = c.dec.dms

    # Truncate RA to two decimal places (i.e., 0.01s)
    ra_sec_trunc = truncate(ra_h.s, 2)
    ra_str = f"{int(ra_h.h):02d}{int(ra_h.m):02d}{int(ra_sec_trunc):02d}.{int((ra_sec_trunc % 1) * 100):02d}"

    # Truncate Dec to one decimal place (i.e., 0.1")
    sign = '+' if dec_d.d >= 0 else '-'
    dec_d_abs = abs(dec_d.d)
    dec_m_abs = abs(dec_d.m)
    dec_s_trunc = truncate(abs(dec_d.s), 1)
    dec_str = f"{sign}{int(dec_d_abs):02d}{int(dec_m_abs):02d}{int(dec_s_trunc):02d}.{int((dec_s_trunc % 1) * 10):1d}"

    return f"J{ra_str}${sign}${dec_str[1:]}"  # replace - with $-$ in LaTeX

# === Output PDF Name Construction ===
output_dir = Path.cwd().parents[1] / 'plots' / 'seds'
sed_dir = output_dir / 'individual'

# Directory to read in the .out file
lephare_out_dir = Path.home() / 'lephare' / 'lephare_dev' / 'test'

det_list = [f for f, t in filters.items() if t['type'] == 'detection']
if det_list:
    base_det = 'det_' + '_'.join(det_list)
else:
    stack_list = [f for f, t in filters.items() if t['type'] == 'stacked-detection']
    stack_filters = stack_list[0].split('+')
    base_det = 'det_' + '_'.join(stack_filters)
    det_list = stack_filters

base_nondet = 'nonDet_' + '_'.join([f for f, t in filters.items() if t['type'] == 'non-detection'])
filename_components = [field_name, base_det, base_nondet]
if object_type:
    filename_components.append(object_type)
if run_type:
    filename_components.append(run_type)

output_pdf = indiv_pdf_name if individual_sed else '_'.join(filename_components) + '.pdf'
print('Saving to: ', output_pdf)

# === SED Input & Output Directories ===
zphot_folder = f"{base_det}_{run_type}_{object_type}" if run_type else f"{base_det}_{object_type}"
zphot_dir = Path.cwd().parents[1] / 'data' / 'sed_fitting' / 'zphot' / field_name / 'best_fits' / zphot_folder
print('Taking SEDs from: ', zphot_dir)

dusty_dir = Path.cwd().parents[1] / 'data' / 'sed_fitting' / 'zphot' / field_name / (f"{base_det}_{run_type}_dusty" if run_type else f"{base_det}_dusty")
bd_dir = Path.cwd().parents[1] / 'data' / 'sed_fitting' / 'zphot' / field_name / (f"{base_det}_{run_type}_bd" if run_type else f"{base_det}_bd")

# === Crossmatch & Input Files ===
crossmatch = Table.read(Path.cwd().parents[1] / 'data' / 'ref_catalogues' / 'all_COSMOS_highz.fits', format='fits')
crossmatch['Redshift'] = crossmatch['Redshift'].astype(float)

input_dir = Path.home().parents[1] / 'hoy' / 'temporaryFilesROHAN' / 'lephare' / 'inputs' / 'euclid'
input_base = base_det.replace('_DR3', '') + (f"_{run_type}" if run_type else '')
input_suffix = '_bd' if bools[0] else '_dusty' if bools[1] else '_lya' if bools[2] else ''
input_name = f"{input_base}{input_suffix}.in"
print('Reading in: ', input_name)
flux_table = Table.read(input_dir / input_name, format='ascii.commented_header')


# === Filters ===
filter_dict = filter_widths()
if not bools[1]:
    for f in ['f444w', 'ch1cds', 'ch2cds']:
        filter_dict.pop(f, None)
if det_list in [['Y', 'J'], ['HSC-Z_DR3']]:
    for f in ['f115w', 'f150w', 'f277w']:
        filter_dict.pop(f, None)
if bools[0]:
    for f in ['HSC-G_DR3', 'HSC-R_DR3']:
        filter_dict.pop(f, None)
if ('no_euclid' in input_name) or (not run_type):
    for f in ['VIS', 'Ye', 'Je', 'He']:
        filter_dict.pop(f, None)
if field_name == 'XMM':
    filter_dict['ch1servs'] = filter_dict.pop('ch1cds', None)
    filter_dict['ch2servs'] = filter_dict.pop('ch2cds', None)

if bools[1]:
    filter_dict.pop('f444w')

n_in = len(filter_dict) * 2
print(filter_dict)
print('Number of input filters:', n_in)

# Define Euclid filters and their colors
euclid_filters = ['VIS', 'Ye', 'Je', 'He']
euclid_colors = ['tab:purple', 'tab:blue', 'tab:green', 'tab:red']
euclid_fwhms = [[0.5413062969825945, 0.9111926674843263],
                [0.9496834246372521, 1.2122804577368327],
                [1.1677419926865804, 1.5668010030527695],
                [1.5215645976366097, 2.020962958853444]]

# Define the VISTA filters and their colors
vista_filters = ['Y', 'J', 'H', 'Ks']
vista_colors = ['tab:orange', 'tab:orange', 'tab:orange', 'tab:orange']
vista_fwhms = [
[0.9740576634297727, 1.0666447029277641],
[1.166227990951164, 1.3387055091259747],
[1.499405612450609, 1.7910553348248728],
[1.991994299052835, 2.3011976928006836]
]

# Multiply the FWHM values by 10000 to convert to Angstroms
euclid_fwhms = [[f * 10000 for f in fwhm] for fwhm in euclid_fwhms]
vista_fwhms = [[f * 10000 for f in fwhm] for fwhm in vista_fwhms]


# Filter FWHMS


# === Parent Catalogue ===
print('Reading in parent catalogue:', generate_selection_name(field_name, filters))
#parent_cat = Table.read(Path.cwd().parents[1] / 'data' / 'catalogues' / generate_selection_name(field_name, filters), format='fits')
parent_cat = Table.read(Path.cwd().parents[1] / 'data' / 'catalogues' / 'COSMOS_5sig_Y_J_nonDet_HSC_G_nonDet_HSC_R_nonDet_HSC_I_5percent_IRACfloor.fits')
if label_crosstalk:
    label_ct(str(Path.cwd().parents[1] / 'data' / 'catalogues' / generate_selection_name(field_name, filters)), fieldName=field_name)

# === Load .spec files ===
spec_files = glob.glob(str(zphot_dir / '*.spec'))

# === Column names for renaming ===
# Column names
names_phot=['phot', 'yerr', 'wlen', 'xerr', 'modelPhot', 'col6', 'col7']
names_sed = ['wlen', 'flux']
names_zpdf = ['z', 'P']
names_param = ['Type', 'Nline', 'Model', 'Library', 'Nband', 'Zphot', 'Zinf', 'Zsup', 'Chi2', 'PDF', 'Extlaw', 'EB-V', 'Lir', 'Age', 'Mass', 'SFR', 'SSFR']

# === Load Sample Catalogue and Sort ===
candidate_path = Path.cwd().parents[1] / 'data' / 'catalogues' / 'candidates'
if field_name == 'COSMOS':
    #sample_cat = Table.read(candidate_path / 'COSMOS_5sig_Y_J_nonDet_HSC_G_nonDet_HSC_R_nonDet_HSC_I_candidates_2025_02_14_with_euclid.fits', format='fits')
    #sample_cat = Table.read(candidate_path / 'COSMOS_5sig_Y_J_nonDet_HSC_G_nonDet_HSC_R_nonDet_HSC_I_best_bd_INTERLOPERS_2024_11_26_with_euclid.fits') #BROWN DWARFS
    sample_cat = Table.read(candidate_path / 'COSMOS_5sig_Y_J_nonDet_HSC_G_nonDet_HSC_R_nonDet_HSC_I_dustyInterlopers_2025_08_14.fits', format='fits') # DUSTY INTERLOPERS
elif field_name == 'XMM':
    sample_cat = Table.read(candidate_path / 'XMM_5sig_HSC_Z_nonDet_HSC_G_nonDet_HSC_R_candidates_2025_05_14.fits', format='fits')

if sort_by_Muv:
    sample_cat.sort('Muv')
    sorted_ids = list(sample_cat['ID'])
    spec_dict = {int(Path(f).stem.replace('Id', '').lstrip('0')): f for f in spec_files}
    spec_files = [spec_dict[ID] for ID in sorted_ids if ID in spec_dict]

if test:
    spec_files = spec_files[:N]

with PdfPages(str(output_dir/output_pdf)) as pdf:
    #! Loop through files
    for i, spec_file in enumerate(spec_files):
        print(f'Object {i+1} of {len(spec_files)}')
        print(spec_file)

        # Read in the .spec file
        file = parse_spec_file(spec_file)

        phot = file.get('phot')
        zpdf = file.get('zpdf')
        sed = file.get('sed')
        params = file.get('model')

        # Rename the column names for each table according to names_{} above
        phot.rename_columns(phot.colnames, names_phot)
        zpdf.rename_columns(zpdf.colnames, names_zpdf)
        params.rename_columns(params.colnames, names_param)
        for table in sed:
            table.rename_columns(table.colnames, names_sed)

        # Get each of the SEDs
        model_wlens = [table['wlen'] for table in sed]
        model_fluxes = [table['flux'] for table in sed]  
        # Get the ID of the object from the file name
        ID = spec_file.split('/')[-1].split('Id')[-1].lstrip('0').split('.spec')[0]
        #print(ID)
        
        if individual_sed:
            if ID != indiv_ID:
                continue
        # testing by only plotting a couple (first is problematic, second is fine)
        # if (ID != '954198') and (ID != '468417'):
        #     continue

        # Index in input file of this ID
        if len(np.where(flux_table['ID'] == int(ID))[0]) == 0:
            raise ValueError('ID not found: check the correct directory is in use.')

        idx = np.where(flux_table['ID'] == int(ID))[0][0]

        # Match ID with the input file to find the object.
        object_flux = flux_table[idx]
        
        # Take values of astropy row and put into array
        object_flux = np.array(object_flux[:])

        #? Print statements to check filters are in correct order.
        # print(object_flux.colnames)
        # exit()

        # Initialize flux and error arrays, to add 329 data.
        flux = []
        error = []

        # In object_FLUX, ignore first and last columns. Get flux and errors, which are in format flux, err, flux, err,...
        for i in range(1, n_in, 2):
            flux.append(object_flux[i])
            error.append(object_flux[i+1])

        flux = np.array(flux)
        error = np.array([np.abs(e) for e in error])

        # Filter centres and widths
        central_wavelengths = [value[0] for value in filter_dict.values()]
        central_wavelengths = [c*10000 for c in central_wavelengths]

        fwhm_values = [value[1] for value in filter_dict.values()]
        fwhm_values = [f*10000 for f in fwhm_values]

        # Get each of the SEDs
        #? Primary solution
        primary_sed = mag_to_flux(model_fluxes[0])
        primary_wlen = model_wlens[0]

        if len(model_fluxes) < 2:
            print('Only one solution')
            continue

        #? Secondary solution
        if np.sum(bools) == 0:
            sec_file = parse_spec_file(dusty_dir / spec_file.split('/')[-1])
            sec_params = sec_file.get('model')
            params.rename_columns(params.colnames, names_param)
            secondary_sed = sec_file.get('sed')
            if secondary_sed != None:
                for table in secondary_sed:
                    table.rename_columns(table.colnames, names_sed)
                # Get each of the SEDs
                model_wlens = [table['wlen'] for table in secondary_sed]
                model_fluxes = [table['flux'] for table in secondary_sed]
                
                if sec_params['Zphot'][1] == -1:
                    secondary_sed = None
                    secondary_wlen = None
                else:
                    secondary_sed = mag_to_flux(model_fluxes[1])
                    secondary_wlen = model_wlens[1]

        else:
            secondary_sed = mag_to_flux(model_fluxes[1])
            secondary_wlen = model_wlens[1]

        #? Stellar solution
        if np.sum(bools) == 0:
            stellar_file = parse_spec_file(bd_dir / spec_file.split('/')[-1])
            stellar_params = stellar_file.get('model')
            #print(stellar_params)
            params.rename_columns(params.colnames, names_param)
            stellar_sed = stellar_file.get('sed')
            for table in stellar_sed:
                table.rename_columns(table.colnames, names_sed)
            # Get each of the SEDs
            model_wlens = [table['wlen'] for table in stellar_sed]
            model_fluxes = [table['flux'] for table in stellar_sed]
            if len(model_fluxes) > 2:
                stellar_sed = mag_to_flux(model_fluxes[2])
                stellar_wlen = model_wlens[2]
            else:
                stellar_sed = mag_to_flux(model_fluxes[1])
                stellar_wlen = model_wlens[1]

        else:
            if len(model_fluxes) > 2:
                stellar_sed = mag_to_flux(model_fluxes[2])
                stellar_wlen = model_wlens[2]
            else:
                stellar_sed = mag_to_flux(model_fluxes[1])
                stellar_wlen = model_wlens[1]

        # Get parameters
        secondary_is_stellar = False

        zphot_1 = round(params['Zphot'][0], 2)
        dz_sup = round(params['Zsup'][0], 2)
        dz_inf = round(params['Zinf'][0], 2)
        chi2_1 = round(params['Chi2'][0], 1)

        print('GAL 1 SOLUTION: ', zphot_1, chi2_1)

        if secondary_sed is None:
            secondary_sed = stellar_sed
            secondary_wlen = stellar_wlen

        #! Read in the secondary solutions from the lephare_out tables
        #? STAR
        if np.sum(bools) == 0:
            #* get from secondary solution
            chi2_star = round(stellar_params['Chi2'][-1], 1)
        else:
            bd_out_name = base_det.replace('_DR3', '') + '_bd.out'
            bd_out = ascii.read(lephare_out_dir / bd_out_name, format='no_header')
            chi2_star = round(bd_out[bd_out['col1'] == int(ID)]['col21'][0], 1)

        #? DSFG
        if np.sum(bools) == 0:
            #* get from secondary solution
            chi2_2 = round(sec_params['Chi2'][1], 1)
            zphot_2 = round(sec_params['Zphot'][1], 2)
        else:
            dusty_out_name = base_det.replace('_DR3', '') + '_dusty.out'
            dusty_out = ascii.read(lephare_out_dir / dusty_out_name, format='no_header')
            chi2_2 = round(dusty_out[dusty_out['col1'] == int(ID)]['col15'][0], 1)
            zphot_2 = round(params['Zphot'][1], 2)

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
        left, bottom, width, height = [0.67, 0.26, 0.2, 0.2]
        ax2 = fig.add_axes([left, bottom, width, height])

        ############! PLOT RESULTS OF FITTING #############
        # Plot best model
        ax1.plot(primary_wlen, primary_sed, color='deepskyblue', alpha=1.0, label= rf'$z = {zphot_1}^{{+{round(dz_sup-zphot_1, 2)}}}_{{-{round(zphot_1-dz_inf, 2)}}}$, ' + r'$\chi^{2} = $' + f'{chi2_1}', linewidth=3, zorder=4)

        # plot secondary galaxy solution and stellar model
        if not secondary_is_stellar:
            ax1.plot(secondary_wlen, secondary_sed, color='darkorange', label=f'z = {zphot_2}, ' + r'$\chi^{2} = $' + f'{chi2_2}', linewidth=2, alpha=0.7) #dusty
            ax1.plot(stellar_wlen, stellar_sed, color='tab:red', label=r'$\chi^{2} = $' + f'{chi2_star}, type = {mlt_type}', linewidth=2, alpha=0.7)            # star

        if secondary_is_stellar:
            ax1.plot(secondary_wlen, secondary_sed, color='tab:red', label=r'$\chi^{2} = $' + f'{chi2_star}, type = {mlt_type}', linewidth=2, alpha=0.7) # only star if no dusty galaxy

        ax2.plot(zpdf['z'], zpdf['P'], color='black', linewidth=2)
        ax2.set_xlim(0, 10)
        ax2.set_xlabel(r'$z_{\mathrm{phot}}$', fontsize=23)
        ax2.set_ylabel(r'$P(z)$', fontsize=23)
        
        # Make tick labels larger
        ax2.tick_params(axis='both', which='major', labelsize=23)

        # ticks on x at 0, 3, 5, 7, 9
        ax2.set_xticks([1, 3, 5, 7, 9])
        ax2.set_yticks([0, 0.5, 1])

        # Show left-side ticks without labels
        ax2.yaxis.set_ticks_position('both')  # Enable ticks on both left and right
        ax2.tick_params(axis='y', which='both', labelright=False)


        ###########! PLOT MODEL AND REAL PHOTOMETRY ##############

        # Model photometry
        model_photometry = mag_to_flux(phot['modelPhot'])
        print(len(model_photometry))
        print(len(central_wavelengths))
        ax1.scatter(central_wavelengths, model_photometry, marker='o', s=100, alpha=0.6, zorder=5, edgecolor='black', facecolor='none', linewidth=2)

        # Plot real photometry, styling Euclid filters separately
        for i, filt in enumerate(filter_dict):
            # Set defaults
            color = 'black'
            marker = 'o'

            # Override for Euclid filters
            if filt in euclid_filters:
                color = euclid_colors[euclid_filters.index(filt)]
                marker = 's'  # square
                markersize = 10
            elif filt in vista_filters:
                color = vista_colors[vista_filters.index(filt)]
                marker = 'D'
                markersize = 10
            else:
                markersize = 8

            # Plot upper limits
            if sigma[i] < 2:
                ax1.scatter(central_wavelengths[i], flux[i] + 2 * error[i], marker='v', color=color, s=100, zorder=6)
            else:
                ax1.errorbar(central_wavelengths[i], flux[i], yerr=error[i], fmt=marker, color=color, 
                             markersize=markersize, zorder=6, elinewidth=3, markeredgecolor='black', ecolor='black')
                
        # Plot the Euclid and VISTA FWHMs as hlines with the correct colours, at 9e-30
        if run_type == 'with_euclid':
            bump = [0, 0, 0.5e-30, 1e-30]
            for i, (fwhm, color) in enumerate(zip(euclid_fwhms, euclid_colors)):
                ax1.hlines(8e-30-bump[i], fwhm[0], fwhm[1], colors=color, linewidth=3.5)

        for i, (fwhm, color) in enumerate(zip(vista_fwhms, vista_colors)):
            ax1.hlines(8.5e-30, fwhm[0], fwhm[1], colors=color, linewidth=3.5)


        ############! CROSSMATCH WITH EXISTING SOURCES ############
        # Get the RA and DEC of the object from the parent catalog
        obj_ra = parent_cat[np.where(parent_cat['ID'] == int(ID))]['RA'][0]
        obj_dec = parent_cat[np.where(parent_cat['ID'] == int(ID))]['DEC'][0]

        if Muv_avail:
            obj_Muv = sample_cat[np.where(sample_cat['ID'] == int(ID))]['Muv'][0]
            print(obj_Muv)

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

            

        #? TITLE STUFF: OLD VERSION USING JUST THE SEXTRACTOR IDS
        #     if Muv_avail:
        #         title_string = f'ID {ID}, ' +  r'$M_{\rm{UV}}=$'+f'{obj_Muv:.2f}, ' + f'{name}'
        #     else:
        #         title_string = f'ID {ID}, z = {redshift}, {name}'
        #     # if ct > 0:
        #     #     title_string += f', POSSIBLE CROSSTALK'
        #     #ax1.set_title(title_string, pad=23)
        #     ax1.set_title(title_string, pad=10)
        # else:
        #     # If no match is found, just set the title with the ID
        #     #title_string = f'ID {ID}'
        #     if Muv_avail:
        #         title_string = f'ID {ID}, ' +  r'$M_{\rm{UV}}=$'+f'{obj_Muv:.2f}'
        #     else:
        #         title_string = f'ID {ID}'
            
        #     # if ct > 0:
        #     #     title_string += f', POSSIBLE CROSSTALK'
        #     #ax1.set_title(title_string, pad=23)
        #     ax1.set_title(title_string, pad=10)

        #? TITLE STUFF: NEW VERSION USING EUCLID NAMING CONVENTIONS
            ID_str = make_euclid_name(obj_ra, obj_dec)
            if Muv_avail:
                title_string = r'$\rm{EUCL}\,$' + ID_str + '\n' + r'$M_{\rm{UV}}=-$' + f'{np.abs(obj_Muv):.2f}, ' + f'{name}'

            else:
                title_string = f'ID {ID_str}, z = {redshift}, {name}'
            # if ct > 0:
            #     title_string += f', POSSIBLE CROSSTALK'
            #ax1.set_title(title_string, pad=23)
            ax1.set_title(title_string, pad=5, fontsize=fontsize)
        else:
            ID_str = make_euclid_name(obj_ra, obj_dec)
            # If no match is found, just set the title with the ID
            #title_string = f'ID {ID}'
            if Muv_avail:
                title_string = r'$\rm{EUCL}\,$' + ID_str + '\n' + r'$M_{\rm{UV}}=-$' + f'{np.abs(obj_Muv):.2f} '
            else:
                title_string = ID_str
            
            # if ct > 0:
            #     title_string += f', POSSIBLE CROSSTALK'
            #ax1.set_title(title_string, pad=23)
            ax1.set_title(title_string, pad=5, fontsize=fontsize)

        # Remove title
        if remove_title:
            ax1.set_title('')

        ax1.set_yscale('log')
        ax1.legend(loc='upper right', fontsize=17)
        ax1.set_ylim(3e-32, 1e-29)
        ax1.set_xlim(3000, 40000)
        if field_name == 'XMM' or field_name == 'CDFS':
            ax1.set_ylim(5e-32, 2e-29)
            ax1.set_xlim(3000, 20000)

        # Convert the x axis into microns by dividing by 10000
        ax1.set_xticks([5000, 10000, 15000, 20000, 25000, 30000, 35000])
        ax1.set_xticklabels([0.5, 1, 1.5, 2, 2.5, 3, 3.5])

        # Set custom y-ticks at powers of 10, and label them as log10(f_nu)
        yticks = [1e-31, 1e-30, 1e-29]
        ax1.set_yticks(yticks)
        ax1.set_yticklabels([r"$-31$", r"$-30$", r"$-29$"])

        ax1.set_ylabel(r'$\log_{10}(f_{\nu}\,/\,\rm{erg}\,\rm{s}^{-1}\,\rm{cm}^{-2}\,\rm{Hz}^{-1})$', fontsize=fontsize) #[erg s$^{-1}$ cm$^{-2}$ Hz$^{-1}$]
        ax1.set_xlabel(r'$\lambda \, [\rm{\mu m}]$', fontsize=fontsize)

        # Add magnitude axis
        secax = ax1.secondary_yaxis('right', functions=(flux_to_mag, mag_to_flux))

        mags = [30, 29, 28, 27, 26, 25, 24, 23, 22]
        secax.set_yticks(mags)

        secax.yaxis.set_major_formatter(mticker.ScalarFormatter())

        ax1.tick_params(which='major', length=7, width=4)
        ax1.tick_params(axis='both', which='minor', length=4, width=3)

        secax.tick_params(which='major', length=7, width=4)
        secax.tick_params(axis='both', which='minor', length=5, width=3)

        secax.set_ylabel(r'$\mathrm{m_{AB}}$', fontsize=fontsize)

        # Increas marker label sizes
        ax1.tick_params(axis='both', which='major', labelsize=fontsize)
        secax.tick_params(axis='both', which='major', labelsize=fontsize)

        ax1.minorticks_on()


        pdf.savefig(fig, bbox_inches='tight')

        # Also save fig as its own pdf
        if save_indiv:
            plt.savefig(str(sed_dir / f'{ID}_SED_EUCLname.pdf'), bbox_inches='tight')

        plt.close(fig)

        ############! CUTOUT ##############
        if individual_sed == False:
            # Get ra, dec of object from parent catalogue
            obj = parent_cat[np.where(parent_cat['ID'] == int(ID))]
            ra = obj['RA'][0]
            dec = obj['DEC'][0]

            contained_in = isCoordInSurveyFootprints(ra, dec)
            if field_name == 'COSMOS':
                print(contained_in)

            # Use smaller cutout size for visual inspection
            cutout_size = 6. if det_list == ['Y', 'J'] else 10. # arcsec

            # If det_list is ['Y', 'J'] and euclid_blind is True, use VistaCutout
            if det_list == ['Y', 'J'] and euclid_blind or (det_list == ['Y', 'J'] and contained_in[0][0] == '0') or (det_list == ['HSC-Z_DR3'] and euclid_blind):
                print('Doing Vista cutout')
                cutout_fig, cutout_axs = VistaCutout(ra, dec, size=cutout_size, save_cutout=False)
                #cutout_fig, cutout_axs = AllCutout(ra, dec, size=cutout_size, save_cutout=False)
            if field_name == 'XMM':
                cutout_fig, cutout_axs = XMMCutout(ra, dec, size=cutout_size, save_cutout=False)
            if not euclid_blind:
                print('Doing all cutouts')
                cutout_fig, cutout_axs = AllCutout(ra, dec, size=6., save_cutout=False)
                #cutout_fig, cutout_axs = Cutout(ra, dec, size=6., save_cutout=False)

            pdf.savefig(cutout_fig)
            plt.close(cutout_fig)


