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
import matplotlib.ticker as mticker

# Backend for running on queue
import matplotlib as mpl
mpl.use('Agg')

cutout_path = Path.cwd().parents[0] / 'cutouts'
sys.path.append(str(cutout_path))
from cutout_codes import *
from vista_cutouts import *
from all_cutouts import AllCutout
from xmm_cutouts import XMMCutout

crosstalk_path = Path.cwd().parents[3] / 'HSC_SSP_DR3' / 'codes'
sys.path.append(str(crosstalk_path))
from full_catalogue_codes import label_ct

sed_path = Path.cwd().parents[0] / 'sed_fitting'
sys.path.append(str(sed_path))
from sed_fitting_codes import parse_spec_file

plt.rcParams['axes.linewidth'] = 2.5
plt.rcParams.update({'font.size': 15})
plt.rcParams['figure.dpi'] = 100

#! field name
field_name = 'XMM'

#! If I want to plot individual things, overwrite the standard output pdf
individual_sed = False
indiv_ID = '661703'
indiv_pdf_name = 'FAINT_BD_SED.pdf'

#! TEST by plotting first N objects
test = False
N = 10

#! Plot final samples sorted by Muv
sort_by_Muv = False

if len(sys.argv) > 1:
    filters_json = sys.argv[1]
    filters = json.loads(filters_json)
    bools_json = sys.argv[2]
    bools = json.loads(bools_json)
    all_filters_json = sys.argv[3]
    all_filters = json.loads(all_filters_json)
    run_type_json = sys.argv[4]
    run_type = json.loads(run_type_json)
    object_type_json = sys.argv[5]
    object_type = json.loads(object_type_json)
    field_name_json = sys.argv[6]
    field_name = json.loads(field_name_json)

if run_type == '':
    euclid_blind = True
else:
    euclid_blind = False

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
        -1: 'No fit'
    }
    return stellar_dict[model]

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

#! Label crosstalk?
label_crosstalk = False

#! Output PDF file
output_dir = Path.cwd().parents[1] / 'plots' / 'seds'

#! Output dir for individual SEDs
sed_dir = Path.cwd().parents[1] / 'plots' / 'seds' / 'individual'

#! Set up the output PDF name from the det and nondet filters
# Set up the detection filter base
det_list = [f for f, t in filters.items() if t['type'] == 'detection']
if len(det_list) != 0:
    base_det = 'det_' + '_'.join(det_list)
else:
    stack_list = [f for f, t in filters.items() if t['type'] == 'stacked-detection']
    stack_filters = stack_list[0].split('+')
    base_det = 'det_' + '_'.join(stack_filters)
    det_list = stack_filters

# Set up the non-detection filter base
base_nondet = 'nonDet_' + '_'.join([f for f, t in filters.items() if t['type'] == 'non-detection'])

# Collect the components for the filename
filename_components = [field_name, base_det, base_nondet]

# Add object type if available
if object_type:
    filename_components.append(object_type)

# Add run_type if available
if run_type != '':
    filename_components.append(run_type)

# Construct the final output PDF name
output_pdf = '_'.join(filename_components) + '.pdf'
if individual_sed:
    output_pdf = indiv_pdf_name

print('Saving to: ', output_pdf)

#! Set up directories
# zphot_dir = Path.cwd().parents[1] / 'data' / 'sed_fitting' / 'zphot' / 'best_fits' / base_det
zphot_folder = base_det + '_' + run_type + f'_{object_type}' if run_type != '' else base_det + f'_{object_type}'

zphot_dir = Path.cwd().parents[1] / 'data' / 'sed_fitting' / 'zphot' / field_name / 'best_fits' / zphot_folder
print('Taking SEDs from: ', zphot_dir)

# Set up the dusty and brown dwarf dir to get the same correspoinding SEDs
dusty_folder = base_det + '_' + run_type + '_dusty' if run_type != '' else base_det + '_dusty'
dusty_dir = Path.cwd().parents[1] / 'data' / 'sed_fitting' / 'zphot' / field_name / dusty_folder

bd_folder = base_det + '_' + run_type + '_bd' if run_type != '' else base_det + '_bd'
bd_dir = Path.cwd().parents[1] / 'data' / 'sed_fitting' / 'zphot' / field_name / bd_folder

# Crossmatched catalogue name to get existing sources
crossmatch_name = 'all_COSMOS_highz.fits'

# Load the crossmatched catalogue of known objects
crossmatch_dir = Path.cwd().parents[1] / 'data' / 'ref_catalogues'
crossmatch = Table.read(crossmatch_dir / crossmatch_name, format='fits')
crossmatch['Redshift'] = crossmatch['Redshift'].astype(float)

# Read in the input .in file
input_dir = Path.home().parents[1] / 'hoy' / 'temporaryFilesROHAN' / 'lephare' / 'inputs' / 'euclid'
if not bools[0]:    
    input_name = base_det.replace('_DR3', '') + (f"_{run_type}" if run_type != '' else '') + '.in'
if bools[0]:
    input_name = base_det.replace('_DR3', '') + (f"_{run_type}" if run_type != '' else '') + '_bd' + '.in'
if bools[1]:
    input_name = base_det.replace('_DR3', '') + (f"_{run_type}" if run_type != '' else '') + '_dusty' + '.in'
if bools[2]:
    input_name = base_det.replace('_DR3', '') + (f"_{run_type}" if run_type != '' else '') + '_lya' + '.in'

print('Reading in: ', input_name)

flux_table = Table.read(input_dir / input_name, format='ascii.commented_header')

# Directory to read in the .out file
lephare_out_dir = Path.home() / 'lephare' / 'lephare_dev' / 'test'

# Column names
names_phot=['phot', 'yerr', 'wlen', 'xerr', 'modelPhot', 'col6', 'col7']
names_sed = ['wlen', 'flux']
names_zpdf = ['z', 'P']
names_param = ['Type', 'Nline', 'Model', 'Library', 'Nband', 'Zphot', 'Zinf', 'Zsup', 'Chi2', 'PDF', 'Extlaw', 'EB-V', 'Lir', 'Age', 'Mass', 'SFR', 'SSFR']

# Get filters
filter_dict = filter_widths()

# Remove f444w
if not bools[1]:
    filter_dict.pop('f444w')
    filter_dict.pop('ch1cds')
    filter_dict.pop('ch2cds')
    # filter_dict.pop('ch1servs')
    # filter_dict.pop('ch2servs')
    #print('Removed f444w, ch1cds, ch2cds')

if det_list == ['Y', 'J'] or det_list == ['HSC-Z_DR3']:
    # filter_dict.pop('VIS')
    # filter_dict.pop('Ye')
    # filter_dict.pop('Je')
    # filter_dict.pop('He')
    # filter_dict.pop('f115w')
    # filter_dict.pop('f150w')
    # filter_dict.pop('f277w')
    # filter_dict.pop('f444w')'
    print('f115w', 'f150w', 'f277w', 'f444w')

if bools[0]:
    filter_dict.pop('HSC-G_DR3')
    filter_dict.pop('HSC-R_DR3')
    #filter_dict.pop('VIS')
    #print('Removed f277w, HSC-G_DR3, HSC-R_DR3')

# If running only VISTA: Remove items with keys VIS, Ye, Je, He
if ('no_euclid' in input_name) or (run_type == ''):
    filter_dict.pop('VIS')
    filter_dict.pop('Ye')
    filter_dict.pop('Je')
    filter_dict.pop('He')

if field_name == 'XMM':
    # Rename ch1cds and ch2cds to ch1servs and ch2servs
    filter_dict['ch1servs'] = filter_dict.pop('ch1cds')
    filter_dict['ch2servs'] = filter_dict.pop('ch2cds')

# Get n_in from the length of filter_dict
n_in = len(filter_dict) * 2

#print('Filters: ', filter_dict)

# Load the parent catalogue to get RA,DEC
# Generate the name of the parent catalogue
parent_cat_dir = Path.cwd().parents[1] / 'data' / 'catalogues'
parent_cat_name = generate_selection_name(field_name, filters)

# Label crosstalk
if label_crosstalk:
    label_ct(str(parent_cat_dir / parent_cat_name), fieldName=field_name)

parent_cat = Table.read(parent_cat_dir / parent_cat_name, format='fits')


#! Get list of .spec files
spec_files = glob.glob(str(zphot_dir / '*.spec'))

#! Catalogue of final sample with Muv
#sample_cat = Table.read(parent_cat_dir / 'candidates' / 'COSMOS_5sig_Y_J_nonDet_HSC_G_nonDet_HSC_R_nonDet_HSC_I_candidates_2025_02_14_with_euclid.fits', format='fits')
sample_cat = Table.read(parent_cat_dir / 'candidates' / 'XMM_5sig_HSC_Z_nonDet_HSC_G_nonDet_HSC_R_candidates_2025_05_14.fits', format='fits')

if sort_by_Muv:
    sample_cat.sort('Muv')
    sorted_ids = list(sample_cat['ID'])

    # Build a dictionary from ID -> filepath (extracting ID from filename)
    spec_dict = {}
    for f in spec_files:
        fname = Path(f).name 
        id_str = fname.replace('Id', '').replace('.spec', '')  
        id_int = int(id_str.lstrip('0'))  
        spec_dict[id_int] = f

        spec_files = [spec_dict[ID] for ID in sorted_ids if ID in spec_dict]

if test:
    # Limit the number of files to process
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

        # Initialize flux and error arrays, to add .in data.
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
                stellar_sed = None
                stellar_wlen = None

        else:
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
        left, bottom, width, height = [0.70, 0.23, 0.18, 0.18]
        ax2 = fig.add_axes([left, bottom, width, height])

        ############! PLOT RESULTS OF FITTING #############
        # Plot best model
        ax1.plot(primary_wlen, primary_sed, color='deepskyblue', alpha=1.0, label= rf'$z = {zphot_1}^{{+{round(dz_sup-zphot_1, 2)}}}_{{-{round(zphot_1-dz_inf, 2)}}}$, ' + r'$\chi^{2} = $' + f'{chi2_1}', linewidth=3, zorder=4)

        # plot secondary galaxy solution and stellar model
        if not secondary_is_stellar:
            ax1.plot(secondary_wlen, secondary_sed, color='darkorange', label=f'z = {zphot_2}, ' + r'$\chi^{2} = $' + f'{chi2_2}', linewidth=2, alpha=0.7) #dusty
            ax1.plot(stellar_wlen, stellar_sed, color='red', label=r'$\chi^{2} = $' + f'{chi2_star}, type = {mlt_type}', linewidth=2, alpha=0.7)            # star

        if secondary_is_stellar:
            ax1.plot(secondary_wlen, secondary_sed, color='red', label=r'$\chi^{2} = $' + f'{chi2_star}, type = {mlt_type}', linewidth=2, alpha=0.7) # only star if no dusty galaxy

        ax2.plot(zpdf['z'], zpdf['P'], color='black', linewidth=2)
        ax2.set_xlim(0, 10)
        ax2.set_xlabel(r'$z_{\mathrm{phot}}$', fontsize=20)
        ax2.set_ylabel(r'$P(z)$', fontsize=20)
        
        # Make tick labels larger
        ax2.tick_params(axis='both', which='major', labelsize=20)

        ###########! PLOT MODEL AND REAL PHOTOMETRY ##############

        # Model photometry
        model_photometry = mag_to_flux(phot['modelPhot'])
        print(len(model_photometry))
        print(len(central_wavelengths))
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

        #obj_Muv = sample_cat[np.where(sample_cat['ID'] == int(ID))]['Muv'][0]
        #print(obj_Muv)

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
            #title_string = f'ID {ID}, z = {redshift}, {name}'
            #title_string = f'ID {ID}, ' +  r'$M_{\rm{UV}}=$'+f'{obj_Muv:.2f}, ' + f'{name}'
            # if ct > 0:
            #     title_string += f', POSSIBLE CROSSTALK'
            #ax1.set_title(title_string, pad=23)
            ax1.set_title(title_string, pad=10)
        else:
            # If no match is found, just set the title with the ID
            title_string = f'ID {ID}'
            #title_string = f'ID {ID}, ' +  r'$M_{\rm{UV}}=$'+f'{obj_Muv:.2f}'
            # if ct > 0:
            #     title_string += f', POSSIBLE CROSSTALK'
            #ax1.set_title(title_string, pad=23)
            ax1.set_title(title_string, pad=10)



        ax1.set_yscale('log')
        ax1.legend(loc='upper right', fontsize=15)
        ax1.set_ylim(3e-32, 1e-29)
        ax1.set_xlim(3000, 40000)
        if field_name == 'XMM' or field_name == 'CDFS':
            ax1.set_ylim(5e-32, 2e-29)
            ax1.set_xlim(3000, 20000)

        fontsize = 21

        # Convert the x axis into microns by dividing by 10000
        ax1.set_xticks([5000, 10000, 15000, 20000, 25000, 30000, 35000])
        ax1.set_xticklabels([0.5, 1, 1.5, 2, 2.5, 3, 3.5])

        ax1.set_ylabel(r'$f_{\nu}$ (erg s$^{-1}$ cm$^{-2}$ Hz$^{-1}$)', fontsize=fontsize)
        ax1.set_xlabel(r'$\lambda (\mu m)$', fontsize=fontsize)

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


        pdf.savefig(fig, bbox_inches='tight')

        # Also save fig as its own pdf
        plt.savefig(str(sed_dir / f'{ID}_SED.pdf'), bbox_inches='tight')

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
            if det_list == ['Y', 'J'] and euclid_blind or (det_list == ['Y', 'J'] and contained_in[0][0] == '0'):
                cutout_fig, cutout_axs = VistaCutout(ra, dec, size=cutout_size, save_cutout=False)
                #cutout_fig, cutout_axs = AllCutout(ra, dec, size=cutout_size, save_cutout=False)
            if field_name == 'XMM':
                cutout_fig, cutout_axs = XMMCutout(ra, dec, size=cutout_size, save_cutout=False)
            else:
                print('Doing all cutouts')
                cutout_fig, cutout_axs = AllCutout(ra, dec, size=6., save_cutout=False)
                #cutout_fig, cutout_axs = Cutout(ra, dec, size=6., save_cutout=False)

            pdf.savefig(cutout_fig)
            plt.close(cutout_fig)


