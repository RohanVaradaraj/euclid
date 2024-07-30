"""
make tables to do the comparison of SED fits with and without Euclid filters.

Created: Monday 29th July 2024.
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

# Define a function which maps the stellar model to its stellar type
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

#! Output directory for comparison catalogue
save_dir = Path.cwd().parents[1] / 'data' / 'sed_fitting' / 'catalogues' / 'comparison'

# Set up directories
zphot1_dir = Path.cwd().parents[1] / 'data' / 'sed_fitting' / 'zphot' / 'det_Ye_y_LBG_no_euclid'
zphot2_dir = Path.cwd().parents[1] / 'data' / 'sed_fitting' / 'zphot' / 'det_Ye_y_LBG'

# Read in the input .in file
input_dir = Path.home().parents[1] / 'hoy' / 'temporaryFilesROHAN' / 'lephare' / 'inputs' / 'euclid'

input1_name = 'det_Ye_Y_no_euclid.in'
input2_name = 'det_Ye_Y.in'

# Save name based on input 2
save_name = input2_name.split('.in')[0] + '_comparison.fits'

flux1_table = Table.read(input_dir / input1_name, format='ascii.commented_header')
flux2_table = Table.read(input_dir / input2_name, format='ascii.commented_header')


# Define start and end points of each section in the LePhare .spec file
#! Without euclid
ds_phot_1 = 9
de_phot_1 = 24

ds_mod_1 = 205
de_mod_1 = -1

ds_pz_1 = 24
de_pz_1 = 205

ds_param_1 = 3
de_param_1 = 9

#! With euclid
ds_phot_2 = 9
de_phot_2 = 28

ds_mod_2 = 209
de_mod_2 = -1

ds_pz_2 = 28
de_pz_2 = 209

ds_param_2 = 3
de_param_2 = 9



# From these values, we can get the number of filters which will be useful for extracting fluxes from the .in file.
n_bands_1 = de_phot_1 - ds_phot_1
n_bands_2 = de_phot_2 - ds_phot_2

# Double the value to correspond to each flux and error
n_in_1 = 2 * n_bands_1
n_in_2 = 2 * n_bands_2

# Column names
names_phot=['phot', 'yerr', 'wlen', 'xerr', 'modelPhot', 'col6', 'col7']
names_sed = ['wlen', 'flux']
names_zpdf = ['z', 'P']
names_param = ['Type', 'Nline', 'Model', 'Library', 'Nband', 'Zphot', 'Zinf', 'Zsup', 'Chi2', 'PDF', 'Extlaw', 'EB-V', 'Lir', 'Age', 'Mass', 'SFR', 'SSFR']

# Get filters
filter_dict_1 = filter_widths()
filter_dict_2 = filter_widths()

# If running only VIISTA: Remove items with keys VIS, Ye, Je, He
filter_dict_1.pop('VIS')
filter_dict_1.pop('Ye')
filter_dict_1.pop('Je')
filter_dict_1.pop('He')

# Load the crossmatched catalogue of known objects
crossmatch_dir = Path.cwd().parents[1] / 'data' / 'catalogues'
crossmatch_name = 'XMATCH_COSMOS_5sig_Ye_2sig_HSC_Y_nonDet_HSC_G_nonDet_HSC_R_nonDet_HSC_I.fits'
crossmatch = Table.read(crossmatch_dir / crossmatch_name, format='fits')
crossmatch['Redshift'] = crossmatch['Redshift'].astype(float)

# Load the parent catalogue to get RA,DEC
parent_cat = Table.read(crossmatch_dir / 'COSMOS_5sig_Ye_2sig_Y_nonDet_HSC_G_nonDet_HSC_R_nonDet_HSC_I.fits', format='fits')

# Collect all .spec files
spec_files_1 = glob.glob(str(zphot1_dir / '*.spec'))
spec_files_2 = glob.glob(str(zphot2_dir / '*.spec'))

#? Make empty lists to get results of fitting

#! id
IDs = []

#! High redshift solution
zphot_highz_1 = []
zphot_highz_2 = []

#! Errors
highz_sup_1 = []
highz_sup_2 = []
highz_inf_1 = []
highz_inf_2 = []

#! Low redshift solution
zphot_lowz_1 = []
zphot_lowz_2 = []

#! Stellar model
stellar_1 = []
stellar_2 = []

stellar_type_1 = []
stellar_type_2 = []

#! Chi-squared values
chi2_highz_1 = []
chi2_highz_2 = []

chi2_lowz_1 = []
chi2_lowz_2 = []

chi2_star_1 = []
chi2_star_2 = []

#! Loop through files
for i, spec_file in enumerate(spec_files_1):
    print(f'Object {i+1} of {len(spec_files_1)}')

    #! without euclid
    phot_1 = ascii.read(spec_file, format='basic', data_start=ds_phot_1, data_end=de_phot_1, delimiter=' ', names=names_phot)
    zpdf_1 = ascii.read(spec_file, format='basic', data_start=ds_pz_1, data_end=de_pz_1, delimiter=' ', names=names_zpdf)
    sed_1 = ascii.read(spec_file, format='basic', data_start=ds_mod_1, data_end=de_mod_1, delimiter=' ', names=names_sed)
    params_1 = ascii.read(spec_file, format='basic', data_start=ds_param_1, data_end=de_param_1, delimiter=' ', names=names_param)

    # High redshift solution
    zphot_highz_1.append(params_1['Zphot'][0])
    highz_sup_1.append(params_1['Zsup'][0])
    highz_inf_1.append(params_1['Zinf'][0])

    # Low redshift solution
    zphot_lowz_1.append(params_1['Zphot'][1])

    # Chi2 values
    chi2_highz_1.append(params_1['Chi2'][0])
    chi2_lowz_1.append(params_1['Chi2'][1])
    chi2_star_1.append(params_1['Chi2'][-1])

    # Stellar model
    stellar_1.append(params_1['Model'][-1])
    stellar_type_1.append(stellar_type(params_1['Model'][-1]))

    # Get the file name
    file_name = spec_file.split('/')[-1]

    #! with euclid
    phot_2 = ascii.read(str(zphot2_dir/file_name), format='basic', data_start=ds_phot_2, data_end=de_phot_2, delimiter=' ', names=names_phot)
    zpdf_2 = ascii.read(str(zphot2_dir/file_name), format='basic', data_start=ds_pz_2, data_end=de_pz_2, delimiter=' ', names=names_zpdf)
    sed_2 = ascii.read(str(zphot2_dir/file_name), format='basic', data_start=ds_mod_2, data_end=de_mod_2, delimiter=' ', names=names_sed)
    params_2 = ascii.read(str(zphot2_dir/file_name), format='basic', data_start=ds_param_2, data_end=de_param_2, delimiter=' ', names=names_param)

    # Get the ID of the object from the file name
    ID_1 = spec_file.split('/')[-1].split('Id')[-1].lstrip('0').split('.spec')[0]
    ID_2 = file_name.split('Id')[-1].lstrip('0').split('.spec')[0]

    if ID_1 == ID_2:
        IDs.append(ID_1)
    else:
        raise ValueError('IDs do not match in the two files.')

    # High redshift solution
    zphot_highz_2.append(params_2['Zphot'][0])
    highz_sup_2.append(params_2['Zsup'][0])
    highz_inf_2.append(params_2['Zinf'][0])

    # Low redshift solution
    zphot_lowz_2.append(params_2['Zphot'][1])

    # Chi2 values
    chi2_highz_2.append(params_2['Chi2'][0])
    chi2_lowz_2.append(params_2['Chi2'][1])
    chi2_star_2.append(params_2['Chi2'][-1])

    # Stellar model
    stellar_2.append(params_2['Model'][-1])
    stellar_type_2.append(stellar_type(params_2['Model'][-1]))

# Create astropy table for the fitting results
# First name the columns
names = ('ID', 'zphot_pri_1', 'dz_sup_1', 'dz_inf_1', 'zphot_sec_1', 'chi2_pri_1', 'chi2_sec_1', 'chi2_star_1', 'star_model_1', 'stellar_type_1', 'zphot_pri_2', 'dz_sup_2', 'dz_inf_2', 'zphot_sec_2', 'chi2_pri_2', 'chi2_sec_2', 'chi2_star_2', 'star_model_2', 'stellar_type_2')
columns = [IDs, zphot_highz_1, highz_sup_1, highz_inf_1, zphot_lowz_1, chi2_highz_1, chi2_lowz_1, chi2_star_1, stellar_1, stellar_type_1, zphot_highz_2, highz_sup_2, highz_inf_2, zphot_lowz_2, chi2_highz_2, chi2_lowz_2, chi2_star_2, stellar_2, stellar_type_2]

t = Table(columns, names=names)
t.write(save_dir / save_name, format='fits', overwrite=True)
