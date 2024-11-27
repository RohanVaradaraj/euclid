"""
Take the folder of objects that made it through the SED fitting.

Make a new astropy catalogue of these. Can then append other properties such as Muv, beta slope, stellar mass, etc.

Created: Friday 8th November 2024.
"""

from astropy.table import Table, Column
from pathlib import Path
import numpy as np
import glob
import datetime
import sys

sed_path = Path.cwd().parents[0] / 'sed_fitting'
sys.path.append(str(sed_path))
from sed_fitting_codes import parse_spec_file, LymanAlphaModel

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

#! Run type - governing the filter set used.
run_type = 'with_euclid'

#! Run flag - options are 'z7', 'BD'/'best_bd', 'dustyInterlopers'
run_flag = 'really_good_LBGs'

# Only get lya if we are looking at the LBG sample
run_lya = (run_flag == 'z7')

# Name of the directory we want to use to make the catalogue
folder = f'det_Y_J_{run_type}_{run_flag}' if run_type != '' else f'det_Y_J_{run_flag}'
lya_folder = f'det_Y_J_{run_type}_lya' if run_type != '' else 'det_Y_J_lya'

# Parent catalogue from which to get fluxes
cat_name = 'COSMOS_5sig_Y_J_nonDet_HSC_G_nonDet_HSC_R_nonDet_HSC_I.fits'

today_date = datetime.datetime.now().strftime('%Y_%m_%d')

# Name of the new catalogue
if run_flag == 'z7':
    new_cat_name = f'COSMOS_5sig_Y_J_nonDet_HSC_G_nonDet_HSC_R_nonDet_HSC_I_candidates_{today_date}_{run_type}.fits' if run_type != '' else f'COSMOS_5sig_Y_J_nonDet_HSC_G_nonDet_HSC_R_nonDet_HSC_I_candidates_{today_date}.fits'
else:
    new_cat_name = f'COSMOS_5sig_Y_J_nonDet_HSC_G_nonDet_HSC_R_nonDet_HSC_I_{run_flag}_INTERLOPERS_{today_date}_{run_type}.fits' if run_type != '' else f'COSMOS_5sig_Y_J_nonDet_HSC_G_nonDet_HSC_R_nonDet_HSC_I_{run_flag}_INTERLOPERS_{today_date}.fits'

# Read in the parent catalogue
cat_dir = Path.cwd().parents[1] / 'data' / 'catalogues'
t = Table.read(cat_dir / cat_name)

# Get the list of objects that made it through the SED fitting
obj_dir = Path.cwd().parents[1] / 'data' / 'sed_fitting' / 'zphot' / 'best_fits'
obj_list = glob.glob(str(obj_dir / folder / '*.spec'))

print(obj_dir / folder)

# And lyman-alpha emitters
if run_lya:
    lya_obj_list = glob.glob(str(obj_dir / lya_folder / '*.spec'))
    lya_IDs = [spec_file.split('/')[-1].split('Id')[-1].lstrip('0').split('.spec')[0] for spec_file in lya_obj_list]
    lya_IDs = [int(ID) for ID in lya_IDs]

# Get the IDs
IDs = [spec_file.split('/')[-1].split('Id')[-1].lstrip('0').split('.spec')[0] for spec_file in obj_list]
IDs = [int(ID) for ID in IDs]


# Take the objects from the parent catalogue with these IDs
mask = np.isin(t['ID'], IDs)
if run_lya:
    mask_lya = np.isin(t['ID'], lya_IDs)

if run_lya:
    t_candidates = t[mask | mask_lya]
else:
    t_candidates = t[mask]

# Add new columns to the table for primary solution
t_candidates['Zphot'] = Column(np.zeros(len(t_candidates)))
t_candidates['Zinf'] = Column(np.zeros(len(t_candidates)))
t_candidates['Zsup'] = Column(np.zeros(len(t_candidates)))
t_candidates['Chi2'] = Column(np.zeros(len(t_candidates)))

# And for secondary solution
t_candidates['Zphot_sec'] = Column(np.zeros(len(t_candidates)))
t_candidates['Chi2_sec'] = Column(np.zeros(len(t_candidates)))

# And for stellar solution
t_candidates['Stellar_model'] = Column(np.zeros(len(t_candidates), dtype='U2'))
t_candidates['Chi2_star'] = Column(np.zeros(len(t_candidates)))

# And for the equivalent width of the Lyman-alpha line galaxies, only if we are running Lya
if run_lya:
    t_candidates['Lyman_alpha_EW'] = Column(np.zeros(len(t_candidates)))


# Loop through the objects in the table and get its SED properties
for i, ID in enumerate(IDs):

    # Find the row index in the table with this ID
    row_index = np.where(t_candidates['ID'] == ID)[0]

    # Open the SED solution file
    file_name = obj_list[i]
    spec_data = parse_spec_file(file_name)

    # Get model table
    model_table = spec_data.get('model')

    # Get the properties of the primary solution
    t_candidates['Zphot'][row_index] = model_table[model_table['ID'] == 'GAL-1']['Zphot'][0]
    t_candidates['Zinf'][row_index] = model_table[model_table['ID'] == 'GAL-1']['Zinf'][0]
    t_candidates['Zsup'][row_index] = model_table[model_table['ID'] == 'GAL-1']['Zsup'][0]
    t_candidates['Chi2'][row_index] = model_table[model_table['ID'] == 'GAL-1']['Chi2'][0]

    # Get the properties of the secondary solution
    t_candidates['Zphot_sec'][row_index] = model_table[model_table['ID'] == 'GAL-2']['Zphot'][0]
    t_candidates['Chi2_sec'][row_index] = model_table[model_table['ID'] == 'GAL-2']['Chi2'][0]

    # Get the properties of the stellar solution
    BD_model = model_table[model_table['ID'] == 'STAR']['Model'][0]
    BD_type = stellar_type(int(BD_model))
    t_candidates['Stellar_model'][row_index] = BD_type
    t_candidates['Chi2_star'][row_index]= model_table[model_table['ID'] == 'STAR']['Chi2'][0]

#! Now loop through Lya objects
if run_lya:
    for i, lya_ID in enumerate(lya_IDs):

        row_index = np.where(t_candidates['ID'] == lya_ID)[0]

        # Open the SED solution file
        file_name = lya_obj_list[i]
        spec_data = parse_spec_file(file_name)

        # Get model table
        model_table = spec_data.get('model')

        model_number = model_table[model_table['ID'] == 'GAL-1']['Model'][0]

        lya_EW = LymanAlphaModel(model_number)

        t_candidates['Lyman_alpha_EW'][row_index] = lya_EW

        # And all the other table params
        # Get the properties of the primary solution
        t_candidates['Zphot'][row_index] = model_table[model_table['ID'] == 'GAL-1']['Zphot'][0]
        t_candidates['Zinf'][row_index] = model_table[model_table['ID'] == 'GAL-1']['Zinf'][0]
        t_candidates['Zsup'][row_index] = model_table[model_table['ID'] == 'GAL-1']['Zsup'][0]
        t_candidates['Chi2'][row_index] = model_table[model_table['ID'] == 'GAL-1']['Chi2'][0]

        # Get the properties of the secondary solution
        t_candidates['Zphot_sec'][row_index] = model_table[model_table['ID'] == 'GAL-2']['Zphot'][0]
        t_candidates['Chi2_sec'][row_index] = model_table[model_table['ID'] == 'GAL-2']['Chi2'][0]

        # Get the properties of the stellar solution
        BD_model = model_table[model_table['ID'] == 'STAR']['Model'][0]
        BD_type = stellar_type(int(BD_model))
        t_candidates['Stellar_model'][row_index] = BD_type
        t_candidates['Chi2_star'][row_index]= model_table[model_table['ID'] == 'STAR']['Chi2'][0]

print(t_candidates)

# Save the new catalogue
t_candidates.write(cat_dir / 'candidates' / new_cat_name, overwrite=True)

 