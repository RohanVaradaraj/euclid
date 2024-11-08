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
from sed_fitting_codes import parse_spec_file


# Name of the directory we want to use to make the catalogue
folder = 'det_Y_J_z7'

# Parent catalogue from which to get fluxes
cat_name = 'COSMOS_5sig_Y_J_nonDet_HSC_G_nonDet_HSC_R_nonDet_HSC_I.fits'

today_date = datetime.datetime.now().strftime('%Y_%m_%d')

# Name of the new catalogue
new_cat_name = f'COSMOS_5sig_Y_J_nonDet_HSC_G_nonDet_HSC_R_nonDet_HSC_I_candidates_{today_date}.fits'

# Read in the parent catalogue
cat_dir = Path.cwd().parents[1] / 'data' / 'catalogues'
t = Table.read(cat_dir / cat_name)

# Get the list of objects that made it through the SED fitting
obj_dir = Path.cwd().parents[1] / 'data' / 'sed_fitting' / 'zphot' / 'best_fits'
obj_list = glob.glob(str(obj_dir / folder / '*.spec'))

# Get the IDs
IDs = [spec_file.split('/')[-1].split('Id')[-1].lstrip('0').split('.spec')[0] for spec_file in obj_list]

# Take the objects from the parent catalogue with these IDs
mask = np.isin(t['ID'], IDs)

t_candidates = t[mask]

file_name = obj_list[0]
parsed_data = parse_spec_file(file_name)

print(parsed_data.keys())

phot_table = parsed_data.get('phot')
model_table = parsed_data.get('model')
zpdf_table = parsed_data.get('zpdf')
sed_table = parsed_data.get('sed')

print(phot_table)
print(model_table)
print(zpdf_table)
print(sed_table)


# Now loop through the objects and get the properties of the SED solutions.

# Save the new catalogue
t_candidates.write(cat_dir / 'candidates' / new_cat_name, overwrite=True)

 