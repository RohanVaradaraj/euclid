"""
Measure UV mags for galaxy candidates.

Created: Friday 8th November 2024.
"""

from astropy.table import Table, Column
from pathlib import Path
import numpy as np
import glob
import datetime

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

# Save the new catalogue
t_candidates.write(cat_dir / 'candidates' / new_cat_name, overwrite=True)
 