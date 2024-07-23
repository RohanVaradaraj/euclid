"""
Run the SED fitting.

Created: Friday 12th July 2024.
"""

from sed_fitting_codes import *
from pathlib import Path
import subprocess

selection_script = Path.cwd().parents[0] / 'catalogues' / 'selection.py'
convert_script = Path.cwd() / 'convert_fits_txt.py'

# #! FIRST RUN SELECTION AND CONVERSION STEPS
# print("Running selection.py...")
# subprocess.run(['python3', str(selection_script)], check=True)

# # Execute convert_fits.txt.py
# print("Running convert_fits_txt.py...")
# subprocess.run(['python3', str(convert_script)], check=True)

#! LBG
zphot_dir = Path.cwd().parents[1] / 'data' / 'sed_fitting' / 'zphot' / 'det_Je_J_LBG'
#buildLePhareLibrary(parameter_file='euclid.para', build_libs=True, build_filters=True, build_mags=True)
if not zphot_dir.exists():
    zphot_dir.mkdir(parents=True)
runPhotometricRedshifts(parameter_file='euclid.para', zphot_dir=zphot_dir, overwrite=True)

#! LAE
# zphot_dir = Path.cwd().parents[1] / 'data' / 'sed_fitting' / 'zphot' / 'det_Ye_y_LAE'
# buildLePhareLibrary(parameter_file='euclid_lae.para', build_libs=True, build_filters=True, build_mags=True)
# if not zphot_dir.exists():
#     zphot_dir.mkdir(parents=True)
# runPhotometricRedshifts(parameter_file='euclid_lae.para', zphot_dir=zphot_dir, overwrite=True)