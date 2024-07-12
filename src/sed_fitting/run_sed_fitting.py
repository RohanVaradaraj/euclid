"""
Run the SED fitting.

Created: Friday 12th July 2024.
"""

from sed_fitting_codes import *
from pathlib import Path

zphot_dir = Path.cwd().parents[1] / 'data' / 'sed_fitting' / 'zphot' / 'test_euclid'

buildLePhareLibrary(parameter_file='euclid.para', build_libs=True, build_filters=True, build_mags=True)
runPhotometricRedshifts(parameter_file='euclid.para', zphot_dir=zphot_dir, overwrite=True)