"""
Measure continuum level of galaxy
"""

from pathlib import Path
import numpy as np
import sys
import matplotlib.pyplot as plt

zphot_dir = Path.cwd().parents[1] / 'data' / 'sed_fitting' / 'zphot' / 'COSMOS' / 'best_fits' / 'det_Y_J_with_euclid_best_lya'

sed_path = Path.cwd().parents[0] / 'sed_fitting'
sys.path.append(str(sed_path))
from sed_fitting_codes import parse_spec_file

def mag_to_flux(mag):
    return 10**(-0.4 * (mag + 48.6))

file_name = 'Id000591643.spec'
spec_file = zphot_dir / file_name

file = parse_spec_file(spec_file)
phot = file.get('phot')
zpdf = file.get('zpdf')
sed = file.get('sed')
params = file.get('model')

names_phot=['phot', 'yerr', 'wlen', 'xerr', 'modelPhot', 'col6', 'col7']
names_sed = ['wlen', 'flux']
names_zpdf = ['z', 'P']
names_param = ['Type', 'Nline', 'Model', 'Library', 'Nband', 'Zphot', 'Zinf', 'Zsup', 'Chi2', 'PDF', 'Extlaw', 'EB-V', 'Lir', 'Age', 'Mass', 'SFR', 'SSFR']

phot.rename_columns(phot.colnames, names_phot)
zpdf.rename_columns(zpdf.colnames, names_zpdf)
params.rename_columns(params.colnames, names_param)
for table in sed:
    table.rename_columns(table.colnames, names_sed)

# Get each of the SEDs
model_wlens = [table['wlen'] for table in sed]
model_fluxes = [table['flux'] for table in sed]  

primary_sed = mag_to_flux(model_fluxes[0])
primary_wlen = model_wlens[0]

plt.plot(primary_wlen, primary_sed, label='Primary SED', color='blue')
plt.yscale('log')
plt.show()

# Select wavelength range: 11000 to 20000 Å
mask = (primary_wlen >= 11000) & (primary_wlen <= 20000)

wlen_subset = primary_wlen[mask]
flux_subset = primary_sed[mask]

# Compute average and median continuum level
mean_continuum = np.mean(flux_subset)
median_continuum = np.median(flux_subset)

print(f'Mean Continuum Level: {mean_continuum}')
print(f'Median Continuum Level: {median_continuum}')

