"""
Measure continuum level of galaxy
"""

from pathlib import Path
import numpy as np
import sys
import matplotlib.pyplot as plt
from astropy.table import Table
from scipy.optimize import least_squares
from astropy.cosmology import Planck18 as cosmo
from scipy.interpolate import interp1d


zphot_dir = Path.cwd().parents[1] / 'data' / 'sed_fitting' / 'zphot' / 'COSMOS' / 'best_fits' / 'det_Y_J_with_euclid_z7'

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

# FLUXES in HSC-y, VISTA-Y, Euclid-Y
fluxes = [1.283925e-30, 2.204734e-30, 1.61798e-30]
errors = [3.832285e-31, 2.756622e-31, 2.611918e-31]

fluxes, errors = np.array(fluxes), np.array(errors)

# Redshift and Lyα properties
z = 7.19
lya_rest = 1215.67
lya_obs = lya_rest * (1 + z)
fwhm_kms = 300.0
c_kms = 2.998e5

# Read filter transmission curves
filter_dir = Path.home() / 'lephare' / 'lephare_dev' / 'filt' / 'myfilters'
hsc_y = Table.read(filter_dir / 'HSC' / 'y_HSC.txt', format='ascii')
vista_Y = Table.read(filter_dir / 'VISTA' / 'VISTA_Y.txt', format='ascii')
euclid_Y = Table.read(filter_dir / 'Euclid' / 'Euclid_Y.txt', format='ascii')
filters = [hsc_y, vista_Y, euclid_Y]


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


# -------------------
# Helper function: Gaussian line
# -------------------
def gaussian_line(wlen, EW, F_cont):
    """Return Gaussian line flux (erg/s/cm^2/Å) normalized to EW and continuum."""
    fwhm_obs = fwhm_kms * (1 + z)  # km/s
    fwhm_ang = lya_obs * fwhm_obs / c_kms
    sigma = fwhm_ang / (2 * np.sqrt(2 * np.log(2)))
    
    # Integrated line flux = EW * continuum
    F_line = EW * F_cont
    norm = F_line / (sigma * np.sqrt(2 * np.pi))
    return norm * np.exp(-0.5 * ((wlen - lya_obs)/sigma)**2)

# -------------------
# Helper: synthetic photometry
# -------------------
def filter_average(wlen, spec, filt):
    lam = filt['col1'].data
    trans = filt['col2'].data
    interp_spec = np.interp(lam, wlen, spec, left=0, right=0)
    num = np.trapz(interp_spec * trans * lam, lam)
    denom = np.trapz(trans * lam, lam)
    return num / denom

def model_fluxes(EW, F_cont):
    """Return fluxes in all filters for a given Lyα EW and continuum."""
    spec = np.copy(primary_sed)
    spec += gaussian_line(primary_wlen, EW, F_cont)
    return np.array([filter_average(primary_wlen, spec, f) for f in filters])

def residuals(EW, F_cont, data, errs):
    return (model_fluxes(EW, F_cont) - data) / errs

# -------------------
# Fit Lyα EW to photometry
# -------------------
F_cont = median_continuum  # use measured continuum from SED
p0 = [240]  # starting guess for EW in Å

res = least_squares(residuals, p0, args=(F_cont, fluxes, errors))
best_EW = res.x[0]
print(f"Best-fit Lyα EW: {best_EW:.1f} Å")

# -------------------
# Monte Carlo for EW uncertainty
# -------------------
nmc = 2000
rng = np.random.default_rng(42)
EW_samples = []

for i in range(nmc):
    perturbed_fluxes = fluxes + rng.normal(0, errors)
    r = least_squares(residuals, p0, args=(F_cont, perturbed_fluxes, errors))
    EW_samples.append(r.x[0])

EW_samples = np.array(EW_samples)
EW_med = np.median(EW_samples)
EW_sigma = np.std(EW_samples)

print(f"Lyα EW = {EW_med:.1f} ± {EW_sigma:.1f} Å (1σ)")
print(f"Lyα EW = {EW_med:.1f} ± {3*EW_sigma:.1f} Å (3σ)")

# Optional: plot histogram
plt.hist(EW_samples, bins=50, color='lightblue', edgecolor='k')
plt.axvline(EW_med, color='r', linestyle='--', label='median')
plt.xlabel("Lyα EW [Å]")
plt.ylabel("Counts")
plt.legend()
plt.title("Monte Carlo distribution of Lyα EW")
plt.show()
