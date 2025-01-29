from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table
from astropy.coordinates import SkyCoord
import astropy.units as u
from astropy.cosmology import FlatLambdaCDM

# --- Set Global Configurations and Constants ---
plt.rcParams['axes.linewidth'] = 2.5
plt.rcParams.update({'font.size': 20})
plt.rcParams['figure.dpi'] = 100

# Cosmology constants
cosmo = FlatLambdaCDM(H0=70, Om0=0.3)

z = 7
DL = cosmo.luminosity_distance(z).to(u.pc).value
M = np.arange(-26, -18, 0.1)

# --- Utility Functions ---
def size_luminosity_relation(r_0, M_0, M, beta):
    """Calculate the size-luminosity relation."""
    return r_0 * 10 ** (0.4 * beta * (M_0 - M))

def flux_to_mag(flux):
    """Convert flux to magnitude."""
    return -2.5 * np.log10(flux) - 48.6

def mAB_to_MUV(mAB, z=z, DL=DL):
    """Convert apparent magnitude to absolute magnitude."""
    return mAB - 5 * np.log10(DL / 10) + 2.5 * np.log10(1 + z)

def MUV_to_mAB(MUV, z=z, DL=DL):
    """Convert absolute magnitude to apparent magnitude."""
    return MUV + 5 * np.log10(DL / 10) - 2.5 * np.log10(1 + z)

def FWHM_arcsec_to_Re_kpc(fwhm, z):
    """Convert FWHM in arcseconds to effective radius (Re) in kpc."""
    as_per_kpc = cosmo.arcsec_per_kpc_proper(z).value
    return fwhm / (2 * np.sqrt(2 * np.log(2))) / as_per_kpc

def Re_kpc_to_FWHM_arcsec(r_e, z, psf_broadening=True):
    """Convert effective radius (Re) in kpc to FWHM in arcseconds."""
    as_per_kpc = cosmo.arcsec_per_kpc_proper(z).value
    fwhm = r_e * (4 * np.sqrt(2 * np.log(2))) / 1.678 * as_per_kpc
    if psf_broadening:
        fwhm = np.sqrt(fwhm**2 + 0.51**2)  # Include PSF broadening
    return fwhm

def calc_stats(fwhms, mags, mag_limit=28):
    """Calculate mean and standard deviation for FWHM and magnitude below a limit."""
    fwhm_mean = np.mean(fwhms[mags < mag_limit])
    fwhm_std = np.std(fwhms[mags < mag_limit])
    mag_mean = np.mean(mags[mags < mag_limit])
    mag_std = np.std(mags[mags < mag_limit])
    return fwhm_mean, fwhm_std, mag_mean, mag_std

# --- Crossmatching ---
def crossmatch_and_get_fwhms(candidates, depth_cat, filter_name, max_distance=1 * u.arcsec):
    """Crossmatch candidates with the depth catalog and return FWHMs and magnitudes."""
    coords_candidates = SkyCoord(ra=candidates['RA'], dec=candidates['DEC'], frame='icrs')
    coords_depth = SkyCoord(ra=depth_cat['ALPHA_J2000'], dec=depth_cat['DELTA_J2000'], frame='icrs')
    idx, d2d, _ = coords_candidates.match_to_catalog_sky(coords_depth)
    matched_idx = idx[d2d < max_distance]
    fwhms = depth_cat['FWHM_IMAGE'][matched_idx] * 0.1  # Convert to arcsec
    flux = candidates[f'flux_{filter_name}e'][d2d < max_distance]
    mag = flux_to_mag(flux)
    return fwhms, mag

# --- Directories and Data Loading ---
filter_name = 'J'
candidate_dir = Path.cwd().parents[1] / 'data' / 'catalogues' / 'candidates'
euclid_depth_dir = Path.cwd().parents[1] / 'data' / 'depths' / 'COSMOS' / 'catalogues'
vista_depth_dir = Path.home().parent.parent / 'vardy' / 'vardygroupshare' / 'data' / 'depths' / 'COSMOS' / 'catalogues'

# Load candidate catalogues
euclid_lbg = Table.read(candidate_dir / 'COSMOS_5sig_Y_J_nonDet_HSC_G_nonDet_HSC_R_nonDet_HSC_I_candidates_2024_11_28_with_euclid.fits')
euclid_bd = Table.read(candidate_dir / 'COSMOS_5sig_Y_J_nonDet_HSC_G_nonDet_HSC_R_nonDet_HSC_I_BD_INTERLOPERS_2024_11_28_with_euclid.fits')
vista_lbg = Table.read(candidate_dir / 'COSMOS_5sig_Y_J_nonDet_HSC_G_nonDet_HSC_R_nonDet_HSC_I_candidates_2024_11_20.fits')
vista_bd = Table.read(candidate_dir / 'COSMOS_5sig_Y_J_nonDet_HSC_G_nonDet_HSC_R_nonDet_HSC_I_BD_INTERLOPERS_2024_11_26.fits')

# Load depth catalogues
euclid_depth_cat = Table.read(euclid_depth_dir / f'd{filter_name}.fits')
vista_depth_cat = Table.read(vista_depth_dir / f'd{filter_name}.fits')

# --- Crossmatch Candidates ---
euclid_lbg_fwhms, euclid_lbg_mag = crossmatch_and_get_fwhms(euclid_lbg, euclid_depth_cat, filter_name)
euclid_bd_fwhms, euclid_bd_mag = crossmatch_and_get_fwhms(euclid_bd, euclid_depth_cat, filter_name)
vista_lbg_fwhms, vista_lbg_mag = crossmatch_and_get_fwhms(vista_lbg, vista_depth_cat, filter_name)
vista_bd_fwhms, vista_bd_mag = crossmatch_and_get_fwhms(vista_bd, vista_depth_cat, filter_name)



#! Begin figure
plt.figure(figsize=(10, 8))

# --- Plotting Function ---
def plot_fwhm_vs_mag(mags, fwhms, color, marker, size, alpha, zorder, label=None):
    """Plot FWHM vs Magnitude."""
    plt.scatter(mags, fwhms, color=color, marker=marker, alpha=alpha, s=size, edgecolor='none', label=label, zorder=zorder)

# --- Plot CWEB PSF sizes ---
fwhm_median = np.load('cweb_psfs_euclid_fwhm_median.npy')
fwhm_std = np.load('cweb_psfs_euclid_fwhm_std.npy')
fwhm_1sigma = np.load('cweb_psfs_euclid_fwhm_1sigma.npy')
fwhm_upper = fwhm_1sigma[1] + fwhm_median
fwhm_lower = fwhm_median - fwhm_1sigma[0]
# fwhm_upper = fwhm_median + fwhm_std
# fwhm_lower = fwhm_median - fwhm_std
bins = np.arange(20, 30, 0.5)
plt.fill_between(bins[:-1], fwhm_upper, fwhm_lower, color='black', alpha=0.2, edgecolor='none', zorder=0)

# --- Euclid PSF FWHM ---
plt.axhline(0.51, color='black', linestyle='dotted', linewidth=3, alpha=0.3, zorder=1)


# --- Plot FWHM vs J-Magnitude for Euclid ---
plot_fwhm_vs_mag(euclid_lbg_mag, euclid_lbg_fwhms, 'blue', 'o', 80, 0.65, 3)
plot_fwhm_vs_mag(euclid_bd_mag, euclid_bd_fwhms, 'red', '*', 230, 0.9, 5)

# --- Size-Luminosity Relations ---
mAB = MUV_to_mAB(M, z, DL)

relations = {
    "Roper+2022": size_luminosity_relation(1.126, -21, M, 0.29),
    "Kawamata+2018": size_luminosity_relation(0.94, -21, M, 0.46),
    "Shibuya+2015": size_luminosity_relation(0.75, -21, M, 0.27),
    "Yang+2022": size_luminosity_relation(10**(-0.09), -21, M, 0.50),
    "Sun+2022": size_luminosity_relation(0.85, -21, M, 0.28),
}

linestyles = ['-', '--', '-.', (0, (1, 1)), (0, (3, 1, 1, 1, 1, 1))]
colours = ['black', 'mediumorchid', 'seagreen', 'orange']
colours = ['orange', 'green', 'purple', 'black', 'turquoise']
alpha=[0.95, 0.95, 0.95, 0.95]

# Convert sizes to FWHM in arcseconds
for key in relations:
    relations[key] = Re_kpc_to_FWHM_arcsec(relations[key], z)

# Add literature relations
for i, (label, values) in enumerate(relations.items()):
    plt.plot(mAB, values, label=label, linestyle=linestyles[i], linewidth=4.5, alpha=0.8, zorder=4, color=colours[i])

# Labels and legends
ax2 = plt.gca().secondary_xaxis('top', functions=(mAB_to_MUV, MUV_to_mAB))
ax2.set_xlabel(r'$M_{\rm{UV}}$')

plt.xlabel(r'$J_{E}$ (mag)')
plt.ylabel('FWHM (arcsec)')

plt.tick_params(axis='both', which='major', width=2.5, length=5)
ax2.tick_params(axis='both', which='major', width=2.5, length=5)

plt.xlim(22.9, 28)
plt.ylim(0.3, 3.1)
plt.legend()
plt.tight_layout()
plt.savefig(Path.cwd().parents[1] / 'plots' / 'sizes' / 'euclid_fwhm_vs_Jmag.pdf')
plt.show()

