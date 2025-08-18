from pathlib import Path
import glob
import matplotlib.pyplot as plt
import numpy as np
import sys

# Set plotting defaults
plt.rcParams.update({'axes.linewidth': 4, 'font.size': 25, 'figure.dpi': 100})

plt.rcParams['axes.linewidth'] = 4
plt.rcParams.update({'font.size': 25})
plt.rcParams['figure.dpi'] = 100

plt.rcParams['xtick.top'] = True
plt.rcParams['xtick.bottom'] = True
plt.rcParams['ytick.left'] = True
plt.rcParams['ytick.right'] = True

plt.rcParams['xtick.direction'] = 'in'
plt.rcParams['ytick.direction'] = 'in'

# Make ticks larger and thicker
plt.rcParams['xtick.major.size'] = 10  # Length of major ticks
plt.rcParams['ytick.major.size'] = 10

plt.rcParams['xtick.major.width'] = 3  # Width of major ticks
plt.rcParams['ytick.major.width'] = 3

# (Optional) control minor ticks too
plt.rcParams['xtick.minor.size'] = 5
plt.rcParams['ytick.minor.size'] = 5
plt.rcParams['xtick.minor.width'] = 2
plt.rcParams['ytick.minor.width'] = 2


# Constants and config
MAG_FLOOR = 30
LW = 5

# Paths
base_path = Path.cwd().parents[0]
sed_fitting_path = base_path / 'sed_fitting'
sys.path.append(str(sed_fitting_path))
from sed_fitting_codes import parse_spec_file

filter_dir = Path.home() / 'lephare' / 'lephare_dev' / 'filt' / 'myfilters'
plot_dir = base_path.parent / 'plots' / 'filters'
sed_dir = base_path.parent / 'data' / 'sed_fitting' / 'zphot' / 'COSMOS' / 'best_fits' / 'det_Y_J_with_euclid_best_highz_good'

# Functions

def mag_to_flux(mag):
    return 10**(-0.4 * (mag + 48.6))

def flux_to_mag(flux):
    '''Convert flux to mag, for the secondary y axis'''
    mag = -2.5*np.log10(flux)-48.6
    return mag

def load_filter_data(filepath):
    with open(filepath, 'r') as f:
        lines = [line.strip() for line in f if line[0] != '#']
    data = np.array([[float(val) for val in line.split()] for line in lines])
    return data[:, 0] * 1e-4, data[:, 1]  # Convert Angstrom to micron

def normalize_transmission(trans, depth):
    return MAG_FLOOR - (trans / np.max(trans)) * (MAG_FLOOR - depth)

def compute_fwhm(wavelength, trans, depth):
    norm = (MAG_FLOOR - trans) / (MAG_FLOOR - depth)
    crossings = np.where(np.diff(np.sign(norm - 0.5)))[0]
    if len(crossings) >= 2:
        interp = lambda x1, x2, y1, y2: x1 + (0.5 - y1) * (x2 - x1) / (y2 - y1)
        left = interp(wavelength[crossings[0]], wavelength[crossings[0]+1], norm[crossings[0]], norm[crossings[0]+1])
        right = interp(wavelength[crossings[1]], wavelength[crossings[1]+1], norm[crossings[1]], norm[crossings[1]+1])
        print(left, right)
        return left, right, (left + right) / 2
    return None, None, None

def process_filter_group(filters, depths, color=None, label=None, labels=None, text=False, colours=None):
    centres = []
    for i, fpath in enumerate(filters):
        wave, trans = load_filter_data(fpath)
        trans_mag = normalize_transmission(trans, depths[i])
        left, right, centre = compute_fwhm(wave, trans_mag, depths[i])

        # Pick color
        this_color = colours[i] if colours is not None else color

        if centre:
            centres.append(centre)
            plt.hlines(depths[i], left, right, colors=this_color, linewidth=6, alpha=1, zorder=15)
            if text and labels:
                plt.text(centre, depths[i] + 0.25, labels[i], color=this_color, fontsize=33, ha='center', va='center', zorder=20, fontweight='bold')

        if i == 0 and label:
            plt.plot([], [], color=this_color, lw=LW, alpha=1, label=label)

    return centres

def compute_mag_eff(sed_flux, sed_wlen, filt_path):
    with open(filt_path, 'r') as f:
        lines = [line for line in f if not line.startswith('#')]

    wavelength = np.array([float(line.split()[0]) for line in lines]) * 1e-4  # μm
    transmission = np.array([float(line.split()[1]) for line in lines])
    transmission /= transmission.max()

    # Convert to frequency (Hz)
    freq = c / (wavelength * 1e-6)
    sed_freq = c / (sed_wlen * 1e-6)

    # Sort arrays for interpolation
    transmission_interp = np.interp(
        np.sort(sed_freq), np.sort(freq), np.sort(transmission), left=0, right=0
    )

    sed_flux_sorted = sed_flux[np.argsort(sed_freq)]
    sed_freq_sorted = np.sort(sed_freq)

    numerator = np.trapz(sed_flux_sorted * transmission_interp, sed_freq_sorted)
    denominator = np.trapz(transmission_interp, sed_freq_sorted)

    flux_eff = numerator / denominator
    return flux_to_mag(flux_eff)


def plot_model_mags(flux_array, wlen_array, filters, centres, colours, marker, label=None, zorder=30):
    for i, filt in enumerate(filters):
        mag_eff = compute_mag_eff(flux_array, wlen_array, filt)
        plt.scatter(
            centres[i], mag_eff,
            color=colours[i], s=300, zorder=zorder,
            marker=marker, edgecolors='black'
        )
    if label:
        plt.plot([], [], color='gray', marker=marker, markersize=20,
                 label=label, alpha=1, lw=LW, markeredgecolor='black')
        

def load_sed_components(file_path, stretch=1.03):
    file = parse_spec_file(file_path)
    sed = file.get('sed')
    for table in sed:
        table.rename_columns(table.colnames, ['wlen', 'flux'])

    model_wlens = [table['wlen'].data / 1e4 for table in sed]  # Convert to microns
    model_fluxes = [table['flux'].data for table in sed]

    primary_flux = model_fluxes[0] / stretch
    secondary_flux = model_fluxes[1] / np.max(model_fluxes[1])
    stellar_flux = model_fluxes[2] / stretch

    return model_wlens, [primary_flux, secondary_flux, stellar_flux]

def compute_effective_mag(wavelength, transmission, sed_wlen, sed_flux):
    # Constants
    c = 2.99792458e8  # Speed of light in m/s

    # Convert to frequency space
    wl_m = wavelength * 1e-6
    freq = c / wl_m

    sed_wl_m = sed_wlen * 1e-6
    sed_freq = c / sed_wl_m

    # Sort for interpolation
    freq_sorted = freq[np.argsort(freq)]
    trans_sorted = transmission[np.argsort(freq)]
    sed_freq_sorted = sed_freq[np.argsort(sed_freq)]
    sed_flux_sorted = sed_flux[np.argsort(sed_freq)]

    # Interpolate filter transmission to SED frequency grid
    trans_interp = np.interp(sed_freq_sorted, freq_sorted, trans_sorted, left=0, right=0)

    # Integrate flux * transmission over frequency
    numerator = np.trapz(sed_flux_sorted * trans_interp, sed_freq_sorted)
    denominator = np.trapz(trans_interp, sed_freq_sorted)

    flux_eff = numerator / denominator
    return flux_to_mag(flux_eff)

def load_filter_curve(filter_path):
    with open(filter_path, 'r') as f:
        lines = [line for line in f if line[0] != '#']
    data = np.array([[float(val) for val in line.split()] for line in lines])
    wl = data[:, 0] * 1e-4  # Angstrom to microns
    trans = data[:, 1]
    trans /= trans.max()  # Normalize
    return wl, trans

# Euclid filters
order = [r'$I_{\rm{E}}$', r'$Y_{\rm{E}}$', r'$J_{\rm{E}}$', r'$H_{\rm{E}}$']
euclid_colours = ['tab:purple', 'tab:blue', 'tab:green', 'tab:red']
euclid_depths = [27.7, 26.4, 26.4, 26.3]
euclid_filters = sorted(
    glob.glob(str(filter_dir / 'Euclid' / 'Euclid_*')),
    key=lambda f: next((i for i, o in enumerate(order) if o in f), len(order))
)

# VISTA filters
vista_filters = glob.glob(str(filter_dir / 'VISTA' / 'VISTA_*'))
vista_depths = [25.7, 25.3, 26.0, 26.2]

# HSC filters
desired_order = ['g_HSC.txt', 'r_HSC.txt', 'i_HSC.txt', 'z_HSC.txt', 'y_HSC.txt', 'nb816_HSC.txt', 'nb921_HSC.txt']
hsc_labels = ['g', 'r', 'i', 'z', 'y', 'nb816', 'nb921']
hsc_filters = sorted(
    glob.glob(str(filter_dir / 'HSC' / '*')),
    key=lambda f: next((i for i, name in enumerate(desired_order) if f.endswith(name)), len(desired_order))
)
hsc_depths = [27.8, 27.4, 27.2, 26.7, 26.0, 26.2, 26.0]

# Plot setup
plt.figure(figsize=(18, 9))

print('Euclid')
fwhm_centres_euclid = process_filter_group(
    euclid_filters, euclid_depths, colours=euclid_colours, labels=order, text=True
)
for i, colour in enumerate(euclid_colours):
    plt.hlines([], [], [], colors=colour, linewidth=6, alpha=1)
print('VISTA')
fwhm_centres_vista = process_filter_group(vista_filters, vista_depths, color='tab:orange', label='VISTA')
print('HSC')
fwhm_centres_hsc = process_filter_group(hsc_filters, hsc_depths, color='black', label='HSC')

# Final plot tweaks
plt.xlabel(r'$\lambda \ [\rm{\mu} \mathrm{m}]$', fontsize=35)
plt.ylabel('Depth [mag]', fontsize=30)
plt.tick_params(axis='both', which='major', size=10, width=4)
plt.xticks(fontsize=30)
plt.yticks(fontsize=30)
plt.legend()
plt.tight_layout()



# Load SED components
model_wlens, model_fluxes = load_sed_components(sed_dir / 'Id000468417.spec')
primary_wlen, secondary_wlen, stellar_wlen = model_wlens
primary_flux, secondary_flux, stellar_flux = model_fluxes

# Plot SEDs
plt.plot(stellar_wlen, stellar_flux, color='gray', lw=LW, zorder=2, alpha=0.7)
plt.plot(primary_wlen, primary_flux, color='gray', lw=LW, zorder=2, linestyle='--', alpha=0.7)
plt.gca().invert_yaxis()

# Combine all filters
all_filters = euclid_filters + vista_filters + hsc_filters
all_colours = (
    ['tab:purple', 'tab:blue', 'tab:green', 'tab:red']
    + ['tab:orange'] * len(vista_filters)
    + ['black'] * len(hsc_filters)
)
fwhm_centres_all = fwhm_centres_euclid + fwhm_centres_vista + fwhm_centres_hsc

# Plot effective magnitudes for both SEDs
for i, fpath in enumerate(all_filters):
    wl, trans = load_filter_curve(fpath)

    mag_stellar = compute_effective_mag(wl, trans, stellar_wlen, mag_to_flux(stellar_flux))
    mag_primary = compute_effective_mag(wl, trans, primary_wlen, mag_to_flux(primary_flux))

    plt.scatter(fwhm_centres_all[i], mag_stellar, color=all_colours[i], s=300, zorder=30,
                marker='*', edgecolors='black')
    plt.scatter(fwhm_centres_all[i], mag_primary, color=all_colours[i], s=300, zorder=29,
                marker='o', edgecolors='black')

    #print(f"Filter: {Path(fpath).name}, mag (primary): {mag_primary:.3f}")

# Dummy legend markers
plt.plot([], [], color='gray', marker='o', markersize=20, linestyle='--',
         label=r'LBG at $z\simeq7$', alpha=1, lw=LW, markeredgecolor='black')
plt.plot([], [], color='gray', marker='*', markersize=20,
         label='Ultra-cool dwarf', alpha=1, lw=LW, markeredgecolor='black')

plt.minorticks_on()

plt.xlim(0.4, 2.5)
plt.ylim(29,23.5)
plt.legend(loc='lower right', fontsize=25, ncols=2)
plt.tight_layout()
plt.savefig(plot_dir / 'depths_and_model_mags.pdf', bbox_inches='tight')
plt.show()