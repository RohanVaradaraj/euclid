"""
filter_transmission_curves.py

Created: Tuesday 16th April 2024.
"""

from pathlib import Path
import glob
import matplotlib.pyplot as plt
import re
import numpy as np
import sys
from scipy.constants import c

sed_fitting_path = Path.cwd().parents[0] / 'sed_fitting'
sys.path.append(str(sed_fitting_path))
from sed_fitting_codes import parse_spec_file

mag_floor = 30  # or whatever faintest mag you want as baseline

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
plt.rcParams['xtick.minor.width'] = 1.5
plt.rcParams['ytick.minor.width'] = 1.5

def mag_to_flux(mag):
    return 10**(-0.4 * (mag + 48.6))

def mag_to_Jy(m):
    '''Convert mags and their errors to Jy'''
    flux = 10**(-0.4*(m+48.6))
    return flux

def flux_to_mag(flux):
    '''Convert flux to mag, for the secondary y axis'''
    mag = -2.5*np.log10(flux)-48.6
    return mag

filter_dir = Path().home() / 'lephare' / 'lephare_dev' / 'filt' / 'myfilters'
plot_dir = Path.cwd().parent.parent / 'plots' / 'filters'
sed_dir = Path.cwd().parents[1] / 'data' / 'sed_fitting' / 'zphot' / 'COSMOS' / 'best_fits' / 'det_Y_J_with_euclid_z7'

euclid_filters = glob.glob(str(filter_dir / 'Euclid' / 'Euclid_*'))
order = [r'$I_{E}$', r'$Y_{E}$', r'$J_{E}$', r'$H_{E}$']

# Labels like order but bold
euclid_labels = [
    r'$\boldsymbol{I_E}$',
    r'$\boldsymbol{Y_E}$',
    r'$\boldsymbol{J_E}$',
    r'$\boldsymbol{H_E}$'
]

cweb_filters = glob.glob(str(filter_dir / 'JWST' / '*'))

# Only take cweb filters in F115W, F150W, F277W, and F444W
cweb_filters = [f for f in cweb_filters if 'f115w' in f or 'f150w' in f or 'f277w' in f or 'f444w' in f]
print(cweb_filters)

lw = 5

# Sort the euclid filter paths by the order of the filters above
def custom_sort(filename):
    for idx, substring in enumerate(order):
        if substring in filename:
            return idx

    return len(order)

euclid_filters = sorted(euclid_filters, key=custom_sort)

euclid_colours = ['tab:purple', 'tab:blue', 'tab:green', 'tab:red']
euclid_depths = [27.7, 26.4, 26.4, 26.3]  # Depths for each filter in magnitudes
#euclid_colours = ['blue', 'blue', 'blue', 'blue']

plt.figure(figsize=(18, 9))

fwhm_centres_euclid = []

for i, euclid_filter in enumerate(euclid_filters[1:2]):

    # Open filter in two column format
    with open(euclid_filter, 'r') as f:
        lines = f.readlines() 

    # Extract wavelength and transmission
    wavelength = []
    transmission = []
    for line in lines:
        if line[0] != '#':
            values = line.split()
            wavelength.append(float(values[0]))
            transmission.append(float(values[1]))

    # Convert wavelength to microns
    wavelength = [w * 1e-4 for w in wavelength]

    # Normalise to self
    transmission = [t / max(transmission) for t in transmission]
  
    # Plot filter
    plt.plot(wavelength, transmission, color='tab:blue', lw=lw, alpha=0.8, zorder=10)

    # # Draw line of FWHM of the filtr


#     # Add text at midpoint of each filter according to order labels
#     mid_wavelength = (max(wavelength) - min(wavelength)) / 2 + min(wavelength)
#    # plt.text(mid_wavelength, 28, order[i], color=euclid_colours[i], fontsize=25, ha='center', va='center', zorder=20)

vista_filters = glob.glob(str(filter_dir / 'VISTA' / 'VISTA_*'))
vista_depths = [25.7, 25.3, 26.0, 26.2]

fwhm_centres_vista = []

for i, vista_filter in enumerate(vista_filters[-1:]):
    
    # Open filter in two column format
    with open(vista_filter, 'r') as f:
        lines = f.readlines() 

    # Extract wavelength and transmission
    wavelength = []
    transmission = []
    for line in lines:
        if line[0] != '#':
            values = line.split()
            wavelength.append(float(values[0]))
            transmission.append(float(values[1]))

    # Convert wavelength to microns
    wavelength = [w * 1e-4 for w in wavelength]

    # Normalise to self
    transmission = [t / max(transmission) for t in transmission]

    # Convert to numpy arrays for easy indexing
    # wavelength_arr = np.array(wavelength)
    # norm_trans_arr = np.array(normalized_transmission)

    # # Find indices where transmission crosses 0.5 (half max)
    # crossings = np.where(np.diff(np.sign(norm_trans_arr - 0.5)))[0]

    # if len(crossings) >= 2:
    #     # Linear interpolation to get exact crossing wavelengths
    #     def interp_x(x1, x2, y1, y2, y_target):
    #         return x1 + (y_target - y1) * (x2 - x1) / (y2 - y1)

    #     left_idx = crossings[0]
    #     right_idx = crossings[1]

    #     left_wavelength = interp_x(wavelength_arr[left_idx], wavelength_arr[left_idx + 1],
    #                             norm_trans_arr[left_idx], norm_trans_arr[left_idx + 1], 0.5)
    #     right_wavelength = interp_x(wavelength_arr[right_idx], wavelength_arr[right_idx + 1],
    #                             norm_trans_arr[right_idx], norm_trans_arr[right_idx + 1], 0.5)
        
    #     # Find fwhm centre
    #     fwhm_centre = (left_wavelength + right_wavelength) / 2
    #     fwhm_centres_vista.append(fwhm_centre)

        # Plot horizontal line at filter depth magnitude
    #     plt.hlines(vista_depths[i], left_wavelength, right_wavelength, colors='tab:orange', 
    #             linewidth=6, alpha=1, zorder=15)
    # else:
    #     print(f"Could not find FWHM crossings for filter {vista_filter}")

    # Plot filter
    plt.plot(wavelength, transmission, color='tab:orange', lw=lw, alpha=0.5)
    if vista_filter == vista_filters[0]:
        # Dummy label
        plt.plot([], [], color='tab:orange', lw=lw, alpha=1, label='VISTA')

hsc_labels = ['g', 'r', 'i', 'z', 'y']
hsc_filters = list( set(glob.glob(str(filter_dir / 'HSC' / '*')))  - set(glob.glob(str(filter_dir / 'HSC' / '*nb*'))))
desired_order = ['g_HSC.txt', 'r_HSC.txt', 'i_HSC.txt', 'z_HSC.txt', 'y_HSC.txt']
def hsc_sort_key(path):
    for i, name in enumerate(desired_order):
        if path.endswith(name):
            return i
    return len(desired_order)  # Put any unmatched filters at the end
hsc_filters = sorted(hsc_filters, key=hsc_sort_key)
print([h.split('/')[-1] for h in hsc_filters])
hsc_depths = [27.8, 27.4, 27.2, 26.7, 26.]  # Depths for each filter in magnitudes

fwhm_centres_hsc = []

for i, hsc_filter in enumerate(hsc_filters[-1:]):
        
        # Open filter in two column format
        with open(hsc_filter, 'r') as f:
            lines = f.readlines() 
    
        # Extract wavelength and transmission
        wavelength = []
        transmission = []
        for line in lines:
            if line[0] != '#':
                values = line.split()
                wavelength.append(float(values[0]))
                transmission.append(float(values[1]))
    
        # Convert wavelength to microns
        wavelength = [w * 1e-4 for w in wavelength]
            
        # Normalise to self
        transmission = [t / max(transmission) for t in transmission]

        # After normalizing transmission as you did:
        # max_trans = max(transmission)
        # half_max = mag_floor - (0.5 * (mag_floor - hsc_depths[i]))

        # # To find FWHM, first compute original normalized transmission (0 to 1 scale)
        # # Since your current 'transmission' is scaled to mags, revert that scale temporarily:
        # normalized_transmission = [(mag_floor - t) / (mag_floor - hsc_depths[i]) for t in transmission]

        # Convert to numpy arrays for easy indexing
        # wavelength_arr = np.array(wavelength)
        # norm_trans_arr = np.array(normalized_transmission)

        # # Find indices where transmission crosses 0.5 (half max)
        # crossings = np.where(np.diff(np.sign(norm_trans_arr - 0.5)))[0]

        # if len(crossings) >= 2:
        #     # Linear interpolation to get exact crossing wavelengths
        #     def interp_x(x1, x2, y1, y2, y_target):
        #         return x1 + (y_target - y1) * (x2 - x1) / (y2 - y1)

        #     left_idx = crossings[0]
        #     right_idx = crossings[1]

        #     left_wavelength = interp_x(wavelength_arr[left_idx], wavelength_arr[left_idx + 1],
        #                             norm_trans_arr[left_idx], norm_trans_arr[left_idx + 1], 0.5)
        #     right_wavelength = interp_x(wavelength_arr[right_idx], wavelength_arr[right_idx + 1],
        #                             norm_trans_arr[right_idx], norm_trans_arr[right_idx + 1], 0.5)
            
        #     # Find fwhm centre
        #     fwhm_centre = (left_wavelength + right_wavelength) / 2
        #     fwhm_centres_hsc.append(fwhm_centre)

        #     # Plot horizontal line at filter depth magnitude
        #     plt.hlines(hsc_depths[i], left_wavelength, right_wavelength, colors='black', 
        #             linewidth=6, alpha=1, zorder=15)
        #     # plt.text((left_wavelength + right_wavelength) / 2, hsc_depths[i] - 0.2, hsc_labels[i],
        #     #         color='black', fontsize=30, ha='center', va='center', zorder=20, fontweight='bold')
        # else:
        #     print(f"Could not find FWHM crossings for filter {hsc_filter}")

        # Plot filter
        plt.plot(wavelength, transmission, color='black', lw=lw, alpha=0.5)
        if hsc_filter == hsc_filters[0]:
            # Dummy label
            plt.plot([], [], color='black', lw=lw, alpha=1, label='HSC')






plt.xlim(0.7, 1.4)
plt.ylim(0, 1.1)
#plt.gca().invert_yaxis()  # Magnitudes: lower is brighter]
# plt.legend(loc='lower right', fontsize=25, ncols=2)
plt.xlabel(r'$\lambda \, [\rm{\mu m}]$', fontsize=30)
plt.ylabel('Transmission', fontsize=30)
plt.tight_layout()
plt.savefig(plot_dir / 'all_the_Y_filters.pdf', bbox_inches='tight')
plt.show()