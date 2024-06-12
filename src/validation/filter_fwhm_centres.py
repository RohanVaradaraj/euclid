"""
filter_transmission_curves.py

Created: Tuesday 16th April 2024.
"""

from pathlib import Path
import glob
import matplotlib.pyplot as plt
import re

plt.rcParams['axes.linewidth'] = 2.5
plt.rcParams.update({'font.size': 15})
plt.rcParams['figure.dpi'] = 100

filter_dir = Path().home() / 'lephare' / 'lephare_dev' / 'filt' / 'myfilters'
plot_dir = Path.cwd().parent.parent / 'plots' / 'filters'

euclid_filters = glob.glob(str(filter_dir / 'Euclid' / 'Euclid_*'))
order = ['VIS', 'Y', 'J', 'H']

# Sort the euclid filter paths by the order of the filters above
def custom_sort(filename):
    for idx, substring in enumerate(order):
        if substring in filename:
            return idx

    return len(order)

euclid_filters = sorted(euclid_filters, key=custom_sort)

euclid_colours = ['purple', 'blue', 'green', 'red']

plt.figure(figsize=(14, 6))

for i, euclid_filter in enumerate(euclid_filters):

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

    # Normalise
    transmission = [t / max(transmission) for t in transmission]

    # Find the midpoint and FWHM of the filter
    midpoint = wavelength[transmission.index(max(transmission))]
    fwhm = wavelength[transmission.index(max(transmission))] - wavelength[transmission.index(max(transmission)) - 1]

    # Plot these on top of the filter
    plt.plot(midpoint, 1.1, 'o', color=euclid_colours[i], zorder=20)
    # Plot fwhm
    plt.plot([midpoint - fwhm/2, midpoint + fwhm/2], [1.1, 1.1], color=euclid_colours[i], lw=2., zorder=20)

    # Plot filter
    plt.plot(wavelength, transmission, color=euclid_colours[i], lw=2., alpha=0.8, zorder=10)

    # Add text at midpoint of each filter according to order labels
    mid_wavelength = (max(wavelength) - min(wavelength)) / 2 + min(wavelength)
    plt.text(mid_wavelength, 0.8, order[i], color=euclid_colours[i], fontsize=15, ha='center', va='center', zorder=20)


plt.xlabel(r'$\lambda \ (\mu \mathrm{m})$')
plt.ylabel('Relative Transmission')
plt.legend()
plt.show()







vista_filters = glob.glob(str(filter_dir / 'VISTA' / 'VISTA_*'))

for vista_filter in vista_filters:
    
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

    # Normalise
    transmission = [t / (max(transmission)*1.3) for t in transmission]

    # Plot filter
    plt.plot(wavelength, transmission, color='darkorange', lw=2., alpha=0.5)
    if vista_filter == vista_filters[0]:
        # Dummy label
        plt.plot([], [], color='darkorange', lw=2., alpha=0.5, label='VISTA')

hsc_labels = ['g', 'r', 'i', 'z', 'y']
hsc_filters = list( set(glob.glob(str(filter_dir / 'HSC' / '*'))) - set(glob.glob(str(filter_dir / 'HSC' / '*nb*'))))
hsc_filters = sorted(hsc_filters, key=custom_sort)
print(hsc_filters)

for hsc_filter in hsc_filters:
        
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
    
        # Normalise
        transmission = [t / (max(transmission)*1.3) for t in transmission]
    
        # Plot filter
        plt.plot(wavelength, transmission, color='steelblue', lw=2., alpha=0.5)
        if hsc_filter == hsc_filters[0]:
            # Dummy label
            plt.plot([], [], color='steelblue', lw=2., alpha=0.5, label='HSC')


jwst_filters = glob.glob(str(filter_dir / 'JWST' / '*angstroms.txt'))
order = ['f115w', 'f150w', 'f277w', 'f444w']
jwst_filters = sorted(jwst_filters, key=custom_sort)


