"""
filter_transmission_curves.py

Created: Tuesday 16th April 2024.
"""

from pathlib import Path
import glob
import matplotlib.pyplot as plt
import re

plt.rcParams['axes.linewidth'] = 4
plt.rcParams.update({'font.size': 25})
plt.rcParams['figure.dpi'] = 100

filter_dir = Path().home() / 'lephare' / 'lephare_dev' / 'filt' / 'myfilters'
plot_dir = Path.cwd().parent.parent / 'plots' / 'filters'

euclid_filters = glob.glob(str(filter_dir / 'Euclid' / 'Euclid_*'))
order = [r'$I_{E}$', r'$Y_{E}$', r'$J_{E}$', r'$H_{E}$']

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

euclid_colours = ['purple', 'blue', 'green', 'red']
euclid_colours = ['blue', 'blue', 'blue', 'blue']

plt.figure(figsize=(16, 9))

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

    # Plot filter
    plt.plot(wavelength, transmission, color=euclid_colours[i], lw=lw, alpha=0.8, zorder=10)

    # Add text at midpoint of each filter according to order labels
    mid_wavelength = (max(wavelength) - min(wavelength)) / 2 + min(wavelength)
    #plt.text(mid_wavelength, 0.85, order[i], color=euclid_colours[i], fontsize=25, ha='center', va='center', zorder=20)

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
    plt.plot(wavelength, transmission, color='darkorange', lw=lw, alpha=0.5)
    if vista_filter == vista_filters[0]:
        # Dummy label
        plt.plot([], [], color='darkorange', lw=lw, alpha=0.5, label='VISTA')

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
        plt.plot(wavelength, transmission, color='steelblue', lw=lw, alpha=0.5)
        if hsc_filter == hsc_filters[0]:
            # Dummy label
            plt.plot([], [], color='steelblue', lw=lw, alpha=0.5, label='HSC')

# Add CWeb filters
for cweb_filter in cweb_filters:
    
    # Open filter in two column format
    with open(cweb_filter, 'r') as f:
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
    #wavelength = [w * 1e-4 for w in wavelength]

    # Normalise
    transmission = [t / (max(transmission)*1.6) for t in transmission]

    # Plot filter
    plt.plot(wavelength, transmission, color='red', lw=7.) #, alpha=0.5, linestyle='--')
    if cweb_filter == cweb_filters[0]:
        # Dummy label
        plt.plot([], [], label='COSMOS-Web', color='red', lw=4.) #, alpha=0.5, linestyle='--')

plt.xlabel(r'$\lambda \ (\mu \mathrm{m})$', fontsize=35)
plt.ylabel('Relative Transmission', fontsize=30)


# Full filter set
plt.xlim(0.2, 2.55)
plt.ylim(-0.09, 1.6)
plt.tight_layout()

# Make ticks markers bigger
plt.tick_params(axis='both', which='major', size=10, width=4)

# And tick labels bigger
plt.xticks(fontsize=30)
plt.yticks(fontsize=30)

# Remove ylabel markers
#plt.yticks([])

# Define the x-range (wavelength in microns) and the height
x_range = [1.2, 3.0]  # Corresponds to 12000 Å - 30000 Å
y_min = 0.0
y_max = 0.8

# Get redshifts of lyman break at z=8, 9, 10
lyman_break_rest = 0.1216
lyman_break_observed = [lyman_break_rest * (1 + z) for z in [8, 9, 10]]
heights = [0.8, 0.6, 0.4]
color = ['lightgrey', 'grey', 'black']

# Plot the Lyman break as filled regions
for i, z in enumerate(lyman_break_observed):
    plt.fill_between([z, z + 1000], y_min, heights[i], color=color[i], alpha=0.7, label=r'$z={}$'.format(8+i))

#plt.fill_between(x_range, y_min, y_max, color='grey', alpha=0.7)
plt.legend(loc='upper right', fontsize=29)
plt.savefig(plot_dir / 'filter_transmission_curves_with_JWST.pdf', bbox_inches='tight')



# Zoomed in around Y
#plt.xlim(0.7, 1.7)
#plt.savefig(plot_dir / 'filter_transmission_curves_around_Ye.png', bbox_inches='tight')

plt.show()