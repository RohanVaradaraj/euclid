"""
brown_dwarf_colours.py

Look at the colours of brown dwarfs when Euclid data is included.

Created: Thursday 18th April 2024.
"""

import numpy as np
import matplotlib.pyplot as plt
from astropy.io import ascii, fits
from pathlib import Path
import glob
from astropy.constants import c
from scipy.integrate import simps

plt.rcParams['axes.linewidth'] = 2.5
plt.rcParams.update({'font.size': 15})
plt.rcParams['figure.dpi'] = 100

# Set up directories
filter_dir = Path().home() / 'lephare' / 'lephare_dev' / 'filt' / 'myfilters'
plot_dir = Path.cwd().parent.parent / 'plots' / 'brown_dwarfs'
dwarf_dir = Path().home() / 'lephare' / 'lephare_dev' / 'sed' / 'STAR' / 'DWARFSTARS'



# Sort the euclid filter paths by the order of the filters above
def custom_sort(filename, order):
    for idx, substring in enumerate(order):
        if substring in filename:
            return idx
        


def mag_to_flux(m):
	'''Convert mags and their errors to flux'''
	flux = 10**(-0.4*(m+48.6))
	return flux



def flux_to_mag(flux):

    '''Convert flux to mag'''
    # Deal with negative flux case
    if flux <= 0:
        mag = np.nan
        print('Negative flux')
    else:
        mag = -2.5*np.log10(flux)-48.6
    return mag




def getFilters(instrument: str, plot: bool = False, plot_kwargs: dict = None) -> dict:
    """
    By specifying the instrument name, retrieve/plot all of its photometric filters.

    Based on src/validation/filter_transmission_curves.py.

    Parameters:
    -----------
    instrument: str
        The name of the instrument. One of 'Euclid', 'VISTA', 'HSC', 'Spitzer'
    plot: bool, optional
        Whether to plot the filters. Default is False.
    plot_kwargs: dict, optional
        Additional keyword arguments to customize the plot, such as color, linewidth, alpha, etc.

    Returns:
    --------
    dict
        A dictionary with the filter names as keys and the wavelength (microns) and transmission as values.
        The naming convention is as follows:
            Euclid: VIS, Ye, Je, He
            VISTA: Y, J, H, Ks
            HSC: g, r, i, nb816, z, nb921, y
            Spitzer: ch1, ch2
    """

    # Initialize plot_kwargs if not provided
    if plot_kwargs is None:
        plot_kwargs = {}

    # Get the filters and define the filter order
    if instrument.lower() == 'euclid':
        filter_names = glob.glob(str(filter_dir / 'Euclid' / 'Euclid_*'))
        order = ['VIS', 'Y', 'J', 'H']
    elif instrument.upper() == 'VISTA':
        filter_names = glob.glob(str(filter_dir / 'VISTA' / 'VISTA_*')) 
        order = ['Y', 'J', 'H', 'Ks']
    elif instrument.upper() == 'HSC':
        filter_names = glob.glob(str(filter_dir / 'HSC' / '*'))
        order = ['g', 'r', 'i', 'nb816', 'z', 'nb921', 'y']
    elif instrument.upper() == 'SPITZER':
        filter_names = glob.glob(str(filter_dir / 'SPITZER' / 'irac_*'))
        order = ['ch1', 'ch2']
    else:
        raise ValueError("Invalid instrument name")

    # Sort the filters
    filter_names = sorted(filter_names, key=custom_sort(filter_names, order))

    # Empty objects for transmission and wavelength
    wavelengths = []
    transmissions = []

    for i, filter_name in enumerate(filter_names):
        # Open filter in two column format
        with open(filter_name, 'r') as f:
            lines = f.readlines() 

        # Extract wavelength and transmission
        wavelength = []
        transmission = []
        for line in lines:
            if line[0] != '#':
                values = line.split()
                wavelength.append(float(values[0]))
                transmission.append(float(values[1]))

        # If vista, divide transmission by 100
        if instrument.lower() == 'vista':
            transmission = [t / 100 for t in transmission] #! Makes no difference in convolution, just for plotting.


        # Add to empty lists
        wavelengths.append(wavelength)
        transmissions.append(transmission)


        if plot:
            # Plot filter with customized plotting attributes
            plt.plot(wavelength, transmission, **plot_kwargs)

    # Modify labels for Euclid filters if necessary
    if instrument.lower() == 'euclid':
        order = ['VIS', 'Ye', 'Je', 'He']

    # Return a dictionary with the order name and the data
    return {order[i]: (wavelengths[i], transmissions[i]) for i, filter_name in enumerate(filter_names)}





def loadBrownDwarfTemplates(spectral_types: list[str]=['M', 'L', 'T'], sub_types: list[int] = None) -> dict:

    """
    Load the brown dwarf templates from the SpeX Prism Library.

    Parameters:
    -----------

    spectral_types: list[str]
        The spectral types of the brown dwarfs to load. Any combo of 'M', 'L', 'T'.
    sub_types: list[str]
        The subtypes of the brown dwarfs to load, if something specific is required. Any combo of 0 through 9 is allowed.
        The available subtypes for each are
            M: 4-9
            L: 0-9
            T: 0-8

    Returns:
    --------
    dict
        A dictionary with the spectral type as the key and the wavelength and flux as the values in a tuple.

    """

    if sub_types == None:
        sub_types = np.arange(0, 10)

    # Load the brown dwarf templates
    brown_dwarfs = {}
    
    for spectral_type in spectral_types:
        for sub_type in sub_types:

            if spectral_type == 'M' and int(sub_type) < 4:
                continue
            if spectral_type == 'T' and int(sub_type) > 8:
                continue

            brown_dwarf_file = dwarf_dir / f'{spectral_type}{sub_type}_reformat2.txt_trim'

            # Read data
            data = ascii.read(brown_dwarf_file, data_start=0)
            wavelength = np.array(data['col1'])
            flux = np.array(data['col2'])

            brown_dwarfs[f'{spectral_type}{sub_type}'] = (wavelength, flux)

    return brown_dwarfs


def convolveFilters(filter_set: list[dict], dwarf_templates: dict) -> dict:
    """
    Given a set of filters and brown dwarf templates, convolve the filters with the templates to get resultant fluxes.

    Parameters:
    -----------
    filter_set: list[dict]
        A list of dictionaries with the filter names as keys and the wavelength and transmission as values.
        Multiple dictionaries can be passed, the function will stack them.
    dwarf_templates: dict
        A dictionary with the spectral type as the key and the wavelength and flux as the values in a tuple.

    Returns:
    --------
    dict{spectral_type_1: {filter1: mag1, filter2: mag2, ...}, spectral_type_2: {filter1: mag1, filter2: mag2, ...}, ...}
        A dictionary with the spectral type as the key and a dictionary with the filter name and the resultant magnitude as the value.
    """

    # If there are multiple dicts in filter_set, turn it into one big dict
    filters = {}
    for filter_dict in filter_set:
        filters.update(filter_dict)

    # Go through the brown dwarf templates and convolve them with the filters
    convolved_fluxes = {}

  # Loop through the dwarf templates
    for spectral_type, (bd_wlen_grid, flux) in dwarf_templates.items():

        # Loop through the filters
        for filter_name, (filter_wlen_grid, transmission) in filters.items():

            # Interpolate the BD SED grid to the filter
            flux_interp = np.interp(filter_wlen_grid, bd_wlen_grid, flux)

            # Normalise the BD sed to a flux of 1e-29 erg/s/cm^2/A, roughly 24-26 AB mag
            flux_interp /= np.max(flux_interp)
            flux_interp *= 1e-29

            # Convert to frquency space (f = c/wlen)
            filter_freq_grid = np.array([c.value / (wlen*1e-10) for wlen in filter_wlen_grid])

            # Calculate the convolved flux, normalised by filter area under curve
            convolved_flux = simps(flux_interp * transmission, filter_freq_grid) / simps(transmission, filter_freq_grid)

            mag = flux_to_mag(convolved_flux)

            # Store with keys for each BD corresponding to a dict with a key for each filter
            convolved_fluxes.setdefault(spectral_type, {})[filter_name] = mag
                
    return convolved_fluxes







bds = loadBrownDwarfTemplates()

# Plot BD templates
for spectral_type, (wavelength, flux) in bds.items():
    #plt.figure(figsize=(10, 6))

    print(spectral_type)

    # Filters
    euclid_filters = getFilters('euclid') #, plot=True, plot_kwargs={'linewidth': 2.5, 'alpha':0.6, 'color': 'blue'})
    vista_filters = getFilters('vista') #, plot=True, plot_kwargs={'linewidth': 2.5, 'alpha':0.6, 'color': 'orange'})

    # Get the magnitudes
    mags = convolveFilters([vista_filters, euclid_filters], {spectral_type: (wavelength, flux)})
    # Get the Y-Ye, J-Je, H-He colours
    Ye = mags[spectral_type]['Ye']
    Y = mags[spectral_type]['Y']
    Je = mags[spectral_type]['Je']
    J = mags[spectral_type]['J']
    He = mags[spectral_type]['He']
    H = mags[spectral_type]['H']

    # Print to 1dp
    # print(f'Y - Ye: {Y-Ye:.1f}')
    # print(f'J - Je: {J-Je:.1f}')
    # print(f'H - He: {H-He:.1f}')
    print(f'Y - J: {Y-J:.1f}')

    # # Dummy plots for filter labels
    # plt.plot([], [], label='Euclid', color='blue', linewidth=2.5, alpha=0.6)
    # plt.plot([], [], label='VISTA', color='orange', linewidth=2.5, alpha=0.6)

    # # Brown dwarf template
    # plt.plot(wavelength, flux/np.max(flux), label=f'{spectral_type} dwarf', color='red', alpha=0.8, linewidth=2.5)

    # plt.xlabel(r'$\lambda \ (\AA)$')
    # plt.ylabel('Relative Transmission/Flux')
    # plt.xlim(5000, 25000)
    # plt.ylim(0, 1.15)
    # plt.legend(loc='upper right')
    # plt.tight_layout()

    # plt.savefig(plot_dir / f'{spectral_type}_template.png')
    # plt.close()

    #plt.show()

exit()

mags = convolveFilters([vista_filters, euclid_filters], loadBrownDwarfTemplates())

# Plot the colours, Ye-Y vs Je - J
plt.figure(figsize=(10, 8))

# Loop through the mags
for spectral_type, mags in mags.items():
        plt.plot(mags['Ye'] - mags['J'], mags['Ye'] - mags['Y'], 'o', label=spectral_type, color='black')

plt.xlabel('Ye - Y')
plt.ylabel('Ye - J')

plt.show()







