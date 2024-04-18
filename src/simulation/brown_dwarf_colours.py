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




def getFilters(instrument: str, plot: bool = False, plot_kwargs: dict = None) -> dict:
    """
    By specifying the instrument name, plot all of its photometric filters.

    Based on validation/filter_transmission_curves.py.

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

        # Convert wavelength to microns
        wavelength = [w * 1e-4 for w in wavelength]

        # Normalise
        transmission = [t / max(transmission) for t in transmission]

    if plot:
        # Plot filter with customized plotting attributes
        plt.plot(wavelength, transmission, **plot_kwargs)

    # Modify labels for Euclid filters if necessary
    if instrument.lower() == 'euclid':
        order = ['VIS', 'Ye', 'Je', 'He']

    # Return a dictionary with the order name and the data
    return {order[i]: (wavelength, transmission) for i, filter_name in enumerate(filter_names)}





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

            # Read uncommented two column file, no header
            data = ascii.read(brown_dwarf_file, data_start=0)
            wavelength = data['col1']
            flux = data['col2']

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
    dict
        A dictionary with the filter names as keys and the flux as values.
    """

    # If there are multiple dicts in filter_set, turn it into one big dict
    filters = {}
    for filter_dict in filter_set:
        filters.update(filter_dict)

    return filters




euclid_filters = getFilters('euclid')
vista_filters = getFilters('vista')
hsc_filters = getFilters('hsc')

filters = convolveFilters([hsc_filters, vista_filters, euclid_filters], loadBrownDwarfTemplates())

print(filters.keys())






