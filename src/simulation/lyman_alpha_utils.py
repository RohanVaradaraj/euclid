#!/usr/bin/env python3

"""
lyman_alpha_utils.py

Functions to simulate Lyman-alpha emitters (LAEs) passing through the Euclid+VISTA filter set.

Created: Friday 10th May 2024.
"""

import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table
from astropy.io import ascii
import sys
from pathlib import Path
from astropy.cosmology import FlatLambdaCDM
from astropy.constants import c
from synphot import etau_madau
from synphot import SourceSpectrum, Empirical1D, units
from scipy.integrate import simps
import astropy.units as u
from matplotlib.animation import FuncAnimation
import matplotlib.animation as animation
import matplotlib as mpl
from typing import Union, Tuple

sys.path.append(str(Path.cwd()))
from brown_dwarf_colours import *

plt.rcParams['axes.linewidth'] = 2.5
plt.rcParams.update({'font.size': 15})
plt.rcParams['figure.dpi'] = 100

# Use backend agg
#mpl.use('Agg')


filter_dir = Path().home() / 'lephare' / 'lephare_dev' / 'filt' / 'myfilters'
plot_dir = Path.cwd().parent.parent / 'plots' / 'LAEs'
table_dir = Path.cwd().parent.parent / 'data' / 'simulations' / 'tables'
rebels_dir = Path.cwd().parent.parent / 'data' / 'ref_catalogues' / 'fluxes'
sed_dir = Path.cwd().parent.parent / 'data' / 'simulations' / 'LAEs'



def makeLBG_vectorised(redshifts: np.array, SFH_component: str, Muv: np.array = None,
                       age: Union[np.array, Tuple[np.array, np.array]] = None,
                       tau: np.array = None, tmax: np.array = None, fwhm: np.array = None,
                       massformed: np.array = None, metallicity: np.array = None,
                       dust_type: str = None, Av: np.array = None,
                       nebular: bool = False, logU: np.array = None,
                       filter_list: Path = Path.cwd() / 'filter_list.txt') -> tuple[np.array, np.array]:
    """
    Create model Lyman-break galaxies using BAGPIPES with vectorized input.

    Parameters:
    -----------
    redshift: np.array
        The redshifts of the galaxies.
    SFH_component: str
        The star formation history component to use. One of:
        'burst', 'constant', 'exponential' (e^-(t/tau)), 'delayed' (t*e^-(t/tau)), 'lognormal'.
    Muv: np.array, optional
        The target absolute rest-UV magnitudes of the galaxies. SEDs are scaled to these if provided.
    age: Union[np.array, Tuple[np.array, np.array]], optional
        The ages of the galaxies in Gyr. If SFH is 'constant', provide a tuple (age_min, age_max) (= time since SF turned off/on).
    tau: np.array, optional
        The e-folding times in Gyr. Required for 'exponential' and 'delayed' SFHs.
    tmax: np.array, optional
        The times of maximum star formation in Gyr. Required for 'lognormal' SFH.
    fwhm: np.array, optional
        The full widths at half maximum of the lognormal distributions in Gyr.
    massformed: np.array, optional
        The masses of the galaxies formed in log_10(M*/M_solar)
    metallicity: np.array, optional
        The metallicities of the galaxies in Z/Z_solar.
    dust: str, optional
        The dust attenuation law to use.
    Av: np.array, optional
        The dust attenuations in magnitudes.
    nebular: bool, optional
        Whether to include nebular emission.
    logU: np.array, optional
        The ionization parameters.
    filter_list: Path, optional
        The path to the filter list file. Default is in this directory, 'filter_list.txt'.

    Returns:
    --------
    tuple[np.array, np.array]
        Arrays of the wavelength grids and fluxes of the model galaxies.
    """

    wlen_list = []
    flux_list = []

    # Divide length of EW by z to get the EW iteration
    
    ew_counter = 0

    for z in redshifts:

        if z == redshifts[0]:
            ew_counter += 1
            print(f'########## ITERATION NUMBER {ew_counter} ##########')



        print(f'Creating LBG at z = {round(z, 2)}')
        wlen, flux = makeLBG(z, SFH_component, Muv, age, tau, tmax, fwhm,
                                         massformed, metallicity, dust_type, Av,
                                         nebular, logU, filter_list)
        wlen_list.append(wlen)
        flux_list.append(flux)

    return np.array(wlen_list), np.array(flux_list)



def set_Muv(z: float, Muv_target: float, wlen: np.array, flux: np.array) -> (np.array, np.array):
    """
    Set the absolute magnitude of the LAEs at a given redshift.
    
    Parameters
    ----------
    z : float
        Redshift of the LAE.
    Muv : float
        Absolute magnitude of the LAE.
    wlen : array
        Wavelength grid of the input SED.
    flux : array
        Flux of the input SED.
    
    Returns
    -------
    wlen : array
        Wavelength grid of the output SED (same as input).
    flux : array
        Flux of the output SED, with the new absolute magnitude.
    """


    # Set up a cosmology
    cosmo = FlatLambdaCDM(H0=70, Om0=0.3)

    # Luminosity distance
    DL = cosmo.luminosity_distance(z).value * 10 ** 6 # convert to pc

    # Create top-hat filter at 1500A in observed frame
    tophat_filter = np.zeros(len(wlen))

    tophat_filter[(wlen > 1450*(1+z)) & (wlen < 1550*(1+z))] = 1.0

    #for i, lam in enumerate(wlen):
        # Set as 1 if within 1450 and 1550 in observed frame
    #    if (lam > 1450*(1+z)) and (lam < 1550*(1+z)):
    #        tophat_filter[i] = 1.0

    # Convolve the filter with the SED
    conv = np.sum(flux * tophat_filter) / np.sum(tophat_filter)

    # Apparent 1500A magnitude
    m1500 = flux_to_mag(conv)

    # Compute absolute magnitude
    M1500 = m1500 - 5*np.log10(DL/10) + 2.5*np.log10(1+z)

    # Target Muv is some factor of the current M1500
    A = Muv_target / M1500

    # Convert this factor A into a prefactor for the SED (see notebook for derivation)
    B = 10 ** (0.4 * M1500 * (1-A))

    # Multiply the SED by this factor
    flux *= B

    # Now compute the new M1500 as a check
    conv = np.sum(flux * tophat_filter) / np.sum(tophat_filter)
    m1500 = flux_to_mag(conv)
    M1500 = m1500 - 5*np.log10(DL/10) + 2.5*np.log10(1+z)

    return wlen, flux



def set_Muv_vectorised(redshifts: np.array, Muv_target: np.array, wlen: np.array, flux: np.array) -> (np.array, np.array):
    """
    Set the absolute magnitude of the LAEs at a given redshift, vectorised for speed.

    Parameters
    ----------
    z : np.array
        Redshifts of the LAEs.
    Muv : np.array
        Absolute magnitudes of the LAEs.
    wlen : np.array

    Returns
    -------
    wlen : np.array
        Wavelength grid of the output SED (same as input).
    flux : np.array
        Flux of the output SED, with the new absolute magnitude.
    """ 

    wlen_list = []
    flux_list = []

    for i, z in enumerate(redshifts):

        lam, sed = set_Muv(z, Muv_target, wlen[i], flux[i])

        wlen_list.append(lam)
        flux_list.append(sed)

    return np.array(wlen_list), np.array(flux_list)



def redshift_SED(z: float, wlen: np.array, flux: np.array) -> (np.array, np.array):
    """
    Redshift the SED to a given redshift, applying Madau (1995) IGM absorption
    
    Parameters
    ----------
    z : float
        Target redshift of the LAE.
    wlen : array
        Wavelength grid of the input SED.
    flux : array
        Flux of the input SED.
    
    Returns
    -------
    wlen : array
        Wavelength grid of the output SED, redshifted to z.
    flux : array
        Flux of the output SED, redshifted to z.
    """

    # Shift the wavelength grid
    wlen = wlen * (1+z)

    # Get luminosity distance
    cosmo = FlatLambdaCDM(H0=70, Om0=0.3)
    DL = cosmo.luminosity_distance(z).value * 10 ** 6 # convert to pc

    # Shift the SED! Convert to magnitude for this step.
    flux = flux_to_mag(flux)
    flux = flux - (5 * np.log10(DL/10)) + 2.5 * np.log10(1+z)
    flux = mag_to_flux(flux)

    # Apply Madau absorption
    igm_abs = etau_madau(wlen, z)

    flux_new = flux * igm_abs.model.lookup_table

    return wlen, flux_new



def add_emission_line(EW: float, z: float, wlen: np.array, flux: np.array) -> tuple[np.array, np.array]:
    """
    Add a Lyman-alpha emission line to the SED.

    Parameters
    ----------
    EW : float
        Equivalent width of the emission line.
    wlen : array
        Wavelength grid of the input SED.
    flux : array
        Flux of the input SED.
    
    Returns
    -------
    wlen : array
        Wavelength grid of the output SED (same as original)
    flux : array
        Flux of the output SED, with the emission line added.
    """

    # EWobs = EW * (1+z), and negative for emission
    EW = EW * (1+z)

    # First measure the continuum flux between 1250-1300A (but observed frame)
    tophat_filter = np.zeros(len(wlen))
    for i, lam in enumerate(wlen):
        if (lam > 1249.0*(1+z)) and (lam < 1301.0*(1+z)):
            tophat_filter[i] = 1.0

    continuum_flux = np.sum(flux * tophat_filter) / np.sum(tophat_filter)

    # Find the index and value where the wavelength is closest to 1215.67A
    idx = np.argmin(np.abs(wlen - 1215.67 * (1+z)))

    # Resolution of the SED
    h = wlen[idx] - wlen[idx-1]

    # Compute the flux density of the SED point at 1215.67A
    line_flux = continuum_flux * ((EW/h)-1)

    # Add the emission line to the SED
    flux[idx] += line_flux
    # Change wlen value to Lyman-alpha wavelength
    wlen[idx] = 1215.67 * (1+z)

    return wlen, flux



def add_emission_line_vectorised(EW: np.array, z: np.array, wlen: np.array, flux: np.array) -> tuple[np.array, np.array]:
    """
    Same as add_emission_line, but vectorised for speed and can take in arrays of EW and z

    Parameters
    ----------
    EW : np.array
        Equivalent width of the emission line.
    z : np.array    
        Redshift of the LAE.
    wlen : np.array
        Wavelength grid of the input SED.
    flux : np.array
        Flux of the input SED.
    
    Returns
    -------
    wlen_new : np.array
        Wavelength grid of the output SED for each combination of z and EW.
    flux_new : np.array
        Flux of the output SED for each combination of z and EW, with the emission line added.
    """
    # Ensure EW and z are 1D arrays
    #EW = np.atleast_1d(EW)
    #z = np.atleast_1d(z)
    
    # Create arrays to store wlen and flux for each combination of z and EW
    wlen_new = np.empty((len(wlen[0, :]), len(z)))
    flux_new = np.empty((len(wlen[0, :]), len(z)))

    for j, redshift in enumerate(z):

        ew = EW[j]

        # Calculate the corresponding wavelength and flux for each combination of z and EW, need correct index of wlen and flux
        wlen_temp, flux_temp = add_emission_line(ew, redshift, wlen[j, :], flux[j, :])
        wlen_new[:, j] = wlen_temp
        flux_new[:, j] = flux_temp
    
    return wlen_new, flux_new



def convolveFilters(filter_set: list[dict], sed: tuple[np.array, np.array], magnitudes=False) -> dict:
    """
    Given a set of filters and SEDs, convolve the SEDs with the filters to get magnitudes.

    Parameters:
    -----------
    filter_set: list[dict]
        A list of dictionaries with the filter names as keys and the wavelength and transmission as values.
        Multiple dictionaries can be passed, the function will stack them.
    sed: tuple[np.array, np.array]
        A tuple with the wavelength grid and flux of the SED.
    magnitudes: bool
        If True, the function will return magnitudes. If False, it will return fluxes.

    Returns:
    --------
    dict{filter_1: flux_1, filter_2: flux_2, ...}
        A dictionary with the filter names as keys and the magnitudes as values.
    """

    # If there are multiple dicts in filter_set, turn it into one big dict
    filters = {}
    for filter_dict in filter_set:
        filters.update(filter_dict)

    convolved_fluxes = {}

    # Loop through the filters
    for filter_name, (filter_wlen_grid, transmission) in filters.items():

        # Interpolate the SED grid to the filter
        flux_interp = np.interp(filter_wlen_grid, sed[0], sed[1])

        # Convert to frquency space (f = c/wlen)
        filter_freq_grid = np.array([c.value / (wlen*1e-10) for wlen in filter_wlen_grid])

        # Calculate the convolved flux, normalised by filter area under curve
        convolved_flux = simps(flux_interp * transmission, filter_freq_grid) / simps(transmission, filter_freq_grid)

        convolved_fluxes[filter_name] = convolved_flux

        if magnitudes:
            convolved_fluxes = flux_to_mag(convolved_flux)
             
    return convolved_fluxes



def convolveFilters_vectorised(filter_set: list[dict], sed: np.array, magnitudes=False) -> dict:
    """
    same as convolveFilters, but vectorised for speed

    Parameters:
    -----------
    filter_set: list[dict]
        A list of dictionaries with the filter names as keys and the wavelength and transmission as values.
        Multiple dictionaries can be passed, the function will stack them.
    sed: np.array[np.array, np.array]
        Array of SEDs. First dimension is all the combos, second is the wavelength grid, third is the flux grid.
    magnitudes: bool
        If True, the function will return magnitudes. If False, it will return fluxes.
    
    Returns:
    --------
    list[dict{filter_1: flux_1, filter_2: flux_2, ...}]
        A list of dictionaries with the filter names as keys and the magnitudes as values.
    """

    dict_list = []

    # sed has shape (2, 1296, 15), iterate through last axis and get the elements of first axis
    for i in range(sed.shape[-1]):

        wlen = sed[0, :, i]
        flux = sed[1, :, i]

        convolved_fluxes = convolveFilters(filter_set, (wlen.flatten(), flux.flatten()), magnitudes)

        dict_list.append(convolved_fluxes)
    
    return dict_list



def simulateDepths(filter_set: list[dict]) -> tuple[dict, dict]:

    """
    Given a list of filters, select a random local depth from its depth table.

    Parameters:
    -----------
    filter_set: list[dict]
        A list of dictionaries with the filter names as keys and the wavelength and transmission as values.
        Multiple dictionaries can be passed, the function will stack them.
    
    Returns:
    --------
    tuple[dict{filter_1: mag_depth_1, filter_2: mag_depth_2, ...}, dict{filter_1: flux_depth_1, filter_2: flux_depth_2, ...}]
        A tuple with two dictionaries, one with the magnitudes and one with the flux depths.
    """

    # Set random seed
    np.random.seed(42)

    # If there are multiple dicts in filter_set, turn it into one big dict
    filters = {}
    for filter_dict in filter_set:
        filters.update(filter_dict)

    mag_depths = {}
    flux_depths = {}

    # Loop through the filters
    for filter_name, (filter_wlen_grid, transmission) in filters.items():

        # Set up different directories based on filter name
        if filter_name in ['Ye', 'Je', 'He']:
            depth_dir = Path.cwd().parent.parent / 'data' / 'depths' / 'COSMOS' / 'phot'
            depth_table = filter_name.split('e')[0] + '_1.2as_gridDepths_300_200.fits'

        elif filter_name == 'VIS':
            depth_dir = Path.cwd().parent.parent / 'data' / 'depths' / 'COSMOS' / 'phot'
            depth_table = filter_name + '_0.3as_gridDepths_300_200.fits'

        elif filter_name in ['Y', 'J', 'H', 'Ks']:
            depth_dir = Path.home().parent.parent/ 'vardy' / 'vardygroupshare' / 'data' / 'depths' / 'COSMOS' / 'phot'
            depth_table = filter_name.split('s')[0] + '_DR6_1.8as_gridDepths_300_200.fits'

        elif filter_name in ['g', 'r', 'i', 'nb816', 'z', 'nb0921', 'y']:

            name_dict = {'g': 'HSC-G_DR3', 'r': 'HSC-R_DR3', 'i': 'HSC-I_DR3', 'nb816': 'HSC-NB0816_DR3', 'z': 'HSC-Z_DR3', 'nb921': 'HSC-NB0921_DR3', 'y': 'HSC-Y_DR3'}

            depth_dir = Path.home().parent.parent / 'vardy' / 'vardygroupshare' / 'data' / 'depths' / 'COSMOS' / 'phot'
            depth_table = name_dict[filter_name] + '_1.8as_gridDepths_300_200.fits'

        # Open the depth table
        t = Table.read(depth_dir / depth_table)

        # Ensure mask == 1
        t = t[t['mask'] == 1]

        # Select a random depth
        idx = np.random.randint(0, len(t))
        mag_depth = t['depths'][idx]

        flux_depth = mag_to_flux(mag_depth)

        mag_depths[filter_name] = mag_depth
        flux_depths[filter_name] = flux_depth

    return mag_depths, flux_depths



def scatterFluxes(fluxes: dict, depths: dict, magnitude=False) -> dict:
    """
    Given a set of fluxes and depths, scatter the fluxes by the depth values.

    Parameters:
    -----------
    fluxes: dict{filter_1: flux_1, filter_2: flux_2, ...}
        A dictionary with the filter names as keys and the fluxes as values.
    depths: dict{filter_1: mag_depth_1, filter_2: mag_depth_2, ...}
        A dictionary with the filter names as keys and the depths as values.
    magnitude: bool
        If True, the function will return magnitudes. If False, it will return fluxes.
    
    Returns:
    --------
    dict{filter_1: scattered_flux_1, filter_2: scattered_flux_2, ...}
        A dictionary with the filter names as keys and the scattered fluxes as values.
    """

    scattered_fluxes = {}

    for filter_name, flux in fluxes.items():
        depth = depths[filter_name]

        # Set random seed
        np.random.seed(42)

        # Scatter the fluxes
        fluxes[filter_name] = np.random.normal(flux, depth)

        sigma = 0.2 * 10 ** (-0.4 * (48.6 + depth)) * 5

        scatter_value = np.random.normal(0, sigma)

        scattered_flux = flux + scatter_value

        scattered_fluxes[filter_name] = scattered_flux

        if magnitude:
            scattered_fluxes[filter_name] = flux_to_mag(scattered_flux)

    return scattered_fluxes



def scatterFluxes_vectorised(fluxes: list[dict], depths: dict, magnitude=False) -> list[dict]:
    """
    Given a set of fluxes and depths, scatter the fluxes by the depth values.

    Parameters:
    -----------
    fluxes: list[dict{filter_1: flux_1, filter_2: flux_2, ...}]
        A list of dictionaries with the filter names as keys and the fluxes as values.
    depths: dict{filter_1: mag_depth_1, filter_2: mag_depth_2, ...}
        A dictionary with the filter names as keys and the depths as values.
    magnitude: bool
        If True, the function will return magnitudes. If False, it will return fluxes.
    
    Returns:
    --------
    list[dict{filter_1: scattered_flux_1, filter_2: scattered_flux_2, ...}]
        A list of dictionaries with the filter names as keys and the scattered fluxes as values.
    """

    dict_list = []

    for i, flux_dict in enumerate(fluxes):

        scattered_fluxes = scatterFluxes(flux_dict, depths, magnitude=magnitude)

        dict_list.append(scattered_fluxes)

    return dict_list



def getErrors(depths: dict) -> dict:
    """
    Given a set of depths, compute the errors on the fluxes.

    Parameters:
    -----------
    depths: dict{filter_1: mag_depth_1, filter_2: mag_depth_2, ...}
        A dictionary with the filter names as keys and the depths as values.

    Returns:
    --------
    dict{filter_1: error_1, filter_2: error_2, ...}
        A dictionary with the filter names as keys and the errors as values.
    """

    errors = {}
    for filter_name, depth in depths.items():

        exponent = - (48.6 + depth) / 2.5
        error = 0.2 * (10 ** exponent)
        errors[filter_name] = error

    return errors



def getSignalToNoise(fluxes: dict, errors: dict) -> dict:
    """
    Given a set of fluxes and errors, compute the signal-to-noise ratio.

    Parameters:
    -----------
    fluxes: dict{filter_1: flux_1, filter_2: flux_2, ...}
        A dictionary with the filter names as keys and the fluxes as values.
    errors: dict{filter_1: error_1, filter_2: error_2, ...}
        A dictionary with the filter names as keys and the errors as values.

    Returns:
    --------
    dict{filter_1: snr_1, filter_2: snr_2, ...}
        A dictionary with the filter names as keys and the signal-to-noise ratios as values.
    """

    signal_to_noise = {}
    for filter_name, flux in fluxes.items():
        error = errors[filter_name]

        signal_to_noise[filter_name] = flux / error

    return signal_to_noise



def filterCentreAndWidth(filter_name: str, instrument: str) -> tuple[float, float]:
    """
    Get the central wavelength and width of a filter.

    Parameters
    ----------
    filter_names : list
        List of filter names.
    instrument : str
        Name of the filter instrument.

    Returns
    -------

    tuple[float, float]
        Tuple of central wavelength and width.
    """

    if instrument.lower() == 'euclid':

        t = ascii.read(Path.cwd().parent / 'validation' / 'euclid_filters.txt')

        centre = t[t['filter'] == filter_name]['centre'][0]
        upper_edge = t[t['filter'] == filter_name]['upper_edge'][0]
        lower_edge = t[t['filter'] == filter_name]['lower_edge'][0]
        width = upper_edge - lower_edge

    if instrument.lower() == 'vista':
            
        t = ascii.read(Path.cwd().parent / 'validation' / 'vista_filters.txt')

        centre = t[t['filter'] == filter_name]['centre'][0]
        width = t[t['filter'] == filter_name]['width'][0]

    if instrument.lower() == 'hsc':
                
            t = ascii.read(Path.cwd().parent / 'validation' / 'hsc_filters.txt')
    
            centre = t[t['filter'] == filter_name]['centre'][0]
            width = t[t['filter'] == filter_name]['width'][0]

    return centre, width



def update_plot(redshift, Muv, EW):

    """
    Update the plot with the new redshift and EW values, in order to make an animation! This is at fixed muv
    """

    euclid_filters = getFilters('Euclid')
    vista_filters = getFilters('VISTA')
    hsc_filters = getFilters('HSC')

    plt.clf()  # Clear the previous plot
    plt.rcParams['figure.dpi'] = 75

    wlen, flux_sed = makeLBG(redshift=redshift, SFH_component='constant', age=(0, 13.8), massformed=11., metallicity=0.2, 
                                dust_type='Calzetti', Av=0.2, nebular=True, logU=-1.)
    wlen, flux_sed = set_Muv(z=redshift, Muv_target=Muv, wlen=wlen, flux=flux_sed)
    wlen, flux_sed = add_emission_line(EW=EW, z=redshift, wlen=wlen, flux=flux_sed)
    fluxes = convolveFilters([euclid_filters, vista_filters, hsc_filters], (wlen, flux_sed), magnitudes=False)
    depths, _ = simulateDepths([euclid_filters, vista_filters, hsc_filters])
    errors = getErrors(depths)
    signal_to_noise = getSignalToNoise(fluxes, errors)
    scattered_fluxes = scatterFluxes(fluxes, depths)
    
    # Plotting
    plt.plot(wlen, flux_sed, color='deepskyblue', lw=2.5, alpha=0.8)

    # Get the centre of the Euclid, VISTA and HSC filters and then put into one big dictionary
    filter_centres = {}
    filter_widths = {}
    for filter_name in euclid_filters.keys():
        filter_key = filter_name.split('e')[0]
        centre, width = filterCentreAndWidth(filter_key, 'Euclid')
        filter_centres[filter_name] = centre
        filter_widths[filter_name] = width

    for filter_name in vista_filters.keys():
        centre, width = filterCentreAndWidth(filter_name, 'VISTA')
        filter_centres[filter_name] = centre
        filter_widths[filter_name] = width

    for filter_name in hsc_filters.keys():
        filter_keys = {'g': 'HSC-G_DR3', 'r': 'HSC-R_DR3', 'i': 'HSC-I_DR3', 'nb816': 'HSC-NB0816_DR3', 'z': 'HSC-Z_DR3', 'nb921': 'HSC-NB0921_DR3', 'y': 'HSC-Y_DR3'}
        centre, width = filterCentreAndWidth(filter_keys[filter_name], 'HSC')
        filter_centres[filter_name] = centre
        filter_widths[filter_name] = width

    # convert filter centres to angstroms from microns
    filter_centres = {k: v*1e4 for k, v in filter_centres.items()}
    filter_widths = {k: v*1e4 for k, v in filter_widths.items()}

    # If signal to noise is larger than 2, plot the fluxes as errorbars
    for filter_name, flux in scattered_fluxes.items():

        # Set colour based on filter name
        if filter_name in ['VIS', 'Ye', 'Je', 'He']:
            colour = 'red'
            fmt = 'o'
        elif filter_name in ['Y', 'J', 'H', 'Ks']:
            colour = 'blue'
            fmt = 's'
        elif filter_name in ['g', 'r', 'i', 'nb816', 'z', 'nb921', 'y']:
            colour = 'black'
            fmt='D'

        if signal_to_noise[filter_name] > 2:
            plt.errorbar(filter_centres[filter_name], flux, yerr=errors[filter_name], xerr=filter_widths[filter_name]/2, fmt=fmt, color=colour)

        # Else, plot a 2sigma upper limit
        else:
            plt.errorbar(filter_centres[filter_name], flux+2*errors[filter_name], fmt='v', uplims=True, color=colour)

    # Dummy plots for labels for each instrument
    plt.errorbar(0, 1, yerr=0, xerr=0, fmt='o', label='Euclid', color='red')
    plt.errorbar(0, 1, yerr=0, xerr=0, fmt='s', label='VISTA', color='blue')
    plt.errorbar(0, 1, yerr=0, xerr=0, fmt='D', label='HSC', color='black')

    plt.xlabel(r'Wavelength ($\AA$)')
    plt.ylabel('Flux')
    plt.legend()
    plt.xlim(5000, 30000)
    plt.yscale('log')
    plt.ylim(3e-32, 1e-27)
    plt.tight_layout()

    plt.title(r'$z = $' + f'{round(redshift, 2)}')



def runLymanAlphaEmitter_vectorised(redshift, Muv, EW) -> dict:

    """
    Run the Lyman-alpha emitter simulation.

    Parameters:
    -----------
    redshift: float
        Redshift of the LAE.
    Muv: float
        Absolute magnitude of the LAE.
    EW: float
        Equivalent width of the emission line.
    
    Returns:
    --------
    dict{filter_1: mag_flux_1, filter_2: mag_flux_2, ...}
        A dictionary with the filter names as keys and the magnitudes as values.

    """

    # Get filter sets
    euclid_filters = getFilters('Euclid')
    vista_filters = getFilters('VISTA')
    hsc_filters = getFilters('HSC')

    # Make a Lyman-break galaxy
    wlen, flux_sed = makeLBG_vectorised(redshifts=redshift, SFH_component='constant', age=(0, 13.8), massformed=11., metallicity=0.2, 
                                dust_type='Calzetti', Av=0.2, nebular=True, logU=-1.)

    # Set the absolute magnitude
    wlen, flux_sed = set_Muv_vectorised(redshifts=redshift, Muv_target=Muv, wlen=wlen, flux=flux_sed)

    # Add the Lyman-alpha emission line
    wlen, flux_sed = add_emission_line_vectorised(EW=EW, z=redshift, wlen=wlen, flux=flux_sed)

    sed = np.array([wlen, flux_sed])

    # Compute the flux through the filters
    fluxes = convolveFilters_vectorised([euclid_filters, vista_filters, hsc_filters], sed, magnitudes=False)

    # Pull random depths from the depth table
    depths, _ = simulateDepths([euclid_filters, vista_filters, hsc_filters])

    # Scatter the photometry within the errors
    scattered_fluxes = scatterFluxes_vectorised(fluxes, depths, magnitude=False)

    # get the errors
    errors = getErrors(depths)

    return scattered_fluxes, errors



def runLymanAlphaEmitter(redshift, Muv, EW) -> dict:

    """
    Run the Lyman-alpha emitter simulation.

    Parameters:
    -----------
    redshift: float
        Redshift of the LAE.
    Muv: float
        Absolute magnitude of the LAE.
    EW: float
        Equivalent width of the emission line.
    
    Returns:
    --------
    dict{filter_1: mag_flux_1, filter_2: mag_flux_2, ...}
        A dictionary with the filter names as keys and the magnitudes as values.

    """

    # Get filter sets
    euclid_filters = getFilters('Euclid')
    vista_filters = getFilters('VISTA')
    hsc_filters = getFilters('HSC')

    # Make a Lyman-break galaxy
    wlen, flux_sed = makeLBG(redshift=redshift, SFH_component='constant', age=(0, 13.8), massformed=11., metallicity=0.2, 
                                dust_type='Calzetti', Av=0.2, nebular=True, logU=-1.)

    # Set the absolute magnitude
    wlen, flux_sed = set_Muv(z=redshift, Muv_target=Muv, wlen=wlen, flux=flux_sed)

    # Add the Lyman-alpha emission line
    wlen, flux_sed = add_emission_line(EW=EW, z=redshift, wlen=wlen, flux=flux_sed)

    # Compute the flux through the filters
    fluxes = convolveFilters([euclid_filters, vista_filters, hsc_filters], (wlen, flux_sed), magnitudes=False)

    # Pull random depths from the depth table
    depths, _ = simulateDepths([euclid_filters, vista_filters, hsc_filters])

    # Scatter the photometry within the errors
    scattered_fluxes = scatterFluxes(fluxes, depths, magnitude=True)

    return scattered_fluxes













