#!/usr/bin/env python3

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
import bagpipes as pipes
from typing import Union, Tuple
from astropy.cosmology import FlatLambdaCDM
import astropy.units as u
import re
from collections import defaultdict
from astropy.table import Table

plt.rcParams['axes.linewidth'] = 2.5
plt.rcParams.update({'font.size': 15})
plt.rcParams['figure.dpi'] = 100

# Set up directories
filter_dir = Path().home() / 'lephare' / 'lephare_dev' / 'filt' / 'myfilters'
plot_dir = Path.cwd().parent.parent / 'plots' / 'brown_dwarfs'
dwarf_dir = Path().home() / 'lephare' / 'lephare_dev' / 'sed' / 'STAR' / 'DWARFSTARS'
spectra_dir = Path.cwd().parent.parent / 'data' / 'SEDs' / 'dwarfs'
table_dir = Path.cwd().parent.parent / 'data' / 'simulations' / 'tables'
rebels_dir = Path.cwd().parent.parent / 'data' / 'ref_catalogues' / 'fluxes'



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
    # # Deal with negative flux case
    # if flux <= 0:
    #     mag = np.nan
    # else:
    mag = -2.5*np.log10(flux)-48.6
    return mag



def flat_spectrum(wlen, target_flux):
    """
    Given a wavelength grid and a flux grid, return a flat spectrum.
    """
    return np.ones(len(wlen)) * target_flux



def powerLawSpectrum(wlen, target_flux, power, redshift):
    """
    Given a wavelength grid, a target flux and a power, return a power law spectrum.
    Then given a redshift, add the Lyman break by setting the flux to zero bluewards of (1+z)*1216 Angstroms
    """

    # Make the power law spectrum
    flux = target_flux * (wlen) ** power

    # Add the Lyman break
    flux[wlen < 1216 * (1 + redshift)] = 0

    return flux



def max_age_at_redshift(redshift: float) -> float:

    cosmo = FlatLambdaCDM(H0=70, Om0=0.3) 
    age_at_z = cosmo.age(redshift).to(u.Gyr).value

    return age_at_z



def getFilters(instrument: str, plot: bool = False, plot_kwargs: dict = None) -> dict:
    """
    By specifying the instrument name, retrieve/plot all of its photometric filters.

    Based on src/validation/filter_transmission_curves.py.

    Parameters:
    -----------
    instrument: str
        The name of the instrument. One of 'Euclid', 'VISTA', 'HSC', 'Spitzer', 'JWST.
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
    elif instrument.lower() == 'jwst':
        filter_names = glob.glob(str(filter_dir / 'JWST' / '*'))
        filter_names = [f for f in filter_names if 'f115w_angstroms' in f or 'f150w_angstroms' in f or 'f277w_angstroms' in f or 'f444w_angstroms' in f]
        order = ['f115w', 'f150w', 'f277w', 'f444w']
    else:
        raise ValueError("Invalid instrument name")

    # Sort the filters
    if instrument.lower() == 'euclid' or instrument.lower() == 'spitzer' or instrument.lower() == 'vista':
        filter_names = sorted(filter_names, key=lambda x: order.index(x.split('/')[-1].split('_')[-1].split('.')[0]))
    if instrument.lower() == 'jwst':
        filter_names = sorted(filter_names, key=lambda x: order.index(x.split('/')[-1].split('_')[0]))
    if instrument.lower() == 'hsc':  
        filter_names = sorted(filter_names, key=lambda x: order.index(x.split('/')[-1].split('_HSC.txt')[0]))

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

            # normalize to ~ 1e-29 erg/s/cm^2/A, roughly 24-26 AB mag
            flux /= np.max(flux)
            flux *= 1e-29

            brown_dwarfs[f'{spectral_type}{sub_type}'] = (wavelength, flux)

    return brown_dwarfs



def loadBrownDwarfSpectra() -> dict:

    """
    Load brown dwarf spectra taken from https://cass.ucsd.edu/~ajb/browndwarfs/spexprism/library.html.

    The commented header contains the reference for the object and the spectral type, which is extracted.

    Parameters:
    -----------
    None

    Returns:
    --------
    dict
        A dictionary with the spectral type as the key and a list of wavelength and flux tuples as the values.
    """

    # Collect the brown dwarf spectra file names
    spectra_files = glob.glob(str(spectra_dir / '*'))

    # Use defaultdict to automatically create a list for each spectral type
    brown_dwarfs = defaultdict(list)

    for file in spectra_files:

        # Get the spectral type from the file name
        spectral_type = file.split('/')[-1].split('_')[0]

        # Read data, commented out lines are skipped
        data = ascii.read(file)

        # Find the line containing "Near infrared spectral type:", and if not, search for "Optical spectral type:"
        for line in data.meta['comments']:
            if 'Near infrared spectral type:' in line:
                spectral_type = line.split(':')[-1].strip()
                break
            elif 'Optical spectral type:' in line:
                spectral_type = line.split(':')[-1].strip()
                break
            else:
                spectral_type = 'Unknown'

        wavelength = np.array(data['col1'])
        flux = np.array(data['col2'])

        # Convert wavelength from microns to angstroms
        wavelength *= 1e4

        # Normalise to ~1e-29 ergs, around 24-26 AB mag. Only interested in relative colours.
        flux /= np.max(flux)
        flux *= 1e-29

        brown_dwarfs[spectral_type].append((wavelength, flux))

    # Convert defaultdict to regular dictionary
    brown_dwarfs = dict(brown_dwarfs)

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

            # Convert to frquency space (f = c/wlen)
            filter_freq_grid = np.array([c.value / (wlen*1e-10) for wlen in filter_wlen_grid])

            # Calculate the convolved flux, normalised by filter area under curve
            convolved_flux = simps(flux_interp * transmission, filter_freq_grid) / simps(transmission, filter_freq_grid)

            mag = flux_to_mag(convolved_flux)

            # Store with keys for each BD corresponding to a dict with a key for each filter
            convolved_fluxes.setdefault(spectral_type, {})[filter_name] = mag
                
    return convolved_fluxes



def makeLBG(redshift: float, SFH_component: str, Muv: float = None,
            age: Union[float, Tuple[float, float]] = None, tau: float = None, 
            tmax: float = None, fwhm: float = None,
            massformed: float = None, metallicity: float = None,
            dust_type: str = None, Av: float = None,
            nebular: bool=False, logU: float=None,
            filter_list: Path = Path.cwd() / 'filter_list.txt') -> pipes.model_galaxy:
    """
    Create a model Lyman-break galaxy using BAGPIPES. Can specify different star formation histories.
    Most arguments default to None, but certain ones are required for certain SFHs. These are:
        burst: age.
        constant: age = (age_min, age_max).
        exponential: age, tau.
        delayed: age, tau.
        lognormal: tmax, fwhm.

    Parameters:
    -----------
    redshift: float
        The redshift of the galaxy.
    SFH_component: str
        The star formation history component to use. One of: 
        'burst', 'constant', 'exponential' (e^-(t/tau)), 'delayed' (t*e^-(t/tau)), 'lognormal'.
    Muv: float, optional
        The target absolute rest-UV magnitude of the galaxy. SED is scaled to this if provided.
    age: Union[float, Tuple[float, float]], optional
        The age of the galaxy in Gyr. If SFH is 'constant', provide a tuple (age_min, age_max) (= time since SF turned off/on).
    tau: float, optional
        The e-folding time in Gyr. Required for 'exponential' and 'delayed' SFHs.
    tmax: float, optional
        The time of maximum star formation in Gyr. Required for 'lognormal' SFH.
    fwhm: float, optional
        The full width at half maximum of the lognormal distribution in Gyr.
    massformed: float, optional
        The mass of the galaxy formed in log_10(M*/M_solar)
    metallicity: float, optional
        The metallicity of the galaxy in Z/Z_solar.
    dust: str, optional
        The dust attenuation law to use.
    Av: float, optional
        The dust attenuation in magnitudes.
    nebular: bool, optional
        Whether to include nebular emission.
    logU: float, optional
        The ionization parameter.
    filter_list: Path, optional
        The path to the filter list file. Default is in this directory, 'filter_list.txt'.

    Returns:
    --------
    pipes.model_galaxy
        The model galaxy object. 
    """

    # Load filter list
    filter_list = np.loadtxt(filter_list, dtype="str")

    # Ensure age or age_max is no more than age of the Universe at the given redshift
    if age is not None:
        if SFH_component == 'constant':
            if age[1] > max_age_at_redshift(redshift):
                age = (age[0], max_age_at_redshift(redshift))
        elif age > max_age_at_redshift(redshift):
            age = max_age_at_redshift(redshift)

    # Create mandatory components of star formation history
    sfh = {}
    if SFH_component == 'burst':
        sfh['age'] = age
    elif SFH_component == 'constant':
        sfh['age_min'] = age[0]
        sfh['age_max'] = age[1]
    elif SFH_component == 'exponential':
        sfh['age'] = age
        sfh['tau'] = tau
    elif SFH_component == 'delayed':
        sfh['age'] = age
        sfh['tau'] = tau
    elif SFH_component == 'lognormal':
        sfh['tmax'] = tmax
        sfh['fwhm'] = fwhm

    # Add optional components
    if massformed is not None:
        sfh['massformed'] = massformed
    if metallicity is not None:
        sfh['metallicity'] = metallicity
    
    # Dust component
    dust = {}
    if dust is not None:
        dust['type'] = dust_type
        dust['Av'] = Av
    
    # Nebular component
    nebular = {}
    if nebular:
        nebular['logU'] = logU

    # Create model component dictionary
    model_components = {}
    model_components['redshift'] = redshift
    if dust is not None:
        model_components['dust'] = dust
    if nebular:
        model_components['nebular'] = nebular

    # Add SFH
    if SFH_component == 'burst':
        model_components['burst'] = sfh
    elif SFH_component == 'constant':
        model_components['constant'] = sfh
    elif SFH_component == 'exponential':
        model_components['exponential'] = sfh
    elif SFH_component == 'delayed':
        model_components['delayed'] = sfh
    elif SFH_component == 'lognormal':
        model_components['lognormal'] = sfh


    model = pipes.model_galaxy(model_components, filt_list=filter_list)

    wlen = model.wavelengths
    flux = model.spectrum_full

    # Shift wavelength to observed frame
    wlen *= (1 + redshift)

    # Scale to Muv
    if Muv is not None:

        # Set up a cosmology
        cosmo = FlatLambdaCDM(H0=70, Om0=0.3)

        # Need to convert flux from f_lambda to f_nu. # Conversion from Nathan, then into cgs
        flux = flux * (wlen**2)/(10**-29*2.9979*10**18) * 1e-19
    
        # Compute location of Lyman break
        lyman_break = 1216 * (1+redshift) # Angstroms

        # Get the distance modulus from the redshift.
        DL = cosmo.luminosity_distance(redshift).value * 10 ** 6 # in pc

        # Compute apparent magnitude from absolute magnitude.
        m = Muv + 5*np.log10(DL/10) - 2.5*np.log10(1+redshift)

        # Compute the observed wavelength of 1500A
        uv_obs = 1500 * (1+redshift) # Angstroms

        # Compute flux corresponding to apparent uv magnitude.
        app_flux = mag_to_flux(m)

        # Compute the ratio between this and the BAGPIPES flux at 1500 A
        idx = np.abs(wlen - uv_obs).argmin()
        ratio = app_flux / flux[idx]

        # Match the model to our desired magnitude at 1500A
        flux = flux * ratio

        #model.line_fluxes["H  1  1215.67A"] = 1e-29


    return wlen, flux



#? --------------------------------------------------------------------------
#? Uncomment this to do some plotting checks on top of filter response curves
#__name__ = '__none__'
#? --------------------------------------------------------------------------



#####################################################################################
if __name__ == '__main__':

    #! ------------------------------------------------------
    #! Generate a bunch of brown dwarf colours, save to table
    #! ------------------------------------------------------

    # #bds = loadBrownDwarfTemplates()
    # bds = loadBrownDwarfSpectra()

    # # Make empty lists for all the VISTA and Euclid filters

    # spec_type = []

    # z_mags_BD = []
    # y_mags_BD = []

    # Y_mags_BD = []
    # J_mags_BD = []
    # H_mags_BD = []
    # Ks_mags_BD = []

    # Ye_mags_BD = []
    # Je_mags_BD = []
    # He_mags_BD = []

    # # Plot BD templates
    # for spectral_type, (spectra_list) in bds.items(): ###! Spectral sample
    # #for spectral_type, (wavelength, flux) in bds.items(): ####! SpeX templates
        
    # #! Uncomment and tab this loop if reading the spectra
    #     for spectrum in spectra_list:
    #         wavelength, flux = spectrum

    #         spec_type.append(spectral_type)

    #         # Filters
    #         euclid_filters = getFilters('euclid')
    #         vista_filters = getFilters('vista')
    #         hsc_filters = getFilters('hsc')

    #         mags = convolveFilters([hsc_filters, vista_filters, euclid_filters], {spectral_type: (wavelength, flux)})

    #         z_mags_BD.append(mags[spectral_type]['z'])
    #         y_mags_BD.append(mags[spectral_type]['y'])

    #         Y_mags_BD.append(mags[spectral_type]['Y'])
    #         J_mags_BD.append(mags[spectral_type]['J'])
    #         H_mags_BD.append(mags[spectral_type]['H'])
    #         Ks_mags_BD.append(mags[spectral_type]['Ks'])

    #         Ye_mags_BD.append(mags[spectral_type]['Ye'])
    #         Je_mags_BD.append(mags[spectral_type]['Je'])
    #         He_mags_BD.append(mags[spectral_type]['He'])

    # # Convert all the mags lists into numpy arrays
    # z_mags_BD = np.array(z_mags_BD)
    # y_mags_BD = np.array(y_mags_BD)

    # Ye_mags_BD = np.array(Ye_mags_BD)
    # Je_mags_BD = np.array(Je_mags_BD)
    # He_mags_BD = np.array(He_mags_BD)

    # Y_mags_BD = np.array(Y_mags_BD)
    # J_mags_BD = np.array(J_mags_BD)
    # H_mags_BD = np.array(H_mags_BD)
    # Ks_mags_BD = np.array(Ks_mags_BD)

    # flux_table = Table([z_mags_BD, y_mags_BD, Y_mags_BD, J_mags_BD, H_mags_BD, Ks_mags_BD, Ye_mags_BD, Je_mags_BD, He_mags_BD], names=['z', 'y', 'Y', 'J', 'H', 'Ks', 'Ye', 'Je', 'He'])
    # flux_table.add_column(Table.Column(name='Spectral Type', data=spec_type), index=0)
    # flux_table.write(table_dir / 'dwarf_spectra_mags.fits', overwrite=True)
    # exit()

    #! --------------------------------------------------
    #!     Generate a bunch of LBG colours, save to table
    #! --------------------------------------------------

    redshifts = np.arange(4., 10., 0.05)
    Av_vals = np.arange(0, 0.6, 0.1)
    ages = np.arange(0.05, 0.51, 0.05)
    

    # Make empty lists for all the VISTA and Euclid filters
    g_mags_LBG = []
    r_mags_LBG = []
    i_mags_LBG = []
    nb816_mags_LBG = []
    z_mags_LBG = []
    nb921_mags_LBG = []
    y_mags_LBG = []

    Y_mags_LBG = []
    J_mags_LBG = []
    H_mags_LBG = []
    Ks_mags_LBG = []

    VIS_mags_LBG = []
    Ye_mags_LBG = []
    Je_mags_LBG = []
    He_mags_LBG = []

    f115w_mags_LBG = []
    f150w_mags_LBG = []
    f277w_mags_LBG = []
    f444w_mags_LBG = []

    redshift_array = []
    Av_array = []
    age_array = []


    euclid_filters = getFilters('euclid')
    vista_filters = getFilters('vista')
    hsc_filters = getFilters('hsc')
    jwst_filters = getFilters('jwst')

    for redshift in redshifts:
        for Av in Av_vals:
            for age in ages:
                print(f'Generating LBG at z={redshift:.2f} for Av={Av:.2f} and age={age:.2f} Gyr.')

                redshift_array.append(redshift)
                Av_array.append(Av)
                age_array.append(age)

                wlen, LBG_flux = makeLBG(redshift=redshift, SFH_component='constant', age=(0, age), massformed=10., metallicity=0.2, 
                                dust_type='Calzetti', Av=Av, nebular=True, logU=-2., Muv=-22.)
                
                #LBG_flux = powerLawSpectrum(wlen, 1e-30, -2., redshift=redshift)
                
                # Make a similar BD dictionary
                lbg_dict = {redshift: (wlen, LBG_flux)}

                # Convolve
                mags = convolveFilters([hsc_filters, vista_filters, euclid_filters, jwst_filters], lbg_dict)

                g_mags_LBG.append(mags[redshift]['g'])
                r_mags_LBG.append(mags[redshift]['r'])
                i_mags_LBG.append(mags[redshift]['i'])
                nb816_mags_LBG.append(mags[redshift]['nb816'])
                z_mags_LBG.append(mags[redshift]['z'])
                nb921_mags_LBG.append(mags[redshift]['nb921'])
                y_mags_LBG.append(mags[redshift]['y'])

                Y_mags_LBG.append(mags[redshift]['Y'])
                J_mags_LBG.append(mags[redshift]['J'])
                H_mags_LBG.append(mags[redshift]['H'])
                Ks_mags_LBG.append(mags[redshift]['Ks'])

                VIS_mags_LBG.append(mags[redshift]['VIS'])
                Je_mags_LBG.append(mags[redshift]['Je'])
                Ye_mags_LBG.append(mags[redshift]['Ye'])
                He_mags_LBG.append(mags[redshift]['He'])

                f115w_mags_LBG.append(mags[redshift]['f115w'])
                f150w_mags_LBG.append(mags[redshift]['f150w'])
                f277w_mags_LBG.append(mags[redshift]['f277w'])
                f444w_mags_LBG.append(mags[redshift]['f444w'])

    # Convert to arrays
    g_mags_LBG = np.array(g_mags_LBG)
    r_mags_LBG = np.array(r_mags_LBG)
    i_mags_LBG = np.array(i_mags_LBG)
    nb816_mags_LBG = np.array(nb816_mags_LBG)
    z_mags_LBG = np.array(z_mags_LBG)
    nb921_mags_LBG = np.array(nb921_mags_LBG)
    y_mags_LBG = np.array(y_mags_LBG)

    Y_mags_LBG = np.array(Y_mags_LBG)
    J_mags_LBG = np.array(J_mags_LBG)
    H_mags_LBG = np.array(H_mags_LBG)
    Ks_mags_LBG = np.array(Ks_mags_LBG)

    VIS_mags_LBG = np.array(VIS_mags_LBG)
    Ye_mags_LBG = np.array(Ye_mags_LBG)
    Je_mags_LBG = np.array(Je_mags_LBG)
    He_mags_LBG = np.array(He_mags_LBG)

    f115w_mags_LBG = np.array(f115w_mags_LBG)
    f150w_mags_LBG = np.array(f150w_mags_LBG)
    f277w_mags_LBG = np.array(f277w_mags_LBG)
    f444w_mags_LBG = np.array(f444w_mags_LBG)

    all_mags_list = [g_mags_LBG, r_mags_LBG, i_mags_LBG, nb816_mags_LBG, z_mags_LBG, nb921_mags_LBG, y_mags_LBG, Y_mags_LBG, J_mags_LBG, H_mags_LBG, Ks_mags_LBG, VIS_mags_LBG, Ye_mags_LBG, Je_mags_LBG, He_mags_LBG, f115w_mags_LBG, f150w_mags_LBG, f277w_mags_LBG, f444w_mags_LBG]
    names = ['g', 'r', 'i', 'nb816', 'z', 'nb921', 'y', 'Y', 'J', 'H', 'Ks', 'VIS', 'Ye', 'Je', 'He', 'f115w', 'f150w', 'f277w', 'f444w']

    flux_table = Table(all_mags_list, names=names)
    flux_table.add_column(Table.Column(name='Redshift', data=redshift_array), index=0)
    flux_table.add_column(Table.Column(name='Av', data=Av_array), index=1)
    flux_table.add_column(Table.Column(name='Age', data=age_array), index=2)
    flux_table.write(table_dir / 'lbg_spectra_mags.fits', overwrite=True)
    exit()

    #! ########################################################################
    #! ------------------------------------------------------------------------

    #! PLOT COLOURS OF BROWN DWARFS AND LBGs

    #! ------------------------------------------------------------------------
    #! ########################################################################

    # Define the list of magnitudes to extract
    magnitudes = ['z', 'y', 'Y', 'J', 'H', 'Ks', 'Ye', 'Je', 'He']

    # Read in the tables
    bd_table = Table.read(table_dir / 'dwarf_spectra_mags.fits') #! Spectra
    spex_table = Table.read(table_dir / 'spex_template_mags.fits') #! Templates
    lbg_table = Table.read(table_dir / 'lbg_spectra_mags.fits')

    # Restrict lbg table to certain redshift range
    #lbg_table = lbg_table[lbg_table['Redshift'] < 7.55]
    #lbg_table = lbg_table[lbg_table['Redshift'] > 6.45]
    lbg_table = lbg_table[lbg_table['Redshift'] > 6.]
    lbg_table = lbg_table[lbg_table['Redshift'] < 7.5]
    # Modify plot dir accordingly
    plot_dir = Path.cwd().parent.parent / 'plots' / 'brown_dwarfs'

    lbg_table = lbg_table[lbg_table['Age'] == 0.5]

    # Create dictionaries to store magnitudes for brown dwarfs and LBGs
    bd_mags = {}
    spex_mags = {}
    lbg_mags = {}

    # Add spectral type to bd_mags
    bd_mags['Spectral Type'] = np.array(bd_table['Spectral Type'])
    spex_mags['Spectral Type'] = np.array(spex_table['Spectral Type'])

    # Add redshift, Av and age to lbg_mags
    lbg_mags['Redshift'] = np.array(lbg_table['Redshift'])
    lbg_mags['Av'] = np.array(lbg_table['Av'])
    lbg_mags['Age'] = np.array(lbg_table['Age'])

    # Extract magnitudes for brown dwarfs
    for mag in magnitudes:
        bd_mags[f"{mag}"] = np.array(bd_table[mag])

    # Extract magnitudes for dwarf templates
    for mag in magnitudes:
        spex_mags[f"{mag}"] = np.array(spex_table[mag])


    # Extract magnitudes for LBGs
    for mag in magnitudes:
        lbg_mags[f"{mag}"] = np.array(lbg_table[mag])


    #! Plot J-Je vs Y-J
    plt.figure(figsize=(8, 8))

    plt.plot(bd_mags['Y'] - bd_mags['J'], bd_mags['J'] - bd_mags['Je'], 'o', label='Brown Dwarfs', color='red', alpha=0.8, marker='*')

    # Do a LBG plot for each value of Av and age
    for Av in np.unique(lbg_mags['Av']):
        for age in np.unique(lbg_mags['Age']):
            idx = np.logical_and(lbg_mags['Av'] == Av, lbg_mags['Age'] == age)
            plt.plot(lbg_mags['Y'][idx] - lbg_mags['J'][idx], lbg_mags['J'][idx] - lbg_mags['Je'][idx], label=f'LBGs, Av={Av}, Age={age}', alpha=0.6, color='gray', marker='none', lw=2.5)

    plt.scatter(spex_mags['Y'] - spex_mags['J'], spex_mags['J'] - spex_mags['Je'], label='SpeX Templates', color='black', alpha=0.8, marker='*', s=20, zorder=10)

    for i, spectral_type in enumerate(spex_mags['Spectral Type']):
        spec = str(spectral_type).split('b')[-1].replace("'", "")
        plt.text(spex_mags['Y'][i] - spex_mags['J'][i]+0.01, spex_mags['J'][i] - spex_mags['Je'][i]+0.01, f'{spec}', fontsize=12, zorder=10)

    #Check to see where the z>7.5 LBGs are
    # for redshift in np.unique(lbg_mags['Redshift']):
    #     if redshift >= 7.:
    #         idx = lbg_mags['Redshift'] == redshift
    #         if round(redshift, 2) == 7.50:
    #             plt.plot(lbg_mags['Y'][idx] - lbg_mags['J'][idx], lbg_mags['J'][idx] - lbg_mags['Je'][idx], label=f'LBGs, z={redshift}', alpha=0.6, marker='none', lw=2.5, c='blue')
    #             #plt.text(lbg_mags['Y'][idx][-1] - lbg_mags['Ye'][idx][-1], lbg_mags['J'][idx][-1] - lbg_mags['Je'][idx][-1], f'z={redshift:.2f}', fontsize=12)
    #         else:
    #             plt.plot(lbg_mags['Y'][idx] - lbg_mags['J'][idx], lbg_mags['J'][idx] - lbg_mags['Je'][idx], label=f'LBGs, z={redshift}', alpha=0.4, marker='none', lw=2.5, c='gray')


    plt.xlabel(r'$Y - J$')
    plt.ylabel(r'$J - J_{\mathrm{e}}$')

    plt.tight_layout()

    plt.xlim(-0.55, 1.76)
    plt.ylim(-0.67, 0.1)

    plt.savefig(plot_dir / 'Y-J_vs_J-Je.png')


    #! Plot z-Y vs Y-Ye
    # plt.figure(figsize=(8, 8))

    # plt.plot(bd_mags['z'] - bd_mags['Y'], bd_mags['Y'] - bd_mags['Ye'], 'o', label='Brown Dwarfs', color='red', alpha=0.8, marker='*')

    # # Do a LBG plot for each value of Av and age
    # for Av in np.unique(lbg_mags['Av']):
    #     for age in np.unique(lbg_mags['Age']):
    #         idx = np.logical_and(lbg_mags['Av'] == Av, lbg_mags['Age'] == age)
    #         plt.plot(lbg_mags['z'][idx] - lbg_mags['Y'][idx], lbg_mags['Y'][idx] - lbg_mags['Ye'][idx], label=f'LBGs, Av={Av}, Age={age}', alpha=0.4, color='gray', marker='none', lw=0.5)

    # plt.scatter(spex_mags['z'] - spex_mags['Y'], spex_mags['Y'] - spex_mags['Ye'], label='SpeX Templates', color='black', alpha=0.8, marker='*', s=20, zorder=10)

    # for i, spectral_type in enumerate(spex_mags['Spectral Type']):
    #     spec = str(spectral_type).split('b')[-1].replace("'", "")
    #     plt.text(spex_mags['z'][i] - spex_mags['Y'][i]+0.01, spex_mags['Y'][i] - spex_mags['Ye'][i]+0.01, f'{spec}', fontsize=12, zorder=10)

    # plt.xlabel(r'$z - Y$')
    # plt.ylabel(r'$Y - Y_{\mathrm{e}}$')

    # plt.tight_layout()

    # # plt.xlim(-1., 1.)
    # # plt.ylim(-0., 5.)

    # plt.savefig(plot_dir / 'z-Y_vs_Y-Ye.png')

    #! Plot Y-J vs z-Y
    # plt.figure(figsize=(8, 8))

    # plt.plot(bd_mags['Y'] - bd_mags['J'], bd_mags['z'] - bd_mags['Y'], 'o', label='Brown Dwarfs', color='red', alpha=0.8, marker='*')

    # # Do a LBG plot for each value of Av and age
    # for Av in np.unique(lbg_mags['Av']):
    #     for age in np.unique(lbg_mags['Age']):
    #         idx = np.logical_and(lbg_mags['Av'] == Av, lbg_mags['Age'] == age)
    #         plt.plot(lbg_mags['Y'][idx] - lbg_mags['J'][idx], lbg_mags['z'][idx] - lbg_mags['Y'][idx], label=f'LBGs, Av={Av}, Age={age}', alpha=0.6, color='gray', marker='none', lw=2.5)

    # plt.scatter(spex_mags['Y'] - spex_mags['J'], spex_mags['z'] - spex_mags['Y'], label='SpeX Templates', color='black', alpha=0.8, marker='*', s=20, zorder=10)

    # for i, spectral_type in enumerate(spex_mags['Spectral Type']):
    #     spec = str(spectral_type).split('b')[-1].replace("'", "")
    #     plt.text(spex_mags['Y'][i] - spex_mags['J'][i]+0.01, spex_mags['z'][i] - spex_mags['Y'][i]+0.01, f'{spec}', fontsize=12, zorder=10)

    # #? Open and plot REBELS
    # rebels_g = Table.read(rebels_dir / 'REBELS_ground_fluxes.fits')
    # rebels_e = Table.read(rebels_dir / 'REBELS_euclid_fluxes.fits')

    # rebels_z = flux_to_mag(rebels_g['flux_HSC-Z_DR3'])
    # rebels_Y = flux_to_mag(rebels_g['flux_Y'])
    # rebels_J = flux_to_mag(rebels_g['flux_J'])

    # # Plot
    # plt.plot(rebels_Y - rebels_J, rebels_z - rebels_Y, 'o', label='REBELS', color='green', alpha=0.8, marker='s')

    # plt.xlabel(r'$Y - J$')
    # plt.ylabel(r'$z - Y$')

    # plt.tight_layout()

    # plt.xlim(-1., 1.)
    # plt.ylim(-0., 5.)

    # plt.savefig(plot_dir / 'Y-J_vs_z-Y.png')

    #! Plot Y-J vs Ye-Je
        
    # plt.figure(figsize=(8, 8))

    # plt.plot(bd_mags['Y'] - bd_mags['J'], bd_mags['Ye'] - bd_mags['Je'], 'o', label='Brown Dwarfs', color='red', alpha=0.8, marker='*')

    # # Do a LBG plot for each value of Av and age
    # for Av in np.unique(lbg_mags['Av']):
    #     for age in np.unique(lbg_mags['Age']):
    #         idx = np.logical_and(lbg_mags['Av'] == Av, lbg_mags['Age'] == age)
    #         plt.plot(lbg_mags['Y'][idx] - lbg_mags['J'][idx], lbg_mags['Ye'][idx] - lbg_mags['Je'][idx], label=f'LBGs, Av={Av}, Age={age}', alpha=0.6, color='gray', marker='none', lw=2.5)

    # plt.scatter(spex_mags['Y'] - spex_mags['J'], spex_mags['Ye'] - spex_mags['Je'], label='SpeX Templates', color='black', alpha=0.8, marker='*', s=20, zorder=10)

    # for i, spectral_type in enumerate(spex_mags['Spectral Type']):
    #     spec = str(spectral_type).split('b')[-1].replace("'", "")
    #     plt.text(spex_mags['Y'][i] - spex_mags['J'][i]+0.01, spex_mags['Ye'][i] - spex_mags['Je'][i]+0.01, f'{spec}', fontsize=12, zorder=10)

    # plt.xlabel(r'$Y - J$')
    # plt.ylabel(r'$Y_{\mathrm{e}} - J_{\mathrm{e}}$')

    # plt.tight_layout()

    # # plt.xlim(-0.55, 1.25)
    # # plt.ylim(-0.6, 0.6)

    # plt.savefig(plot_dir / 'Y-J_vs_Ye-Je.png')

    #! Plot Ye-Je vs Je-He

    # plt.figure(figsize=(8, 8))

    # plt.plot(bd_mags['Ye'] - bd_mags['Je'], bd_mags['Je'] - bd_mags['He'], 'o', label='Brown Dwarfs', color='red', alpha=0.8, marker='*')

    # # Do a LBG plot for each value of Av and age
    # for Av in np.unique(lbg_mags['Av']):
    #     for age in np.unique(lbg_mags['Age']):
    #         idx = np.logical_and(lbg_mags['Av'] == Av, lbg_mags['Age'] == age)
    #         plt.plot(lbg_mags['Ye'][idx] - lbg_mags['Je'][idx], lbg_mags['Je'][idx] - lbg_mags['He'][idx], label=f'LBGs, Av={Av}, Age={age}', alpha=0.6, color='gray', marker='none', lw=2.5)

    # plt.scatter(spex_mags['Ye'] - spex_mags['Je'], spex_mags['Je'] - spex_mags['He'], label='SpeX Templates', color='black', alpha=0.8, marker='*', s=20, zorder=10)

    # for i, spectral_type in enumerate(spex_mags['Spectral Type']):
    #     spec = str(spectral_type).split('b')[-1].replace("'", "")
    #     plt.text(spex_mags['Ye'][i] - spex_mags['Je'][i]+0.01, spex_mags['Je'][i] - spex_mags['He'][i]+0.01, f'{spec}', fontsize=12, zorder=10)

    # plt.xlabel(r'$Y_{\mathrm{e}} - J_{\mathrm{e}}$')
    # plt.ylabel(r'$J_{\mathrm{e}} - H_{\mathrm{e}}$')

    # plt.tight_layout()

    # plt.xlim(-0.6, 0.6)
    # plt.ylim(-1.2, 0.4)

    # plt.savefig(plot_dir / 'Ye-Je_vs_Je-He.png')

    #! Plot Y-Ye vs J-Je

    # plt.figure(figsize=(8, 8))

    # plt.plot(bd_mags['Y'] - bd_mags['Ye'], bd_mags['J'] - bd_mags['Je'], 'o', label='Brown Dwarfs', color='red', alpha=0.8, marker='*')

    # # Do a LBG plot for each value of Av and age
    # for Av in np.unique(lbg_mags['Av']):
    #     for age in np.unique(lbg_mags['Age']):
    #         idx = np.logical_and(lbg_mags['Av'] == Av, lbg_mags['Age'] == age)
    #         plt.plot(lbg_mags['Y'][idx] - lbg_mags['Ye'][idx], lbg_mags['J'][idx] - lbg_mags['Je'][idx], label=f'LBGs, Av={Av}, Age={age}', alpha=0.6, color='gray', marker='none', lw=2.5)

    # #Check to see where the z>7.5 LBGs are
    # for redshift in np.unique(lbg_mags['Redshift']):
    #     if redshift >= 7.:
    #         idx = lbg_mags['Redshift'] == redshift
    #         if round(redshift, 2) == 7.50:
    #             plt.plot(lbg_mags['Y'][idx] - lbg_mags['Ye'][idx], lbg_mags['J'][idx] - lbg_mags['Je'][idx], label=f'LBGs, z={redshift}', alpha=0.6, marker='none', lw=2.5, c='blue')
    #             #plt.text(lbg_mags['Y'][idx][-1] - lbg_mags['Ye'][idx][-1], lbg_mags['J'][idx][-1] - lbg_mags['Je'][idx][-1], f'z={redshift:.2f}', fontsize=12)
    #         else:
    #             plt.plot(lbg_mags['Y'][idx] - lbg_mags['Ye'][idx], lbg_mags['J'][idx] - lbg_mags['Je'][idx], label=f'LBGs, z={redshift}', alpha=0.4, marker='none', lw=2.5, c='gray')


    # plt.scatter(spex_mags['Y'] - spex_mags['Ye'], spex_mags['J'] - spex_mags['Je'], label='SpeX Templates', color='black', alpha=0.8, marker='*', s=20, zorder=10)

    # for i, spectral_type in enumerate(spex_mags['Spectral Type']):
    #     spec = str(spectral_type).split('b')[-1].replace("'", "")
    #     plt.text(spex_mags['Y'][i] - spex_mags['Ye'][i]+0.01, spex_mags['J'][i] - spex_mags['Je'][i]+0.01, f'{spec}', fontsize=12, zorder=10)



    # plt.xlabel(r'$Y - Y_{\mathrm{e}}$')
    # plt.ylabel(r'$J - J_{\mathrm{e}}$')

    # plt.tight_layout()

    # # plt.xlim(-0.5, 0.8)
    # plt.ylim(-0.7, 0.0)

    # plt.savefig(plot_dir / 'Y-Je_vs_J-Je.png')

    #! Plot J-Je vs H-He

    # plt.figure(figsize=(8, 8))

    # plt.plot(bd_mags['J'] - bd_mags['Je'], bd_mags['H'] - bd_mags['He'], 'o', label='Brown Dwarfs', color='red', alpha=0.8, marker='*')

    # # Do a LBG plot for each value of Av and age
    # for Av in np.unique(lbg_mags['Av']):
    #     for age in np.unique(lbg_mags['Age']):
    #         idx = np.logical_and(lbg_mags['Av'] == Av, lbg_mags['Age'] == age)
    #         plt.plot(lbg_mags['J'][idx] - lbg_mags['Je'][idx], lbg_mags['H'][idx] - lbg_mags['He'][idx], label=f'LBGs, Av={Av}, Age={age}', alpha=0.6, color='gray', marker='none', lw=2.5)

    # # Plot a line at each LBG redshift value
    # for redshift in np.unique(lbg_mags['Redshift']):
    #     idx = lbg_mags['Redshift'] == redshift
    #     plt.plot(lbg_mags['J'][idx] - lbg_mags['Je'][idx], lbg_mags['H'][idx] - lbg_mags['He'][idx], label=f'LBGs, z={redshift}', alpha=0.6, marker='none', lw=2.5, c='gray')

    # plt.scatter(spex_mags['J'] - spex_mags['Je'], spex_mags['H'] - spex_mags['He'], label='SpeX Templates', color='black', alpha=0.8, marker='*', s=20, zorder=10)

    # for i, spectral_type in enumerate(spex_mags['Spectral Type']):
    #     spec = str(spectral_type).split('b')[-1].replace("'", "")
    #     plt.text(spex_mags['J'][i] - spex_mags['Je'][i]+0.01, spex_mags['H'][i] - spex_mags['He'][i]+0.01, f'{spec}', fontsize=12, zorder=10)

    # plt.xlabel(r'$J - J_{\mathrm{e}}$')
    # plt.ylabel(r'$H - H_{\mathrm{e}}$')

    # plt.tight_layout()

    # #plt.xlim(-0.6, 0.6)
    # #plt.ylim(-1.2, 0.4)

    # plt.savefig(plot_dir / 'J-Je_vs_H-He.png')

    #! Y-Ye vs H-He
    # plt.figure(figsize=(8, 8))

    # plt.plot(bd_mags['Y'] - bd_mags['Ye'], bd_mags['H'] - bd_mags['He'], 'o', label='Brown Dwarfs', color='red', alpha=0.8, marker='*')

    # # Do a LBG plot for each value of Av and age
    # for Av in np.unique(lbg_mags['Av']):
    #     for age in np.unique(lbg_mags['Age']):
    #         idx = np.logical_and(lbg_mags['Av'] == Av, lbg_mags['Age'] == age)
    #         plt.plot(lbg_mags['Y'][idx] - lbg_mags['Ye'][idx], lbg_mags['H'][idx] - lbg_mags['He'][idx], label=f'LBGs, Av={Av}, Age={age}', alpha=0.6, color='gray', marker='none', lw=2.5)

    # plt.scatter(spex_mags['Y'] - spex_mags['Ye'], spex_mags['H'] - spex_mags['He'], label='SpeX Templates', color='black', alpha=0.8, marker='*', s=20, zorder=10)

    # for i, spectral_type in enumerate(spex_mags['Spectral Type']):
    #     spec = str(spectral_type).split('b')[-1].replace("'", "")
    #     plt.text(spex_mags['Y'][i] - spex_mags['Ye'][i]+0.01, spex_mags['H'][i] - spex_mags['He'][i]+0.01, f'{spec}', fontsize=12, zorder=10)

    # plt.xlabel(r'$Y - Y_{\mathrm{e}}$')
    # plt.ylabel(r'$H - H_{\mathrm{e}}$')

    # plt.tight_layout()

    #plt.xlim(-0.6, 0.6)
    #plt.ylim(-1.2, 0.4)

    # plt.savefig(plot_dir / 'Y-Ye_vs_H-He.png') 

    #! J-Je vs Je - He
    # plt.figure(figsize=(8, 8))

    # plt.plot(bd_mags['J'] - bd_mags['Je'], bd_mags['Je'] - bd_mags['He'], 'o', label='Brown Dwarfs', color='red', alpha=0.8, marker='*')

    # # Do a LBG plot for each value of Av and age
    # for Av in np.unique(lbg_mags['Av']):
    #     for age in np.unique(lbg_mags['Age']):
    #         idx = np.logical_and(lbg_mags['Av'] == Av, lbg_mags['Age'] == age)
    #         plt.plot(lbg_mags['J'][idx] - lbg_mags['Je'][idx], lbg_mags['Je'][idx] - lbg_mags['He'][idx], label=f'LBGs, Av={Av}, Age={age}', alpha=0.6, color='gray', marker='none', lw=2.5)

    # plt.scatter(spex_mags['J'] - spex_mags['Je'], spex_mags['Je'] - spex_mags['He'], label='SpeX Templates', color='black', alpha=0.8, marker='*', s=20, zorder=10)

    # for i, spectral_type in enumerate(spex_mags['Spectral Type']):
    #     spec = str(spectral_type).split('b')[-1].replace("'", "")
    #     plt.text(spex_mags['J'][i] - spex_mags['Je'][i]+0.01, spex_mags['Je'][i] - spex_mags['He'][i]+0.01, f'{spec}', fontsize=12, zorder=10)

    # plt.xlabel(r'$J - J_{\mathrm{e}}$')
    # plt.ylabel(r'$J{\mathrm{e}} - H_{\mathrm{e}}$')

    # plt.tight_layout()

    # # plt.xlim(-0.6, 0.6)
    # # plt.ylim(-1.2, 0.4)

    # plt.savefig(plot_dir / 'J-Je_vs_Je-He.png') 

    #! y-Y vs Y-Ye
    # plt.figure(figsize=(8, 8))

    # plt.plot(bd_mags['y'] - bd_mags['Y'], bd_mags['Y'] - bd_mags['Ye'], 'o', label='Brown Dwarfs', color='red', alpha=0.8, marker='*')

    # # Do a LBG plot for each value of Av and age
    # for Av in np.unique(lbg_mags['Av']):
    #     for age in np.unique(lbg_mags['Age']):
    #         idx = np.logical_and(lbg_mags['Av'] == Av, lbg_mags['Age'] == age)
    #         plt.plot(lbg_mags['y'][idx] - lbg_mags['Y'][idx], lbg_mags['Y'][idx] - lbg_mags['Ye'][idx], label=f'LBGs, Av={Av}, Age={age}', alpha=0.6, color='gray', marker='none', lw=2.5)

    # plt.scatter(spex_mags['y'] - spex_mags['Y'], spex_mags['Y'] - spex_mags['Ye'], label='SpeX Templates', color='black', alpha=0.8, marker='*', s=20, zorder=10)

    # for i, spectral_type in enumerate(spex_mags['Spectral Type']):
    #     spec = str(spectral_type).split('b')[-1].replace("'", "")
    #     plt.text(spex_mags['y'][i] - spex_mags['Y'][i]+0.01, spex_mags['Y'][i] - spex_mags['Ye'][i]+0.01, f'{spec}', fontsize=12, zorder=10)

    # plt.xlabel(r'$y_{\mathrm{HSC}} - Y_{\mathrm{VISTA}}$')
    # plt.ylabel(r'$Y_{\mathrm{VISTA}} - Y_{\mathrm{e}}$')

    # plt.tight_layout()

    # # plt.xlim(-0.6, 0.6)
    # # plt.ylim(-1.2, 0.4)

    # #? Open and plot REBELS
    # rebels_g = Table.read(rebels_dir / 'REBELS_ground_fluxes.fits')
    # rebels_e = Table.read(rebels_dir / 'REBELS_euclid_fluxes.fits')

    # rebels_y = flux_to_mag(rebels_g['flux_HSC-Y_DR3'])
    # rebels_Y = flux_to_mag(rebels_g['flux_Y'])
    # rebels_Ye = flux_to_mag(rebels_e['flux_Y'])

    # # Plot
    # plt.plot(rebels_y - rebels_Y, rebels_Y - rebels_Ye, 'o', label='REBELS', color='green', alpha=0.8, marker='s')

    # plt.savefig(plot_dir / 'y-Y_vs_Y-Ye.png')
        
    #! J - H vs Je - He
    # plt.figure(figsize=(8, 8))

    # plt.plot(bd_mags['J'] - bd_mags['H'], bd_mags['Je'] - bd_mags['He'], 'o', label='Brown Dwarfs', color='red', alpha=0.8, marker='*')

    # # Do a LBG plot for each value of Av and age
    # for Av in np.unique(lbg_mags['Av']):
    #     for age in np.unique(lbg_mags['Age']):
    #         idx = np.logical_and(lbg_mags['Av'] == Av, lbg_mags['Age'] == age)
    #         plt.plot(lbg_mags['J'][idx] - lbg_mags['H'][idx], lbg_mags['Je'][idx] - lbg_mags['He'][idx], label=f'LBGs, Av={Av}, Age={age}', alpha=0.6, color='gray', marker='none', lw=2.5)

    # # Plot a line at each LBG redshift value
    # for redshift in np.unique(lbg_mags['Redshift']):
    #     idx = lbg_mags['Redshift'] == redshift
    #     plt.plot(lbg_mags['J'][idx] - lbg_mags['H'][idx], lbg_mags['Je'][idx] - lbg_mags['He'][idx], label=f'LBGs, z={redshift}', alpha=0.6, marker='none', lw=2.5, c='gray')

    # plt.scatter(spex_mags['J'] - spex_mags['H'], spex_mags['Je'] - spex_mags['He'], label='SpeX Templates', color='black', alpha=0.8, marker='*', s=20, zorder=10)

    # for i, spectral_type in enumerate(spex_mags['Spectral Type']):
    #     spec = str(spectral_type).split('b')[-1].replace("'", "")
    #     plt.text(spex_mags['J'][i] - spex_mags['H'][i]+0.01, spex_mags['Je'][i] - spex_mags['He'][i]+0.01, f'{spec}', fontsize=12, zorder=10)

    # plt.xlabel(r'$J - J_{\mathrm{e}}$')
    # plt.ylabel(r'$H - H_{\mathrm{e}}$')

    # plt.tight_layout()

    # #plt.xlim(-0.6, 0.6)
    # #plt.ylim(-1.2, 0.4)

    # plt.savefig(plot_dir / 'J-H_vs_Je-He.png')

    # ! Add text
    # On the LBGs add text of the redshift and Av
    #for i, redshift in enumerate(redshifts):
    #    for j, Av in enumerate(Av_vals):
    #        plt.text(Y_mags_LBG[i*len(Av_vals) + j] - Ye_mags_LBG[i*len(Av_vals) + j], J_mags_LBG[i*len(Av_vals) + j] - Je_mags_LBG[i*len(Av_vals) + j], f'z={redshift:.2f}, Av={Av:.2f}', fontsize=8)
    
    #On the BDs add text of the spectral type
    # for i, spectral_type in enumerate(bd_mags.keys):
    #    plt.text(J_mags_BD[i] - Je_mags_BD[i], H_mags_BD[i] - He_mags_BD[i], f'{spectral_type}', fontsize=8)

    

    plt.show()
    plt.close()



if __name__ == '__none__':
#! ----------------------------------------
#! -------- PLOTTING ON FILTERS -----------
#! ----------------------------------------
    bds = loadBrownDwarfTemplates(spectral_types=['M'], sub_types=[5])

    # Plot BD templates
    for spectral_type, (wavelength, flux) in bds.items():

        # Make figure
        plt.figure(figsize=(10, 6))

        for redshift in np.arange(6.5, 7.5, 0.1):

            print(f'Making figure at z={redshift:.2f}')

            # Make LBG model
            wlen, LBG_flux = makeLBG(redshift=redshift, SFH_component='constant', age=(0, 13.8), massformed=10.5, metallicity=0.2, 
                dust_type='Calzetti', Av=0.2, nebular=True, logU=-2.5, Muv=-21.)
            
            LBG_flux = powerLawSpectrum(wlen, 1e-30, -2., redshift=redshift)

            wlen = np.array(wlen)
            LBG_flux = np.array(LBG_flux)            

            # Filters
            euclid_filters = getFilters('euclid', plot=True, plot_kwargs={'linewidth': 2.5, 'alpha':0.6, 'color': 'blue'})
            vista_filters = getFilters('vista', plot=True, plot_kwargs={'linewidth': 2.5, 'alpha':0.6, 'color': 'orange'})

            # Get the magnitudes
            mags = convolveFilters([vista_filters, euclid_filters], {redshift: (wlen, LBG_flux)})

            # # Dummy plots for filter labels
            plt.plot([], [], label='Euclid', color='blue', linewidth=2.5, alpha=0.6)
            plt.plot([], [], label='VISTA', color='orange', linewidth=2.5, alpha=0.6)

            # # Brown dwarf template
            #plt.plot(wavelength, flux/np.max(flux), label=f'{spectral_type} dwarf', color='red', alpha=0.8, linewidth=2.5)
            
            # LBG model
            plt.plot(wlen, LBG_flux/np.max(LBG_flux), label='LBG', color='red', alpha=0.8, linewidth=2.5)

            plt.xlabel(r'$\lambda \ (\AA)$')
            plt.ylabel('Relative Transmission/Flux')
            plt.xlim(5000, 25000)
            plt.ylim(0, 1.15)
            #plt.legend(loc='upper right')
            plt.tight_layout()

            #plt.savefig(plot_dir.parent / 'LBG_models' / f'z{redshift:.2f}_LBG.png')
            #plt.close()

        plt.show()








