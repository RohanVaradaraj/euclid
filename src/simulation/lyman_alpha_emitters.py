"""
lyman_alpha_emitters.py

Simulate Lyman-alpha emitters (LAEs) passing through the Euclid+VISTA filter set.

Created: Friday 10th May 2024.
"""

import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table
import sys
from pathlib import Path
from astropy.cosmology import FlatLambdaCDM
from synphot import etau_madau
from synphot import SourceSpectrum, Empirical1D, units

import astropy.units as u


sys.path.append(str(Path.cwd()))
from brown_dwarf_colours import mag_to_flux, flux_to_mag, max_age_at_redshift, getFilters, convolveFilters, makeLBG

plt.rcParams['axes.linewidth'] = 2.5
plt.rcParams.update({'font.size': 15})
plt.rcParams['figure.dpi'] = 100

filter_dir = Path().home() / 'lephare' / 'lephare_dev' / 'filt' / 'myfilters'
plot_dir = Path.cwd().parent.parent / 'plots' / 'LAEs'
table_dir = Path.cwd().parent.parent / 'data' / 'simulations' / 'tables'
rebels_dir = Path.cwd().parent.parent / 'data' / 'ref_catalogues' / 'fluxes'
sed_dir = Path.cwd().parent.parent / 'data' / 'simulations' / 'LAEs'

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
    for i, lam in enumerate(wlen):
        if (lam > 1450.0*(1+z)) and (lam < 1550.0*(1+z)):
            tophat_filter[i] = 1.0

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

    # Plot flux with and without absorption
    # plt.figure()
    # plt.plot(wlen, flux, label='Original', alpha=0.6, lw=2.5, color='deepskyblue')
    # plt.plot(wlen, flux_new, label='Madau Absorption', alpha=0.6, lw=2.5, color='red')
    # plt.xlabel('Wavelength (A)')
    # plt.ylabel('Flux')
    # plt.legend()
    # plt.xlim(5000, 70000)
    # plt.yscale('log')
    # plt.show()

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

    # First measure the continuum flux between 1250-1300A (but observed frame)

    # Make top-hat filter between these values
    tophat_filter = np.zeros(len(wlen))
    for i, lam in enumerate(wlen):
        if (lam > 1249.0) and (lam < 1301.0):
            tophat_filter[i] = 1.0

    # Measure the flux through this top-hat filter
    continuum_flux = np.sum(flux * tophat_filter) / np.sum(tophat_filter)

    print('Continuum: ', continuum_flux)

    # Find the index and value where the wavelength is closest to 1215.67A
    idx = np.argmin(np.abs(wlen - 1215.67 * (1+z)))

    # Compute the flux density of the SED point at 1215.67A

    line_flux =  np.abs(3 * (1 - (EW / (wlen[1] - wlen[0])))*continuum_flux - (flux[idx-1] + flux[idx+1]))
    line_flux_2 = np.abs((1 - (EW / (wlen[1] - wlen[0])))*continuum_flux)

    print('#########')
    print(line_flux)
    print(line_flux_2)
    print(flux[idx])


    return wlen, flux





# wlen, flux = makeLBG(redshift=0.0, SFH_component='constant', age=(0, 13.8), massformed=11., metallicity=0.2, 
#                                 dust_type='Calzetti', Av=0.2, nebular=True, logU=-1.)


t = Table.read(sed_dir / 'bc2003_lr_m62_chab_tau_0_05_age_0_05_EW240A.ascii', format='ascii')

wlen = t['col1']
flux = t['col2']


wlen, flux = redshift_SED(z=7, wlen=wlen, flux=flux)
wlen, flux = set_Muv(z=7, Muv_target=-24.0, wlen=wlen, flux=flux)

print(wlen[1] - wlen[0])

wlen, flux = add_emission_line(EW=240, z=7, wlen=wlen, flux=flux)
exit()

#plt.figure(figsize=(10, 6))
plt.plot(wlen, flux, lw=2.5, color='deepskyblue')
plt.xlabel(r'Wavelength ($\AA$)')
plt.ylabel('Flux')
plt.legend()
plt.xlim(5000, 40000)
plt.yscale('log')
plt.ylim(3e-31, 1e-29)
plt.show()





