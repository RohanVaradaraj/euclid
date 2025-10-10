"""
Save the SEDs of some interesting galaxies into a format we can put into the JWST exposure time calc (calc is slang for calculator btw).

Created: Wednesday 10th September 2025.
"""

from pathlib import Path
import numpy as np
import sys
from scipy.signal import savgol_filter
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt

sed_fitting_path = Path.cwd().parent / 'sed_fitting'
sys.path.append(str(sed_fitting_path))
from sed_fitting_codes import parse_spec_file

def mag_to_flux(mag):
    return 10**(-0.4 * (np.array(mag) + 48.6))

def smooth_continuum(flux, window=51, polyorder=3):
    """Estimate smooth continuum using Savitzky-Golay filter."""
    return savgol_filter(flux, window_length=min(window, len(flux)//2*2+1), polyorder=polyorder)

def add_emission_lines(wlen, flux, z, lines, R=3000, resample_factor=5):
    """
    Add emission lines to the SED.

    Parameters
    ----------
    wlen : array
        Wavelength in microns (observed frame).
    flux : array
        Flux density in mJy.
    z : float
        Redshift of the galaxy.
    lines : list of dicts
        Each dict: {'name': 'CIII]', 'lam_rest': 1909, 'EW_rest': 20} (lam in Å, EW in Å)
    R : float
        Spectral resolution for line width (lambda/dlambda).
    resample_factor : int
        Finer sampling factor around lines.

    Returns
    -------
    wlen_new, flux_new : arrays with added lines
    """
    # Interpolate SED to finer grid
    wl_fine = np.linspace(wlen.min(), wlen.max(), len(wlen) * resample_factor)
    f_interp = interp1d(wlen, flux, kind='linear', bounds_error=False, fill_value="extrapolate")
    f_fine = f_interp(wl_fine)

    f_cont = smooth_continuum(f_fine, window=11, polyorder=3)

    for line in lines:
        lam_obs = line['lam_rest'] * 1e-4 * (1 + z)  # Convert Å to microns
        print(f"Adding line {line['name']} at observed λ = {lam_obs:.4f} microns")
        EW_obs = line['EW_rest'] * (1 + z) * 1e-4    # Convert Å to microns

        # Find continuum at this λ
        f_cont_line = np.interp(lam_obs, wl_fine, f_cont)

        # Gaussian line profile
        sigma = lam_obs / R / (2*np.sqrt(2*np.log(2)))  # convert FWHM to sigma
        amplitude = (f_cont_line * EW_obs) / (sigma * np.sqrt(2*np.pi))
        f_line = amplitude * np.exp(-0.5 * ((wl_fine - lam_obs) / sigma)**2)

        f_fine += f_line

    return wl_fine, f_fine

#! Define some paths
sed_path = Path.cwd().parents[1] / 'data' / 'sed_fitting' / 'zphot' / 'COSMOS' / 'best_fits' / 'det_Y_J_with_euclid_z7'
output_path = Path.cwd().parents[1] / 'data' / 'jwst_proposal' / 'etc'

#interesting_IDs = ['178396', '869246', '510801', '640903', '878786', '499162']
#interesting_IDs = ['591643', '869246', '510801', '640903', '878786', '499162']

#! Rest-UV lines
emission_lines = [
    {'name': 'CIII]', 'lam_rest': 1909, 'EW_rest': 20},
    {'name': 'OIII]', 'lam_rest': 1666, 'EW_rest': 15},
    {'name': 'Lyα',   'lam_rest': 1216, 'EW_rest': 10},
    {'name': 'HeII',  'lam_rest': 1640, 'EW_rest': 5},
    {'name': 'CIV',   'lam_rest': 1549, 'EW_rest': 5},
    {'name': 'NV',    'lam_rest': 1240, 'EW_rest': 5},
    {'name': '[OIII]', 'lam_rest': 5007, 'EW_rest': 300},
    {'name': 'Hβ',    'lam_rest': 4861, 'EW_rest': 100},
    {'name': 'Hα',    'lam_rest': 6563, 'EW_rest': 300},
]

# emission_lines = [
#     {'name': 'CIII]', 'lam_rest': 1909, 'EW_rest': 100},
#     {'name': 'OIII]', 'lam_rest': 1666, 'EW_rest': 100},
#     {'name': 'Lyα',   'lam_rest': 1216, 'EW_rest': 100},
#     {'name': 'HeII',  'lam_rest': 1640, 'EW_rest': 100},
#     {'name': 'CIV',   'lam_rest': 1549, 'EW_rest': 100},
#     {'name': 'NV',    'lam_rest': 1240, 'EW_rest': 100},
# ]

#! Things which are clumpy, and have a grism redshift in agreement with zphot.
interesting_IDs = ['387777', '615679', '878786', '591643']
redshifts = {'387777': 6.7614, '615679': 6.8484, '878786': 7.0975, '591643': 7.1105}

for ID in interesting_IDs:

    spec_file_name = f'Id000{ID}.spec'

    file = parse_spec_file(sed_path / spec_file_name)

    sed = file.get('sed')

    wlen = [table['lambda'] for table in sed]
    mag = [table['flux'] for table in sed] # In mag
    print(mag[0])

    flux = np.array(mag_to_flux(mag[0]))  # In cgs

    # Convert flux from cgs to millijanskys
    flux *= 1e26  # 1 mJy = 1e-26 erg/s/cm^2/Hz

    # Convert wlen from Angstroms to microns
    wlen = np.array(wlen[0]) * 1e-4  # 1 micron = 10,000 Angstroms

    z = redshifts[ID]
    wl_new, flux_new = add_emission_lines(wlen, flux, z, emission_lines, R=1000)

    plt.plot(wl_new, flux_new, label='With Emission Lines')
    plt.title(f'ID {ID} at z={z}')

#
    # Plot the savgol filtered continuum
    #plt.plot(wlen, flux, label='Original SED')
    #plt.plot(wl_new, smooth_continuum(flux_new, window=11, polyorder=3), label='Smoothed Continuum')
    # exit()
    plt.show()
    # Write to text file
    # output_file = output_path / f'ID_{ID}.txt'
    # with open(output_file, 'w') as f:
    #     f.write('# Wavelength (microns)    Flux (mJy)\n')
    #     for wl, fl in zip(wlen, flux):
    #         f.write(f'{wl:.6f}    {fl:.6f}\n')
    


