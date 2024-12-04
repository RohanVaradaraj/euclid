"""
Creates source injection class.

Created: Wednesday 4th December 2024.
"""

import numpy as np
from scipy.interpolate import interp1d
from astropy.cosmology import FlatLambdaCDM
from astropy import units as u
from utils import load_config, filter_files
import numpy as np
import matplotlib.pyplot as plt
from scipy.constants import physical_constants
from pathlib import Path
from scipy.integrate import simps



class SourceInjector:
    def __init__(self, samples, params):
        """
        Initialize the SourceInjector with the necessary parameters.
        
        :param samples: Number of samples to draw from the luminosity function.
        :param z_range: Range of redshift values to draw from.
        :param beta_mean: Mean value of the beta slope.
        :param beta_std: Standard deviation of the beta slope.
        :param filters: Dictionary of filter responses.
        """
        self.samples = samples
        self.zmin = params['zmin']
        self.zmax = params['zmax']
        self.beta_mean = params['beta_mean']
        self.beta_std = params['beta_std']
        self.filters = params['filters']

        filter_dir = Path.home() / 'lephare' / 'lephare_dev' / 'filt'

        #! Read in the filter response curves
        filter_response_curves = {}

        filter_dict = filter_files()
        
        # Get filter file paths from self.filter keys
        filters = {band: filter_dict[band] for band in self.filters}
        print(filters)
        
        for filter_name, filter_path in filters.items():

            # Open the filter response curve
            filter_data = np.genfromtxt(filter_dir / filter_path)

            # Add the filter response curve to the dictionary
            filter_response_curves[filter_name] = filter_data
        
        self.filter_response_curves = filter_response_curves

        #! Cosmology
        self.cosmo = FlatLambdaCDM(H0=70, Om0=0.3)
    


    def draw_parameters(self):
        """
        Draw random redshift and beta slope values for all samples.
        
        :return: Arrays of redshifts and beta slopes.
        """
        n = len(self.samples)
        z = np.random.uniform(self.zmin, self.zmax, n)
        beta = np.random.normal(self.beta_mean, self.beta_std, n)
        return z, beta



    def generate_seds(self, z, beta):
        """
        Generate un-normalized SEDs for all samples in a vectorized way.
        
        :param z: Array of redshifts.
        :param beta: Array of beta slopes.
        :return: 2D array of wavelengths and 2D array of fluxes (one row per sample).
        """
        # Define the wavelength grid for all SEDs
        wavelengths = np.linspace(1000, 30000, 1000)  # Angstroms

        lyman_breaks = 1216 * (1 + z[:, None])  # Lyman break for each redshift

        # Calculate the flux based on the power-law SED
        flux = wavelengths[None, :] ** beta[:, None]

        # Apply redshift and Lyman break
        flux /= (1 + z[:, None])**2
        flux[wavelengths < lyman_breaks] = 0

        # Convert flux to f_nu
        c = physical_constants['speed of light in vacuum'][0]
        frequencies = c / wavelengths  # Frequencies in Hz
        flux = flux * c / (frequencies[None, :]**2)

        print(flux.shape)
        print(wavelengths.shape)

        return wavelengths, flux



    def scale_seds_to_muv(self, sed_wavelengths, sed_fluxes, Muv, z):
        """
        Scale the SEDs to match given Muv magnitudes for all samples.
        
        :param sed_wavelengths: 2D array of wavelengths.
        :param sed_fluxes: 2D array of fluxes.
        :param Muv: Array of desired Muv magnitudes.
        :param z: Array of redshifts.
        :return: Scaled fluxes (2D array).
        """
        # Tophat filter for Muv
        tophat_filter = np.zeros_like(sed_fluxes)
        filter_mask = (sed_wavelengths >= 1450 * (1 + z[:, None])) & (sed_wavelengths <= 1550 * (1 + z[:, None]))
        tophat_filter[filter_mask] = 1

        # Flux in the tophat filter
        flux_tophat = np.trapz(sed_fluxes * tophat_filter, sed_wavelengths, axis=1) / 100  # Width of the filter

        # Luminosity distance
        DL = self.cosmo.luminosity_distance(z).value * 10**6  # in pc

        # Current Muv
        current_Muv = -2.5 * np.log10(flux_tophat) - 48.6 - 5 * np.log10(DL / 10) + 2.5 * np.log10(1 + z)

        # Scaling factors
        scaling_factors = 10 ** ((Muv - current_Muv) / -2.5)

        # Scale SEDs
        scaled_fluxes = sed_fluxes * scaling_factors[:, None]

        return scaled_fluxes


    def calculate_fluxes(self, sed_wavelengths, sed_fluxes):
        """
        Calculate the fluxes in the filters for all samples in a vectorized way.
        
        :param sed_wavelengths: 2D array of wavelengths.
        :param sed_fluxes: 2D array of scaled fluxes.
        :return: Dictionary of fluxes for each filter (key: filter name, value: 1D array).
        """
        fluxes = {}

        for filter_name, filter_data in self.filter_response_curves.items():
            filter_wlen = filter_data[:, 0]
            filter_trans = filter_data[:, 1]

            # Normalize the transmission curve
            filter_trans /= np.max(filter_trans)

            # Interpolate SEDs to filter wavelength grid
            interpolated_seds = np.array([np.interp(filter_wlen, sed_wavelengths, sed_fluxes[i, :])
                                        for i in range(len(sed_fluxes))])

            # Calculate fluxes in the filter for all samples
            filter_area = simps(filter_trans, filter_wlen)
            flux = np.array([simps(interpolated_seds[i, :] * filter_trans, filter_wlen) / filter_area
                            for i in range(len(interpolated_seds))])

            fluxes[filter_name] = flux

        return fluxes


    def inject_sources(self, image, n_sources, Muv_range):
        """
        Inject n_sources into the image based on the randomly generated parameters.
        
        :param image: The 2D image where sources will be injected.
        :param n_sources: The number of sources to inject.
        :param Muv_range: The range of Muv values (for scaling).
        :return: The image with the injected sources.
        """
        for _ in range(n_sources):
            z, beta = self.draw_parameters()
            Muv = np.random.uniform(Muv_range[0], Muv_range[1])  # Random Muv for each source
            
            # Generate and scale the SED
            sed_wavelengths, sed_flux = self.generate_sed(z, beta)
            scaled_flux = self.scale_sed_to_muv(sed_wavelengths, sed_flux, Muv)
            
            # Calculate magnitudes in the Y and J bands
            magnitudes = self.calculate_magnitudes(sed_wavelengths, scaled_flux)
            
            # Here you would inject the source into the image (e.g., adding a point source)
            # Example: you could inject this source into a random location
            x_pos, y_pos = np.random.randint(0, image.shape[1]), np.random.randint(0, image.shape[0])
            self._inject_single_source(image, x_pos, y_pos, magnitudes, sed_wavelengths, scaled_flux)
        
        return image

    def _inject_single_source(self, image, x, y, magnitudes, sed_wavelengths, scaled_flux):
        """
        Helper function to inject a single source into the image.
        
        :param image: The image to inject the source into.
        :param x, y: The position of the source.
        :param magnitudes: Magnitudes of the source in Y, J bands.
        :param sed_wavelengths: Wavelengths of the scaled SED.
        :param scaled_flux: Scaled fluxes for the source.
        """
        # This could be a simple Gaussian, point source, or more complex PSF
        # Add a point source at (x, y) with the calculated magnitude
        pass


#! Testing
if __name__ == '__main__':

    config = load_config("config.yaml")
    injection_config = config['source_injection']

    Muv_sample = np.linspace(-22, -20, 5)
    source_injector = SourceInjector(samples=Muv_sample, params=injection_config)
    z, beta = source_injector.draw_parameters()
    # plt.hist(z)
    # plt.hist(beta)
    # plt.show()
    wavelengths, fluxes = source_injector.generate_seds(z, beta)
    scaled_fluxes = source_injector.scale_seds_to_muv(wavelengths, fluxes, Muv_sample, z)
    filter_fluxes = source_injector.calculate_fluxes(wavelengths, scaled_fluxes)

    # x_Y = np.full(len(Muv_sample), 10214)
    # x_J = np.full(len(Muv_sample), 12544)
    # plt.plot(wavelengths, scaled_fluxes.T)
    # plt.scatter(x_Y, filter_fluxes['Y'], color='red')
    # plt.scatter(x_J, filter_fluxes['J'], color='red')
    
    # for i, Muv in enumerate(samples):
    #     wavelengths, flux = source_injector.generate_sed(z[i], beta[i])
    #     scaled_flux = source_injector.scale_sed_to_muv(wavelengths, flux, Muv, z[i])
    #     fluxes = source_injector.calculate_fluxes(wavelengths, scaled_flux)

        # plt.plot(wavelengths, scaled_flux, label=f'Muv = {Muv:.2f}')
        # plt.scatter([10214], fluxes['Y'], label='Y', color='red')
        # plt.scatter([12544], fluxes['J'], label='J', color='red')

    # plt.yscale('log')
    # plt.ylim(3e-32, 1e-29)
    # plt.xlim(3000, 40000)
    #plt.show()







