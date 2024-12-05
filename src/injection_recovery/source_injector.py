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
from photutils.aperture import CircularAperture
from astropy.io import fits


class SourceInjector:
    def __init__(self, samples, params):
        """
        Initialize the SourceInjector with the necessary parameters.
        
        :param samples: Number of samples to draw from the luminosity function.
        :param zmin: Minimum redshift value.
        :param zmax: Maximum redshift value.
        :param beta_mean: Mean value of the beta slope.
        :param beta_std: Standard deviation of the beta slope.
        :param filters: Dictionary of filter responses.
        :param image_size: Size of the images (in arcminutes).
        :param n_images: Number of images to generate.
        :param base_image: Path to the base image.
        """
        self.samples = samples
        self.zmin = params['zmin']
        self.zmax = params['zmax']
        self.dz = params['dz']
        self.beta_mean = params['beta_mean']
        self.beta_std = params['beta_std']
        self.filters = params['filters']
        self.image_size = params['image_size_arcmin']
        self.n_images = params['n_images']
        self.base_image = params['base_image']

        filter_dir = Path.home() / 'lephare' / 'lephare_dev' / 'filt'

        # Read in the filter response curves
        filter_response_curves = {}
        filter_dict = filter_files()  # Assuming this function exists
        filters = {band: filter_dict[band] for band in self.filters}
        
        for filter_name, filter_path in filters.items():
            filter_data = np.genfromtxt(filter_dir / filter_path)
            filter_response_curves[filter_name] = filter_data
        
        self.filter_response_curves = filter_response_curves

        # Cosmology
        self.cosmo = FlatLambdaCDM(H0=70, Om0=0.3)



    def draw_parameters(self):
        """
        Draw random redshift and beta slope values for all samples.
        
        :return: Arrays of redshifts and beta slopes.
        """
        n = len(self.samples)
        
        # Define a grid for redshift with dz = 0.05
        z_grid = np.arange(self.zmin, self.zmax + self.dz, self.dz)

        z_prob = np.ones_like(z_grid)  # Uniform distribution over z_grid
        
        # Normalize the probability distribution to form a valid CDF
        z_cdf = np.cumsum(z_prob) / np.sum(z_prob)
        
        # Sample redshifts from the grid using inverse transform sampling
        uniform_randoms = np.random.uniform(0, 1, n)
        z_samples = np.interp(uniform_randoms, z_cdf, z_grid)

        # Draw beta slopes from a normal distribution
        beta_samples = np.random.normal(self.beta_mean, self.beta_std, n)
        
        return z_samples, beta_samples



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
        # Area of filter
        tophat_area = simps(tophat_filter, sed_wavelengths, axis=1)

        # Flux in the tophat filter
        flux_tophat = np.trapz(sed_fluxes * tophat_filter, sed_wavelengths, axis=1) / tophat_area # Width of the filter

        # Luminosity distance
        DL = self.cosmo.luminosity_distance(z).value * 10**6  # in pc

        # Current Muv
        current_Muv = -2.5 * np.log10(flux_tophat) - 48.6 - 5 * np.log10(DL / 10) + 2.5 * np.log10(1 + z)

        # Scaling factors
        scaling_factors = 10 ** ((Muv - current_Muv) / -2.5)

        # Scale SEDs
        scaled_fluxes = sed_fluxes * scaling_factors[:, None]

        return scaled_fluxes



    def calculate_fluxes(self, sed_wavelengths, sed_fluxes, average=True):
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

        # Get average fluxes if required
        if average:
            # Get the combined filter name from the keys - append them
            combined_filter = ''.join(filter_name for filter_name, _ in self.filter_response_curves.items())
            fluxes[combined_filter] = np.sum([fluxes[filter_name] for filter_name in self.filters], axis=0)
            self.combined_filter = combined_filter
        else:
            self.combined_filter = None

        return fluxes



    def generate_random_positions(self, image_size):
        """
        Generate random positions for the sources within the image.
        
        :param image_size: Size of the image in pixels.
        :return: Arrays of x and y positions.
        """

        # Convert image size from arcmin to pixels
        image_size = image_size * 60 / 0.15

        x = np.random.randint(0, image_size, len(self.samples))
        y = np.random.randint(0, image_size, len(self.samples))

        return x, y



    def get_psf(self):
        """
        Add PSF to the class.
        """

        psf_dir = Path.cwd().parents[1] / 'data' / 'psf' / 'COSMOS' / 'ref_psf'
        hdu = fits.open(psf_dir / 'Y_DR6_psf.fits')
        psf = hdu[0].data

        self.psf = psf

        return None



    def measure_psf_flux(self):

        # Place 1.8 arcsec diameter aperture at centre of image
        aperture = CircularAperture((self.psf.shape[1] / 2, self.psf.shape[0] / 2), r=1.8 / 0.15)

        # Measure the flux in the aperture
        flux, _ = aperture.do_photometry(self.psf)

        # PSF correction
        flux /= 0.69

        return flux

    
    
    def scale_psf_to_Muv(self, fluxes):
        """
        Scale the N = len(samples) psfs we will inject into the image to the desired Y+J band magnitude.

        :param fluxes: Array of fluxes in Y+J.
        :return: Array of scaled PSFs.
        """

        # Get the PSF flux
        psf_flux = self.measure_psf_flux()

        # Check if there exists a combined filter
        if self.combined_filter is None:
            fluxes = fluxes[self.filters]
        else:
            fluxes = fluxes[self.combined_filter]

        vista_zpt = 30.

        # Convert the cgs fluxes to counts, which are the units of the image
        fluxes *= 10** ((48.6 + vista_zpt)/2.5)

        # Calculate the scaling factor
        scaling_factors = fluxes / psf_flux

        # Make N copies of the PSF
        scaled_psfs = np.array([scaling_factors[i] * self.psf for i in range(len(scaling_factors))])

        return scaled_psfs



    def inject_sources(self, image_name, x_positions, y_positions, scaled_psfs, overwrite=False):
        """
        Inject the sources into the image at the given positions.
        
        :param image_name: Path to the FITS file containing the image to inject the sources into.
        :param x_positions: Array of x positions.
        :param y_positions: Array of y positions.
        :param scaled_psfs: List or array of scaled PSFs.
        :return: None. Saves the injected image to disk.
        """

        image_dir = Path.cwd() / 'images' / 'cutouts' 
        injected_dir = Path.cwd() / 'images' / 'injected'

        if overwrite:
            for file in injected_dir.glob('*'):
                file.unlink()

        # Open the image
        with fits.open(image_dir / image_name) as hdul:
            image = hdul[0].data
            header = hdul[0].header

        image_height, image_width = image.shape
        n_psfs = len(scaled_psfs)

        # Create a blank overlay to hold the injections
        overlay = np.zeros_like(image)

        # Iterate over each PSF (vectorization limited by differing injection positions)
        for i in range(n_psfs):
            psf = scaled_psfs[i]
            x = x_positions[i]
            y = y_positions[i]

            # Get PSF size
            psf_height, psf_width = psf.shape

            # Ensure injection is within bounds, and calculate slices for both PSF and image
            x_start, x_end = max(0, x - psf_width // 2), min(image_width, x + psf_width // 2)
            y_start, y_end = max(0, y - psf_height // 2), min(image_height, y + psf_height // 2)

            # Adjust psf_x and psf_y slices based on the image slices
            psf_x_start = max(0, psf_width // 2 - x)
            psf_x_end = psf_x_start + (x_end - x_start)

            psf_y_start = max(0, psf_height // 2 - y)
            psf_y_end = psf_y_start + (y_end - y_start)

            # Ensure the slices match in size
            overlay[y_start:y_end, x_start:x_end] += psf[psf_y_start:psf_y_end, psf_x_start:psf_x_end]

        # Add the overlay to the image
        image += overlay

        # Save the injected image
        hdu = fits.PrimaryHDU(image, header=header)
        hdu.writeto(injected_dir / image_name, overwrite=True)

        return None



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


    #! Testing injection






