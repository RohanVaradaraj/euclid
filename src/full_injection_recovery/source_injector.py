#!/usr/bin/env python3

"""
Creates source injection class.

Created: Wednesday 4th December 2024.
"""

from scipy.interpolate import interp1d
from astropy.cosmology import FlatLambdaCDM
from astropy import units as u
import numpy as np
from utils import filter_files
import matplotlib.pyplot as plt
from scipy.constants import physical_constants
from pathlib import Path
from scipy.integrate import simpson
from photutils.aperture import CircularAperture
from astropy.io import fits
# import bagpipes as pipes
from synthesizer.emission_models.attenuation import Asada25
from unyt import angstrom

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
        self.zeropoints = params['zeropoints']
        self.image_size = params['image_size_arcmin']
        self.n_images = params['n_images']
        self.base_image = params['base_image']
        self.ucd_ref_filter = params['ucd_ref_filter']
        self.ucd_mAB_range = params['ucd_mAB_range']
        self.ucd_types = params['ucd_types']
        self.filter_zeropoints = dict(zip(self.filters, self.zeropoints))
        self.filter_flux_zeropoints = {
            filter_name: 10 ** (-0.4 * (zeropoint + 48.6))
            for filter_name, zeropoint in self.filter_zeropoints.items()
        }
        self.ucd_ref_flux_zeropoint = self.filter_flux_zeropoints[self.ucd_ref_filter]

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
        self._psf_flux_cache = {}

    
    def load_ucd_templates(self, template_dir=Path.home() / 'lephare' / 'lephare_dev' / 'sed' / 'STAR' / 'DWARFSTARS'):

        # Get all dwarfs by stellar type, starting with M and ending with .txt_trim
        template_files_M = list(template_dir.glob('M*.txt_trim'))
        template_files_L = list(template_dir.glob('L*.txt_trim'))   
        template_files_T = list(template_dir.glob('T*.txt_trim'))

        # Read the templates into a dictionary, and get type by splitting filename by '_' and taking the first part
        templates = {}
        for template_file in template_files_M + template_files_L + template_files_T:
            template_type = template_file.stem.split('_')[0]
            data = np.genfromtxt(template_file)
            templates[template_type] = data
        self.ucd_templates = templates

        # Sort the dictionary by stellar type, in order M, L, T
        templates = dict(sorted(templates.items()))
        print('Loaded UCD templates:', list(templates.keys()))

        return templates
    
    def sample_ucd_mags(self):
        """
        Sample magnitudes for ultra-cool dwarfs from a uniform distribution within the specified range.
        
        :param n_samples: Number of samples to draw for each UCD type.
        :param mAB_range: Tuple specifying the (min, max) magnitude range.
        :return: Dictionary of sampled magnitudes for each UCD type.
        """
        n = len(self.samples)

        ucd_mags = {}
        for ucd_type in self.ucd_types:
            mags = np.random.uniform(self.ucd_mAB_range[0], self.ucd_mAB_range[1], n)
            ucd_mags[ucd_type] = mags
        return ucd_mags
    
    def draw_ucd_types(self):

        n = len(self.samples)

        ucd_sample = {}

        for type in self.ucd_types:
            # Draw n_samples of this type, at random (so e.g. can be any of M0-M9 for M type
            types = np.random.choice([template_type for template_type in self.ucd_templates.keys() if template_type.startswith(type)], n)

            ucd_sample[type] = types

        return ucd_sample
    
    
    def scale_ucds_to_mags(self, ucd_mags, ucd_types, pix_scale):
        """
        Scale the UCD templates to the desired magnitudes in the reference filter.
        
        :param ucd_mags: Dictionary of magnitudes for each UCD type.
        :param ucd_types: Dictionary of sampled UCD types for each source.
        :param pix_scale: Pixel scale in arcsec/pixel.
        :return: Wavelength arrays and scaled SEDs in reference-filter count units.
        """
        scaled_seds = []
        wlens = []

        for ucd_type in self.ucd_types:
            mags = ucd_mags[ucd_type]
            types = ucd_types[ucd_type]

            for mag, template_type in zip(mags, types):
                template_data = self.ucd_templates[template_type]
                wavelengths = template_data[:, 0]  # Angstroms
                fluxes = template_data[:, 1]  # Arbitrary units

                # Interpolate the template to the filter wavelength grid
                filter_data = self.filter_response_curves[self.ucd_ref_filter]
                filter_wlen = filter_data[:, 0]

                interp_func = interp1d(wavelengths, fluxes, bounds_error=False, fill_value=0)
                interpolated_fluxes = interp_func(filter_wlen)

                # Calculate the flux in the reference filter
                filter_trans = filter_data[:, 1] 
                filter_trans /= np.max(filter_trans)  # Normalize the transmission curve

                filter_area = simpson(filter_trans, filter_wlen)    
                flux_in_filter = simpson(interpolated_fluxes * filter_trans, filter_wlen) / filter_area 

                # Calculate the scaling factor to match the desired magnitude
                target_flux = 10 ** (-0.4 * (mag + 48.6))  # Convert mag to flux

                scaling_factor = target_flux / flux_in_filter
                scaled_fluxes = fluxes * scaling_factor
                scaled_counts = scaled_fluxes / self.ucd_ref_flux_zeropoint
                scaled_seds.append(scaled_counts)
                wlens.append(wavelengths)

        return wlens, scaled_seds



    def draw_parameters(self):
        """
        Draw random redshift and beta slope values for all samples.
        
        :return: Arrays of redshifts and beta slopes.
        """
        n = len(self.samples)
        
        # Define a grid for redshift with dz = 0.05
        z_grid = np.arange(self.zmin, self.zmax + self.dz, self.dz)

        # Draw redshifts from a uniform distribution
        z_samples = np.random.uniform(self.zmin, self.zmax, n)

        # Then snap the redshifts to the grid
        z_samples = np.round(z_samples / self.dz) * self.dz

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
        # Define the wavelength grid for all SEDs. This is in the observed frame.
        wavelengths = np.linspace(1000, 30000, 1000)  # Angstroms

        lyman_breaks = 1216 * (1 + z[:, None])  # Lyman break for each redshift

        # Calculate the flux based on the power-law SED
        flux = wavelengths[None, :] ** beta[:, None]

        # Apply redshift and Lyman break
        flux /= (1 + z[:, None])**2
        
        #flux[wavelengths < lyman_breaks] = 0
        asada = Asada25()

        # For every SED, apply the IGM transmission
        for i in range(len(flux)):
            flux[i, :] *= asada.get_transmission(lam_obs=wavelengths*angstrom, redshift=z[i])

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
        tophat_area = simpson(tophat_filter, sed_wavelengths, axis=1)

        # Flux in the tophat filter
        flux_tophat = simpson(sed_fluxes * tophat_filter, sed_wavelengths) / tophat_area # Width of the filter

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
            filter_area = simpson(filter_trans, filter_wlen)
            flux = np.array([simpson(interpolated_seds[i, :] * filter_trans, filter_wlen) / filter_area
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

        # Convert these cgs fluxes back into counts
        fluxes = {
            filter_name: flux / self.filter_flux_zeropoints[filter_name]
            for filter_name, flux in fluxes.items()
            if filter_name in self.filter_flux_zeropoints
        }


        return fluxes

    def calculate_ucd_fluxes(self, sed_wavelengths, sed_counts, average=False):
        """
        Calculate the image counts in each filter for scaled UCD SEDs.

        The input SEDs are expected to be in the reference-filter count-density
        units returned by ``scale_ucds_to_mags``.

        :param sed_wavelengths: 2D array/list of wavelength grids.
        :param sed_counts: 2D array/list of scaled count-density SEDs.
        :param average: If True, also return the summed counts across all filters.
        :return: Dictionary of fluxes for each filter (key: filter name, value: 1D array).
        """
        fluxes = {}

        if len(sed_wavelengths) != len(sed_counts):
            raise ValueError("sed_wavelengths and sed_counts must have matching lengths.")

        for filter_name, filter_data in self.filter_response_curves.items():
            filter_wlen = filter_data[:, 0]
            filter_trans = filter_data[:, 1].copy()
            filter_trans /= np.max(filter_trans)

            interpolated_seds = np.array([
                np.interp(filter_wlen, sed_wavelengths[i], sed_counts[i], left=0.0, right=0.0)
                for i in range(len(sed_counts))
            ])

            filter_area = simpson(filter_trans, filter_wlen)
            ref_filter_counts = np.array([
                simpson(interpolated_seds[i, :] * filter_trans, filter_wlen) / filter_area
                for i in range(len(interpolated_seds))
            ])

            fluxes[filter_name] = (
                ref_filter_counts
                * self.ucd_ref_flux_zeropoint
                / self.filter_flux_zeropoints[filter_name]
            )

        if average:
            combined_filter = ''.join(filter_name for filter_name, _ in self.filter_response_curves.items())
            fluxes[combined_filter] = np.sum([fluxes[filter_name] for filter_name in self.filters], axis=0)

        return fluxes


    def get_filter_centre_and_fwhm(self):
       
        filter_centre = []
        filter_fwhm = []

        for filter_name, filter_data in self.filter_response_curves.items():
            filter_wlen = filter_data[:, 0]
            filter_trans = filter_data[:, 1]

            # Normalize the transmission curve
            filter_trans /= np.max(filter_trans)

            # Calculate the filter centre
            centre = simpson(filter_wlen * filter_trans, filter_wlen) / simpson(filter_trans, filter_wlen)

            # Calculate the FWHM
            cumulative_trans = np.cumsum(filter_trans)

            fwhm = np.interp(0.5 * np.max(filter_trans), cumulative_trans, filter_wlen)
        
            filter_centre.append(centre)
            filter_fwhm.append(fwhm)
        
        return filter_centre, filter_fwhm


    def generate_random_positions(self, image_size, pix_scale, ucd=False):
        """
        Generate random positions for the sources within the image.
        
        :param image_size: Size of the image in pixels.
        :return: Arrays of x and y positions.
        """

        # Convert image size from arcmin to pixels
        image_size = image_size * 60 / pix_scale

        # Buffer for PSF size, plus some extra pixels
        psf_buffer = (75 // 2) + 20

        n_types = len(self.ucd_types)

        if ucd == False:
            x = np.random.randint(psf_buffer, image_size-psf_buffer, len(self.samples))
            y = np.random.randint(psf_buffer, image_size-psf_buffer, len(self.samples))
        else:
            x = np.random.randint(psf_buffer, image_size-psf_buffer, n_types*len(self.samples))
            y = np.random.randint(psf_buffer, image_size-psf_buffer, n_types*len(self.samples))            

        return x, y



    def get_psf(self, field_name):
        """
        Add PSF to the class.
        :param field_name: Name of the field.

        If the PSF file does not exist, it extracts it from the .psf cube file.
        """

        # Store PSFs in an array
        psf_array = []

        if 'XMM' in field_name:
            base_dir = Path.home().parents[1] / 'hoy' / 'temporaryFilesROHAN' / 'psf'
        else:
            base_dir = Path.cwd().parents[3] / 'data' / 'psf'

        for filter_name in self.filters:
            psf_path = Path.cwd().parents[1] / 'data' / 'psf' / field_name / 'ref_psf' / f'{filter_name}_psf.fits'
            if not psf_path.is_file():

                # Open the .psf file
                original_psf_path = base_dir / field_name / 'results' / f'{filter_name}.psf'

                hdu = fits.open(original_psf_path)
                data = hdu[1].data[0][0]
                header = hdu[1].header

                # PSF slice
                yz_slice = np.array(data[0, :, :])
                #  Save the slice as a fits file with the same header
                hdu = fits.PrimaryHDU(yz_slice, header=header)
                hdu.writeto(psf_path, overwrite=True)

            # Open the PSF file
            with fits.open(psf_path) as hdu:
                psf = hdu[0].data
            
            psf_array.append(psf)

        psf_array = np.array(psf_array)

        self.psf = psf_array

        return None



    def measure_psf_flux(self, pix_scale, plot=False):
        cache_key = (float(pix_scale), bool(plot))
        if cache_key in self._psf_flux_cache:
            return self._psf_flux_cache[cache_key]

        fluxes = []

        for psf in self.psf:

            # Place 1.8 arcsec diameter aperture at centre of image
            aperture = CircularAperture((psf.shape[1] / 2, psf.shape[0] / 2), r=1.0 / pix_scale)

            # Measure the flux in the aperture
            flux, _ = aperture.do_photometry(psf)

            fluxes.append(flux[0])

            # Draw aperture on PSF
            if plot:
                plt.imshow(self.psf, origin='lower')
                aperture.plot(color='red')
                plt.show()

        fluxes = np.array(fluxes)
        self._psf_flux_cache[cache_key] = fluxes
        return fluxes

    
    
    def scale_psf_to_Muv(self, fluxes, Muv, z, pix_scale):
        """
        Scale the N = len(samples) psfs we will inject into the image to the desired input image magnitudes.

        :param fluxes: Array of fluxes in Y+J.
        :return: Array of scaled PSFs.
        """

        # Get the PSF fluxes, to be scaled
        flux_counts = self.measure_psf_flux(pix_scale)

        # Check if there exists a combined filter
        if self.combined_filter is None:
            fluxes = fluxes[self.filters]
            target_mags = -2.5 * np.log10(fluxes) - 48.6
        else:
            fluxes = fluxes[self.combined_filter]
            target_mags = -2.5 * np.log10(fluxes) - 48.6

        # Convert flux count to flux
        fluxes = flux_counts  * 10 **(-0.4 * (self.zeropoints + 48.6))

        # Compute apparent magnitude, to be rescaled to mags
        app_mag = -2.5 * np.log10(fluxes) - 48.6

        # Get luminosity distance
        DL = self.cosmo.luminosity_distance(z).value * 10**6  # in pc

        # # Compute absolute magnitude
        abs_mag = app_mag - 5 * np.log10(DL / 10) + 2.5 * np.log10(1 + z)

        # Calculate the scaling factor
        scaling_factors = 10 ** (0.4 * (abs_mag - Muv))

        # Make N copies of the PSF
        scaled_psfs = np.array([scaling_factors[i] * self.psf for i in range(len(scaling_factors))])

        return scaled_psfs

    def scale_PSFs_to_model_fluxes(self, target_fluxes, pix_scale):
        """
        Scale the PSFs in each band according to the model fluxes for each source.

        Parameters
        ----------
        target_fluxes : dict
            Dictionary keyed by filter name. Each value should be either a scalar
            or a 1D array of length n_sources.
        pix_scale : float
            Pixel scale in arcsec/pixel.

        Returns
        -------
        np.ndarray
            Array with shape (n_sources, n_filters, psf_y, psf_x).
        """
        psf_fluxes = self.measure_psf_flux(pix_scale)

        flux_arrays = []
        n_sources = None
        for filter_name in self.filters:
            if filter_name not in target_fluxes:
                raise KeyError(f"Missing target fluxes for filter: {filter_name}")

            flux_array = np.atleast_1d(np.asarray(target_fluxes[filter_name], dtype=float))
            if n_sources is None:
                n_sources = len(flux_array)
            elif len(flux_array) != n_sources:
                raise ValueError("All target_fluxes arrays must have the same length.")

            flux_arrays.append(flux_array)

        flux_matrix = np.stack(flux_arrays, axis=1)
        scaling_factors = flux_matrix / psf_fluxes[np.newaxis, :]
        self.scaled_psfs = self.psf[np.newaxis, :, :, :] * scaling_factors[:, :, np.newaxis, np.newaxis]
        return self.scaled_psfs

    def inject_sources(
        self,
        image_name,
        x_positions,
        y_positions,
        Muv_array,
        z_array,
        scaled_psfs,
        overwrite=False,
        plot_each_source=False,
        filter_name=None,
        image=None,
        header=None,
    ):
        """
        Inject sources into the image.

        :param image_name: Path to the FITS file containing the image to inject the sources into.
        :param x_positions: Array of x positions.
        :param y_positions: Array of y positions.
        :param Muv_array: Array of Muv magnitudes.
        :param z_array: Array of redshifts.
        :param scaled_psfs: List or array of scaled PSFs.
        :param overwrite: Whether to overwrite existing files.
        :return: wcs of new imahe. Saves the injected image to disk.
        """
        injected_dir = Path.cwd() / 'images' / 'injected'
        injected_dir.mkdir(parents=True, exist_ok=True)

        if image is None or header is None:
            image_dir = Path.cwd() / 'images' / 'cutouts'
            with fits.open(image_dir / image_name) as hdul:
                image = np.array(hdul[0].data, copy=True)
                header = hdul[0].header.copy()
        else:
            image = np.array(image, copy=True)
            header = header.copy()

        x_positions = np.atleast_1d(np.asarray(x_positions, dtype=int))
        y_positions = np.atleast_1d(np.asarray(y_positions, dtype=int))
        Muv_array = np.atleast_1d(np.asarray(Muv_array))
        z_array = np.atleast_1d(np.asarray(z_array))

        n_sources = len(x_positions)
        if len(y_positions) != n_sources or len(Muv_array) != n_sources or len(z_array) != n_sources:
            raise ValueError(
                "x_positions, y_positions, Muv_array, and z_array must all have the same length."
            )

        filter_idx = None
        resolved_filter_name = filter_name
        for i, current_filter_name in enumerate(self.filters):
            if resolved_filter_name == current_filter_name:
                filter_idx = i
                break
            if resolved_filter_name is None and f"{current_filter_name}_cutout" in image_name:
                filter_idx = i
                resolved_filter_name = current_filter_name
                break

        if filter_idx is None:
            raise ValueError(f"Could not determine filter from image name: {image_name}")

        scaled_psfs = np.asarray(scaled_psfs)
        if scaled_psfs.ndim == 4:
            if scaled_psfs.shape[0] != n_sources:
                raise ValueError(
                    "For per-source PSFs, scaled_psfs must have shape "
                    "(n_sources, n_filters, psf_y, psf_x)."
                )
        elif scaled_psfs.ndim != 3:
            raise ValueError(
                "scaled_psfs must have shape (n_filters, psf_y, psf_x) or "
                "(n_sources, n_filters, psf_y, psf_x)."
            )

        image_height, image_width = image.shape
        psf_stack = scaled_psfs[:, filter_idx] if scaled_psfs.ndim == 4 else np.broadcast_to(
            scaled_psfs[filter_idx],
            (n_sources, *scaled_psfs[filter_idx].shape),
        )

        # Inject each source using the PSF that matches this image's filter.
        for i in range(n_sources):
            psf = psf_stack[i]

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

            if plot_each_source:
                # Print useful information for checking plots
                print('Image number:' , i+1)
                print('Inject at:', x, y)
                print(x_start, x_end, y_start, y_end)

                # Check if image is empty at these coords
                image_sum = np.sum(image[y_start:y_end, x_start:x_end])
                print('Image sum:', image_sum)
            

            # Inject directly into the image to avoid an extra full-size overlay allocation.
            image[y_start:y_end, x_start:x_end] += psf[psf_y_start:psf_y_end, psf_x_start:psf_x_end]

        #! Plot the injected sources
        if plot_each_source:
            cutout_size = 500
            for i in range(len(x_positions)):

                x_center, y_center = x_positions[i], y_positions[i]
                cutout = image[
                    y_center - cutout_size // 2 : y_center + cutout_size // 2,
                    x_center - cutout_size // 2 : x_center + cutout_size // 2,
                ]
                # cutout_weight = weight_image[
                #     y_center - cutout_size // 2 : y_center + cutout_size // 2,
                #     x_center - cutout_size // 2 : x_center + cutout_size // 2,  
                # ]                

                #* SCI
                plt.figure(figsize=(6, 6))
                plt.imshow(cutout, origin='lower', vmin=0, vmax=2, cmap='viridis')
                
                # Add crosshair
                center = cutout_size // 2
                cross_length = 7
                gap = 3  # Gap around the center

                plt.plot(np.arange(center + gap, center + cross_length), [center] * (cross_length - gap),  color='red', linewidth=1.5)
                plt.plot([center] * (cross_length - gap), np.arange(center + gap, center + cross_length), color='red', linewidth=1.5)
                
                plt.title(f'Muv: {Muv_array[i]:.2f}, z: {z_array[i]:.2f}', fontsize=12)
                plt.show()
                plt.close()

        # Plot the full image and draw big crosshairs at each source
        if plot_each_source:
            plt.figure(figsize=(10, 10))
            plt.imshow(image, origin='lower', vmin=0, vmax=2, cmap='viridis')

            for i in range(len(x_positions)):
                x_center, y_center = x_positions[i], y_positions[i]
                plt.plot(x_center, y_center, 'rx', markersize=10)
            
            plt.show()

        # Save the injected image
        hdu_image = fits.PrimaryHDU(image, header=header)
        hdu_image.writeto(injected_dir / image_name, overwrite=True)

        return None

    # secax = ax.secondary_yaxis('right', functions=(flux_to_mab, mab_to_flux))
    # secax.set_ylabel(r'$m_{\rm AB}$')
    # secax.set_yticks(np.arange(20, 30, 1))

    # plt.show()
