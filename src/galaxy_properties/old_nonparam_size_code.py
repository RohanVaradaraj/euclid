65#!/usr/bin/env python3

'''
nonparametric_sizes.py

Following Roper et al. 2022, use a pixel-based approach to measure the sizes.

We can then do a like-for-like comparison of sizes and the size-luminosity relations.

Created: Tuesday 3rd October 2023 (2nd day of third year)
'''

import os
import numpy as np
from jwst_classes import Catalogue, MultibandGalaxyFitter
from astropy.io import fits
from astropy.table import Table
from astropy.cosmology import FlatLambdaCDM
import matplotlib.pyplot as plt
from pathlib import Path
import autogalaxy as ag
from photutils import CircularAperture, aperture_photometry
from photutils import detect_threshold, detect_sources
from photutils.background import Background2D, MedianBackground

def aperture_30pkpc(z, cosmo):
    '''
    Computes the aperture radius for an object at redshift z such that 30kpc is enclosed.
    '''
    
    ASperkpc = cosmo.arcsec_per_kpc_proper(z)
    #print(ASperkpc)
    aperture_size = (30 * ASperkpc.value) / 2
    return aperture_size

# Visualise objects and apertures
vis = False

# Run testing?
test = False

##### Setting up some testing objects #####

if test:
    # Size in arcseconds for the first object
    size_arcsec1 = 1.0

    # Size in arcseconds for the second object
    size_arcsec2 = 0.5

    # Field size in arcseconds
    field_size_arcsec = 2.0

    # Pixel scale in arcsec per pixel
    pixScale = 0.03

    # Calculate the size in pixels for both objects
    size_pix1 = size_arcsec1 / pixScale
    size_pix2 = size_arcsec2 / pixScale

    field_size_pix = field_size_arcsec / pixScale

    # Create a blank field for the first object
    field1 = np.zeros((int(field_size_pix), int(field_size_pix)))

    # Calculate the center of the field for the first object
    center_x1 = field1.shape[0] // 2
    center_y1 = field1.shape[1] // 2

    # Create a blank field for the second object
    field2 = np.zeros((int(field_size_pix), int(field_size_pix)))

    # Calculate the center of the field for the second object
    center_x2 = field2.shape[0] // 2
    center_y2 = field2.shape[1] // 2

    # Create artificial disks with a gradient for the first object
    for i in range(field1.shape[0]):
        for j in range(field1.shape[1]):
            # Calculate the distance from the center
            distance1 = np.sqrt((i - center_x1)**2 + (j - center_y1)**2)
            if distance1 <= (size_pix1 / 2):
                # Intensity value with a gradient from 1 at the center to 0.1 at the edge
                intensity1 = 1.0 - 0.9 * (distance1 / (size_pix1 / 2))
                # Ensure intensity is within the range [0.1, 1.0]
                intensity1 = max(0.1, min(1.0, intensity1))
                field1[i, j] = intensity1

    # Create artificial disks with a gradient for the second object
    for i in range(field2.shape[0]):
        for j in range(field2.shape[1]):
            # Calculate the distance from the center
            distance2 = np.sqrt((i - center_x2)**2 + (j - center_y2)**2)
            if distance2 <= (size_pix2 / 2):
                # Intensity value with a gradient from 1 at the center to 0.1 at the edge
                intensity2 = 1.0 - 0.9 * (distance2 / (size_pix2 / 2))
                # Ensure intensity is within the range [0.1, 1.0]
                intensity2 = max(0.1, min(1.0, intensity2))
                field2[i, j] = intensity2

# We're 'bout to flex some mad science, peep this lit code 🔬🔥

# Initialize the cosmology, just for the vibes 🌌
cosmo = FlatLambdaCDM(H0=70, Om0=0.3)

# Pixel scale, you know, how we pixelatin' it 📸
pixScale = 0.03  # Arcsec per pixel

# We need that perfect size, in arcsec, it's a vibe 📏
sizeAS = 2

# Calculate the size in pixels, keeping it pixel-perfect 📐
sizePix = (sizeAS / 2) / pixScale

# Calculate the area of a pixel, it's small but mighty 🌟
pixArea = pixScale * pixScale

# Extent of our plot, just keeping it within bounds 📊
min_pix = 333 - int(sizePix)
max_pix = 333 + int(sizePix)

# Catalogue name, where we find all the cool stuff 📦
#catName = 'Z3_crossmatched.fits' # COSMOS
catName = 'Z3_overlap.fits' # UDS

# Redshift string, gotta check that cosmic zip code 🌌📭
zString = catName.split('_')[0]

# Filter names, we got all the shades 🌈
#filterNames = ['f090w', 'f115w', 'f150w', 'f200w', 'f277w', 'f356w', 'f410m', 'f444w']
filterNames = ['f115w', 'f200w', 'f356w', 'f444w']
filterNames = ['f356w']

# Directory vibes, where the stamps live 🏠💌
#baseDir = Path.home() / 'JWST' / 'stamps' / f'{zString}'
baseDir = Path.home().parent.parent / 'extraspace' / 'varadaraj' / 'JWST' / 'UDS' / 'stamps' / f'{zString}'
psfDir = Path.home() / 'JWST' / 'psf'

# And a catalogue name to save the outputs to
outputDir = Path.cwd().parents[1] / 'ref_catalogues' / 'nathan' / 'sizes'
output_table = Table(names=('ID', 'ZPhot') + tuple(f'R_{filterName}' for filterName in filterNames) + tuple(f'dR_{filterName}' for filterName in filterNames))

# Read the catalogue, you know, data source unlocked 📚💽
# cat = Catalogue(Catalogue.jwstDir + 'UDS/', catName)
# cat.sort_column('stellar_mass')
# cat.cut_catalogue('stellar_mass', 9)
# cat.cat = cat.cat[cat.cat['ZPhot'] > 2.75]
# IDs = cat.cat['UID']
# ZPhots = cat.cat['ZPhot']

# #####! CATALOGUE FOR CORRECTIONS: COSMOLOGICAL DIMMING TEST. #######
# cat = Table.read(Path.cwd() / 'corrections' / 'results' / 'cosmo_dim_input_sizes.fits')
# IDs = cat['ID']
# ZPhots = cat['ZPhot']
# #! And update outputDir and baseDir
# outputDir = Path.cwd() / 'corrections' / 'results'
# baseDir = Path.home() / 'JWST' / 'corrections'

ap_sizes = []
# Loop through the IDs, one by one, it's a vibe 🔄
for j, ID in enumerate(IDs):

    #print(aperture_30pkpc(ZPhots[j], cosmo))
    ap_sizes.append(aperture_30pkpc(ZPhots[j], cosmo))

    print('###############################')
    print(f'Object {ID}, {j} of {len(IDs)}')
    print('###############################')


    #if os.path.exists(str(baseDir / f'{filterNames[-1]}' / f'ID_{ID}_{filterNames[-1]}_bgsub.fits')):
    if os.path.exists(str(baseDir / f'ID_{ID}_f356w_bgsub.fits')):

        #! Normal running
        images = [
            ag.Imaging.from_fits(image_path=baseDir / f'{filterName}' / f'ID_{ID}_{filterName}_bgsub.fits',
                                noise_map_path=baseDir / f'{filterName}' / f'ID_{ID}_{filterName}_err.fits',
                                #psf_path=psfDir / f'{filterName}_empirical_psf.fits',
                                psf_path=psfDir / f'{filterName}_UDS_empirical_psf.fits',
                                pixel_scales=pixScale)
            for filterName in filterNames
        ]
        # #! COSMOLOGICAL DIMMING SOURCES
        # images = [
        #     ag.Imaging.from_fits(data_path=baseDir / f'ID_{ID}_{filterName}_bgsub.fits',
        #                         noise_map_path=baseDir / f'ID_{ID}_{filterName}_errNozeros.fits',
        #                         psf_path=psfDir / f'{filterName}_UDS_empirical_psf.fits',
        #                         pixel_scales=pixScale)
        #     for filterName in filterNames
        # ]


        # Masks for that aperture vibe, keepin' it focused 🔍🎯
        masks = [
            ag.Mask2D.circular(
                shape_native=imaging.shape_native,
                pixel_scales=imaging.pixel_scales,
                radius=1.0
            )
            for imaging in images
        ]

        # Apply those masks, it's like wearing shades 😎🕶️
        masked_images = [
            image.apply_mask(mask=mask)
            for image, mask in zip(images, masks)
        ]

        # Initialize lists to store the sizes and size errors for each filter, to save to a table
        sizes = []
        size_errors = []

        # Loop through each filter image, 'cause we got options 🎥🔍
        for k, image in enumerate(masked_images):

            # Extract the data in the central vibe, it's that core sample 🧪🔍
            data = masked_images[k].image.native #[min_pix:max_pix, min_pix:max_pix]

            if test:
                # Modify the image with the test data?
                data = field1 if j == 0 else field2

            # #! Modified for COSMODIM
            # bkg_estimator = MedianBackground()
            # bkg = Background2D(data, (50, 50), filter_size=(3, 3),
            #                 bkg_estimator=bkg_estimator)
            # threshold = 1.5 * bkg.background_rms

            threshold = detect_threshold(data, nsigma=2.0, background=0.) #! ORIGINAL

            # Detect sources in the data using the threshold
            seg = detect_sources(data, threshold, npixels=5)

            # Now you have a segmentation map 'segm' with labeled objects

            # Get the background-subtracted object data by setting background to zero
            #print('Before/after segmentation')
            #print(np.sum(data))
            data_obj = data.copy()

            if seg is not None:
                # Get the background-subtracted object data by setting background to zero
                data_obj = data.copy()
                data_obj[seg.data == 0] = 0
            else:
                # If no sources are detected, set size to 0
                size = 0
                size_error = 0

            # Get the psf too
            psf = masked_images[k].psf

            #print(psf.native.shape)

            ###### Scaling PSF to object photometry ######

            # Create a circular aperture for the data
    #        data_aperture = CircularAperture((data.native.shape[0] / 2, data.native.shape[1] / 2), r=0.5)        
            data_aperture = CircularAperture((masked_images[k].image.native[min_pix:max_pix, min_pix:max_pix].shape[0] / 2, masked_images[k].image.native[min_pix:max_pix, min_pix:max_pix].shape[1] / 2), r=0.5)
            #data_aperture = CircularAperture((masked_images[k].image.native.shape[0] / 2, masked_images[k].image.native.shape[1] / 2), r=0.5) #! COSMODIM

            # Create a circular aperture for the PSF
            psf_aperture = CircularAperture((psf.native.shape[0] / 2, psf.native.shape[1] / 2), r=0.5)

            # Measure the total flux within the data aperture
            data_flux = aperture_photometry(data, data_aperture)['aperture_sum'][0]

            # Measure the total flux within the PSF aperture
            psf_flux = aperture_photometry(psf.native, psf_aperture)['aperture_sum'][0]

            # Scale the PSF to match the luminosity of the data
            scaled_psf = psf * (data_flux / psf_flux)

            # Aperture photometry on the scaled psf
            scaled_aperture = CircularAperture((scaled_psf.native.shape[0] / 2, scaled_psf.native.shape[1] / 2), r=0.5)
            scaled_flux = aperture_photometry(scaled_psf.native, scaled_aperture)['aperture_sum'][0]

            print('SCALED FLUX: ', scaled_flux)
            print('OG FLUX: ', data_flux)

            ############################################


            ####### Size analysis on PSF ##########

            scaled_psf = scaled_psf.flatten()

            scaled_psf = np.sort(scaled_psf)[::-1]

            total_psf = np.sum(scaled_psf)

            cumsum = 0
            num_pix_psf = 0

            for lum in scaled_psf:
                cumsum += lum
                num_pix_psf += 1

                if cumsum >= total_psf / 2:
                    break
            
            area_psf = num_pix_psf * pixArea

            rAng_psf = np.sqrt(area_psf / np.pi)

            ASperkpc = cosmo.arcsec_per_kpc_proper(ZPhots[j])

            psf_size = rAng_psf / ASperkpc.value
            print('PSF size: ', psf_size)

            ########################################

            # Flatten the data, gotta keep it flat for analysis 🥞📉

            data_fl = data.flatten()
            original_order = np.arange(data_fl.size)  # Keep track of original pixel positions

            # Sort that data, high to low, it's how we roll 📊📈
            sorted_indices = np.argsort(data_fl)[::-1]
            invert = np.argsort(sorted_indices)
            data_fl = data_fl[sorted_indices]

            #print(f'Total pixels: {len(data_fl)}')

            # Calculate the total light, all that glow ✨💡
            total_light = np.sum(data_fl)

            # Analyze the data array, pixel by pixel, the grind never stops 🔄🔍
            cumulative_sum = 0
            num_pixels = 0

            highlighted_pixels = np.zeros_like(data_fl)  # Create a binary mask for highlighting

            # Keep adding until we hit that halfway mark, hustle hard 📈💪
            for i, luminosity in enumerate(data_fl):
                cumulative_sum += luminosity
                num_pixels += 1

                highlighted_pixels[i] = 1

                # Check if we've crossed that halfway mark, goals achieved 🏆🔥
                if cumulative_sum >= total_light / 2:
                    break

            print(f'Number of pixels for half-light: {num_pixels}')
            print(f'Number of pixels in mask: {np.sum(highlighted_pixels)}')

    #        highlighted_pixels = highlighted_pixels[invert].reshape(data.native.shape[0], data.native.shape[1])
            highlighted_pixels = highlighted_pixels[invert].reshape(masked_images[k].image.native[min_pix:max_pix, min_pix:max_pix].shape[0], masked_images[k].image.native[min_pix:max_pix, min_pix:max_pix].shape[1])       
            #highlighted_pixels = highlighted_pixels[invert].reshape(masked_images[k].image.native.shape[0], masked_images[k].image.native.shape[1]) #! COSMODIM
            # Calculate the area, in square arcsec, it's all about space 🌌🔲
            area = num_pixels * pixArea

            # Radius in arcsec, we're measuring that vibe 🌟🔍
            rAng = np.sqrt(area / np.pi)

            print(f'ANGULAR SIZE: {rAng} arcsec')

            # Convert angular size to physical size, the science flex 📏🧪
            ASperkpc = cosmo.arcsec_per_kpc_proper(ZPhots[j])
            print('arcsec per kpc at this redshift: ', ASperkpc)

            # Convert angular size to physical size, it's all about the numbers 🔢🔬
            size = rAng / ASperkpc.value
            print(f'Physical size: {size} kpc')

            ################# SUBTRACT PSF SIZE IN QUADRATURE ##################
            size = np.sqrt(size**2 - psf_size**2)
            print(f'PSF unbroadened: {size} kpc')

            # We're 'bout to flex our bootstraps, you know the drill 🥾💪

            # Number of times we flexin' it, no cap 🔄💯
            num_iterations = 1000  # Tweak this if you wanna flex more or less, no cap

            # Stash our bootstrapped sizes, it's all about that data drip 💧📊
            bootstrap_sizes = []

            # Let's dive into that flex zone (aka bootstrapping), we bootstrappin' it 🥾📈
            for _ in range(num_iterations):

                # Get a fresh flex by randomly sampling with replacements, sus vibes 🎲🔄
                bootstrap_sample = np.random.choice(data_fl, size=len(data_fl), replace=True)

                # Calculate our cumulative flex and find how many pixels we need to stay bussin' 📊🔢
                cumulative_flex = np.cumsum(bootstrap_sample)
                half_luminosity_pixels = np.argmax(cumulative_flex >= np.sum(bootstrap_sample) / 2)

                # Calculate the size 'cause size matters, W rizz 📏📈
                size_bootstrap = half_luminosity_pixels * pixArea

                # Add it to our stash of flexed sizes, it's lit 🔥📏
                bootstrap_sizes.append(size_bootstrap)

            # Calculate the average and swag factor of our flexed sizes, on god 📈📊
            average_size = np.mean(bootstrap_sizes)
            size_error = np.std(bootstrap_sizes)

            # Append the size and size error to the lists
            sizes.append(size)
            size_errors.append(size_error)

            # Create a circular aperture with the calculated radius in pixels
            aperture_radius_pix = ap_sizes[j] 
    #        aperture = CircularAperture((data.native.shape[0] / 2, data.native.shape[1] / 2), r=ap_sizes[j])
            aperture = CircularAperture((masked_images[k].image.native[min_pix:max_pix, min_pix:max_pix].shape[0] / 2, masked_images[k].image.native[min_pix:max_pix, min_pix:max_pix].shape[1] / 2), r=ap_sizes[j])
            #aperture = CircularAperture((masked_images[k].image.native.shape[0] / 2, masked_images[k].image.native.shape[1] / 2), r=ap_sizes[j]) #! COSMODIM

            #print('Ap size: ', ap_sizes[j])

            if vis:
                # Visualize the image
                plt.figure(figsize=(8, 8))

    #            plt.imshow(data.native, cmap='viridis', origin='lower', vmin=-0.1, vmax=0.1, zorder=-1)
                plt.imshow(data, cmap='viridis', origin='lower', vmin=-0.1, vmax=1, zorder=-1)

                # Overlay pixels containing half-light
                #plt.imshow(highlighted_pixels, cmap='gray', origin='lower', alpha=0.5) 

                # Overlay the aperture on the image

                aperture.plot(color='red', lw=2)
                plt.xlabel('X (pixels)')
                plt.title(f'Object {ID} - Filter {filterNames[k]}')
                plt.ylabel('Y (pixels)')
                plt.colorbar(label='Flux')

                # Save or show the plot (optional)
                #plt.savefig(f'object_{ID}_filter_{filterNames[k]}.png')
                plt.show()

                plt.clf()
                plt.close()

            # Print or use average_size and swag_factor as your size estimate and its swagger level, on god 💯🔥
            print(f'Size = {round(size, 2)} +/- {round(size_error, 2)} kpc')

        # Create a row of data for the current ID
        data_row = [ID, ZPhots[j]] + sizes + size_errors

        # Add the data row to the output table
        output_table.add_row(data_row)


# Save the output table to a FITS file
output_table.write(outputDir / f'cosmo_dim_output_sizes.fits', overwrite=True)

