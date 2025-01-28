#!/usr/bin/env python3

"""
Combines all components into one high-level workflow.

Created: Wednesday 4th December 2024.
"""

from luminosity_function import LuminosityFunction
from source_injector import SourceInjector
from source_extractor import SourceExtractor
from utils import load_config, cutout_subimage
from pathlib import Path
from astropy.io import fits
from astropy.wcs import WCS
import matplotlib.pyplot as plt
from photutils.aperture import CircularAperture
import numpy as np
import glob
from astropy.table import Table
import os

plt.rcParams['axes.linewidth'] = 2.5
plt.rcParams.update({'font.size': 20})
plt.rcParams['figure.dpi'] = 100

overwrite_cats = True
overwrite_images = True

#! Run in steps
make_cutouts = False
inject_sources = False
run_se = True

def RunFullInjectionRecoveryPipeline(base_image, overwrite=True):

    #! Read in config file
    config = load_config("config.yaml")
    lf_config = config['luminosity_function']
    injection_config = config['source_injection']
    n_images = injection_config['n_images']
    image_size = injection_config['image_size_arcmin']
    se_config = config['source_extraction']
    batch_size = se_config['batch_size']

    injected_dir = Path.cwd() / 'images' / 'injected'
    image_dir = Path.cwd() / 'images' / 'cutouts'
    input_cat_dir = Path.cwd() / 'catalogues' / 'input'
    output_cat_dir = Path.cwd() / 'catalogues' / 'output'

    if make_cutouts:

        print('Generating cutouts of base image')
        cutout_subimage(base_image, image_size, n_images, random=True, overwrite=overwrite)

    if inject_sources:
        #! Get all images to inject into

        images = glob.glob(str(image_dir / '*.fits'))


        #! Catalogue dir for saving the input values of mock sources
        

        if overwrite_cats:
            for file in glob.glob(str(input_cat_dir / '*.fits')):
                os.remove(file)
        if overwrite_images:
            for file in glob.glob(str(injected_dir / '*.fits')):
                os.remove(file)

        for i, image in enumerate(images):

            print(f"Injecting sources into {image.split('/')[-1]}, which is image {i+1} of {len(images)}")

            image_name = image.split('/')[-1]

            #! Draw sample from luminosity function
            luminosity_function = LuminosityFunction(lf_config)
            Muv_sample = luminosity_function.sample_luminosities(uniform=True)

            #! Initiate sample for injection
            source_injector = SourceInjector(samples=Muv_sample, params=injection_config)
            z, beta = source_injector.draw_parameters()
    
            wavelengths, fluxes = source_injector.generate_seds(z, beta)
            scaled_fluxes = source_injector.scale_seds_to_muv(wavelengths, fluxes, Muv_sample, z)
            filter_fluxes = source_injector.calculate_fluxes(wavelengths, scaled_fluxes)

            #! Get random positions
            x, y = source_injector.generate_random_positions(image_size)

            #! Get PSF fluxes corresponding to input Muv
            source_injector.get_psf()
            scaled_psfs = source_injector.scale_psf_to_Muv(filter_fluxes, Muv_sample, z)

            #! Inject sources
            wcs = source_injector.inject_sources(image_name, x, y, Muv_sample, z, scaled_psfs)
            
            #! Convert x,y to RA, Dec
            ra, dec = wcs.all_pix2world(x, y, 0)

            #! Save the input values as an astropy table
            t = Table([x, y, ra, dec, Muv_sample, z, beta, filter_fluxes['YJ']], names=('x', 'y', 'RA', 'DEC', 'Muv', 'z', 'beta_slope', 'flux_YJ'))
            table_name = image_name.replace('.fits', '_input_values.fits')
            t.write(str(cat_dir / table_name), overwrite=overwrite)

    if run_se:
        #! Run Source Extractor!

        injected_images = glob.glob(str(injected_dir / '*.fits'))
        source_extractor = SourceExtractor(injected_images)

        # # If overwrite, clear the output catalogues
        if overwrite_cats:
            for file in glob.glob(str(output_cat_dir / '*.fits')):
                os.remove(file)

        #? Batch 
        batches = source_extractor.batch_image_list(batch_size)

        #! Run SE batches on queue!
        source_extractor.execute_se_batches(batch_size,'YJ', 1.8, queue='normal', overwrite=True, check_interval=15)


    








 
