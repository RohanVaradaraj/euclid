#!/usr/bin/env python3

"""
Contains source extraction class.

Created: Wednesday 4th December 2024.
"""

import numpy as np
from astropy.io import fits
from pathlib import Path
import time
import os
import subprocess


class SourceExtractor:
    def __init__(self, image_list):
        self.image_list = image_list

    
    def batch_image_list(self, batch_size):
        """
        Batch the image list to run SE in bathes
        """

        batches = [self.image_list[i:i+batch_size] for i in range(0, len(self.image_list), batch_size)]
        return batches

    

    def run_se(self, filter_name, image_name, apDiameterAS, overwrite=True):
        """
        Generate a shell script to run Source Extractor for a given detection and measurement filter.
        """
        image_dir = Path.cwd() / 'images' / 'injected'
        catalogue_dir = Path.cwd() / 'catalogues' / 'output'
        shell_dir = Path('./shell_scripts')  # Directory for shell scripts
        shell_dir.mkdir(exist_ok=True)  # Create it if it doesn't exist

        input_sex = Path.cwd() / 'video_mine.sex'
        os.environ['EXTRACTOR_DIR'] = '/mnt/users/videouser/sextractor/share/sextractor'

        with fits.open(image_dir / image_name) as hdulist:
            pixScale = abs(hdulist[0].header['CD2_2']) * 3600.0  # Convert degrees to arcseconds/pixel

        apStringPix = f'{apDiameterAS / pixScale:.2f}'
        output_cat = image_name.split('/')[-1].split('.fits')[0] + '_cat.fits'

        if os.path.isfile(catalogue_dir / output_cat) and not overwrite:
            print(f"Catalog {output_cat} exists. Skipping.")
            return None

        weight_dir = Path.cwd() / 'images' / 'cutouts' / 'weights'
        weight_name = image_name.split('/')[-1].split('.fits')[0] + '_wht.fits'

        command = (
            f"/mnt/users/videouser/sextractor/bin/sex  {str(image_dir / image_name)} "
            f"-c {str(input_sex)} "
            f"-CATALOG_NAME {str(catalogue_dir / output_cat)} "
            f"-MAG_ZEROPOINT 30.0 "
            f"-PHOT_APERTURES {apStringPix} "
            f"-CHECKIMAGE_TYPE SEGMENTATION "
            f"-CHECKIMAGE_NAME seg.fits "
            f"-WEIGHT_TYPE MAP_WEIGHT "
            f"-WEIGHT_IMAGE {weight_dir / weight_name} "
        )

        shellFile = shell_dir / f"run_se_{image_name.split('/')[-1].split('.fits')[0]}.sh"
        with open(shellFile, 'w') as f:
            f.write("#!/bin/bash\n")
            f.write(command + "\n")

        os.system(f"chmod u+x {shellFile}")
        #print(f"Shell script generated: {shellFile}")

        return shellFile


    def run_sextractor_dual_image(
        self,
        detection_image_name,
        measurement_image_name,
        apDiameterAS,
        output_cat_name=None,
        overwrite=True,
    ):
        """
        Run SExtractor in dual-image mode using an injected detection image and
        an injected measurement image. The corresponding weight maps are taken
        from the cutout weight directory.
        """
        image_dir = Path.cwd() / 'images' / 'injected'
        weight_dir = Path.cwd() / 'images' / 'cutouts' / 'weights'
        catalogue_dir = Path.cwd() / 'catalogues' / 'output'
        catalogue_dir.mkdir(parents=True, exist_ok=True)

        input_sex = Path.cwd() / 'video_mine.sex'
        os.environ['EXTRACTOR_DIR'] = '/mnt/users/videouser/sextractor/share/sextractor'

        detection_image = image_dir / detection_image_name
        measurement_image = image_dir / measurement_image_name
        detection_weight = weight_dir / detection_image_name.replace('.fits', '_wht.fits')
        measurement_weight = weight_dir / measurement_image_name.replace('.fits', '_wht.fits')

        if output_cat_name is None:
            output_cat_name = measurement_image_name.replace('.fits', '_cat.fits')

        output_cat = catalogue_dir / output_cat_name
        if output_cat.exists() and not overwrite:
            return output_cat

        with fits.open(measurement_image) as hdulist:
            pixScale = abs(hdulist[0].header['CD2_2']) * 3600.0

        apStringPix = f'{apDiameterAS / pixScale:.2f}'

        command = [
            '/mnt/users/videouser/sextractor/bin/sex',
            f'{detection_image},{measurement_image}',
            '-c', str(input_sex),
            '-CATALOG_NAME', str(output_cat),
            '-MAG_ZEROPOINT', '30.0',
            '-PHOT_APERTURES', apStringPix,
            '-CHECKIMAGE_TYPE', 'NONE',
            '-WEIGHT_TYPE', 'MAP_WEIGHT,MAP_WEIGHT',
            '-WEIGHT_IMAGE', f'{detection_weight},{measurement_weight}',
        ]

        subprocess.run(command, check=True)
        return output_cat


    def run_dual_mode_for_set(
        self,
        detection_image_name,
        measurement_image_names,
        apDiameterAS,
        overwrite=True,
    ):
        """
        Run dual-image mode for every measurement image in a multi-band set.
        """
        output_catalogues = []
        detection_stem = detection_image_name.replace('.fits', '')

        for measurement_image_name in measurement_image_names:
            measurement_stem = measurement_image_name.replace('.fits', '')
            output_cat_name = f'{measurement_stem}__det_{detection_stem}_cat.fits'
            output_cat = self.run_sextractor_dual_image(
                detection_image_name=detection_image_name,
                measurement_image_name=measurement_image_name,
                apDiameterAS=apDiameterAS,
                output_cat_name=output_cat_name,
                overwrite=overwrite,
            )
            output_catalogues.append(output_cat)

        return output_catalogues



    def execute_se_batches(self, batch_size, filter_name, apDiametersAS, queue='normal', overwrite=True, check_interval=10):
        """
        Loop through image batches, create shell scripts, and launch jobs. Wait for all jobs in a batch to complete before queuing the next.

        Parameters:
            batch_size (int): Number of images per batch.
            filter_name (str): Filter name for Source Extractor.
            apDiametersAS (list): List of aperture diameters in arcseconds.
            inputSex (str): Path to the SExtractor configuration file.
            outputDir (str): Output directory for the catalog.
            shellDir (str): Directory to save shell scripts.
            queue (str): Queue to use for job submission.
            overwrite (bool): Overwrite existing catalogs if True.
            check_interval (int): Time interval in seconds to check for completion.
        """

        image_dir = Path.cwd() / 'images' / 'injected'
        catalogue_dir = Path.cwd() / 'catalogues' / 'output'

        batches = self.batch_image_list(batch_size)

        output_dir = catalogue_dir

        for batch_num, batch in enumerate(batches):
            print(f"Processing batch {batch_num + 1}/{len(batches)}")
            
            # Submit jobs for all images in the batch
            shell_files = []
            for image_name in batch:
                print(image_name)
                shell_file = self.run_se(filter_name, image_name, apDiametersAS, overwrite)
                if shell_file:
                    shell_files.append(shell_file)
                    os.system(f'addqueue -c {image_name} -m 2 -q {queue} ./{shell_file}')

            # Wait for all catalogs in the batch to be created
            print(f"Waiting for batch {batch_num + 1} to complete...")
            incomplete = True
            while incomplete:
                incomplete = False
                for image_name in batch:
                    expected_cat = output_dir / (image_name.split('/')[-1].split('.fits')[0] + '_cat.fits')
                    if not expected_cat.is_file():
                        incomplete = True
                        break
                if incomplete:
                    print(f"Some catalogs for batch {batch_num + 1} are not ready. Checking again in {check_interval} seconds...")
                    time.sleep(check_interval)

            print(f"Batch {batch_num + 1} completed. Proceeding to the next batch.")

        










