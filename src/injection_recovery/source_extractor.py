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

        











