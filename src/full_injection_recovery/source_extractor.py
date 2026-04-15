#!/usr/bin/env python3

"""
Source Extractor helpers for the multi-band injection pipeline.
"""

from pathlib import Path
import os
import subprocess


class SourceExtractor:
    def __init__(self, image_list):
        self.image_list = image_list

    def run_sextractor_dual_image(
        self,
        detection_image_name,
        measurement_image_name,
        apDiameterAS,
        pix_scale_arcsec,
        output_cat_name=None,
        overwrite=True,
    ):
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

        ap_string_pix = f'{apDiameterAS / pix_scale_arcsec:.2f}'
        command = [
            '/mnt/users/videouser/sextractor/bin/sex',
            f'{detection_image},{measurement_image}',
            '-c', str(input_sex),
            '-CATALOG_NAME', str(output_cat),
            '-MAG_ZEROPOINT', '30.0',
            '-PHOT_APERTURES', ap_string_pix,
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
        pix_scale_arcsec,
        overwrite=True,
    ):
        output_catalogues = []
        detection_stem = detection_image_name.replace('.fits', '')

        for measurement_image_name in measurement_image_names:
            measurement_stem = measurement_image_name.replace('.fits', '')
            output_cat_name = f'{measurement_stem}__det_{detection_stem}_cat.fits'
            output_cat = self.run_sextractor_dual_image(
                detection_image_name=detection_image_name,
                measurement_image_name=measurement_image_name,
                apDiameterAS=apDiameterAS,
                pix_scale_arcsec=pix_scale_arcsec,
                output_cat_name=output_cat_name,
                overwrite=overwrite,
            )
            output_catalogues.append(output_cat)

        return output_catalogues
