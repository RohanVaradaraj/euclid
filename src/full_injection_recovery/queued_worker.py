#!/usr/bin/env python3

"""
Queue worker: one images.lis row -> cutout -> inject -> source extraction.
"""

from __future__ import annotations

import argparse
import json
import os
import traceback
from pathlib import Path

import numpy as np
from astropy.io import fits
from astropy.nddata import Cutout2D
from astropy.table import Table
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord

from luminosity_function import LuminosityFunction
from source_injector import SourceInjector
from source_extractor import SourceExtractor
from utils import load_config, filter_files


def _resolve_image_directory(field_name: str, directory: str) -> Path:
    data_root = Path.cwd().parents[3] / 'data'
    if directory == 'here':
        return data_root / field_name

    candidate = Path(directory)
    if candidate.is_absolute():
        return candidate

    return (data_root / field_name / directory).resolve()


def _make_single_cutout(
    image_dir: Path,
    image_name: str,
    weight_name: str,
    image_size_arcmin: float,
    pix_scale: float,
):
    cutout_path = Path.cwd() / 'images' / 'cutouts'
    cutout_path.mkdir(parents=True, exist_ok=True)
    weight_path = cutout_path / 'weights'
    weight_path.mkdir(parents=True, exist_ok=True)

    with fits.open(image_dir / image_name) as hdu:
        data = hdu[0].data
        header = hdu[0].header
        wcs = WCS(header)

    with fits.open(image_dir / weight_name) as hdu_w:
        weight = hdu_w[0].data
        wcs_weight = WCS(hdu_w[0].header)

    image_size_pix = image_size_arcmin * 60.0 / pix_scale

    x = np.random.randint(image_size_pix / 2 + 100, data.shape[1] - (image_size_pix / 2 + 100))
    y = np.random.randint(image_size_pix / 2 + 200, data.shape[0] - (image_size_pix / 2 + 100))

    coord = SkyCoord.from_pixel(x, y, wcs)

    cutout = Cutout2D(data, coord, (image_size_pix, image_size_pix), wcs=wcs)
    weight_cutout = Cutout2D(weight, coord, (image_size_pix, image_size_pix), wcs=wcs_weight)

    cutout_header = cutout.wcs.to_header()
    cutout_header['EXPTIME'] = header.get('EXPTIME', 1.0)
    cutout_header['GAIN'] = header.get('GAIN', 1.0)
    cutout_header['SATURATE'] = header.get('SATURATE', 1.0)

    for bad_key in ['LONPOLE', 'LATPOLE', 'MJDREF']:
        if bad_key in cutout_header:
            cutout_header.remove(bad_key)

    if 'PC1_1' in cutout_header and 'PC2_2' in cutout_header:
        cutout_header['CD1_1'] = cutout_header['PC1_1']
        cutout_header['CD1_2'] = 0.0
        cutout_header['CD2_1'] = 0.0
        cutout_header['CD2_2'] = cutout_header['PC2_2']
        cutout_header.remove('PC1_1')
        cutout_header.remove('PC2_2')

    stem = image_name.replace('.fits', '')
    cutout_name = f'{stem}_cutout_{int(x)}_{int(y)}_{int(image_size_pix)}_pix_{int(image_size_arcmin)}_arcmin.fits'
    weight_cutout_name = cutout_name.replace('.fits', '_wht.fits')

    fits.PrimaryHDU(cutout.data, header=cutout_header).writeto(cutout_path / cutout_name, overwrite=True)
    fits.PrimaryHDU(weight_cutout.data, header=cutout_header).writeto(weight_path / weight_cutout_name, overwrite=True)

    return cutout_name, weight_cutout_name


def _cleanup_batch_image_products(cutout_name: str, weight_cutout_name: str):
    cutout_dir = Path.cwd() / 'images' / 'cutouts'
    injected_dir = Path.cwd() / 'images' / 'injected'

    for fp in [
        cutout_dir / cutout_name,
        cutout_dir / 'weights' / weight_cutout_name,
        injected_dir / cutout_name,
    ]:
        if fp.exists():
            fp.unlink()


def run_task(task: dict):
    config = load_config(task['config_file'])

    lf_config = config['luminosity_function']
    injection_config = dict(config['source_injection'])
    field_name = config['field']['name']
    pix_scale = float(config['field']['pix_scale'])

    image_name = task['image']
    weight_name = task['weight']
    directory = task['directory']
    source_name = task['name']

    image_dir = _resolve_image_directory(field_name, directory)

    # If the image Name is available in filter files, use it as detection/filter key.
    known_filters = filter_files()
    if source_name in known_filters:
        injection_config['filters'] = [source_name]

    image_size_arcmin = float(injection_config['image_size_arcmin'])

    cutout_name, weight_cutout_name = _make_single_cutout(
        image_dir=image_dir,
        image_name=image_name,
        weight_name=weight_name,
        image_size_arcmin=image_size_arcmin,
        pix_scale=pix_scale,
    )

    luminosity_function = LuminosityFunction(lf_config)
    Muv_sample = luminosity_function.sample_luminosities(uniform=True)

    source_injector = SourceInjector(samples=Muv_sample, params=injection_config)
    z, beta = source_injector.draw_parameters()

    wavelengths, fluxes = source_injector.generate_seds(z, beta)
    scaled_fluxes = source_injector.scale_seds_to_muv(wavelengths, fluxes, Muv_sample, z)
    filter_fluxes = source_injector.calculate_fluxes(wavelengths, scaled_fluxes)

    x, y = source_injector.generate_random_positions(image_size_arcmin, pix_scale)

    source_injector.get_psf()
    scaled_psfs = source_injector.scale_psf_to_Muv(filter_fluxes, Muv_sample, z, pix_scale)
    wcs = source_injector.inject_sources(cutout_name, x, y, Muv_sample, z, scaled_psfs)

    ra, dec = wcs.all_pix2world(x, y, 0)

    flux_colname = injection_config['filters'][0]
    t = Table(
        [x, y, ra, dec, Muv_sample, z, beta, filter_fluxes[flux_colname]],
        names=('x', 'y', 'RA', 'DEC', 'Muv', 'z', 'beta_slope', f'flux_{flux_colname}'),
    )

    input_cat_dir = Path.cwd() / 'catalogues' / 'input'
    input_cat_dir.mkdir(parents=True, exist_ok=True)
    table_name = cutout_name.replace('.fits', '_input_values.fits')
    t.write(str(input_cat_dir / table_name), overwrite=True)

    se = SourceExtractor([cutout_name])
    shell_file = se.run_se(
        filter_name=source_name,
        image_name=cutout_name,
        apDiameterAS=float(task.get('se_aperture_diameter_arcsec', 2.0)),
        overwrite=bool(task.get('se_overwrite', True)),
    )

    if shell_file is None:
        raise RuntimeError(f'No SExtractor shell file generated for {cutout_name}')

    rc = os.system(f'bash {shell_file}')
    if rc != 0:
        raise RuntimeError(f'SExtractor failed for {cutout_name}, return code {rc}')

    output_cat = Path.cwd() / 'catalogues' / 'output' / cutout_name.replace('.fits', '_cat.fits')
    if not output_cat.exists():
        raise RuntimeError(f'Expected output catalogue missing: {output_cat}')

    if bool(config.get('processing', {}).get('cleanup_batch_images', True)):
        _cleanup_batch_image_products(cutout_name, weight_cutout_name)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--task-file', required=True)
    parser.add_argument('--done-file', required=True)
    args = parser.parse_args()

    done_file = Path(args.done_file)
    fail_file = done_file.with_suffix('.failed')

    with open(args.task_file, 'r') as f:
        task = json.load(f)

    try:
        run_task(task)
        done_file.write_text('ok\n')
    except Exception:
        fail_file.write_text(traceback.format_exc())
        raise


if __name__ == '__main__':
    main()
