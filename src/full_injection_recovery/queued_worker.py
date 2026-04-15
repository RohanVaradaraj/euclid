#!/usr/bin/env python3

"""
Queue worker: one multi-band image set -> shared cutouts -> injection -> truth catalogue.
"""

from __future__ import annotations

import argparse
from concurrent.futures import ThreadPoolExecutor
import json
import traceback
from pathlib import Path

import numpy as np
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.table import Table
from astropy.wcs import WCS

from luminosity_function import LuminosityFunction
from source_extractor import SourceExtractor
from source_injector import SourceInjector
from utils import load_config


def _resolve_image_directory(field_name: str, directory: str) -> Path:
    """
    Resolve the image directory based on the provided field name and images.lis directory entry.
    """
    data_root = Path.cwd().parents[3] / 'data'
    if directory == 'here':
        return data_root / field_name

    candidate = Path(directory)
    if candidate.is_absolute():
        return candidate

    return (data_root / field_name / directory).resolve()


def _open_image_and_wcs(image_dir: Path, image_name: str):
    with fits.open(image_dir / image_name, memmap=True) as hdu:
        data = hdu[0].data
        header = hdu[0].header.copy()
        wcs = WCS(header)

    return data, header, wcs


def _choose_shared_cutout_coord(
    row: dict,
    field_name: str,
    image_size_arcmin: float,
    pix_scale: float,
):
    image_dir = _resolve_image_directory(field_name, row['directory'])
    data, _, wcs = _open_image_and_wcs(image_dir, row['image'])

    image_size_pix = int(round(image_size_arcmin * 60.0 / pix_scale))
    x_margin = image_size_pix // 2 + 100
    y_margin = image_size_pix // 2 + 200

    if data.shape[1] <= 2 * x_margin or data.shape[0] <= 2 * y_margin:
        raise ValueError(
            f"Image {row['image']} is too small for a {image_size_arcmin} arcmin cutout."
        )

    x = np.random.randint(x_margin, data.shape[1] - x_margin)
    y = np.random.randint(y_margin, data.shape[0] - y_margin)
    coord = SkyCoord.from_pixel(x, y, wcs)

    return coord, int(x), int(y)


def _normalise_image_stem(image_name: str) -> str:
    return image_name.replace('.fits', '')


def _resolve_detection_filter_name(base_image_name: str, configured_filters):
    base_stem = _normalise_image_stem(Path(base_image_name).name)
    for filter_name in configured_filters:
        if base_stem == filter_name or base_stem.startswith(f'{filter_name}_'):
            return filter_name
    raise ValueError(
        f"Could not match source_injection.base_image={base_image_name} to configured filters {configured_filters}"
    )


def _build_cutout_header(parent_header, cutout_wcs: WCS, coord: SkyCoord):
    cutout_header = cutout_wcs.to_header()
    cutout_header['EXPTIME'] = parent_header.get('EXPTIME', 1.0)
    cutout_header['GAIN'] = parent_header.get('GAIN', 1.0)
    cutout_header['SATURATE'] = parent_header.get('SATURATE', 1.0)
    cutout_header['CUTRA'] = float(coord.ra.deg)
    cutout_header['CUTDEC'] = float(coord.dec.deg)

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

    return cutout_header


def _make_single_cutout(
    image_dir: Path,
    image_name: str,
    weight_name: str,
    image_size_arcmin: float,
    pix_scale: float,
    coord: SkyCoord,
):
    weight_path = Path.cwd() / 'images' / 'cutouts' / 'weights'
    weight_path.mkdir(parents=True, exist_ok=True)

    data, header, wcs = _open_image_and_wcs(image_dir, image_name)

    image_size_pix = int(round(image_size_arcmin * 60.0 / pix_scale))
    x_parent, y_parent = wcs.world_to_pixel(coord)
    x_parent = int(round(float(x_parent)))
    y_parent = int(round(float(y_parent)))
    half_size = image_size_pix // 2
    x_start = x_parent - half_size
    y_start = y_parent - half_size
    x_end = x_start + image_size_pix
    y_end = y_start + image_size_pix
    if x_start < 0 or y_start < 0 or x_end > data.shape[1] or y_end > data.shape[0]:
        raise ValueError(
            f"Cutout for {image_name} at ({x_parent}, {y_parent}) with size {image_size_pix} exceeds image bounds."
        )

    cutout_data = np.array(data[y_start:y_end, x_start:x_end], copy=True)
    cutout_wcs = wcs.slice((slice(y_start, y_end), slice(x_start, x_end)))
    cutout_header = _build_cutout_header(header, cutout_wcs, coord)

    with fits.open(image_dir / weight_name, memmap=True) as hdu_w:
        weight_data = hdu_w[0].data[y_start:y_end, x_start:x_end]
        weight_cutout = np.array(weight_data, copy=True)

    stem = image_name.replace('.fits', '')
    cutout_name = (
        f'{stem}_cutout_{x_parent}_{y_parent}_'
        f'{int(image_size_pix)}_pix_{int(image_size_arcmin)}_arcmin.fits'
    )
    weight_cutout_name = cutout_name.replace('.fits', '_wht.fits')

    fits.PrimaryHDU(weight_cutout, header=cutout_header).writeto(weight_path / weight_cutout_name, overwrite=True)

    return {
        'cutout_name': cutout_name,
        'weight_cutout_name': weight_cutout_name,
        'cutout_wcs': cutout_wcs,
        'cutout_data': cutout_data,
        'cutout_header': cutout_header.copy(),
    }


def _build_cutout_for_row(row, field_name, image_size_arcmin, pix_scale, coord):
    image_dir = _resolve_image_directory(field_name, row['directory'])
    cutout_product = _make_single_cutout(
        image_dir=image_dir,
        image_name=row['image'],
        weight_name=row['weight'],
        image_size_arcmin=image_size_arcmin,
        pix_scale=pix_scale,
        coord=coord,
    )
    return {
        'filter_name': row['name'],
        'cutout_name': cutout_product['cutout_name'],
        'weight_cutout_name': cutout_product['weight_cutout_name'],
        'cutout_data': cutout_product['cutout_data'],
        'cutout_header': cutout_product['cutout_header'],
        'cutout_wcs': cutout_product['cutout_wcs'],
    }


def _build_truth_catalog(
    x_positions,
    y_positions,
    wcs,
    source_types,
    Muv_sample,
    z,
    beta,
    filter_fluxes,
    filters,
):
    ra, dec = wcs.all_pix2world(x_positions, y_positions, 0)
    columns = [x_positions, y_positions, ra, dec, source_types, Muv_sample, z, beta]
    names = ['x', 'y', 'RA', 'DEC', 'type', 'Muv', 'z', 'beta_slope']

    for filter_name in filters:
        columns.append(np.asarray(filter_fluxes[filter_name]))
        names.append(f'model_flux_{filter_name}')

    return Table(columns, names=names)


def run_task(task: dict):
    config = load_config(task['config_file'])
    lf_config = config['luminosity_function']
    injection_config = dict(config['source_injection'])
    field_name = config['field']['name']
    pix_scale = float(config['field']['pix_scale'])
    processing_cfg = config.get('processing', {})
    source_extraction_cfg = config.get('source_extraction', {})
    cutout_workers = max(1, int(processing_cfg.get('cutout_workers', 1)))

    batch_rows = list(task['rows'])
    if len(batch_rows) == 0:
        raise ValueError('Task contained no images to process.')

    row_by_name = {str(row['name']): row for row in batch_rows}
    all_configured_filters = list(injection_config['filters'])
    all_configured_zeropoints = list(injection_config['zeropoints'])
    configured_filters = [flt for flt in all_configured_filters if flt in row_by_name]
    if not configured_filters:
        raise ValueError('No overlap between configured filters and task rows.')

    batch_rows = [row_by_name[flt] for flt in configured_filters]
    injection_config['filters'] = configured_filters
    injection_config['zeropoints'] = [
        all_configured_zeropoints[all_configured_filters.index(flt)]
        for flt in configured_filters
    ]
    detection_filter_name = _resolve_detection_filter_name(
        injection_config['base_image'],
        configured_filters,
    )

    image_size_arcmin = float(injection_config['image_size_arcmin'])
    shared_coord, shared_x, shared_y = _choose_shared_cutout_coord(
        batch_rows[0],
        field_name,
        image_size_arcmin,
        pix_scale,
    )

    if cutout_workers == 1:
        cutout_results = [
            _build_cutout_for_row(row, field_name, image_size_arcmin, pix_scale, shared_coord)
            for row in batch_rows
        ]
    else:
        max_workers = min(cutout_workers, len(batch_rows))
        with ThreadPoolExecutor(max_workers=max_workers) as pool:
            cutout_results = list(
                pool.map(
                    lambda row: _build_cutout_for_row(
                        row,
                        field_name,
                        image_size_arcmin,
                        pix_scale,
                        shared_coord,
                    ),
                    batch_rows,
                )
            )

    cutout_info = [
        {
            'filter_name': row['filter_name'],
            'cutout_name': row['cutout_name'],
            'weight_cutout_name': row['weight_cutout_name'],
            'cutout_data': row['cutout_data'],
            'cutout_header': row['cutout_header'],
        }
        for row in cutout_results
    ]
    detection_result = next(row for row in cutout_results if row['filter_name'] == detection_filter_name)
    reference_wcs = detection_result['cutout_wcs']

    luminosity_function = LuminosityFunction(lf_config)
    Muv_sample = luminosity_function.sample_luminosities(uniform=True)
    ucd_cfg = config.get('ultra_cool_dwarfs', {})
    n_ucd_per_class = int(ucd_cfg.get('n_samples', len(Muv_sample)))
    ucd_seed_samples = np.arange(n_ucd_per_class)

    source_injector = SourceInjector(samples=Muv_sample, params=injection_config)
    z, beta = source_injector.draw_parameters()

    wavelengths, fluxes = source_injector.generate_seds(z, beta)
    scaled_fluxes = source_injector.scale_seds_to_muv(wavelengths, fluxes, Muv_sample, z)
    lbg_filter_fluxes = source_injector.calculate_fluxes(wavelengths, scaled_fluxes, average=False)

    ucd_injector = SourceInjector(samples=ucd_seed_samples, params=injection_config)
    ucd_injector.load_ucd_templates()
    ucd_mags = ucd_injector.sample_ucd_mags()
    ucd_types_by_class = ucd_injector.draw_ucd_types()
    ucd_wavelengths, ucd_scaled_seds = ucd_injector.scale_ucds_to_mags(
        ucd_mags,
        ucd_types_by_class,
        pix_scale,
    )
    ucd_filter_fluxes = ucd_injector.calculate_ucd_fluxes(
        ucd_wavelengths,
        ucd_scaled_seds,
        average=False,
    )

    ucd_template_types = np.asarray([
        template_type
        for ucd_type in ucd_injector.ucd_types
        for template_type in ucd_types_by_class[ucd_type]
    ])
    n_ucd_total = len(ucd_template_types)
    ucd_placeholder = np.full(n_ucd_total, -99.0)

    filter_fluxes = {
        filter_name: np.concatenate(
            [
                np.asarray(lbg_filter_fluxes[filter_name], dtype=float),
                np.asarray(ucd_filter_fluxes[filter_name], dtype=float),
            ]
        )
        for filter_name in configured_filters
    }

    x_lbg, y_lbg = source_injector.generate_random_positions(image_size_arcmin, pix_scale)
    ucd_position_injector = SourceInjector(samples=np.arange(n_ucd_total), params=injection_config)
    x_ucd, y_ucd = ucd_position_injector.generate_random_positions(image_size_arcmin, pix_scale)
    x_positions = np.concatenate([x_lbg, x_ucd])
    y_positions = np.concatenate([y_lbg, y_ucd])

    source_types = np.asarray(
        np.concatenate([
            np.full(len(Muv_sample), 'LBG', dtype='U3'),
            ucd_template_types.astype(str),
        ]),
        dtype='U8',
    )
    catalog_Muv = np.concatenate([np.asarray(Muv_sample, dtype=float), ucd_placeholder])
    catalog_z = np.concatenate([np.asarray(z, dtype=float), ucd_placeholder])
    catalog_beta = np.concatenate([np.asarray(beta, dtype=float), ucd_placeholder])

    source_injector.get_psf(field_name)
    scaled_psfs = source_injector.scale_PSFs_to_model_fluxes(filter_fluxes, pix_scale)

    for row in cutout_info:
        source_injector.inject_sources(
            image_name=row['cutout_name'],
            x_positions=x_positions,
            y_positions=y_positions,
            Muv_array=catalog_Muv,
            z_array=catalog_z,
            scaled_psfs=scaled_psfs,
            image=row['cutout_data'],
            header=row['cutout_header'],
            filter_name=row['filter_name'],
        )

    detection_row = next(row for row in cutout_info if row['filter_name'] == detection_filter_name)
    measurement_image_names = [row['cutout_name'] for row in cutout_info]
    source_extractor = SourceExtractor(measurement_image_names)
    source_extractor.run_dual_mode_for_set(
        detection_image_name=detection_row['cutout_name'],
        measurement_image_names=measurement_image_names,
        apDiameterAS=float(source_extraction_cfg.get('aperture_diameter_arcsec', 2.0)),
        pix_scale_arcsec=pix_scale,
        overwrite=bool(source_extraction_cfg.get('overwrite', True)),
    )

    truth_table = _build_truth_catalog(
        x_positions=x_positions,
        y_positions=y_positions,
        wcs=reference_wcs,
        source_types=source_types,
        Muv_sample=catalog_Muv,
        z=catalog_z,
        beta=catalog_beta,
        filter_fluxes=filter_fluxes,
        filters=configured_filters,
    )

    input_cat_dir = Path.cwd() / 'catalogues' / 'input'
    input_cat_dir.mkdir(parents=True, exist_ok=True)
    table_name = f"{task['set_id']}_truth_catalog.fits"
    truth_table.write(str(input_cat_dir / table_name), overwrite=True)

    summary = {
        'set_id': task['set_id'],
        'shared_cutout_center_ra_deg': float(shared_coord.ra.deg),
        'shared_cutout_center_dec_deg': float(shared_coord.dec.deg),
        'reference_image_x': shared_x,
        'reference_image_y': shared_y,
        'detection_cutout': detection_row['cutout_name'],
        'n_filters': len(configured_filters),
        'n_sources': len(x_positions),
        'n_lbg_sources': len(Muv_sample),
        'n_ucd_sources': n_ucd_total,
        'truth_catalog': table_name,
        'cutouts': [
            {
                'filter_name': row['filter_name'],
                'cutout_name': row['cutout_name'],
                'weight_cutout_name': row['weight_cutout_name'],
            }
            for row in cutout_info
        ],
    }

    output_dir = Path.cwd() / 'catalogues' / 'input'
    with open(output_dir / f"{task['set_id']}_summary.json", 'w') as f:
        json.dump(summary, f, indent=2)

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
