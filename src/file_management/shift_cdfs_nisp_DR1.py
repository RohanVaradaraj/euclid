#!/usr/bin/env python3

"""
Modify Euclid CDFS DR1 NISP headers to align with VISTA.

Created: Thursday 19th March 2026.
"""

from astropy.io import fits
from pathlib import Path
import glob
import os

euclid_dir = Path.home().parents[1] / 'extraspace' / 'varadaraj' / 'euclid' / 'DR1'

fields = ['CDFS1', 'CDFS2', 'CDFS3']
filter_names = ['Y', 'J', 'H']

for field in fields:

    tile_dir = euclid_dir / field

    for filter_name in filter_names:

        sci_file = f'EDFF_Euclid-DR1_NISP-{filter_name}_200mas_{field}.drz.fits'
        wht_file = f'EDFF_Euclid-DR1_NISP-{filter_name}_200mas_{field}.drz.inverse_var.fits'

        sci_path = tile_dir / sci_file
        wht_path = tile_dir / wht_file

        # Read full HDUs into memory so the file handle can close safely
        with fits.open(sci_path, memmap=False) as sci_hdul:
            sci_data = sci_hdul[0].data
            sci_hdr = sci_hdul[0].header.copy()

            # actual file size (bytes)
            actual_size = os.path.getsize(sci_path)

            # expected data size (bytes)
            bitpix = abs(sci_hdr['BITPIX']) // 8
            naxis = sci_hdr['NAXIS']
            shape = [sci_hdr[f'NAXIS{i+1}'] for i in range(naxis)]
            expected_data_size = bitpix
            for dim in shape:
                expected_data_size *= dim

            print(f"{sci_path}")
            print(f"Actual file size:   {actual_size}")
            print(f"Expected data size: {expected_data_size}")

        with fits.open(wht_path, memmap=False) as wht_hdul:
            wht_data = wht_hdul[0].data
            wht_hdr = wht_hdul[0].header.copy()

            # actual file size (bytes)
            actual_size = os.path.getsize(wht_path)

            # expected data size (bytes)
            bitpix = abs(wht_hdr['BITPIX']) // 8
            naxis = wht_hdr['NAXIS']
            shape = [wht_hdr[f'NAXIS{i+1}'] for i in range(naxis)]
            expected_data_size = bitpix
            for dim in shape:
                expected_data_size *= dim
            
            print(f"{wht_path}")
            print(f"Actual file size:   {actual_size}")
            print(f"Expected data size: {expected_data_size}")

        print(f'~~~~~~~~Processing {field} {filter_name}...~~~~~~~~~')

        print('### VISTA REF ###')
        with fits.open(f'/mnt/vardy/vardygroupshare/data/{field}/{field}_YJ.fits') as vista_hdul:
            vista_hdr = vista_hdul[0].header
        print('NAXIS1, NAXIS2, CRPIX1, CRPIX2')
        print(vista_hdr['NAXIS1'], vista_hdr['NAXIS2'], vista_hdr['CRPIX1'], vista_hdr['CRPIX2'])
        print('CD1_1, CD1_2, CD2_1, CD2_2')
        print(vista_hdr['CD1_1'], vista_hdr['CD1_2'], vista_hdr['CD2_1'], vista_hdr['CD2_2'])

        # In the sci and wht header, modify CD1_1, CD1_2, CD2_1, CD2_2 to match VISTA
        sci_hdr['CD1_1'] = vista_hdr['CD1_1']
        sci_hdr['CD1_2'] = vista_hdr['CD1_2']
        sci_hdr['CD2_1'] = vista_hdr['CD2_1']
        sci_hdr['CD2_2'] = vista_hdr['CD2_2']

        wht_hdr['CD1_1'] = vista_hdr['CD1_1']
        wht_hdr['CD1_2'] = vista_hdr['CD1_2']
        wht_hdr['CD2_1'] = vista_hdr['CD2_1']
        wht_hdr['CD2_2'] = vista_hdr['CD2_2']

        print('### SCI HEADER REF ###')
        print('NAXIS1, NAXIS2, CRPIX1, CRPIX2')
        print(sci_hdr['NAXIS1'], sci_hdr['NAXIS2'], sci_hdr['CRPIX1'], sci_hdr['CRPIX2'])
        print('CD1_1, CD1_2, CD2_1, CD2_2')
        print(sci_hdr['CD1_1'], sci_hdr['CD1_2'], sci_hdr['CD2_1'], sci_hdr['CD2_2'])

        print('### WHT HEADER REF ###')
        print('NAXIS1, NAXIS2, CRPIX1, CRPIX2')
        print(wht_hdr['NAXIS1'], wht_hdr['NAXIS2'], wht_hdr['CRPIX1'], wht_hdr['CRPIX2'])
        print('CD1_1, CD1_2, CD2_1, CD2_2')
        print(wht_hdr['CD1_1'], wht_hdr['CD1_2'], wht_hdr['CD2_1'], wht_hdr['CD2_2'])

        fits.PrimaryHDU(data=sci_data, header=sci_hdr).writeto(sci_path, overwrite=True)
        fits.PrimaryHDU(data=wht_data, header=wht_hdr).writeto(wht_path, overwrite=True)

        print(f'Saved modified headers to files for {field} {filter_name}')