#!/usr/bin/env python3
"""
Mosaic together the SWarped Euclid DR1 chunked tiles.
Recreates full reprojected images and weight maps aligned to the VISTA reference grid.

Created: Sunday 10th November 2025
"""

from pathlib import Path
import os
import glob
import subprocess

instrument = 'euclid'
fields = ['COSMOS']
filter_names = ['VIS', 'Y', 'J', 'H']   # Euclid filters

# Paths
data_dir = Path.cwd().parents[3] / 'data'
euclid_dir = Path.home() / 'euclid' / 'COSMOS'
swarp_base = '~/swarp/bin/swarp'
swarp_conf = 'regrid.swarp'

# SWarp settings
opts_sci = (
    '-COMBINE_TYPE WEIGHTED '
    '-VMEM_MAX 16384 -VMEM_DIR . '
    '-MEM_MAX 2048 -COMBINE_BUFSIZE 2048 '
)
opts_wht = (
    '-WEIGHT_TYPE NONE '
    '-COMBINE_TYPE AVERAGE '
    '-VMEM_MAX 16384 -VMEM_DIR . '
    '-MEM_MAX 2048 -COMBINE_BUFSIZE 2048 '
)

for field in fields:
    for filt in filter_names:
        print(f'\n--- Mosaicking {field} {filt} tiles ---')

        # Find all reprojected tile files
        tile_pattern = str(euclid_dir / 'junk' / f'{field}_{filt}_resamp_tile_*.fits')
        wht_pattern  = str(euclid_dir / 'junk' / f'{field}_{filt}_resamp_wht_tile_*.fits')

        tile_list = sorted(glob.glob(tile_pattern))
        wht_list  = sorted(glob.glob(wht_pattern))

        if not tile_list:
            print(f'No science tiles found for {field} {filt}. Skipping...')
            continue
        if not wht_list:
            print(f'No weight tiles found for {field} {filt}. Skipping...')
            continue

        # Output filenames
        out_sci = euclid_dir / f'{field}_{filt}_DR1_resamp_full.fits'
        out_wht = euclid_dir / f'{field}_{filt}_DR1_resamp_full_wht.fits'

        # --- Combine science tiles ---
        print(f'Combining {len(tile_list)} science tiles...')
        cmd_sci = (
            f"{swarp_base} {' '.join(tile_list)} "
            f"-WEIGHT_IMAGE {' '.join(wht_list)} "
            f"-c {swarp_conf} "
            f"-IMAGEOUT_NAME {out_sci} "
            f"{opts_sci}"
        )
        subprocess.run(cmd_sci, shell=True, check=True)

        # --- Combine weight tiles ---
        print(f'Combining {len(wht_list)} weight tiles...')
        cmd_wht = (
            f"{swarp_base} {' '.join(wht_list)} "
            f"-c {swarp_conf} "
            f"-IMAGEOUT_NAME {out_wht} "
            f"{opts_wht}"
        )
        subprocess.run(cmd_wht, shell=True, check=True)

        print(f'Finished mosaicking {field} {filt}')
        print(f' → Science: {out_sci}')
        print(f' → Weight : {out_wht}')

print('\nAll mosaics completed successfully.')
