#!/usr/bin/env python3
"""
Euclid DR1 images are too big. So match to VISTA in chunks/tiles.
"""

#!/usr/bin/env python3
"""
Chunk up Euclid DR1 mosaics to match VISTA in manageable tiles.
Reprojects each tile using SWarp.
"""

from pathlib import Path
from astropy.io import fits
from astropy.nddata import Cutout2D
from astropy.wcs import WCS
import numpy as np
import os
import subprocess

# Parameters
fields = ['COSMOS']
filter_names = ['VIS', 'Y', 'J', 'H']
CHUNK_SIZE = 20000  # pixels per tile

# Directories
data_dir = Path.cwd().parents[3] / 'data'
euclid_dir = Path.home() / 'euclid' / 'COSMOS'
junk_dir = euclid_dir / 'junk'
junk_dir.mkdir(exist_ok=True)

# SWarp configuration
swarp_base = '~/swarp/bin/swarp'
swarp_conf = 'regrid.swarp'

extra_opts_sci = (
    '-COMBINE_BUFSIZE 2048 -COMBINE_TYPE WEIGHTED '
    '-VMEM_MAX 16384 -VMEM_DIR . -MEM_MAX 2048'
)
extra_opts_wht = (
    '-WEIGHT_TYPE NONE -COMBINE_BUFSIZE 2048 -COMBINE_TYPE AVERAGE '
    '-VMEM_MAX 16384 -VMEM_DIR . -MEM_MAX 2048'
)

# Loop over fields and filters
for field in fields:
    for filt in filter_names:
        sci_file = euclid_dir / f'{field}_{filt}_DR1.fits'
        wht_file = euclid_dir / f'{field}_{filt}_DR1_WHT.fits'

        print(f'\n--- Processing {sci_file.name} ---')

        with fits.open(sci_file, memmap=True) as hdul:
            data = hdul[0].data
            wcs = WCS(hdul[0].header)
            ny, nx = data.shape

        nx_tiles = int(np.ceil(nx / CHUNK_SIZE))
        ny_tiles = int(np.ceil(ny / CHUNK_SIZE))

        for iy in range(ny_tiles):
            for ix in range(nx_tiles):
                x0, x1 = ix * CHUNK_SIZE, min((ix + 1) * CHUNK_SIZE, nx)
                y0, y1 = iy * CHUNK_SIZE, min((iy + 1) * CHUNK_SIZE, ny)

                # --- Cut out science tile ---
                cut_sci = Cutout2D(data, position=((x0 + x1)/2, (y0 + y1)/2),
                                   size=(y1 - y0, x1 - x0), wcs=wcs)
                tile_sci = junk_dir / f'{field}_{filt}_DR1_tile_{iy}_{ix}.fits'
                fits.writeto(tile_sci, cut_sci.data.astype('float32'),
                             cut_sci.wcs.to_header(), overwrite=True)

                # --- Cut out weight tile ---
                with fits.open(wht_file, memmap=True) as hdul_wht:
                    wht_data = hdul_wht[0].data
                    cut_wht = Cutout2D(wht_data,
                                       position=((x0 + x1)/2, (y0 + y1)/2),
                                       size=(y1 - y0, x1 - x0), wcs=wcs)
                tile_wht = junk_dir / f'{field}_{filt}_DR1_wht_tile_{iy}_{ix}.fits'
                fits.writeto(tile_wht, cut_wht.data.astype('float32'),
                             cut_sci.wcs.to_header(), overwrite=True)

                # --- Reproject science tile ---
                out_sci = junk_dir / f'{field}_{filt}_resamp_tile_{iy}_{ix}.fits'
                cmd_sci = (
                    f"{swarp_base} {tile_sci} "
                    f"-WEIGHT_IMAGE {tile_wht} "
                    f"-c {swarp_conf} "
                    f"-IMAGEOUT_NAME {out_sci} "
                    f"{extra_opts_sci}"
                )
                print(f"Running SWarp science tile {iy},{ix} ...")
                subprocess.run(cmd_sci, shell=True, check=True)

                # --- Reproject weight tile ---
                out_wht = junk_dir / f'{field}_{filt}_resamp_wht_tile_{iy}_{ix}.fits'
                cmd_wht = (
                    f"{swarp_base} {tile_wht} "
                    f"-c {swarp_conf} "
                    f"-IMAGEOUT_NAME {out_wht} "
                    f"{extra_opts_wht}"
                )
                print(f"Running SWarp weight tile {iy},{ix} ...")
                subprocess.run(cmd_wht, shell=True, check=True)

                # --- Clean up temporary cutouts ---
                print('Deleting temporary input tiles...')
                os.remove(tile_sci)
                os.remove(tile_wht)

        print(f"Finished tiling {field}_{filt}")

print("All tiles reprojected. Mosaic with SWarp using -COMBINE_TYPE WEIGHTED.")
