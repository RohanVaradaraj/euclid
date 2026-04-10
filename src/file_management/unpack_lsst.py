"""
Extract the science and variance layers from the LSST images.

Created: Tuesday 17th March 2026.
"""

from pathlib import Path
import numpy as np
from astropy.io import fits
import os
import glob
import sys
import traceback

# adjust base path as needed
lsst_dir = Path.home().parents[1] / 'extraspace' / 'varadaraj' / 'lsst'
filters = ['u', 'g', 'r', 'i', 'z', 'y']

def safe_inv_variance(var_array):
    # produce inverse variance as float32; invalid/zero/negative => 0
    var = np.array(var_array, dtype=np.float64)  # work in float64 for safety
    inv = np.zeros_like(var, dtype=np.float32)
    mask = np.isfinite(var) & (var > 0.0)
    inv[mask] = (1.0 / var[mask]).astype(np.float32)
    return inv

for filt in filters:
    folder = lsst_dir / filt
    if not folder.exists():
        print(f"Skipping missing folder: {folder}", file=sys.stderr)
        continue

    for file_path in glob.glob(str(folder / '*.fits')):
        try:
            p = Path(file_path)
            tile_id = p.stem
            print(f"Processing {tile_id}")

            with fits.open(p, mode='readonly') as hdul:
                sci_hdu = hdul[1]
                var_hdu = hdul[3]

                sci_data = sci_hdu.data
                sci_header = sci_hdu.header.copy()

                var_data = var_hdu.data
                var_header = var_hdu.header.copy()

            # compute inverse variance safely
            invvar = safe_inv_variance(var_data)

            # prepare output filenames
            out_sci = p.with_name(f"{tile_id}_SCI.fits")
            out_wht = p.with_name(f"{tile_id}_WHT.fits")

            # write SCI file (PrimaryHDU with original sci header)
            pri_sci = fits.PrimaryHDU(data=np.array(sci_data, dtype=np.float32), header=sci_header)
            # ensure extension name or other keys not conflicting; keep header mostly intact
            pri_sci.header['EXTNAME'] = 'SCI'

            pri_sci.writeto(out_sci, overwrite=True)

            # write WHT file (inverse variance) using VAR header as base, mark EXTNAME
            pri_wht = fits.PrimaryHDU(data=invvar, header=var_header)
            pri_wht.header['EXTNAME'] = 'WHT'
            pri_wht.writeto(out_wht, overwrite=True)

            # if both files exist and are non-empty, delete original
            if out_sci.exists() and out_wht.exists() and out_sci.stat().st_size > 0 and out_wht.stat().st_size > 0:
                os.remove(p)
                print(f"Written {out_sci.name}, {out_wht.name} and removed original {p.name}")
            else:
                print(f"Write failed for {tile_id} — original retained.", file=sys.stderr)

        except Exception:
            print(f"Error processing {file_path}:", file=sys.stderr)
            traceback.print_exc()
            # continue with next file
            continue
