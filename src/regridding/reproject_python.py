#!/usr/bin/env python3
"""
Reproject Euclid images onto the UltraVISTA reference grid
using the `reproject` package, since swarp struggled with memory issues.

Created: 10 November 2025
"""

from pathlib import Path
from astropy.io import fits
from astropy.wcs import WCS
from reproject import reproject_interp
import numpy as np

#! Config
instrument = 'euclid'
fields = ['COSMOS']
#filter_names = ['VIS', 'Y', 'J', 'H']
filter_names = ['H']

#! Directories
data_dir = Path.cwd().parents[3] / 'data'
euclid_dir = Path.home() / 'euclid' / 'COSMOS'
reference_dir = data_dir / 'COSMOS'
reference_file = reference_dir / 'UVISTA_Y_DR6.fits'  # reference wcs file to regrid to
output_suffix = '_DR1_resamp.fits' # new name for resamped files

# ! load UVISTA header and wcs
ref_hdu = fits.open(reference_file)[0]
ref_header = ref_hdu.header
ref_wcs = WCS(ref_header)
ref_shape = (ref_header['NAXIS2'], ref_header['NAXIS1'])
print(f"Reference grid loaded: {reference_file}, size {ref_shape} pixels")

#! Loop through fields and filters
for field in fields:
    for filt in filter_names:
        print(f"\nProcessing {field} {filt}...")

        # Input files
        sci_in = euclid_dir / f'{field}_{filt}_DR1.fits'
        wht_in = euclid_dir / f'{field}_{filt}_DR1_WHT.fits'

        if not sci_in.exists():
            print(f"Missing science image {sci_in.name}, skipping.")
            continue

        sci_out = euclid_dir / f'{field}_{filt}{output_suffix}'
        wht_out = euclid_dir / f'{field}_{filt}_WHT{output_suffix}'

        # Load science image
        with fits.open(sci_in, memmap=True) as hdul:
            sci_data = hdul[0].data
            sci_header = hdul[0].header
            sci_wcs = WCS(sci_header)

        print(f"  Reprojecting {sci_in.name} to {sci_out.name}")

        # Reproject science frame
        #? Block size does the resampling in chunks. size is some factor of the full UVISTA image size.
        reproj_sci, footprint_sci = reproject_interp((sci_data, sci_wcs),
                                                     ref_wcs,
                                                     shape_out=ref_shape, block_size=(459, 1871),
                                                     return_footprint=True)

        # Write reprojected science image
        fits.writeto(sci_out, reproj_sci.astype('float32'), ref_header, overwrite=True)

        # Reproject weight map if available
        if wht_in.exists():
            print(f"  Reprojecting {wht_in.name} to {wht_out.name}")
            with fits.open(wht_in, memmap=True) as hdul:
                wht_data = hdul[0].data
                wht_wcs = WCS(hdul[0].header)

            reproj_wht, footprint_wht = reproject_interp((wht_data, wht_wcs),
                                                         ref_wcs,
                                                         shape_out=ref_shape, block_size=(459, 1871),
                                                         return_footprint=True)

            fits.writeto(wht_out, reproj_wht.astype('float32'), ref_header, overwrite=True)
        else:
            print("  No weight file found, skipping weight reprojection.")

        print(f"  Done: resampled image located at {sci_out}")

print("\nAll Euclid images reprojected onto the VISTA grid.")

