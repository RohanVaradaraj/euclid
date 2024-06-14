"""
compute_psf_kernel.py

Find the kernel to go from the Euclid/JWST psfs to the VISTA psf.

Created: Friday 14th June 2024.
"""

from pathlib import Path
from astropy.io import fits
import numpy as np
import os

#! Switches and setup
instrument = 'euclid'

fields = ['COSMOS']

filter_names_euclid = ['VIS', 'Y', 'J', 'H']
filter_names_jwst = ['F115W', 'F150W', 'F277W', 'F444W']


for field in fields:

    if instrument.lower() == 'euclid':

        filter_names = filter_names_euclid
        euclid_psf_dir = Path.cwd().parents[1] / 'data' / 'psf' / f'{field}' / 'results'

        output_dir = Path.cwd().parents[1] / 'data' / 'psf' / f'{field}' / 'kernel'

        # Get reference VISTA PSF
        vista_psf_dir = Path.cwd().parents[1] / 'data' / 'psf' / f'{field}' / 'ref_psf'
        vista_psf = vista_psf_dir / 'Y_DR6_psf.fits'

        for filter_name in filter_names:

            if not Path.exists(euclid_psf_dir / f'{filter_name}_psf.fits'):

                hdu = fits.open(euclid_psf_dir / f'{filter_name}.psf')
                data = hdu[1].data[0][0]
                header = hdu[1].header

                # PSF slice
                yz_slice = np.array(data[0, :, :])

                # Save the slice as a fits file with the same header
                hdu = fits.PrimaryHDU(yz_slice, header=header)
                hdu.writeto(euclid_psf_dir / f'{filter_name}_psf.fits', overwrite=True)

            euclid_psf = euclid_psf_dir / f'{filter_name}_psf.fits'

            output_name = output_dir / f'{filter_name}_to_VISTA_kernel.fits'
            os.system(f'addpixscl {str(euclid_psf)} 0.1')
            os.system(f'addpixscl {str(vista_psf)} 0.15')
            os.system(f'pypher {str(euclid_psf)} {str(vista_psf)} {str(output_name)}')

    if instrument.lower() == 'jwst':

        print('Do JWST stuff here.')
        



