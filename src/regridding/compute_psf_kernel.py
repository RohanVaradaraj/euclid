"""
compute_psf_kernel.py

Find the kernel to go from the Euclid/JWST psfs to the VISTA psf.

Created: Friday 14th June 2024.
"""

from pathlib import Path
from astropy.io import fits
import numpy as np
import os
import matplotlib.pyplot as plt

plt.rcParams['axes.linewidth'] = 2.5
plt.rcParams.update({'font.size': 15})
plt.rcParams['figure.dpi'] = 100

#! Switches and setup
instrument = 'jwst'

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

                # plt.imshow(yz_slice, origin='lower')
                # plt.show()
                # continue

                # Save the slice as a fits file with the same header
                hdu = fits.PrimaryHDU(yz_slice, header=header)
                hdu.writeto(euclid_psf_dir / f'{filter_name}_psf.fits', overwrite=True)

            euclid_psf = euclid_psf_dir / f'{filter_name}_psf.fits'

            output_name = output_dir / f'{filter_name}_to_VISTA_kernel.fits'
            os.system(f'addpixscl {str(euclid_psf)} 0.1')
            os.system(f'addpixscl {str(vista_psf)} 0.15')
            os.system(f'pypher {str(euclid_psf)} {str(vista_psf)} {str(output_name)}')

    if instrument.lower() == 'jwst':

        filter_names = filter_names_jwst
        jwst_psf_dir = Path.cwd().parents[1] / 'data' / 'psf' / f'{field}' / 'results'

        output_dir = Path.cwd().parents[1] / 'data' / 'psf' / f'{field}' / 'kernel'

        # Get reference VISTA PSF
        vista_psf_dir = Path.cwd().parents[1] / 'data' / 'psf' / f'{field}' / 'ref_psf'
        vista_psf = vista_psf_dir / 'Y_DR6_psf.fits'

        for filter_name in filter_names:

            if not Path.exists(jwst_psf_dir / f'{filter_name.lower()}_psf.fits'):

                hdu = fits.open(jwst_psf_dir / f'{filter_name.lower()}.psf')
                data = hdu[1].data[0][0]
                header = hdu[1].header

                # PSF slice
                yz_slice = np.array(data[0, :, :])

                # plt.imshow(yz_slice, origin='lower', vmin=0, vmax=0.0001)
                # plt.show()
                # continue

                # Save the slice as a fits file with the same header
                hdu = fits.PrimaryHDU(yz_slice, header=header)
                hdu.writeto(jwst_psf_dir / f'{filter_name.lower()}_psf.fits', overwrite=True)

            jwst_psf = jwst_psf_dir / f'{filter_name.lower()}_psf.fits'

            output_name = output_dir / f'{filter_name.lower()}_to_VISTA_kernel.fits'
            os.system(f'addpixscl {str(jwst_psf)} 0.15')
            os.system(f'addpixscl {str(vista_psf)} 0.15')
            os.system(f'pypher {str(jwst_psf)} {str(vista_psf)} {str(output_name)}')




