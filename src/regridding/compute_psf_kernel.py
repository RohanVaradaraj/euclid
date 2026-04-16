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
instrument = 'euclid'

fields = ['COSMOS']

filter_names_euclid = ['YE_DR1', 'JE_DR1', 'HE_DR1', 'VIS_DR1']
#filter_names_jwst = ['F115W', 'F150W', 'F277W', 'F444W']

dr = 'DR1' # DR1

def findPlotLimits(data: np.ndarray) -> tuple:

    mean = np.mean(data)
    std_dev = np.std(data)

    # Sigma clip the data
    data = data[(data < mean + 3 * std_dev)]

    # Recalculate mean and std_dev
    mean = np.mean(data)
    std_dev = np.std(data)

    lower = mean - 2 * std_dev
    upper = mean + 5 * std_dev

    return lower, upper


for field in fields:

    if instrument.lower() == 'euclid':

        filter_names = filter_names_euclid
        # euclid_psf_dir = Path.cwd().parents[1] / 'data' / 'psf' / f'{field}' / 'results'
        euclid_psf_dir = Path.cwd().parents[3] / 'data' / 'psf' / f'{field}' / 'results' #! read in Euclid PSFs from here.
        # Make output directory if it doesn't exist
        euclid_psf_dir.mkdir(parents=True, exist_ok=True)

        # psf_out_dir = Path.cwd().parents[1] / 'data' / 'psf' / f'{field}' / 'results'
        # psf_out_dir.mkdir(parents=True, exist_ok=True)

        output_dir = Path.cwd().parents[1] / 'data' / 'psf' / f'{field}' / 'kernel' #! Output the kernels to here
        output_dir.mkdir(parents=True, exist_ok=True)

        # Get reference VISTA PSF
        vista_psf_dir = Path.cwd().parents[1] / 'data' / 'psf' / f'{field}' / 'ref_psf' #! Read in VISTA PSF from here, and convert .psf to .fits
        # vista_psf_dir.mkdir(parents=True, exist_ok=True)
        vista_psf = vista_psf_dir / 'HSC-Z_DR3_psf.fits'

        y_psf = vista_psf_dir / 'HSC-Z_DR3.psf'
        y_psf_fits = vista_psf_dir / 'HSC-Z_DR3_psf.fits'

        if not y_psf_fits.exists():
            print('Extracting HSC PSF from ', y_psf)

            hdu = fits.open(y_psf)
            data = hdu[1].data[0][0]
            header = hdu[1].header

            yz_slice = np.array(data[0, :, :])

            hdu = fits.PrimaryHDU(yz_slice, header=header)
            hdu.writeto(y_psf_fits, overwrite=True)
            print('Wrote to', y_psf_fits)

        for filter_name in filter_names:

            print(f'Processing filter: {filter_name}')

            #! If the Euclid PSF slice doesn't exist, extract it from the .psf file and save as a fits file with the same header. If it does exist, just read in the fits file.
            if not Path.exists(euclid_psf_dir / f'{filter_name}_psf.fits'):

                print('Reading ', euclid_psf_dir / f'{filter_name}.psf', ' and converting to .fits')
                hdu = fits.open(euclid_psf_dir / f'{filter_name}.psf')
                data = hdu[1].data[0][0]
                header = hdu[1].header

                # PSF slice
                yz_slice = np.array(data[0, :, :])

                vmin, vmax = findPlotLimits(yz_slice)

                plt.imshow(yz_slice, origin='lower', vmin=vmin, vmax=vmax)
                plt.show()

                #  Save the slice as a fits file with the same header
                hdu = fits.PrimaryHDU(yz_slice, header=header)
                hdu.writeto(euclid_psf_dir / f'{filter_name}_psf.fits', overwrite=True)

            euclid_psf = euclid_psf_dir / f'{filter_name}_psf.fits'

            # output_name = output_dir / f'{filter_name}_to_VISTA_kernel_{dr}.fits'
            #output_name = output_dir / f'{filter_name}_to_VISTA_kernel.fits'
            output_name = output_dir / f'{filter_name}_to_HSC_kernel.fits'
            # If field is cosmos, pixel scale is 0.15
            if field == 'COSMOS':
                os.system(f'addpixscl {str(euclid_psf)} 0.15')
                os.system(f'addpixscl {str(vista_psf)} 0.15')
            else:
                os.system(f'addpixscl {str(euclid_psf)} 0.2')
                os.system(f'addpixscl {str(vista_psf)} 0.2')
            os.system(f'pypher {str(euclid_psf)} {str(vista_psf)} {str(output_name)}')

    if instrument.lower() == 'jwst':

        filter_names = filter_names_jwst
        jwst_psf_dir = Path.cwd().parents[1] / 'data' / 'psf' / f'{field}' / 'results'

        output_dir = Path.cwd().parents[1] / 'data' / 'psf' / f'{field}' / 'kernel'

        # Get reference VISTA PSF
        vista_psf_dir = Path.cwd().parents[1] / 'data' / 'psf' / f'{field}' / 'ref_psf'
        vista_psf = vista_psf_dir / 'Y_DR6_psf.fits'

        for filter_name in filter_names:

            #if not Path.exists(jwst_psf_dir / f'{filter_name.lower()}_psf.fits'):
            if not Path.exists(jwst_psf_dir / f'webbpsf_{filter_name.upper()}.fits'):

                print('Manually extracting the PSF slice.')
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

            output_name = output_dir / f'webbpsf_{filter_name.lower()}_to_VISTA_kernel.fits'
            os.system(f'addpixscl {str(jwst_psf)} 0.03')
            os.system(f'addpixscl {str(vista_psf)} 0.15')
            os.system(f'pypher {str(jwst_psf)} {str(vista_psf)} {str(output_name)}')




