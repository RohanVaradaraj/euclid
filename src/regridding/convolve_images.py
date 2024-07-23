#!/usr/bin/env python3

"""
convolve_images.py

Convolve the Euclid/JWST imaging with the VISTA PSF kernel.

Allocate 50GB of memory.

Created: Friday 14th June 2024
"""

from astropy.io import fits
from astropy.convolution import convolve
from pathlib import Path
import time
from datetime import datetime, timedelta

euclid_dir =  Path.cwd().parents[3] / 'data' / 'euclid' / 'images'
save_dir = Path.cwd().parents[3] / 'data' / 'COSMOS'

telescope = 'jwst'

filter_names = ['Y', 'J', 'H', 'VIS']
filter_names = ['f115w', 'f150w', 'f277w', 'f444w']

if telescope == 'euclid':
    for filter_name in filter_names:

        #! Science layer
        print(f'Convolving {filter_name} image with VISTA PSF kernel')

        with fits.open(euclid_dir /f'COSMOS_{filter_name}_resamp.fits') as hdu:
            image = hdu[0].data
            hdr = hdu[0].header

        hdu = fits.open(Path.cwd().parents[1] / 'data' /'psf' / 'COSMOS' / 'kernel' / f'{filter_name}_to_VISTA_kernel.fits')
        kernel = hdu[0].data

        print('Opened image. Now convolving...')

        # Start a clock
        start = time.time()
        start_time = datetime.now()
        print('Convolution started at: ', start_time.strftime('%Y-%m-%d %H:%M:%S'))

        psf_homo = convolve(image, kernel)

        # Print the time the convolution started and finished, and how long it took.
        end = time.time()
        end_time = datetime.now()
        duration = end - start
        duration_hours = duration / 3600

        print('Convolution finished at: ', end_time.strftime('%Y-%m-%d %H:%M:%S'))
        print(f'Total convolution time: {duration_hours:.2f} hours')

        hdu = fits.PrimaryHDU(psf_homo, header=hdr)
        hdu.writeto(save_dir / f'Euclid_{filter_name}_vista_matched.fits', overwrite=True)

        print(f'PSF homogenised image written to', save_dir / f'Euclid_{filter_name}_vista_matched.fits')

        #! Weight layer
        print('Now doing weight layer.')

        with fits.open(euclid_dir /f'COSMOS_{filter_name}_resamp_wht.fits') as hdu:
            image = hdu[0].data
            hdr = hdu[0].header

        hdu = fits.open(Path.cwd().parents[1] / 'data' /'psf' / 'COSMOS' / 'kernel' / f'{filter_name}_to_VISTA_kernel.fits')
        kernel = hdu[0].data

        print('Opened weight image. Now convolving...')

        psf_homo = convolve(image, kernel)

        hdu = fits.PrimaryHDU(psf_homo, header=hdr)
        hdu.writeto(save_dir / f'Euclid_{filter_name}_vista_matched_WHT.fits', overwrite=True)

        print(f'PSF homogenised weight written to', save_dir / f'Euclid_{filter_name}_vista_matched_wht.fits')
        print(f'Finished running on filter {filter_name}')

if telescope == 'jwst':
    jwst_dir = Path.cwd().parents[3] / 'data' / 'CWEB'

    for filter_name in filter_names:

        #! Science layer
        print(f'Convolving {filter_name} image with VISTA PSF kernel')

        with fits.open(jwst_dir /f'CWEB-{filter_name.upper()}_resamp.fits') as hdu:
            image = hdu[0].data
            hdr = hdu[0].header

        hdu = fits.open(Path.cwd().parents[1] / 'data' /'psf' / 'COSMOS' / 'kernel' / f'{filter_name}_to_VISTA_kernel.fits')
        kernel = hdu[0].data

        print('Opened image. Now convolving...')

        # Start a clock
        start = time.time()
        start_time = datetime.now()
        print('Convolution started at: ', start_time.strftime('%Y-%m-%d %H:%M:%S'))

        psf_homo = convolve(image, kernel)

        # Print the time the convolution started and finished, and how long it took.
        end = time.time()
        end_time = datetime.now()
        duration = end - start
        duration_hours = duration / 3600

        print('Convolution finished at: ', end_time.strftime('%Y-%m-%d %H:%M:%S'))
        print(f'Total convolution time: {duration_hours:.2f} hours')

        hdu = fits.PrimaryHDU(psf_homo, header=hdr)
        hdu.writeto(save_dir / f'CWEB_{filter_name}_vista_matched.fits', overwrite=True)

        print(f'PSF homogenised image written to', save_dir / f'CWEB_{filter_name}_vista_matched.fits')

        #! Weight layer
        print('Now doing weight layer.')

        with fits.open(jwst_dir /f'CWEB-{filter_name.upper()}_resamp_rms.fits') as hdu:
            image = hdu[0].data
            hdr = hdu[0].header

        hdu = fits.open(Path.cwd().parents[1] / 'data' /'psf' / 'COSMOS' / 'kernel' / f'{filter_name}_to_VISTA_kernel.fits')
        kernel = hdu[0].data

        print('Opened weight image. Now convolving...')

        psf_homo = convolve(image, kernel)

        hdu = fits.PrimaryHDU(psf_homo, header=hdr)
        hdu.writeto(save_dir / f'CWEB_{filter_name}_vista_matched_RMS.fits', overwrite=True)


