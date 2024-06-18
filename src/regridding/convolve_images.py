#!/usr/bin/env python3

"""
convolve_images.py

Convolve the Euclid/JWST imaging with the VISTA PSF kernel.

Created: Friday 14th June 2024
"""

from astropy.io import fits
from astropy.convolution import convolve
from pathlib import Path

euclid_dir =  Path.cwd().parents[3] / 'data' / 'euclid' / 'images'
save_dir = Path.cwd().parents[3] / 'data' / 'COSMOS'

filter_names = ['Y', 'J', 'H', 'VIS']

for filter_name in filter_names:

    #! Science layer
    print(f'Convolving {filter_name} image with VISTA PSF kernel')

    with fits.open(euclid_dir /f'COSMOS_{filter_name}_resamp.fits') as hdu:
        image = hdu[0].data
        hdr = hdu[0].header

    hdu = fits.open(Path.cwd().parents[1] / 'data' /'psf' / 'COSMOS' / 'kernel' / f'{filter_name}_to_VISTA_kernel.fits')
    kernel = hdu[0].data

    print('Opened image. Now convolving...')

    psf_homo = convolve(image, kernel)

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

    print(f'PSF homogenised weight written to', save_dir / f'Euclid_{filter_name}_vista_matched_WHT.fits')
    print(f'Finished running on filter {filter_name}')

