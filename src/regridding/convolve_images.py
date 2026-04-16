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

telescope = 'euclid'
field_name = 'COSMOS'
data_release = 'DR1'

if field_name == 'COSMOS':
    #euclid_dir =  Path.cwd().parents[3] / 'data' / 'euclid' / 'images' #! Q1
    # if data_release == 'Q1':
    euclid_dir = Path.home() / 'euclid' / 'COSMOS'
    save_dir = Path.cwd().parents[3] / 'data' / 'COSMOS'
    # if data_release == 'DR1':
    #     euclid_dir = Path.home().parents[1] / 'extraspace' / 'varadaraj' / 'euclid' / 'DR1' / 'COSMOS'
    #     save_dir = euclid_dir

if field_name == 'CDFS':
    euclid_dir =  Path.cwd().parents[3] / 'data' / 'euclid' / 'euclid_deep_field_fornax' / 'tmp'
    save_dir = Path.cwd().parents[3] / 'data'


filter_names = ['YE_DR1', 'JE_DR1', 'HE_DR1', 'VIS_DR1'] #? Euclid filters to convolve
# filter_names = ['f115w', 'f150w', 'f277w', 'f444w']

# cdfs_filter_dict = {'YE': 'Y', 'JE': 'J', 'HE': 'H', 'VIS': 'VIS'}

if telescope == 'euclid':
    for filter_name in filter_names:

        if field_name == 'CDFS':
            image_dir = Path.home().parents[1] / 'extraspace' / 'varadaraj' / 'euclid'
            kernel_dir = Path.cwd().parents[1] / 'data' /'psf'

            for field in ['CDFS1', 'CDFS2', 'CDFS3']:

                tile_dir = image_dir / field

                #! Science layer
                print(f'Convolving {filter_name} image with VISTA PSF kernel')
                with fits.open(tile_dir / f'{field}_{cdfs_filter_dict[filter_name]}_MOSAIC.fits') as hdu:
                    image = hdu[0].data
                    hdr = hdu[0].header

                hdu = fits.open(kernel_dir / field / 'kernel' / f'{filter_name}_to_VISTA_kernel_Q1.fits')
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
                output_name = image_dir / field / f'{field}_Euclid_{filter_name}_Q1_psfhom.fits'
                hdu.writeto(output_name)

                print(f'PSF homogenized file for {filter_name} {field} written to ', output_name)

                #! Weight layer
                print('Now doing weight layer.')

                with fits.open(tile_dir / f'{field}_{cdfs_filter_dict[filter_name]}_MOSAIC_WHT.fits') as hdu:
                    image = hdu[0].data
                    hdr = hdu[0].header

                hdu = fits.open(kernel_dir / field / 'kernel' / f'{filter_name}_to_VISTA_kernel_Q1.fits')
                kernel = hdu[0].data

                print('Opened weight image. Now convolving...')

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
                output_name = image_dir / field / f'{field}_Euclid_{filter_name}_Q1_psfhom_WHT.fits'
                hdu.writeto(output_name)

                print(f'PSF homogenized file for {filter_name} {field} written to ', output_name)


        if field_name == 'COSMOS':
            #! Science layer
            print(f'Convolving {filter_name} image with VISTA PSF kernel')

            if filter_name != 'VIS_DR1':
                filter_stem = filter_name.split('E_DR1')[0]
                sci_name = f'COSMOS_{filter_stem}_DR1_resamp.fits'
            else:
                filter_stem = 'VIS'
                sci_name = f'COSMOS_{filter_stem}_DR1_resamp.fits'

            print(euclid_dir / sci_name)

            # with fits.open(euclid_dir /f'COSMOS_{filter_name}_resamp.fits') as hdu:
            with fits.open(euclid_dir / sci_name) as hdu:
                image = hdu[0].data
                hdr = hdu[0].header

            #hdu = fits.open(Path.cwd().parents[1] / 'data' /'psf' / 'COSMOS' / 'kernel' / f'{filter_name}_to_VISTA_kernel.fits')
            hdu = fits.open(Path.cwd().parents[1] / 'data' /'psf' / 'COSMOS' / 'kernel' / f'{filter_name}_to_HSC_kernel.fits')
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
            hdu.writeto(save_dir / f'Euclid_{filter_name}_psfhom_v2.fits', overwrite=True)

            print(f'PSF homogenised image written to', save_dir / f'Euclid_{filter_name}_psfhom_v2.fits')

            #! Weight layer
            print('Now doing weight layer.')

            if filter_name != 'VIS_DR1':
                filter_stem = filter_name.split('E_DR1')[0]
                weight_name = f'COSMOS_{filter_stem}_WHT_DR1_resamp.fits'
            else:
                filter_stem = 'VIS'
                weight_name = f'COSMOS_{filter_stem}_WHT_DR1_resamp.fits'

            print(euclid_dir / weight_name)

            with fits.open(euclid_dir / weight_name) as hdu:
                image = hdu[0].data
                hdr = hdu[0].header

            #hdu = fits.open(Path.cwd().parents[1] / 'data' /'psf' / 'COSMOS' / 'kernel' / f'{filter_name}_to_VISTA_kernel.fits')
            hdu = fits.open(Path.cwd().parents[1] / 'data' /'psf' / 'COSMOS' / 'kernel' / f'{filter_name}_to_HSC_kernel.fits')
            kernel = hdu[0].data

            print('Opened weight image. Now convolving...')

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
            hdu.writeto(save_dir / f'Euclid_{filter_name}_psfhom_WHT_v2.fits', overwrite=True)

            print(f'PSF homogenised weight written to', save_dir / f'Euclid_{filter_name}_psfhom_WHT_v2.fits')
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


