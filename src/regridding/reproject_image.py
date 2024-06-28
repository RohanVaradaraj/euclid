#!/usr/bin/env python3

"""
reproject the euclid images to the VISTA images

Created: Friday 14th June 2024.
"""

from pathlib import Path
import os
from astropy.io import fits
from astropy.wcs import WCS
from reproject import reproject_interp
import glob

#! Instrument can take one of three values:
#! euclid, jwst, vista

instrument ='jwst'

fields = ['COSMOS']

#filter_names = ['VIS', 'Y', 'J', 'H'] # Euclid
#filter_names = ['YJH'] # Euclid stack

filter_names = ['F115W', 'F150W', 'F277W', 'F444W'] # CWEB

#filter_names = ['Y', 'J', 'H', 'K'] # VISTA




data_dir = Path.cwd().parents[3] / 'data'
calib_dir = data_dir / 'bertin_config'

input_swarp = Path.cwd() / 'regrid.swarp'

for field in fields:

    reference_dir = data_dir / f'{field}'
    reference_file = reference_dir / 'UVISTA_Y_DR6.fits'

    if instrument.lower() == 'euclid':

        euclid_dir = data_dir / 'euclid' / 'images'

        for filter_name in filter_names:

            if filter_name != 'YJH':
                euclid_image = euclid_dir / f'{field}_{filter_name}_MOSAIC.fits' #! Original images
                weight_image = euclid_dir / f'{field}_{filter_name}_MOSAIC_WHT.fits'

            if filter_name == 'YJH':
                euclid_image = euclid_dir / f'{field}_{filter_name}_STACK.fits' #! YJH stack
                weight_image = euclid_dir / f'{field}_{filter_name}_STACK_WHT.fits'

            save_image = euclid_dir / f'{field}_{filter_name}_resamp.fits'
            save_weight = euclid_dir / f'{field}_{filter_name}_resamp_wht.fits'

            #? FIRST THE SCI IMAGE
            # Extra memory keywords
            keywords = '-COMBINE_BUFSIZE 2048 -COMBINE_TYPE WEIGHTED -VMEM_MAX 16384 -VMEM_DIR . \
                        -MEM_MAX 2048 -COMBINE_BUFSIZE 2048'

            # Now reproject the image
            os.system(f'~/swarp/bin/swarp {euclid_image} -WEIGHT_IMAGE {weight_image} -c regrid.swarp -IMAGEOUT_NAME {str(save_image)} {keywords}')

            #? NOW THE RMS IMAGE
            # Extra memory keywords
            keywords = '-WEIGHT_TYPE NONE -COMBINE_BUFSIZE 2048 -COMBINE_TYPE AVERAGE -VMEM_MAX 16384 -VMEM_DIR . \
                        -MEM_MAX 2048 -COMBINE_BUFSIZE 2048'

            # Now reproject the weight image
            os.system(f'~/swarp/bin/swarp {weight_image} -c regrid.swarp -IMAGEOUT_NAME {str(save_weight)} {keywords}')

    # If JWST, take all the tiles and resample and combine them.
    if instrument.lower() == 'jwst':

        jwst_dir = data_dir / 'CWEB'
        output_dir = data_dir.parent / 'rohan' / 'images'
        mosaics = ['0A', '0B', '1A', '1B', '2A', '2B', '3A', '3B', '4A', '4B', '5A', '5B', '6A', '6B', '7A', '7B']

        #! Go through mosaics and unpack the image and error layers
        # for mosaic in mosaics:
        #     mosaic_dir = jwst_dir / f'mosaic_{mosaic}'


        #     print('Unpacking mosaic: ', mosaic)

        #     for filter_name in filter_names:
        #         print('Running filter: ', filter_name)
                
        #         # First unpack the image and error layers and temporarily save them.
        #         jwst_image = mosaic_dir / f'CWEB-{filter_name}-{mosaic}_i2dnobg_small.fits'
        #         with fits.open(jwst_image) as hdu:
        #             image = hdu[1].data
        #             img_header = hdu[1].header
        #             err = hdu[2].data
        #             err_header = hdu[2].header
        #         print('Unpacked image and error layers')
                
        #         # Save the image and error layers separately
        #         save_dir = output_dir / 'CWEB' / f'mosaic_{mosaic}'
        #         if Path.exists(save_dir) == False:
        #             os.mkdir(save_dir)

        #         image_layer = save_dir / f'CWEB-{filter_name}-{mosaic}.fits'
        #         err_layer = save_dir / f'CWEB-{filter_name}-{mosaic}_err.fits'

        #         hdu = fits.PrimaryHDU(image, header=img_header)
        #         hdu.writeto(image_layer, overwrite=True)

        #         hdu = fits.PrimaryHDU(err, header=err_header)
        #         hdu.writeto(err_layer, overwrite=True)
        #         print('Saved image and error layers')

        # !Now glob all the image and error layers, looping through filters
        print('Reprojecting the CWEB tiles')
        for filter_name in filter_names:

            print('Running filter: ', filter_name)

            image_list = glob.glob(str(output_dir / 'CWEB' / '*' / f'CWEB-{filter_name}-*[0-7][AB].fits'))
            err_list = glob.glob(str(output_dir / 'CWEB' / '*' / f'CWEB-{filter_name}-*err.fits'))
            print('Globbed image and error layers')

            # Make string of image names to pass to Swarp
            images = image_list[0]
            for i in range(1, len(image_list)):
                images = images + ',' + image_list[i]

            # Same for rms images
            errors = err_list[0]
            for i in range(1, len(err_list)):
                errors = errors + ',' + err_list[i]

            print('Created string of image and error layers')

            save_image = jwst_dir / f'CWEB-{filter_name}_resamp.fits'
            save_error = jwst_dir / f'CWEB-{filter_name}_resamp_rms.fits'

            #? FIRST THE SCI IMAGE
            # Extra memory keywords
            keywords = '-WEIGHT_TYPE MAP_RMS -COMBINE_BUFSIZE 2048 -COMBINE_TYPE WEIGHTED -VMEM_MAX 16384 -VMEM_DIR . \
                        -MEM_MAX 2048 -COMBINE_BUFSIZE 2048'

            # Now reproject the image
            os.system(f'~/swarp/bin/swarp {images} -WEIGHT_IMAGE {errors} -c regrid.swarp -IMAGEOUT_NAME {str(save_image)} {keywords}')

            #? NOW THE RMS IMAGE
            # Extra memory keywords
            keywords = '-WEIGHT_TYPE NONE -COMBINE_BUFSIZE 2048 -COMBINE_TYPE AVERAGE -VMEM_MAX 16384 -VMEM_DIR . \
                        -MEM_MAX 2048 -COMBINE_BUFSIZE 2048'

            # Now reproject the image and error layers
            os.system(f'~/swarp/bin/swarp {errors} -c regrid.swarp -IMAGEOUT_NAME {str(save_error)} {keywords}')

            #!# Remove the temporary image and error layers
            # for image in image_list:
            #     os.remove(image)
            # for error in err_list:
            #     os.remove(error)
            # print('Removed temporary image and error layers')

    if instrument.lower() == 'vista':

        vista_dir = data_dir / f'{field}'

        for filter_name in filter_names:

            vista_image = vista_dir / 'tmp' / f'UVISTA_{filter_name}_DR6_cropped.fits'
            vista_weight = vista_dir / 'tmp' / f'UVISTA_{filter_name}_DR6_weight_cropped.fits'

            save_image = vista_dir / f'UVISTA_{filter_name}_DR6.fits'
            save_weight = vista_dir / f'UVISTA_{filter_name}_DR6_wht.fits'

            os.system(f'~/swarp/bin/swarp {str(vista_image)} -WEIGHT_IMAGE {str(vista_weight)} -c regrid.swarp -IMAGEOUT_NAME {str(save_image)} -WEIGHTOUT_NAME {str(save_weight)}')
