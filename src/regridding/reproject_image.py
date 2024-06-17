#!/usr/bin/env python3

"""
reproject the euclid images to the VISTA images

Created: Friday 14th June 2024.
"""

from pathlib import Path
import os
from astropy.io import fits

instrument ='jwst'

fields = ['COSMOS']

#filter_names = ['VIS', 'Y', 'J', 'H']
#filter_names = ['YJH']

filter_names = ['F115W', 'F150W', 'F277W', 'F444W']



data_dir = Path.cwd().parents[3] / 'data'
calib_dir = data_dir / 'bertin_config'

input_swarp = Path.cwd() / 'regrid.swarp'

for field in fields:

    reference_dir = data_dir / f'{field}'
    reference_file = reference_dir / 'UVISTA_Y_DR6_cropped.fits'

    if instrument.lower() == 'euclid':

        euclid_dir = data_dir / 'euclid' / 'images'

        for filter_name in filter_names:

            euclid_image = euclid_dir / f'{field}_{filter_name}_MOSAIC.fits'
            weight_image = euclid_dir / f'{field}_{filter_name}_MOSAIC_WHT.fits'

            euclid_image = euclid_dir / f'{field}_{filter_name}_STACK.fits'
            weight_image = euclid_dir / f'{field}_{filter_name}_STACK_WHT.fits'

            save_image = euclid_dir / f'{field}_{filter_name}_resamp.fits'
            save_weight = euclid_dir / f'{field}_{filter_name}_resamp_wht.fits'

            #os.system(f'fitsheader {str(reference_file)} > ref_header.txt')
            os.system(f'~/swarp/bin/swarp {str(euclid_image)} -WEIGHT_IMAGE {str(weight_image)} -c regrid.swarp -IMAGEOUT_NAME {str(save_image)} -WEIGHTOUT_NAME {str(save_weight)}')

    if instrument.lower() == 'jwst':

        jwst_dir = data_dir / 'CWEB'
        output_dir = data_dir.parent / 'rohan' / 'images'
        mosaics = ['0A', '0B', '1A', '1B', '2A', '2B', '3A', '3B', '4A', '4B', '5A', '5B', '6A', '6B', '7A', '7B']

        for mosaic in mosaics:
            mosaic_dir = jwst_dir / f'mosaic_{mosaic}'
            print('Running mosaic: ', mosaic)

            for filter_name in filter_names:
                print('Running filter: ', filter_name)
                
                # First unpack the image and error layers and temporarily save them.
                jwst_image = mosaic_dir / f'CWEB-{filter_name}-{mosaic}_i2dnobg_small.fits'
                with fits.open(jwst_image) as hdu:
                    image = hdu[1].data
                    img_header = hdu[1].header
                    err = hdu[2].data
                    err_header = hdu[2].header
                print('Unpacked image and error layers')
                
                # Save the image and error layers separately
                save_dir = output_dir / 'CWEB' / f'mosaic_{mosaic}'
                if Path.exists(save_dir) == False:
                    os.mkdir(save_dir)

                image_layer = save_dir / f'CWEB-{filter_name}-{mosaic}.fits'
                err_layer = save_dir / f'CWEB-{filter_name}-{mosaic}_err.fits'

                hdu = fits.PrimaryHDU(image, header=img_header)
                hdu.writeto(image_layer, overwrite=True)

                hdu = fits.PrimaryHDU(err, header=err_header)
                hdu.writeto(err_layer, overwrite=True)
                print('Saved image and error layers')

                save_image = save_dir / f'CWEB-{filter_name}-{mosaic}_resamp.fits'
                save_error = save_dir / f'CWEB-{filter_name}-{mosaic}_resamp_rms.fits'

                # Now reproject the image and error layers
                os.system(f'~/swarp/bin/swarp {str(image_layer)} -WEIGHT_IMAGE {str(err_layer)} -c regrid.swarp -IMAGEOUT_NAME {str(save_image)} -WEIGHTOUT_NAME {str(save_error)} -WEIGHT_TYPE MAP_RMS')

                # Remove the temporary image and error layers
                os.remove(image_layer)
                os.remove(err_layer)
                print('Deleted temporary image and error layers')
        