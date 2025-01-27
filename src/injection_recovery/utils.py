"""
Functions shared by all components.

Created: Wednesday 4th December 2024.
"""

import yaml
from astropy.io import fits
from pathlib import Path
from astropy.wcs import WCS
from astropy.nddata import Cutout2D
import numpy as np
import matplotlib.pyplot as plt
from astropy.visualization.wcsaxes import WCSAxes
from matplotlib.patches import Rectangle

data_dir = Path.cwd().parents[3] / 'data' / 'COSMOS'

def load_config(config_file):
    """
    Load configuration from a YAML file.
    """
    with open(config_file, 'r') as f:
        return yaml.safe_load(f)



def filter_files():
    """ 
    Returns the location of filter transmission curve files for input into the lephare config file.
    """
    filt_files = {
        'CFHT-u':'cfht/megacam/up.pb',
        'CFHT-g':'cfht/megacam/gp.pb',
        'CFHT-r':'cfht/megacam/rp.pb',
        'CFHT-z':'cfht/megacam/zp.pb',
        'HSC-G_DR3':'myfilters/HSC/g_HSC.txt',
        'HSC-R_DR3':'myfilters/HSC/r_HSC.txt',
        'HSC-I_DR3':'myfilters/HSC/i_HSC.txt',
        'HSC-NB0816_DR3':'myfilters/HSC/nb816_HSC.txt',
        'HSC-Z_DR3':'myfilters/HSC/z_HSC.txt',
        'HSC-NB0921_DR3':'myfilters/HSC/nb921_HSC.txt',
        'HSC-Y_DR3':'myfilters/HSC/y_HSC.txt',
        'Y':'myfilters/VISTA/VISTA_Y.txt',
        'J':'myfilters/VISTA/VISTA_J.txt',
        'H':'myfilters/VISTA/VISTA_H.txt',
        'Ks':'myfilters/VISTA/VISTA_Ks.txt',
        'f115w':'myfilters/JWST/f115w_angstroms.txt',
        'f150w':'myfilters/JWST/f150w_angstroms.txt',
        'f277w':'myfilters/JWST/f277w_angstroms.txt',
        'f444w':'myfilters/JWST/f444w_angstroms.txt',
        'VIS':'myfilters/Euclid/Euclid_VIS.txt',
        'Ye':'myfilters/Euclid/Euclid_Y.txt',
        'Je':'myfilters/Euclid/Euclid_J.txt',
        'He':'myfilters/Euclid/Euclid_H.txt',
        'ch1cds':'myfilters/SPITZER/irac_ch1.txt',
        'ch2cds':'myfilters/SPITZER/irac_ch2.txt',
    }

    return filt_files


def cutout_subimage(image, image_size, n_images, random=True, x=0, y=0, overwrite=False):
    """
    Extract a subimage from a larger image and plot the footprint of the original and cutouts in RA, Dec.
    
    :param image: The image to extract the subimage from.
    :param x: The x-coordinate of the subimage center.
    :param y: The y-coordinate of the subimage center.
    :param size: The size of the subimage, in arcmin
    """

    cutout_path = Path.cwd() / 'images' / 'cutouts'
    cutout_path.mkdir(parents=True, exist_ok=True)
    weight_path = Path.cwd() / 'images' / 'cutouts' / 'weights'
    weight_path.mkdir(parents=True, exist_ok=True)

    # Delete files in cutout_path if overwrite is True
    if overwrite:
        for file in cutout_path.glob('*.fits'):
            file.unlink()
        for file in weight_path.glob('*.fits'):
            file.unlink()

    image_dir = data_dir / image
    weight_dir = data_dir / (image.split('.fits')[0] + '_wht.fits')

    with fits.open(image_dir) as hdu:
        data = hdu[0].data
        header = hdu[0].header
        wcs = WCS(header)

    with fits.open(weight_dir) as hdu_w:
        weight = hdu_w[0].data
        weight_header = hdu_w[0].header
        wcs_weight = WCS(weight_header)

    pix_size = 0.15  # arcsec / pix
    image_size = image_size * 60 / pix_size  # convert to pixels
    print(f"Image size in pixels: {image_size}")

    # Get reference pixel from image
    crpix_x = header['CRPIX1']
    crpix_y = header['CRPIX2']

    crval_ra, crval_dec = wcs.all_pix2world(crpix_x, crpix_y, 0)
    print(f"CRVAL from wcs: {crval_ra}, {crval_dec}")

    crval1 = header['CRVAL1']
    crval2 = header['CRVAL2']
    print(f"CRVAL from header: {crval1}, {crval2}")

    # Prepare for plotting
    # fig = plt.figure(figsize=(10, 10))
    # ax = fig.add_subplot(1, 1, 1, projection=wcs)
    # ax.imshow(data, origin='lower', cmap='gray', alpha=0.5)
    # ax.set_xlabel('RA')
    # ax.set_ylabel('Dec')

    # # Plot original image footprint
    # ax.set_title('Original Image and Cutout Footprints')
    # ax.grid(color='white', ls='dotted')

    for i in range(n_images):
        print(f"Generating image {i + 1}/{n_images}")
        if random:
            x = np.random.randint(image_size / 2 + 100, data.shape[1] - (image_size / 2 + 100))
            y = np.random.randint(image_size / 2 + 200, data.shape[0] - (image_size / 2 + 100))  # avoid edges

        # Convert x, y to RA, DEC
        ra, dec = wcs.all_pix2world(x, y, 0)
        print(f"RA, Dec: {ra}, {dec}")

        cutout = Cutout2D(data, (x, y), (image_size, image_size), wcs=wcs)
        weight_cutout = Cutout2D(weight, (x, y), (image_size, image_size), wcs=wcs_weight)

        # Save the cutout to a new FITS file
        image_name = image.split('.fits')[0] + f'_cutout_{int(x)}_{int(y)}_{int(image_size)}_pix_{int(image_size * pix_size / 60)}_arcmin.fits'
        weight_name = image.split('.fits')[0] + f'_cutout_{int(x)}_{int(y)}_{int(image_size)}_pix_{int(image_size * pix_size / 60)}_arcmin_wht.fits'

        # Reset reference pixel
        cutout.wcs.wcs.crpix = [image_size / 2, image_size / 2]

        # Get CRVALs from original image at this position
        cutout.wcs.wcs.crval = [ra, dec]

        print(cutout.wcs.to_header())

        # Save
        hdu = fits.PrimaryHDU(cutout.data, header=cutout.wcs.to_header())
        hdu.writeto(cutout_path / image_name, overwrite=True)

        # And weight
        hdu = fits.PrimaryHDU(weight_cutout.data, header=weight_cutout.wcs.to_header())
        hdu.writeto(weight_path / weight_name, overwrite=True)

        # Plot the cutout footprint as a rectangle
        ra_min, dec_min = wcs.all_pix2world(x - image_size / 2, y - image_size / 2, 0)
        ra_max, dec_max = wcs.all_pix2world(x + image_size / 2, y + image_size / 2, 0)

    #     ax.add_patch(Rectangle(
    #         (ra_min, dec_min),
    #         ra_max - ra_min,
    #         dec_max - dec_min,
    #         edgecolor='red',
    #         facecolor='none',
    #         lw=2,
    #         label=f'Cutout {i + 1}'
    #     ))

    # plt.legend()
    # plt.show()

    return None

# def cutout_subimage(image, image_size, n_images, random=True, x=0, y=0, overwrite=False):
#     """
#     Extract a subimage from a larger image.
    
#     :param image: The image to extract the subimage from.
#     :param x: The x-coordinate of the subimage center.
#     :param y: The y-coordinate of the subimage center.
#     :param size: The size of the subimage, in arcmin
#     """

#     cutout_path = Path.cwd() / 'images' / 'cutouts'
#     cutout_path.mkdir(parents=True, exist_ok=True)
#     weight_path = Path.cwd() / 'images' / 'cutouts' / 'weights'
#     weight_path.mkdir(parents=True, exist_ok=True)

#     # Delete files in cutout_path if overwrite is True
#     if overwrite:
#         for file in cutout_path.glob('*.fits'):
#             file.unlink()
#         for file in weight_path.glob('*.fits'):
#             file.unlink()

#     image_dir = data_dir / image
#     weight_dir = data_dir / (image.split('.fits')[0] + '_wht.fits')

#     with fits.open(image_dir) as hdu:
#         data = hdu[0].data
#         header = hdu[0].header
#         wcs = WCS(header)

#     with fits.open(weight_dir) as hdu_w:
#         weight = hdu_w[0].data
#         weight_header = hdu_w[0].header
#         wcs_weight = WCS(weight_header)

#     pix_size = 0.15 # arcsec / pix

#     image_size = image_size * 60 / pix_size # convert to pixels
#     print(f"Image size in pixels: {image_size}")

#     # Get reference pixel from image
#     print('Original image size:', data.shape)
#     crpix_x = header['CRPIX1']
#     crpix_y = header['CRPIX2']

#     print(f"CRPIX: {crpix_x}, {crpix_y}")

#     for i in range(n_images):
#         print(f"Generating image {i+1}/{n_images}")
#         if random:
#             x = np.random.randint(image_size/2 +100, data.shape[1]-(image_size/2+100))
#             y = np.random.randint(image_size/2 +200, data.shape[0]-(image_size/2+100))   # avoid edges

#         # Convert x,y to RA,DEC
#         ra, dec = wcs.all_pix2world(x, y, 0)
#         print(f"RA, Dec: {ra}, {dec}")
    
#         cutout = Cutout2D(data, (x, y), (image_size, image_size), wcs=wcs)
#         weight_cutout = Cutout2D(weight, (x, y), (image_size, image_size), wcs=wcs_weight)

#         # Save the cutout to a new FITS file
#         image_name = image.split('.fits')[0] + f'_cutout_{int(x)}_{int(y)}_{int(image_size)}_pix_{int(image_size * pix_size / 60)}_arcmin.fits'
#         weight_name = image.split('.fits')[0] + f'_cutout_{int(x)}_{int(y)}_{int(image_size)}_pix_{int(image_size * pix_size / 60)}_arcmin_wht.fits'

#         # Reset reference pixel
#         cutout.wcs.wcs.crpix = [image_size/2, image_size/2]

#         # Get CRVALs from original image at this position
#         cutout.wcs.wcs.crval = [ra, dec]

#         print(cutout.wcs.to_header())


#         # Save
#         hdu = fits.PrimaryHDU(cutout.data, header=cutout.wcs.to_header())
#         hdu.writeto(cutout_path / image_name, overwrite=True)

#         # And weight
#         hdu = fits.PrimaryHDU(weight_cutout.data, header=weight_cutout.wcs.to_header())
#         hdu.writeto(weight_path / weight_name, overwrite=True)


#     return None



if __name__ == "__main__":

    cutout_subimage('UVISTA_YJ_DR6.fits', 5000, 5000, 10, random=True)