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
    Extract a subimage from a larger image.
    
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

    with fits.open(image_dir) as hdul:
        data = hdul[0].data
        header = hdul[0].header

    with fits.open(weight_dir) as hdul:
        weight = hdul[0].data
        weight_header = hdul[0].header

    wcs = WCS(header)
    wcs_weight = WCS(weight_header)

    pix_size = 0.15 # arcsec / pix

    image_size = image_size * 60 / pix_size # convert to pixels

    for i in range(n_images):
        print(f"Generating image {i+1}/{n_images}")
        if random:
            x = np.random.randint(image_size/2 +100, data.shape[1]-(image_size/2+100))
            y = np.random.randint(image_size/2 +200, data.shape[0]-(image_size/2+100))   # avoid edges
    
        cutout = Cutout2D(data, (x, y), (image_size, image_size), wcs=wcs)
        weight_cutout = Cutout2D(weight, (x, y), (image_size, image_size), wcs=wcs_weight)

        # Save the cutout to a new FITS file
        image_name = image.split('.fits')[0] + f'_cutout_{int(x)}_{int(y)}_{int(image_size)}_pix_{int(image_size * pix_size / 60)}_arcmin.fits'
        weight_name = image.split('.fits')[0] + f'_cutout_{int(x)}_{int(y)}_{int(image_size)}_pix_{int(image_size * pix_size / 60)}_arcmin_wht.fits'
        
        # Save
        hdu = fits.PrimaryHDU(cutout.data, header=cutout.wcs.to_header())
        hdu.writeto(cutout_path / image_name, overwrite=True)

        # And weight
        hdu = fits.PrimaryHDU(weight_cutout.data, header=weight_cutout.wcs.to_header())
        hdu.writeto(weight_path / weight_name, overwrite=True)

    return None



if __name__ == "__main__":

    cutout_subimage('UVISTA_YJ_DR6.fits', 5000, 5000, 10, random=True)