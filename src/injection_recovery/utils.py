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
from astropy.coordinates import SkyCoord
import struct
from scipy.io import FortranFile

data_dir = Path.cwd().parents[3] / 'data'

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



def cutout_subimage(field, image, image_size, pix_scale, n_images, random=True, x=0, y=0, overwrite=False):
    """
    Extract a subimage from a larger image.

    :param field: The VISTA tile name (e.g., 'COSMOS', 'XMM1', 'CDFS2').
    :param image: The image to extract the subimage from.
    :param image_size: The size of the subimage, in arcmin.
    :param pix_scale: The pixel scale of the image, in arcsec/pix.
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

    image_dir = data_dir / field / image
    weight_dir = data_dir / field / (image.split('.fits')[0] + '_wht.fits')

    with fits.open(image_dir) as hdu:
        data = hdu[0].data
        header = hdu[0].header
        wcs = WCS(header)

    with fits.open(weight_dir) as hdu_w:
        weight = hdu_w[0].data
        weight_header = hdu_w[0].header
        wcs_weight = WCS(weight_header)

    pix_size = pix_scale

    image_size = image_size * 60 / pix_size # convert to pixels

    # Get some key information from the original header that doesn't make it into the cutout header
    equinox = header['EQUINOX']
    exptime = header['EXPTIME']
    gain = header['GAIN']
    saturate = header['SATURATE']

    for i in range(n_images):
        print(f"Generating image {i+1}/{n_images}")
        if random:
            #! Note: these x,y are in terms of the original image pixel basis
            x = np.random.randint(image_size/2 +100, data.shape[1]-(image_size/2+100))
            y = np.random.randint(image_size/2 +200, data.shape[0]-(image_size/2+100))   # avoid edges

        # Convert x,y to RA,DEC as a SkyCoord object
        coord = SkyCoord.from_pixel(x, y, wcs)
        ra, dec = coord.to_string('hmsdms').split(' ')

        cutout = Cutout2D(data, coord, (image_size, image_size), wcs=wcs)
        weight_cutout = Cutout2D(weight, coord, (image_size, image_size), wcs=wcs_weight)

        # Save the cutout to a new FITS file
        image_name = image.split('.fits')[0] + f'_cutout_{int(x)}_{int(y)}_{int(image_size)}_pix_{int(image_size * pix_size / 60)}_arcmin.fits'
        weight_name = image.split('.fits')[0] + f'_cutout_{int(x)}_{int(y)}_{int(image_size)}_pix_{int(image_size * pix_size / 60)}_arcmin_wht.fits'

        cutout_header = cutout.wcs.to_header()

        # Add some key information to the header
        #cutout_header['EQUINOX'] = equinox
        cutout_header['EXPTIME'] = exptime
        cutout_header['GAIN'] = gain
        cutout_header['SATURATE'] = saturate

        # Remove some unnecessary header information
        cutout_header.remove('LONPOLE')
        cutout_header.remove('LATPOLE')
        cutout_header.remove('MJDREF')

        # rename PC1_1 to CD1_1 etc
        cutout_header['CD1_1'] = cutout_header['PC1_1']
        cutout_header['CD1_2'] = 0.0
        cutout_header['CD2_1'] = 0.0
        cutout_header['CD2_2'] = cutout_header['PC2_2']

        # And remove the original PCi_j
        cutout_header.remove('PC1_1')
        cutout_header.remove('PC2_2')

        # Save
        hdu = fits.PrimaryHDU(cutout.data, header=cutout_header)
        hdu.writeto(cutout_path / image_name, overwrite=True)

        # And weight
        hdu = fits.PrimaryHDU(weight_cutout.data, header=cutout_header)
        hdu.writeto(weight_path / weight_name, overwrite=True)


    return None



def read_binary_SED_file(file_name, file_dir=Path.home() / 'lephare' / 'lephare_work' / 'lib_mag', verbose=False):
    """
    Read in the .bin file produced by LePhare after running the sedtolib step.

    :param file_dir: The directory where the file is located.
    :param file_name: The name of the file to read.
    """

    #? First read in the .doc file
    doc_info = {}
    
    doc_file = file_dir / (file_name.split('.bin')[0] + '.doc')
    with open(doc_file, "r") as f:
        for line in f:
            # Skip empty lines
            if not line.strip():
                continue
            parts = line.strip().split()
            key = parts[0]
            value = parts[1:]  # rest of the line
            doc_info[key] = value
    
    if verbose:
        print(".doc info:")
        for key, value in doc_info.items():
            print(f"{key}: {' '.join(value)}")

    n_seds = int(doc_info["NUMBER_SED"][0])
    n_rows = int(doc_info["NUMBER_ROWS"][0])
    record_length = int(doc_info["RECORD_LENGTH"][0])  # number of float32 per SED

    n_rows = 153_850
    n_seds = 48

    data = np.fromfile(file_dir / file_name, dtype=np.uint8)

    # assert data.size == n_rows * n_seds, "Unexpected file size!"

    # Reshape: each row is 48 SED fluxes at one wavelength
    # data_matrix = data.reshape((n_rows, 384))

    # To get SED_i as a function of λ:
    # each column is one SED
    # seds = data_matrix.T  # Shape: (48, 153850)

    print("Loaded SEDs:", seds.shape)  # (n_seds, n_wavelengths)

    # plt.plot(seds[24])
    # plt.show()


    # data = np.fromfile(file_dir / file_name, dtype=np.float32)
    # data = data.reshape((n_rows, int(record_length/4)))


    # f = open(file_dir / file_name, "rb")
    # print(f.read()[0:200])

    # 384 bytes / 8 = 48 float64s per row
#     dtype = np.dtype([('values', 'f8', 48)])  # structured row: 48 float64s

#     # Read entire file as structured array
#     data = np.fromfile(file_dir / file_name, dtype=dtype)
#     print(data.shape)
#    #data = data.T

#     plt.plot(data[0:1000])
#     plt.show()

#     print(data)
    # for row in data:
    #     print(row)
    #     plt.plot(row)
    #     plt.show()





if __name__ == "__main__":

    #cutout_subimage('UVISTA_YJ_DR6.fits', 5000, 5000, 10, random=True)
    read_binary_SED_file('GAL_COS.bin', verbose=True)