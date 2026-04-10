"""
Get 10x10 arcmin cutouts for sources in XMM (cosmos handled by my existing cutout code).

Created: Friday 10th October 2025 (Lincare gig day).
"""

import os
import numpy as np
from astropy.io import fits
from pathlib import Path
from astropy.wcs import WCS
from astropy.nddata import Cutout2D
from astropy.coordinates import SkyCoord
import astropy.units as u

image_dir = Path.cwd().parents[3] / 'data'
save_dir = Path.cwd().parents[1] / 'data' / 'proposals' / 'p2_cutouts'

names = ['G1', 'G6', 'G12', 'G13']
fields = ['XMM1', 'XMM2', 'XMM2', 'XMM1']

RAs = [34.553059850109165, 35.81621225464279, 35.29274817965676, 34.364773943187075]
DECs = [-4.491258211612222, -5.458344980602408, -4.669179573052326, -5.141045511552535]

for name, field, RA, DEC in zip(names, fields, RAs, DECs):

    tile_dir = image_dir / field

    image_name = f'{field}_YJ.fits'

    image_path = tile_dir / image_name

    with fits.open(image_path) as hdul:
        image_data = hdul[0].data
        image_header = hdul[0].header
        wcs = WCS(image_header)
        pix_scale = np.abs(image_header['CD1_1']) * 3600

    position = SkyCoord(RA, DEC, unit=(u.deg, u.deg))

    size = 10 * 60 / pix_scale  # 10 arcmin in pixels
    size = (size, size)

    cutout = Cutout2D(image_data, position, size, wcs=wcs)

    # Save cutout
    cutout_header = cutout.wcs.to_header()
    hdu = fits.PrimaryHDU(data=cutout.data, header=cutout_header)

    file_name = f'{name}_cutout_VISTA_YJ_{field}_10x10arcmin.fits'
    save_path = save_dir / file_name
    hdu.writeto(save_path, overwrite=True)



