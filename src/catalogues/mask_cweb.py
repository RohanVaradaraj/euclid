"""
Mask the CWEB bad regions in the catalogue.

Created: Tuesday 28th October 2025.
"""

from astropy.table import Table
import numpy as np
from astropy.io import fits
from pathlib import Path
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
from astropy import units as u
from regions import Regions
#! First read in one of the CWEB images
image_dir = Path.cwd().parents[3] / 'data' / 'COSMOS'
image_name = 'CWEB_f115w_vista_matched.fits'

with fits.open(image_dir / image_name) as hdul:
    cweb_data = hdul[0].data
    cweb_header = hdul[0].header
    cweb_wcs = WCS(cweb_header)

#! Read in the catalogue
cat_dir = Path.cwd().parents[3] / 'data' / 'catalogues' / 'finalCOSMOS' / 'other'
cat_name = 'COSMOS_det_YJHK_masked_1.8as_Euclid_CWEB_2024_10_16.fits'

t = Table.read(cat_dir / cat_name)

#! Read in the CWEB regions file
regions_dir = Path.cwd().parents[3] / 'data' / 'masks' / 'COSMOS' 
regions_name = 'CWEB.reg'
region_list = Regions.read(regions_dir / regions_name)
print(region_list)

#!Convert catalogue RA/Dec to pixel coordinates
cat_coords = SkyCoord(t['RA'], t['DEC'], unit='deg')

#! Column names to mask
columns_to_mask = ['flux_f115w', 'err_f115w', 'flux_f150w', 'err_f150w',
                   'flux_f277w', 'err_f277w', 'flux_f444w', 'err_f444w']

for reg in region_list:
    contained = reg.contains(cat_coords, cweb_wcs)
    print(np.sum(contained))
    # Convert boolean array to -99 where true and original value where false
    for col in columns_to_mask:
        t[col] = np.where(contained, -99, t[col])

#! Save the masked catalogue
output_name = 'COSMOS_det_YJHK_masked_properly.fits'
t.write(output_name, overwrite=True)
