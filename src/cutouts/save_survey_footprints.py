"""
save_survey_footprints.py

Get the footprints of Euclid, CWEB and PRIMER tiles to use in a function that checks if a point lies in the footprint of a tile.

Created: Thursday 21st March 2024.
"""

from astropy.io import fits
from pathlib import Path
import glob
from astropy.wcs import WCS
import numpy as np

#! Save Euclid footprint
euclid_dir = Path.home() / 'euclid' / 'Y' / 'COSMOS'
images = glob.glob(str(euclid_dir / '*BGSUB-*'))

euclid_footprint = []
euclid_labels = []

for i, image in enumerate(images):

    label = image.split('_')[-1].split('.')[0]

    with fits.open(image) as hdul:

        wcs = WCS(hdul[0].header)
        footprint = wcs.calc_footprint()
        euclid_footprint.append(footprint)
        euclid_labels.append(label)

euclid_footprint = np.array(euclid_footprint)
euclid_labels = np.array(euclid_labels)
print(euclid_footprint.shape)
np.save(Path.cwd().parent.parent / 'data' / 'mosaic' / 'euclid_footprint.npy', euclid_footprint)
np.save(Path.cwd().parent.parent / 'data' / 'mosaic' / 'euclid_labels.npy', euclid_labels)

#! Save PRIMER footprint
primer = Path.home() / 'JWST' / 'primer_cosmos_nircam_v0.5_f444w_30mas_sci.fits'
with fits.open(primer) as hdu_primer:
    wcs_primer = WCS(hdu_primer[0].header)
    primer_footprint = wcs_primer.calc_footprint()

np.save(Path.cwd().parent.parent / 'data' / 'mosaic' / 'primer_footprint.npy', primer_footprint)


#! Save the COSMOS-Web footprint
jwst_dir = Path.home().parent.parent / 'extraspace' / 'varadaraj' / 'CWEB'
jwst_files = glob.glob(str(jwst_dir / '*F444W*'))

jwst_footprints = []
jwst_labels = []

for jwst_file in jwst_files:

    label = jwst_file.split('-')[-1].split('_')[0] 
    
    with fits.open(jwst_file) as hdu_jwst:

        wcs_jwst = WCS(hdu_jwst[1].header)
        jwst_footprint = wcs_jwst.calc_footprint()
        jwst_footprints.append(jwst_footprint)
        jwst_labels.append(label)

jwst_footprints = np.array(jwst_footprints)
jwst_labels = np.array(jwst_labels)
np.save(Path.cwd().parent.parent / 'data' / 'mosaic' / 'cweb_footprint.npy', jwst_footprints)
np.save(Path.cwd().parent.parent / 'data' / 'mosaic' / 'cweb_labels.npy', jwst_labels)


