"""
For the selections involving Euclid, mask the sample to the Euclid footprint.

Created: Wednesday 20th November 2024.
"""

from astropy.io import ascii
from astropy.table import Table
import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path
import glob
from sed_fitting_codes import *
import sys
from selection import generate_selection_name
from astropy.coordinates import SkyCoord
import json
from regions import Regions
import shutil
from astropy.wcs import WCS

cutout_path = Path.cwd().parents[0] / 'cutouts'
sys.path.append(str(cutout_path))
from cutout_codes import *
from vista_cutouts import *

sed_path = Path.cwd().parents[0] / 'sed_fitting'
sys.path.append(str(sed_path))
from sed_fitting_codes import parse_spec_file

plt.rcParams['axes.linewidth'] = 2.5
plt.rcParams.update({'font.size': 15})
plt.rcParams['figure.dpi'] = 100

if len(sys.argv) > 1:
    filters_json = sys.argv[1]
    filters = json.loads(filters_json)
    bools_json = sys.argv[2]
    bools = json.loads(bools_json)
    all_filters_json = sys.argv[3]
    all_filters = json.loads(all_filters_json)
    run_type_json = sys.argv[4]
    run_type = json.loads(run_type_json)

#! Set up the output PDF name from the det and nondet filters
# Set up the detection filter base
det_list = [f for f, t in filters.items() if t['type'] == 'detection']
if len(det_list) != 0:
    base_det = 'det_' + '_'.join(det_list)
else:
    stack_list = [f for f, t in filters.items() if t['type'] == 'stacked-detection']
    stack_filters = stack_list[0].split('+')
    base_det = 'det_' + '_'.join(stack_filters)
    det_list = stack_filters

# Set up the non-detection filter base
base_nondet = 'nonDet_' + '_'.join([f for f, t in filters.items() if t['type'] == 'non-detection'])

# Collect the components for the filename
filename_components = [base_det, base_nondet]

# Add run_type if available
if run_type != '':
    filename_components.append(run_type)

# Construct the final output PDF name
output_file = '_'.join(filename_components) + '.fits'

print('Saving to: ', output_file)

#! Set up directories
if run_type != '':
    zphot_folder = base_det + '_' + run_type
else:
    zphot_folder = base_det
zphot_dir = Path.cwd().parents[1] / 'data' / 'sed_fitting' / 'zphot' / zphot_folder
print('Taking SEDs from: ', zphot_dir)

# Get filters
filter_dict = filter_widths()

# Load the parent catalogue to get RA,DEC
# Generate the name of the parent catalogue
parent_cat_dir = Path.cwd().parents[1] / 'data' / 'catalogues'
parent_cat_name = generate_selection_name('COSMOS', filters)

# Load the parent catalogue
t = Table.read(parent_cat_dir / parent_cat_name)

# Load the Euclid mask (DS9 region file)
mask_dir = Path.cwd().parents[3] / 'data' / 'masks'
euclid_mask_file = mask_dir / 'Euclid_square_COSMOS.reg'

# Open a random image to get a WCS
image_dir = Path.cwd().parents[3] / 'data' / 'COSMOS'
image = image_dir / 'UVISTA_Y_DR6.fits'
with fits.open(str(image)) as hdul:
    wcs = WCS(hdul[0].header)

# Read the region
regions = Regions.read(str(euclid_mask_file))

# Assuming it's a single rectangular region
euclid_region = regions[0]

# Extract RA, DEC from the catalogue
coords = SkyCoord(t['RA'], t['DEC'], unit='deg')

# Use the `contains` method of the region, passing `wcs=None`
mask = [euclid_region.contains(coord, wcs) for coord in coords]

# Unpack the mask
mask = np.array([m[0] for m in mask])   

# Apply the mask to the catalogue
t_filtered = t[mask]

# Save or print the filtered catalogue
print(f"Number of objects within the Euclid mask: {len(t_filtered)}")
print(f"Number of objects outside the Euclid mask: {len(t) - len(t_filtered)}")
t_filtered.write(parent_cat_dir / output_file, overwrite=True)

# Plot the filtered and unfiltered objects in RA,DEC
# plt.figure(figsize=(8, 6))
# plt.scatter(t['RA'], t['DEC'], s=3, label='All UVISTA objects')
# plt.scatter(t_filtered['RA'], t_filtered['DEC'], s=5, label='Objects within Euclid pointing')
# plt.xlabel('RA [deg]')
# plt.ylabel('DEC [deg]')
# plt.gca().invert_xaxis()
# plt.legend()
# plt.show()

#! Loop through the sed directory and move those outside the mask to a new directory
# New directory to place things outside footprint
outside_euclid_dir = zphot_dir.parents[0] / (zphot_folder + '_outside_footprint')

# Corresponding directories for the other selections
dusty_folder = base_det + '_' + run_type + '_dusty'
dusty_dir = Path.cwd().parents[1] / 'data' / 'sed_fitting' / 'zphot' / dusty_folder
bd_folder = base_det + '_' + run_type + '_bd'
bd_dir = Path.cwd().parents[1] / 'data' / 'sed_fitting' / 'zphot' / bd_folder
lya_folder = base_det + '_' + run_type + '_lya'
lya_dir = Path.cwd().parents[1] / 'data' / 'sed_fitting' / 'zphot' / lya_folder

# Corresponding not-in-footprint directories for the above
outside_dusty_dir = dusty_dir.parents[0] / (dusty_folder + '_outside_footprint')
outside_bd_dir = bd_dir.parents[0] / (bd_folder + '_outside_footprint')
outside_lya_dir = lya_dir.parents[0] / (lya_folder + '_outside_footprint')

# If it doesn't exist, make it
if not outside_euclid_dir.exists():
    outside_euclid_dir.mkdir(parents=True)
if not outside_dusty_dir.exists():
    outside_dusty_dir.mkdir(parents=True)
if not outside_bd_dir.exists():
    outside_bd_dir.mkdir(parents=True)
if not outside_lya_dir.exists():
    outside_lya_dir.mkdir(parents=True)

# Get the .spec files
spec_files = glob.glob(str(zphot_dir / '*.spec'))

# Loop through the files
for spec_file in spec_files:
    ID = int(spec_file.split('/')[-1].split('Id')[-1].lstrip('0').split('.spec')[0])
    if ID not in t_filtered['ID']:
        # Move the file to the new directory
        shutil.move(spec_file, outside_euclid_dir)

        # Also move the corresponding file in the dusty directory
        if (dusty_dir / spec_file.split('/')[-1]).exists():
            shutil.move(dusty_dir / spec_file.split('/')[-1], outside_dusty_dir)
        
        # Also move the corresponding file in the bd directory
        if (bd_dir / spec_file.split('/')[-1]).exists():
            shutil.move(bd_dir / spec_file.split('/')[-1], outside_bd_dir)
        
        # Also move the corresponding file in the lya directory
        if (lya_dir / spec_file.split('/')[-1]).exists():
            shutil.move(lya_dir / spec_file.split('/')[-1], outside_lya_dir)


