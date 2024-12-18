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

overwrite = True

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

output_file = '_'.join(filename_components) + '.fits'
print('Output file: ', output_file)



#! Set up directories
if run_type != '':
    zphot_folder = base_det + '_' + run_type
else:
    zphot_folder = base_det

zphot_dir = Path.cwd().parents[1] / 'data' / 'sed_fitting' / 'zphot' / zphot_folder
print('Taking SEDs from: ', zphot_dir)

# Get corresponding dusty, bd and lya zphot dirs too
if run_type != '':
    make_zphot_dir = lambda folder: Path.cwd().parents[1] / 'data' / 'sed_fitting' / 'zphot' / (base_det + '_' + run_type + '_' + folder) 
else:
    make_zphot_dir = lambda folder: Path.cwd().parents[1] / 'data' / 'sed_fitting' / 'zphot' / (base_det + '_' + folder) 

dusty_zphot_dir = make_zphot_dir('dusty')
bd_zphot_dir = make_zphot_dir('bd')
lya_zphot_dir = make_zphot_dir('lya')

#! Load the parent catalogue to get RA,DEC
# Generate the name of the parent catalogue
parent_cat_dir = Path.cwd().parents[1] / 'data' / 'catalogues'
parent_cat_name = generate_selection_name('COSMOS', filters)
print('Loading parent catalogue: ', parent_cat_name)

# Load the parent catalogue
t = Table.read(parent_cat_dir / parent_cat_name)

#! Load the Euclid mask (DS9 region file)
mask_dir = Path.cwd().parents[3] / 'data' / 'masks'
euclid_mask_file = mask_dir / 'Euclid_square_COSMOS.reg'

#! Open a random image to get a WCS
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

# Use the `contains` method of the region
mask = [euclid_region.contains(coord, wcs) for coord in coords]

# Unpack the mask
mask = np.array([m[0] for m in mask])   

# Apply the mask to the catalogue
t_filtered = t[mask]

# Save or print the filtered catalogue
print(f"Number of objects within the Euclid mask: {len(t_filtered)}")
print(f"Number of objects outside the Euclid mask: {len(t) - len(t_filtered)}")
t_filtered.write(parent_cat_dir / output_file, overwrite=True)

#? We now have the objects in the Euclid mask.
#? Now, of these, we need those that made the VISTA SED cut.
vista_dir = Path.cwd().parents[1] / 'data' / 'sed_fitting' / 'zphot' / 'best_fits' / (base_det + '_best_highz')

#? If masking the VISTA only, take from the initial run.
if run_type == '':
    vista_dir = Path.cwd().parents[1] / 'data' / 'sed_fitting' / 'zphot' / base_det

# Plot the filtered and unfiltered objects in RA,DEC
# plt.figure(figsize=(8, 6))
# plt.scatter(t['RA'], t['DEC'], s=3, label='All UVISTA objects')
# plt.scatter(t_filtered['RA'], t_filtered['DEC'], s=5, label='Objects within Euclid pointing')
# plt.xlabel('RA [deg]')
# plt.ylabel('DEC [deg]')
# plt.gca().invert_xaxis()
# plt.legend()
# plt.show()
# exit()

#! Loop through the sed directory and move those outside the mask to a new directory
# New directory to place things outside footprint
outside_euclid_dir = zphot_dir.parents[0] / (zphot_folder + '_outside_footprint')

# Corresponding directories for the other selections, to copy into
base_dir = Path.cwd().parents[1] / 'data' / 'sed_fitting' / 'zphot' / 'best_fits'
if run_type != '':
    folder_creator = lambda folder: base_dir / f"{base_det}_{run_type}_best_{folder}"
else:
    folder_creator = lambda folder: base_dir / f"{base_det}_best_{folder}"

highz_dir = folder_creator('highz')
dusty_dir = folder_creator('dusty')
bd_dir = folder_creator('bd')
lya_dir = folder_creator('lya')

print('Copying SEDs to: ', highz_dir)
print('Copying SEDs to: ', dusty_dir)
print('Copying SEDs to: ', bd_dir)
print('Copying SEDs to: ', lya_dir)

# If overwrite, empty the above directories
if overwrite:
    for file in highz_dir.glob('*.spec'):
        file.unlink()
    for file in dusty_dir.glob('*.spec'):
        file.unlink()
    for file in bd_dir.glob('*.spec'):
        file.unlink()
    for file in lya_dir.glob('*.spec'):
        file.unlink()
    print('Deleted all previous .spec files in the best fit highz, dusty, bd and lya directories.')

# Corresponding not-in-footprint directories for the above
if run_type != '':
    outside_folder_creator = lambda folder, dir: dir.parents[0] / f"{base_det}_{run_type}_{folder}_outside_footprint"
else:
    outside_folder_creator = lambda folder, dir: dir.parents[0] / f"{base_det}_{folder}_outside_footprint"

outside_dusty_dir = outside_folder_creator('dusty', dusty_dir)
outside_bd_dir = outside_folder_creator('bd', bd_dir)
outside_lya_dir = outside_folder_creator('lya', lya_dir)

# If it doesn't exist, make it
dirs_to_create = [outside_euclid_dir, outside_dusty_dir, outside_bd_dir, outside_lya_dir]

for dir in dirs_to_create:
    dir.mkdir(parents=True, exist_ok=True)


# Get the .spec files
spec_files = glob.glob(str(zphot_dir / '*.spec'))

# If overwrite, delete all files in the outside directories
if overwrite:
    dirs_to_clear = [outside_euclid_dir, outside_dusty_dir, outside_bd_dir, outside_lya_dir]

    for dir in dirs_to_clear:
        for file in dir.glob('*.spec'):
            file.unlink()
    print('Deleted all previous .spec files in outside directories.')

# Loop through the files
for spec_file in spec_files:
    ID = int(spec_file.split('/')[-1].split('Id')[-1].lstrip('0').split('.spec')[0])

    # Check if the spec file is in vista_dir AND in the filtered catalogue
    if (ID in t_filtered['ID']) and (vista_dir / spec_file.split('/')[-1]).exists():

        # Copy from zphot_dir to the new directory
        shutil.copy(spec_file, highz_dir)

        # Copy from dusty_dir to the new directory
        shutil.copy(dusty_zphot_dir / spec_file.split('/')[-1], dusty_dir)

        # Copy from bd_dir to the new directory
        shutil.copy(bd_zphot_dir / spec_file.split('/')[-1], bd_dir)

        # Copy from lya_dir to the new directory
        shutil.copy(lya_zphot_dir / spec_file.split('/')[-1], lya_dir)


# Print number of files in each directory
print(f'Number of files in highz_dir: {len(list(highz_dir.glob("*.spec")))}')
print(f'Number of files in dusty_dir: {len(list(dusty_dir.glob("*.spec")))}')
print(f'Number of files in bd_dir: {len(list(bd_dir.glob("*.spec")))}')
print(f'Number of files in lya_dir: {len(list(lya_dir.glob("*.spec")))}')