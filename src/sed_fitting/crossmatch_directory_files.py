"""
This code copies files from directory A, which also exist in directory B, into directory C.

Useful for when I have run a selection step, e.g. with Euclid, and want to copy the same files during a separate VISTA selection.

Created: Wednesday 18th December 2024.
"""

import os
from pathlib import Path
from astropy.table import Table
import shutil
import glob
import numpy as np

dir_to_dir = False
cat_to_dir = False
dir_to_cat = True

overwrite = True

#! Original functionality of dir to dir
if dir_to_dir:

    # base sed directory
    base_dir = Path.cwd().parents[1] / 'data' / 'sed_fitting' / 'zphot'

    #! Input directories here

    # The directory to copy files from
    dir_to_copy_from = base_dir / 'det_Y_J_with_euclid_bd'

    # The directory which files will be matched with
    dir_to_match_with = base_dir / 'best_fits' / 'det_Y_J_BD'

    # The directory to copy files to
    dir_to_copy_to = base_dir / 'best_fits' / 'det_Y_J_BD_PLUS_EUCLID_PHOT'

    # If overwrite, empty
    if overwrite:
        for file in dir_to_copy_to.glob('*.spec'):
            file.unlink()
        print('Deleted all previous .spec files in the target directory.')

    # Get the files in the directories
    files_to_copy_from = os.listdir(dir_to_copy_from)
    files_to_match_with = os.listdir(dir_to_match_with)

    # Copy the files
    for file in files_to_copy_from:
        if file in files_to_match_with:
            os.system(f"cp {dir_to_copy_from}/{file} {dir_to_copy_to}/{file}")


#! Copy files from a directory corresponding to those in a catalogue, into another directory
if cat_to_dir:

    # Base SED directory
    base_dir = Path.cwd().parents[1] / 'data' / 'sed_fitting' / 'zphot'

    # Directory to copy files from
    dir_to_copy_from = base_dir / 'det_Y_J'

    # Catalogue to read IDs from
    cat_dir = Path.cwd().parents[1] / 'data' / 'catalogues' / 'candidates'
    cat_name = 'COSMOS_5sig_Y_J_nonDet_HSC_G_nonDet_HSC_R_nonDet_HSC_I_candidates_2024_11_20.fits'
    t = Table.read(cat_dir / cat_name)
    IDs = t['ID']

    # Directory to copy files to
    dir_to_copy_to = base_dir / 'best_fits' / 'det_Y_J_z7'

    for ID in IDs:

        ID = str(ID)

        # Use IDs to generate file names
        file_name = f'Id{ID.zfill(9)}.spec'

        # Copy the file safely
        if (dir_to_copy_from / file_name).exists():
            shutil.copy(dir_to_copy_from / file_name, dir_to_copy_to)
        else:
            print(f'File {file_name} does not exist in {dir_to_copy_from}.')


#! Match files in a directory to a catalogue
if dir_to_cat:

    # Base SED directory
    base_dir = Path.cwd().parents[1] / 'data' / 'sed_fitting' / 'zphot' / 'COSMOS'

    # Directory to copy files from
    dir_to_copy_from = base_dir / 'best_fits' / 'det_Y_J_with_euclid_dustyInterlopers'

    # file names
    file_names = glob.glob(str(dir_to_copy_from / '*.spec'))

    # From file_names, extract the IDs
    file_IDs = [np.int32(file_name.split('/')[-1].split('Id')[1].split('.spec')[0].lstrip('0')) for file_name in file_names]

    # Catalogue to match IDs with
    cat_dir = Path.cwd().parents[1] / 'data' / 'catalogues'
    base_cat_name = 'COSMOS_5sig_Y_J_nonDet_HSC_G_nonDet_HSC_R_nonDet_HSC_I_5percent_IRACfloor.fits'
    t = Table.read(cat_dir / base_cat_name)
    IDs = t['ID']

    # Restrict the catalogue to the IDs in file_IDs
    t = t[np.isin(IDs, file_IDs)]

    print(t)

    # Save the catalogue with some new desired name
    new_cat_name = 'COSMOS_5sig_Y_J_nonDet_HSC_G_nonDet_HSC_R_nonDet_HSC_I_dustyInterlopers_2025_08_14.fits'
    t.write(cat_dir / 'candidates' / new_cat_name, overwrite=True)

