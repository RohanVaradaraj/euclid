"""
This code copies files from directory A, which also exist in directory B, into directory C.

Useful for when I have run a selection step, e.g. with Euclid, and want to copy the same files during a separate VISTA selection.

Created: Wednesday 18th December 2024.
"""

import os
from pathlib import Path

overwrite = True

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
