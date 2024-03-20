"""
format_image_filenames.py

Takes the Euclid image names and formats them to something more readable.

Created: Wednesday 20th March 2024.
"""

from pathlib import Path
import glob
import re
import os

euclid_dir = Path.home() / 'euclid'
filter_names = ['VIS', 'Y', 'J', 'H']

# do bgsub images?
bgsub = False

# Do the rms files?
rms = False

# Do the PSF files?
psf = True

for filter_name in filter_names:
    filter_dir = euclid_dir / filter_name

    #! Do the bgsub images
    if bgsub:

        bgsub_images = glob.glob(str(filter_dir / '*BGSUB-*'))
        bgsub_images = [file_path.split('/')[-1] for file_path in bgsub_images]

        # Pattern to search for: 'TILE' followed by numeric string, up to the next dash
        pattern = r'(TILE\w+)'

        for i, filename in enumerate(bgsub_images):

            # Search for the pattern
            match = re.search(pattern, filename)
            matched_text = match.group(1)
        
            # Add an underscore after 'TILE'
            new_filename = filename[:match.start(1) + len("TILE")] + '_' + filename[match.start(1) + len("TILE"):match.end(1)] + '.fits'
            os.system(f'mv {filter_dir}/{filename} {filter_dir}/{new_filename}')

    #! Do the bgsub images
    if rms:

        rms_images = glob.glob(str(filter_dir / '*-RMS*'))
        rms_images = [file_path.split('/')[-1] for file_path in rms_images]

        # Pattern to search for: 'TILE' followed by numeric string, up to the next dash
        pattern = r'(TILE\w+)'

        for i, filename in enumerate(rms_images):

            # Search for the pattern
            match = re.search(pattern, filename)
            matched_text = match.group(1)
        
            # Add an underscore after 'TILE'
            new_filename = filename[:match.start(1) + len("TILE")] + '_' + filename[match.start(1) + len("TILE"):match.end(1)] + '.fits'
            os.system(f'mv {filter_dir}/{filename} {filter_dir}/{new_filename}')

    #! Do the psf images
    if psf:

        psf_images = glob.glob(str(filter_dir / '*-PSF*'))
        psf_images = [file_path.split('/')[-1] for file_path in psf_images]

        # Pattern to search for: 'TILE' followed by numeric string, up to the next dash
        pattern = r'(TILE\w+)'

        for i, filename in enumerate(psf_images):

            # Search for the pattern
            match = re.search(pattern, filename)
            matched_text = match.group(1)
        
            # Add an underscore after 'TILE'
            new_filename = filename[:match.start(1) + len("TILE")] + '_' + filename[match.start(1) + len("TILE"):match.end(1)] + '.fits'
            os.system(f'mv {filter_dir}/{filename} {filter_dir}/{new_filename}')

