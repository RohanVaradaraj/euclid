"""
rename_with_tile_index.py

Created: Wednesday 20th March 2024.
"""

import numpy as np
from pathlib import Path
import glob
import os

euclid_dir = Path.home() / 'euclid'

tile_dict = np.load(Path.cwd().parent.parent / 'data' / 'mosaic' / 'tile_index_map.npy', allow_pickle=True)
print(tile_dict)

filter_names = ['VIS', 'Y', 'J', 'H']

field_name = 'COSMOS'

for filter_name in filter_names:

    image_dir = euclid_dir / filter_name / 'COSMOS'

    for tile, index in tile_dict.item().items():
        print(tile, index)

        image_stub = glob.glob(str(image_dir/ f'*BGSUB*{tile}*'))
        weight_stub = glob.glob(str(image_dir/ f'*RMS*{tile}*'))

        image = image_stub[0].split('/')[-1]
        new_image = image[:-5] + f'_{index}.fits' # Add tile index to file name


        weight = weight_stub[0].split('/')[-1]
        new_weight = weight[:-5] + f'_{index}.fits'

        os.system(f'mv {image_dir/image} {image_dir/new_image}')
        os.system(f'mv {image_dir/weight} {image_dir/new_weight}')


