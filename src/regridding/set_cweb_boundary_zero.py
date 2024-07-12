#!/usr/bin/env python3

"""
The edge around the cweb mosaic is non-zero, so set everything outside my CWEB square to zero.

Created: Thursday 11th July 2024
"""

from astropy.io import fits
from astropy.wcs import WCS
import numpy as np
from pathlib import Path

filter_names = ['f115w', 'f150w','f277w', 'f444w']

image_dir = Path.cwd().parents[3] / 'data' / 