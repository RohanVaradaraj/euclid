"""
Quick script to just plot the PSFs of some filters in some directory.

Created: Thursday 5th March 2026.
"""

from pathlib import Path
from astropy.io import fits
import numpy as np
import os
import matplotlib.pyplot as plt

plt.rcParams['axes.linewidth'] = 2.5
plt.rcParams.update({'font.size': 15})
plt.rcParams['figure.dpi'] = 100

filter_names = ['YE', 'JE', 'HE', 'VIS']

psf_dir = Path.cwd().parents[3] / 'data' / 'psf' / 'CDFS1' / 'results'
