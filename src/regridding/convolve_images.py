"""
convolve_images.py

Convolve the Euclid/JWST imaging with the VISTA PSF kernel.

Created: Friday 14th June 2024
"""

from astropy.io import fits
from astropy.convolution import convolve
