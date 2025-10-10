"""
Export the BAGPIPES fits to a format the JWST ETC likes.

Created: Wednesday 8th October 2025.
"""

from astropy.table import Table
from pathlib import Path
import numpy as np

output_path = Path.cwd().parents[1] / 'data' / 'jwst_proposal' / 'etc'
pipes_path = Path.cwd().parent / 'sed_fitting' / 'pipes'

IDs = [178396, 387777, 615679, 878786, 879369]
IDs = [178396]

for ID in IDs:
    t = Table.read(pipes_path / f'{ID}_microJy.fits')

    # Remove the restwave column
    t.remove_column('restwave')

    # Convert wave from Angstroms to microns
    t['wave'] *= 1e-4

    # Convert flux from microJy to milliJy
    t['flux'] *= 1e-3

    wlen = t['wave']
    flux = t['flux']

    # Write to a .txt file with column names WAVELENGTH and FLUX
    output_file = output_path / f'ID_{ID}.txt'
    with open(output_file, 'w') as f:
        f.write('# Wavelength (microns)    Flux (mJy)\n')
        for wl, fl in zip(wlen, flux):
            f.write(f'{wl:.6f}    {fl:.6f}\n')
