"""
stack_tables.py
"""

from astropy.table import Table, vstack, unique
import glob

filter_names = ['VIS', 'Y', 'J', 'H']

for filter_name in filter_names:

    files = glob.glob(f'd{filter_name}_*.fits')
    print(files)

    tables = []

    for filename in files:

        t = Table.read(filename)
        tables.append(t)

    stacked_table = vstack(tables)
    stacked_table = unique(stacked_table)
    stacked_table.write(f'd{filter_name}.fits', overwrite=True)

