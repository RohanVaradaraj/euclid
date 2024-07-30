"""
Initial selection of high-z Lyman break candidates.

Created: Friday 12th July 2024.

"""

import numpy as np
from astropy.table import Table
from pathlib import Path
import matplotlib.pyplot as plt
import sys
import os

plt.rcParams['axes.linewidth'] = 2.5
plt.rcParams.update({'font.size': 15})
plt.rcParams['figure.dpi'] = 100

def apply_filters(table, filters):
    for filter_name, threshold in filters.items():
        if threshold['type'] == 'detection':
            table = table[table[f'flux_{filter_name}']/table[f'err_{filter_name}'] > threshold['value']]
            print(f"{threshold['value']}sigma detection in {filter_name}: ", len(table))
        elif threshold['type'] == 'non-detection':
            table = table[table[f'flux_{filter_name}']/table[f'err_{filter_name}'] < threshold['value']]
            print(f"{threshold['value']}sigma non-detection in {filter_name}: ", len(table))
    return table

def generate_selection_name(base_name, filters):
    parts = [base_name]
    for filter_name, threshold in filters.items():
        clean_filter_name = filter_name.replace('_DR3', '')
        clean_filter_name = clean_filter_name.replace('HSC-', 'HSC_')
        if threshold['type'] == 'detection':
            parts.append(f"{threshold['value']}sig_{clean_filter_name}")
        else:
            parts.append(f"nonDet_{clean_filter_name}")
    return '_'.join(parts) + '.fits'

def generate_input_name(base_name, filters):
    parts = ['det']
    for filter_name, threshold in filters.items():
        clean_filter_name = filter_name.replace('_DR3', '')
        if threshold['type'] == 'detection':
            parts.append(clean_filter_name)
    return '_'.join(parts) + '.in'



def main(input_cat_dir, input_cat_name, output_save_dir, base_output_name, filters):
    cat_dir = Path(input_cat_dir)
    cat_name = input_cat_name

    t = Table.read(cat_dir / cat_name)

    print('INITIAL LENGTH: ', len(t))

    # Convert t['RA'] column from string to float
    t['RA'] = t['RA'].astype(float)

    # Apply filters
    t = apply_filters(t, filters)

    # Generate output file name
    output_save_name = generate_selection_name(base_output_name, filters)

    # Set output_save_name to environment variable
    os.environ['SELECTION_CATALOGUE'] = output_save_name


    # Save table
    save_dir = Path(output_save_dir)
    t.write(save_dir / output_save_name, overwrite=True)
    print(f"Filtered catalogue saved as {save_dir / output_save_name}")

# Example usage
if __name__ == "__main__":

    input_cat_dir = Path.cwd().parents[3] / 'data' / 'catalogues' / 'finalCOSMOS' / 'other'
    input_cat_name = "COSMOS_detYJH_masked_1.8as_Euclid_CWEB_2024_07_12.fits"
    output_save_dir = Path.cwd().parents[1] / 'data' / 'catalogues'
    base_output_name = "COSMOS"
    
    filters = {
        'Ye': {'type': 'detection', 'value': 5},
        'Y': {'type': 'detection', 'value': 2},
        'HSC-G_DR3': {'type': 'non-detection', 'value': 2},
        'HSC-R_DR3': {'type': 'non-detection', 'value': 2},
        'HSC-I_DR3': {'type': 'non-detection', 'value': 2},
    }

    # filters = {
    #     'Je': {'type': 'detection', 'value': 5},
    #     'J': {'type': 'detection', 'value': 3},
    #     'HSC-G_DR3': {'type': 'non-detection', 'value': 2},
    #     'HSC-R_DR3': {'type': 'non-detection', 'value': 2},
    #     'HSC-I_DR3': {'type': 'non-detection', 'value': 2},
    #     'HSC-Z_DR3': {'type': 'non-detection', 'value': 2},
    #     'HSC-Y_DR3': {'type': 'non-detection', 'value': 2},
    #     'Y': {'type': 'non-detection', 'value': 2},
    # }

    # Imposing detection in narrowband
    # filters = {
    #     'HSC-NB0921_DR3': {'type': 'detection', 'value': 5},
    #     'Ye': {'type': 'detection', 'value': 3},
    #     'HSC-G_DR3': {'type': 'non-detection', 'value': 2},
    #     'HSC-R_DR3': {'type': 'non-detection', 'value': 2},
    #     'HSC-I_DR3': {'type': 'non-detection', 'value': 2},
    # }


    input_name = generate_input_name(base_output_name, filters)
    os.environ['LEPHARE_INPUT'] = input_name

    main(input_cat_dir, input_cat_name, output_save_dir, base_output_name, filters)
