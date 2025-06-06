#!/usr/bin/env python3

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
import json

plt.rcParams['axes.linewidth'] = 2.5
plt.rcParams.update({'font.size': 15})
plt.rcParams['figure.dpi'] = 100

# Retrieve filters from sys.argv
if len(sys.argv) > 1:
    run_type = sys.argv[1]
    #run_type = json.loads(run_type_json)
    filters_json = sys.argv[2]
    filters = json.loads(filters_json)
    field_name_json = sys.argv[3]
    field_name = json.loads(field_name_json)

else:
    raise ValueError("No filters passed to selection.py")

def apply_filters(table, filters):
    for filter_name, threshold in filters.items():
        if threshold['type'] == 'detection':
            table = table[table[f'flux_{filter_name}']/table[f'err_{filter_name}'] > threshold['value']]
            print(f"{threshold['value']}sigma detection in {filter_name}: ", len(table))
        elif threshold['type'] == 'non-detection':
            table = table[table[f'flux_{filter_name}']/table[f'err_{filter_name}'] < threshold['value']]
            print(f"{threshold['value']}sigma non-detection in {filter_name}: ", len(table))
        elif threshold['type'] == 'stacked-detection':
            stack_filters = filter_name.split('+')
            stacked_flux = table[f'flux_{stack_filters[0]}'] + table[f'flux_{stack_filters[1]}']
            stacked_err = np.sqrt(table[f'err_{stack_filters[0]}']**2 + table[f'err_{stack_filters[1]}']**2)
            table = table[stacked_flux/stacked_err > threshold['value']]
            print(f"{threshold['value']}sigma stacked detection in {stack_filters}: ", len(table))

    return table

def generate_selection_name(base_name, filters):
    parts = [base_name]
    for filter_name, threshold in filters.items():
        clean_filter_name = filter_name.replace('_DR3', '')
        clean_filter_name = clean_filter_name.replace('HSC-', 'HSC_')
        if threshold['type'] == 'detection':
            parts.append(f"{threshold['value']}sig_{clean_filter_name}")
        if threshold['type'] == 'stacked-detection':
            stack_filters = filter_name.split('+')
            parts.append(f"{threshold['value']}sig_{'_'.join(stack_filters)}")
        if threshold['type'] == 'non-detection':
            parts.append(f"nonDet_{clean_filter_name}")
    return '_'.join(parts) + '.fits'

def generate_input_name(filters, run_type, run_brown_dwarfs=False, run_dusty=False, run_lya=False):
    parts = ['det']
    for filter_name, threshold in filters.items():
        clean_filter_name = filter_name.replace('_DR3', '')
        if threshold['type'] == 'detection':
            parts.append(clean_filter_name)
        if threshold['type'] == 'stacked-detection':
            stack_filters = filter_name.split('+')
            for stack_filter in stack_filters:
                parts.append(stack_filter)
    # Add the run_type to the input name
    if run_type != '':
        parts.append(run_type)
    if run_brown_dwarfs:
        return '_'.join(parts) + '_bd.in'
    if run_dusty:
        return '_'.join(parts) + '_dusty.in'
    if run_lya:
        return '_'.join(parts) + '_lya.in'
    if not run_brown_dwarfs and not run_lya and not run_dusty:
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

    #! COSMOS
    input_cat_dir = Path.cwd().parents[3] / 'data' / 'catalogues' / 'finalCOSMOS' / 'other'
    #input_cat_name = "COSMOS_detYJH_masked_1.8as_Euclid_CWEB_2024_07_12.fits"
    #input_cat_name = 'COSMOS_detYJH_masked_1.8as_Euclid_CWEB_2024_10_07.fits' # Newest with JWST
    #input_cat_name = 'COSMOS_det_YJHK_masked_1.8as_Euclid_CWEB_2024_10_16.fits' # VISTA-Selected
    input_cat_name = 'COSMOSFULL_DR3_MASKVISTADET_HSC-Z_DR3_2025_06_05_1.8as_IRAC_2.8as_ALL.fits' #! HSC-Z selected
    output_save_dir = Path.cwd().parents[1] / 'data' / 'catalogues'
    base_output_name = "COSMOS"

    #! XMM
    # input_cat_dir = Path.cwd().parents[3] / 'data' / 'catalogues' / 'XMMFULL'
    # input_cat_name = 'XMMFULL_DR3_MASKVISTADET_HSC-Z_DR3_2.0as_IRAC2.8as_2024_01_18.fits'
    # output_save_dir = Path.cwd().parents[1] / 'data' / 'catalogues'
    # base_output_name = "XMM"
    
    input_name = generate_input_name(filters, run_type)
    os.environ['LEPHARE_INPUT'] = input_name

    main(input_cat_dir, input_cat_name, output_save_dir, base_output_name, filters)
