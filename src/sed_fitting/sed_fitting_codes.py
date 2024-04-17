"""
sed_fitting_codes.py

Functions used in the SED fitting pipeline.

Created: Wednesday 17th April 2024.
"""

from astropy.table import Table, Column
from astropy.io import fits, ascii
import numpy as np
from typing import Optional
from pathlib import Path
import os

vardy_dir = Path.home().parent.parent / 'mnt' / 'vardy' / 'vardygroupshare' / 'data'
euclid_dir = Path.home().parent.parent / 'euclid'


def convertFitsToText(cat_name: str, ground_filters: list[str], euclid_filters: list[str], 
                      cat_dir: Optional[str] = Path.cwd().parent.parent / 'data' / 'catalogues',
                      field: Optional[str] = 'COSMOS', append_unused: bool = False,
                      add_zspec: bool = False, zspec_colname: Optional[str] = None) -> None:
    
    """
    LePhare requires a .txt file of the fluxes and errors in a particular order.
    This function takes in a catalogue and processes it into the required format, saving it in the correct LePhare directory.

    Parameters
    ----------
    cat_name : str
        The name of the catalogue to be converted.
    ground_filters : list[str]
        The names of the ground-based filters in the catalogue.
    euclid_filters : list[str]
        The names of the Euclid filters in the catalogue.3
    cat_dir : Optional[str], optional
        The directory containing the catalogue, by default is 'euclid' / 'data' / 'catalogues'
    field : Optional[str], optional
        The field of the catalogue, by default is 'COSMOS'.
        For now, we only have Euclid PV data in the COSMOS field, but this keeps the function general for the future!
    append_unused : bool, optional
        If True, the function will append unused filters to the end of the .txt file, by default False.

    Returns
    -------
    None
        Directly writes the text file to ~/mnt/hoy/temporaryFilesROHAN/lephare/inputs/.
    """

    output_dir = Path.home().parent.parent / 'mnt' / 'hoy' / 'temporaryFilesROHAN' / 'lephare' / 'inputs'

    # Replace dashes with underscores in the filter names
    ground_filters = [filt.replace('-', '_') for filt in ground_filters]

    # Read in images.lis
    ground_info = Table.read(vardy_dir / field / 'images.lis', format='ascii.commented_header')
    euclid_info = Table.read(euclid_dir / 'data' / 'images.lis', format='ascii.commented_header')

    avail_ground = ground_info['Name']
    avail_euclid = euclid_info['Name']

    # Open the catalogue
    t = Table.read(cat_dir / cat_name)

    #! Create the table in the correct format. First: ground based. Then: Euclid.

    # Make base lists for the column names. We will use this to rearrange the columns later.
    column_names = ['ID']

    # If there are redshifts, get these
    if add_zspec:
        redshifts = t[zspec_colname]

    # Loop through ground-based filters and add to start of list
    for filter_name in ground_filters:

        column_names = column_names + ['flux_{0}'.format(filter_name)]
        column_names = column_names + ['err_{0}'.format(filter_name)]

    # Now do the same for the Euclid filters
    for filter_name in euclid_filters:

        column_names = column_names + ['flux_{0}'.format(filter_name)]
        column_names = column_names + ['err_{0}'.format(filter_name)]

    #! Deprecated
    if append_unused:
        remainder_ground = list(set(avail_ground) - set(ground_filters))
        remainder_euclid = list(set(avail_euclid) - set(euclid_filters))

        for filter_name in remainder_ground:

            column_names = column_names + ['flux_{0}'.format(filter_name)]
            column_names = column_names + ['err_{0}'.format(filter_name)]

        for filter_name in remainder_euclid:

            column_names = column_names + ['flux_{0}'.format(filter_name)]
            column_names = column_names + ['err_{0}'.format(filter_name)]

    n = len(ground_filters) + len(euclid_filters)

    # Compute the LePhare context
    sum = 0
    for i in range(n):
        sum += 2**1
    context = sum
    print(context)

    t['Context'] = context
    t['ID'] = t['ID'].astype(int)
    column_names = column_names + ['Context']

    # Rearrange table
    t = t[column_names]

    # Add spetroscopic redshifts if they exist
    if add_zspec:
        t['z_spec'] = redshifts

    # Save the table
    ascii.write(t, output_dir / cat_name.replace('.fits', '.txt'), format='commented_header', overwrite=True)
    print('Saved to: ', output_dir / cat_name.replace('.fits', '.txt'))
    print(t)

    return None



def buildLePhareLibrary(parameter_file: str, 
                        parameter_dir: Optional[Path] = Path.home().parent.parent / 'lephare' / 'lephare_dev' / 'config',
                        build_libs: Optional[bool] = False, build_filters: Optional[bool] = False, build_mags: Optional[bool] = False) -> None:

    """
    Build the LePhare library.

    Parameters
    ----------
    parameter_file : str
        The name of the parameter file to use.
    parameter_dir : Optional[Path], optional
        The directory containing the parameter file, by default is ~/ 'lephare' / 'lephare_dev' / 'config'.
    build_libs : bool, optional
        If True, build the libraries, by default False.
    build_filters : bool, optional
        If True, build the filters, by default False.
    build_mags : bool, optional
        If True, build the theoretical magnitudes, by default False.

    Returns
    -------
    None
    """

    # Move into lephare directory
    os.chdir(str(parameter_dir))

    if build_libs:

        # Build stars
        os.system(f'$LEPHAREDIR/source/sedtolib -t S -c $LEPHAREDIR/config/{parameter_file}')

        # Build quasars
        os.system(f'$LEPHAREDIR/source/sedtolib -t Q -c $LEPHAREDIR/config/{parameter_file}')

        # Build galaxies
        os.system(f'$LEPHAREDIR/source/sedtolib -t G -c $LEPHAREDIR/config/{parameter_file}')

    if build_filters:

        os.system(f'$LEPHAREDIR/source/filter  -c $LEPHAREDIR/config/{parameter_file}')

    if build_mags:

        # Stars
        os.system(f'$LEPHAREDIR/source/mag_star -c  $LEPHAREDIR/config/{parameter_file}')

        # Quasars
        os.system(f'$LEPHAREDIR/source/mag_gal  -t Q -c $LEPHAREDIR/config/{parameter_file}')

        # Galaxies
        os.system(f'$LEPHAREDIR/source/mag_gal  -t G -c $LEPHAREDIR/config/{parameter_file}')

    return None



def runPhotometricRedshifts(parameter_file: str, zphot_dir: Path,
                        parameter_dir: Optional[Path] = Path.home().parent.parent / 'lephare' / 'lephare_dev' / 'config', 
                        overwrite: bool = False) -> None:

    """
    Run LePhare photometric redshifts

    Parameters
    ----------
    parameter_file : str
        The name of the parameter file to use.
    zphot_dir : Path
        The directory to output the .spec files.
    parameter_dir : Optional[Path], optional
        The directory containing the parameter file, by default is ~/ 'lephare' / 'lephare_dev' / 'config'.
    overwrite : bool, optional
        If True, remove previous files, by default False.
    

    Returns
    -------
    None
    """
    
    # Move into the directory where we will output the .spec files
    os.chdir(zphot_dir)

    # Remove previous files, if desired.
    if overwrite:
        os.system('rm ./*')

    # Run the command
    os.system(f'$LEPHAREDIR/source/zphota -c $LEPHAREDIR/config/{parameter_file}')

    return None
