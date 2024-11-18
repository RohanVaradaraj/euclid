#!/usr/bin/env python3
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
import re

vardy_dir = Path.home().parent.parent / 'mnt' / 'vardy' / 'vardygroupshare' / 'data'
euclid_dir = Path.home().parent.parent / 'euclid'


def skip_initial_lines(lines):
    """Skip the first three uncommented lines."""
    count = 0
    for idx, line in enumerate(lines):
        if not line.startswith("#"):
            count += 1
        if count == 3:
            return lines[idx + 1:]  # Skip three uncommented lines and return the rest
    return lines



def detect_section_changes(lines):
    """Detect the indices where sections change based on line length and abrupt value change."""
    line_lengths = [len(re.split(r'\s+', line.strip())) for line in lines]
    section_indices = []
    
    # Detect changes in line length
    for i in range(1, len(line_lengths)):
        if line_lengths[i] != line_lengths[i - 1]:
            section_indices.append(i)

    # Detect abrupt change in the first column value (from Redshift PDF to SED)
    for i in range(1, len(lines)):
        try:
            prev_val = float(re.split(r'\s+', lines[i - 1].strip())[0])
            curr_val = float(re.split(r'\s+', lines[i].strip())[0])

            # Check for a jump in value (e.g., from ~10 to ~1000)
            if (curr_val - prev_val > 500):
                section_indices.append(i)
                break
        except ValueError:
            continue

    return sorted(set(section_indices))



def split_sed_section(sed_lines):
    """Split the 'sed' section into multiple SEDs based on jumps in the wavelength."""
    sed_sections = []
    current_sed = []
    
    for i, line in enumerate(sed_lines):
        try:
            wavelength = float(re.split(r'\s+', line.strip())[0])
            
            # Detect the jump in wavelength (e.g., from 500 to 200)
            if current_sed and wavelength < float(re.split(r'\s+', current_sed[-1].strip())[0]):
                # New SED detected, save the current one and start a new one
                sed_sections.append(current_sed)
                current_sed = []

            # Add the current line to the ongoing SED
            current_sed.append(line)
        
        except ValueError:
            continue  # Skip lines that don't have a valid wavelength
            
    # Add the last SED section
    if current_sed:
        sed_sections.append(current_sed)
    
    return sed_sections

    

def parse_spec_file(filename):
    """Parses the SPEC file into multiple sections, handling headers appropriately."""
    
    # Read all lines from the file
    with open(filename, 'r') as file:
        lines = file.readlines()

    # Skip the first three uncommented lines
    data_lines = skip_initial_lines(lines)
    data_lines = data_lines[1:]  # Skip the first line of column names (if it's generic)

    # Detect section changes
    section_indices = detect_section_changes(data_lines)
    section_indices.append(len(data_lines))  # Add end of file as the last index

    # Split data into sections based on detected indices
    sections = []
    start_idx = 0
    for end_idx in section_indices:
        sections.append(data_lines[start_idx:end_idx])
        start_idx = end_idx

    # Parse sections into Astropy tables
    tables = {}
    section_names = ['model', 'phot', 'zpdf', 'sed']
    
    for idx, section in enumerate(sections):
        try:
            section_name = section_names[idx] if idx < len(section_names) else f"unknown_{idx}"

            # Define column names for each section
            if section_name == 'model':
                column_names = [
                    'ID', 'Nline', 'Model', 'Library', 'Nband', 'Zphot', 'Zinf', 'Zsup', 
                    'Chi2', 'PDF', 'Extlaw', 'EB-V', 'Lir', 'Age', 'Mass', 'SFR', 'SSFR'
                ]
            elif section_name == 'phot':
                column_names = [f'col{i+1}' for i in range(len(section[0].split()))]
            elif section_name == 'zpdf':
                column_names = ['z', 'P(z)']
            elif section_name == 'sed':
                # Split the 'sed' section into sub-sections
                sed_sections = split_sed_section(section)
                sed_tables = []

                for sed_idx, sed_section in enumerate(sed_sections):
                    # Define columns for SED section
                    column_names = ['lambda', 'flux']
                    data_rows = [re.split(r'\s+', line.strip()) for line in sed_section if line.strip()]
                    
                    # Create the table for this SED section
                    sed_table = Table(rows=data_rows, names=column_names)
                    sed_tables.append(sed_table)
                
                # Store the list of SED tables in the main table under 'sed'
                tables[section_name] = sed_tables
                continue  # Skip the default handling for 'sed' to avoid overwriting

            else:
                # For unknown sections, use generic column names
                column_names = [f'col{i+1}' for i in range(len(section[0].split()))]

            # Read the section into data rows (for non-SED sections)
            data_rows = [re.split(r'\s+', line.strip()) for line in section if line.strip()]

            # Create the table for the section
            tables[section_name] = Table(rows=data_rows, names=column_names)

        except Exception as e:
            print(f"Error reading section {section_name}: {e}")

    return tables



def remove_items(master_list, items_to_remove):
    for item in items_to_remove:
        if item in master_list:
            master_list.remove(item)
    return master_list



def buildLePhareLibrary(parameter_file: str, 
                        parameter_dir: Optional[Path] = Path.home() / 'lephare' / 'lephare_dev' / 'config',
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

    # Get code directory
    code_dir = Path.home().parents[1] / 'vardy' / 'vardygroupshare' / 'rohan' / 'euclid' /'src' / 'sed_fitting'

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

        ##################! $LEPHAREDIR NOT WORKING ###################

    #     # Build stars
    #     os.system(f'~/lephare/lephare_dev/source/sedtolib -t S -c ~/lephare/lephare_dev/config/{parameter_file}')

    #     # Build quasars
    #     os.system(f'~/lephare/l.'ephare_dev/source/sedtolib -t Q -c ~/lephare/lephare_dev/config/{parameter_file}')

    #     # Build galaxies
    #     os.system(f'~/lephare/lephare_dev/source/sedtolib -t G -c ~/lephare/lephare_dev/config/{parameter_file}')

    # if build_filters:

    #     os.system(f'~/lephare/lephare_dev/source/filter  -c ~/lephare/lephare_dev/config/{parameter_file}')

    # if build_mags:

    #     # Stars
    #     os.system(f'~/lephare/lephare_dev/source/mag_star -c  ~/lephare/lephare_dev/config/{parameter_file}')

    #     # Quasars
    #     os.system(f'~/lephare/lephare_dev/source/mag_gal  -t Q -c ~/lephare/lephare_dev/config/{parameter_file}')

    #     # Galaxies
    #     os.system(f'~/lephare/lephare_dev/source/mag_gal  -t G -c ~/lephare/lephare_dev/config/{parameter_file}')

    # Move back to code directory
    os.chdir(code_dir)

    return None



def runPhotometricRedshifts(parameter_file: str, zphot_dir: Path,
                        parameter_dir: Optional[Path] = Path.home().parent.parent / 'lephare' / 'lephare_dev' / 'config') -> None:

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

    Returns
    -------
    None
    """

    code_dir = Path.home().parents[1] / 'vardy' / 'vardygroupshare' / 'rohan' / 'euclid' /'src' / 'sed_fitting'
    
    # Move into the directory where we will output the .spec files
    os.chdir(zphot_dir)

    # Run the command
    os.system(f'$LEPHAREDIR/source/zphota -c $LEPHAREDIR/config/{parameter_file}')

    # Move back to the original directory
    os.chdir(code_dir)

    #################! $LEPHAREDIR NOT WORKING ####################
    #os.system(f'~/lephare/lephare_dev/source/zphota -c ~/lephare/lephare_dev/config/{parameter_file}')

    return None



def filter_widths():
    """
    Returns dictionary of filter central wavelengths and FWHMs.
    """

    filt_dict = {
        'HSC-G_DR3': (0.4816, 0.1386),
        'HSC-R_DR3': (0.6234, 0.1504),
        'HSC-I_DR3': (0.7741, 0.1552),
        'HSC-NB0816_DR3': (0.8177, 0.0113),
        'HSC-Z_DR3': (0.8912, 0.0773),
        'HSC-NB0921_DR3': (0.9214, 0.0134),
        'HSC-Y_DR3': (0.9780, 0.0783),
        'Y': (1.0214, 0.0926),
        'J': (1.2544, 0.1725),
        'H': (1.6465, 0.2916),
        'Ks': (2.1484, 0.3092),
        'f115w': (1.154, 0.225),
        'f150w': (1.501, 0.317),
        'f277w': (2.776, 0.673),
        'f444w': (4.401, 1.023),
        'VIS': (0.7180, 0.3699),
        'Ye': (1.0812, 0.2626),
        'Je': (1.3670, 0.3991),
        'He': (1.7708, 0.4994),
        'ch1cds': (3.551, 0.750),
        'ch2cds': (4.493, 1.010),
    }

    return filt_dict



def filter_files():
    """ 
    Returns the location of filter transmission curve files for input into the lephare config file.
    """
    filt_files = {
        'CFHT-u':'cfht/megacam/up.pb',
        'CFHT-g':'cfht/megacam/gp.pb',
        'CFHT-r':'cfht/megacam/rp.pb',
        'CFHT-z':'cfht/megacam/zp.pb',
        'HSC-G_DR3':'myfilters/HSC/g_HSC.txt',
        'HSC-R_DR3':'myfilters/HSC/r_HSC.txt',
        'HSC-I_DR3':'myfilters/HSC/i_HSC.txt',
        'HSC-NB0816_DR3':'myfilters/HSC/nb816_HSC.txt',
        'HSC-Z_DR3':'myfilters/HSC/z_HSC.txt',
        'HSC-NB0921_DR3':'myfilters/HSC/nb921_HSC.txt',
        'HSC-Y_DR3':'myfilters/HSC/y_HSC.txt',
        'Y':'myfilters/VISTA/VISTA_Y.txt',
        'J':'myfilters/VISTA/VISTA_J.txt',
        'H':'myfilters/VISTA/VISTA_H.txt',
        'Ks':'myfilters/VISTA/VISTA_Ks.txt',
        'f115w':'myfilters/JWST/f115w_angstroms.txt',
        'f150w':'myfilters/JWST/f150w_angstroms.txt',
        'f277w':'myfilters/JWST/f277w_angstroms.txt',
        'f444w':'myfilters/JWST/f444w_angstroms.txt',
        'VIS':'myfilters/Euclid/Euclid_VIS.txt',
        'Ye':'myfilters/Euclid/Euclid_Y.txt',
        'Je':'myfilters/Euclid/Euclid_J.txt',
        'He':'myfilters/Euclid/Euclid_H.txt',
        'ch1cds':'myfilters/SPITZER/irac_ch1.txt',
        'ch2cds':'myfilters/SPITZER/irac_ch2.txt',
    }

    return filt_files



def GenerateLePhareConfig(all_filters: list, det_filters: list, run_brown_dwarfs: bool, run_dusty: bool, run_lya: bool,
    file_name=Path.home() / 'lephare' / 'lephare_dev' / 'config' / 'euclid.para', filter_file_name='FILTERS.filt', z_step=[0.05, 10.0, 0.05]) -> None:
    """
    Generate a LePhare configuration file.

    Parameters
    ----------
    all_filters : list
        List of all filters to include in the config file. Usually all filters available in the cat
    det_filters : list
        List of filters to use for detection. Used to generate the input catalogue name.
    run_brown_dwarfs : bool
        If True, run the brown dwarf SED fitting without the bluest filters.
    run_lya : bool
        If True, run the Lyman-alpha template set (modified BC03)
    file_name : Path
        The name of the parameter file to save. Default is 'euclid.para'.

    Returns
    -------
    None : saves the .para file in the LePhare directory.
    """

    expanded_filename = os.path.expanduser(file_name)

    with open(expanded_filename, 'w') as f:

        f.write('##############################################################################\n')
        f.write('#                CREATION OF LIBRARIES FROM SEDs List                        #\n')
        f.write('#$LEPHAREDIR/source/sedtolib -t (S/Q/G) -c $LEPHAREDIR/config/rohan_test.para#\n')
        f.write('# help : $LEPHAREDIR/source/sedtolib -h (or -help)                           #\n')
        f.write('##############################################################################\n')
        f.write('#\n')
        f.write('#------      STELLAR LIBRARY (ASCII SEDs)\n')
        f.write('STAR_SED	$LEPHAREDIR/sed/STAR/DWARFSTARS/DWARFSTARS_MOD.list	# STAR list (full path)\n')
        f.write('STAR_FSCALE	1. ##3.432E-09				# Arbitrary Flux Scale\n')
        f.write('STAR_LIB	LIB_STAR				# Bin. STAR LIBRARY ->\n')
        f.write('							# $LEPHAREWORK/lib_bin\n')
        f.write('#\n')
        f.write('#------      QSO LIBRARY (ASCII SEDs) \n')
        f.write('QSO_SED		$LEPHAREDIR/sed/QSO/QSO_MOD.list        # QSO list (full path)\n')
        f.write('QSO_FSCALE	1. 	# Arbitrary Flux Scale\n')
        f.write('QSO_LIB   LIB_QSO				# Bin. QSO LIBRARY ->\n')
        f.write('							# $LEPHAREWORK/lib_bin\n')
        f.write('#\n')
        f.write('#------      GALAXY LIBRARY (ASCII or BINARY SEDs)\n')
        #! Lyman alpha/LBG template set
        if not run_lya:
            f.write('GAL_SED	          $LEPHAREDIR/sed/GAL/BC03/ASCII/tmp_BC03_MOD.lis  # GAL list (full path) # REGULAR BC03\n')
        else:
            f.write('GAL_SED	   $LEPHAREDIR/sed/GAL/BC03/ASCII/BC03_MOD.lis  # GAL list (full path) # LYMAN ALPHA\n')
        f.write('GAL_FSCALE	1. 	# Arbitrary Flux Scale\n')
        f.write('GAL_LIB   LIB_GAL				# Bin. GAL LIBRARY ->\n')
        f.write('							# $LEPHAREWORK/lib_bin\n')
        f.write('SEL_AGE   $LEPHAREDIR/sed/GAL/DEFAULT_BC03_CHAB/BC03_AGE.list # Age list(full path)\n')
        f.write('#SEL_AGE   NONE                                          # Age list(full path)\n')
        f.write('							# (def=NONE)	\n')
        f.write('AGE_RANGE  1.0e7,13.8e9                                     # Age Min-Max in yr\n')
        f.write('#\n')
        f.write('#############################################################################\n')
        f.write('#                       FILTERS                                             #\n')
        f.write('#  $LEPHAREDIR/source/filter  -c $LEPHAREDIR/config/rohan_test.para         #\n')
        f.write('#  help: $LEPHAREDIR/source/filter  -h (or -help)                           #\n')
        f.write('#############################################################################\n')
        f.write('\n')

        filter_names = all_filters

        if run_brown_dwarfs:
            filters_to_remove = ['CFHT-u', 'CFHT-g', 'CFHT-r', 'HSC-G_DR3', 'HSC-R_DR3', 'f277w', 'f444w', 'ch1cds', 'ch2cds']
            filter_names = remove_items(filter_names, filters_to_remove)
            print('Running brown dwarfs: blue filters and long-wavelength filters removed in config file.')

        if run_dusty == False:
            filters_to_remove = ['f444w', 'ch1cds', 'ch2cds']
            filter_names = remove_items(filter_names, filters_to_remove)
            print('Not running dusty galaxies: reddest filters removed in config file.')

        filter_dict = filter_files()

        filter_list = 'FILTER_LIST '
        for filter_name in filter_names:
            filter_list += f'{filter_dict[filter_name]},'
        filter_list = filter_list[:-1] # remove last comma

        f.write(filter_list + '\n')
        f.write('					# (in $LEPHAREDIR/filt/*)\n')
        f.write('TRANS_TYPE	0			# TRANSMISSION TYPE\n')
        f.write('                                        # 0[-def]: Energy, 1: Nb of photons\n')
        f.write('FILTER_CALIB    0                       # 0[-def]:  fnu=ctt \n')
        f.write('                                        # 1      :  nu.fnu=ctt\n')
        f.write('                                        # 2      :  fnu=nu\n')
        f.write('                                        # 3      :  fnu=Black Body @ T=10000K\n')
        f.write('                                        # 4      :  for MIPS (leff with nu fnu=ctt and flux with BB @ 10000K)\n')
        f.write(f'FILTER_FILE	{filter_file_name}		# output name of filters file  ->\n')
        f.write('                                        # $LEPHAREWORK/filt/\n')
        f.write('#\n')
        f.write('############################################################################\n')
        f.write('#                 THEORETICAL  MAGNITUDES                                  #\n')
        f.write('#  $LEPHAREDIR/source/mag_star -c $LEPHAREDIR/config/rohan_test.para       #\n')
        f.write('#  help: $LEPHAREDIR/source/mag_star -h (or -help)                         #\n')
        f.write('# $LEPHAREDIR/source/mag_gal  -t (Q or G) -c $LEPHAREDIR/config/rohan_test.para #\n')
        f.write('#                                                         (for gal. & QSO) #\n')
        f.write('# help: $LEPHAREDIR/source/mag_gal  -h (or -help)                          #\n')
        f.write('############################################################################\n')
        f.write('#\n')
        f.write('#-------     From STELLAR LIBRARY  \n')
        f.write('STAR_LIB_IN	LIB_STAR      # Input  STELLAR LIBRARY in $LEPHAREWORK/lib_bin/\n')
        f.write('STAR_LIB_OUT	STAR_EUC      # Output STELLAR MAGN    -> $LEPHAREWORK/lib_mag/\n')
        f.write('#\n')
        f.write('#-------     From QSO     LIBRARY  \n')
        f.write('QSO_LIB_IN	LIB_QSO	      # Input  QSO LIBRARY  in $LEPHAREWORK/lib_bin/\n')
        f.write('QSO_LIB_OUT	QSO_EUC      # Output QSO MAGN     -> $LEPHAREWORK/lib_mag/\n')
        f.write('#\n')
        f.write('#-------     From GALAXY  LIBRARY  \n')
        f.write('GAL_LIB_IN	LIB_GAL	      # Input  GAL LIBRARY  in $LEPHAREWORK/lib_bin/\n')
        f.write('GAL_LIB_OUT	GAL_EUC      # Output GAL MAGN     -> $LEPHAREWORK/lib_mag/\n')
        f.write('#\n')
        f.write('#-------   MAG + Z_STEP + EXTINCTION + COSMOLOGY \n')
        f.write('MAGTYPE         AB		     # Magnitude type (AB or VEGA)\n')

        z_string = str(z_step[0]) + ',' + str(z_step[1]) + ',' + str(z_step[2])

        f.write(f'Z_STEP 		{z_string} 	     # dz, zmax, dzsup(if zmax>6) \n')
        f.write('COSMOLOGY	70,0.3,0.7	     # H0,om0,lbd0    (if lb0>0->om0+lbd0=1)\n')
        f.write('MOD_EXTINC 	0,12		     # model range for extinction \n')
        f.write('EXTINC_LAW	calzetti.dat	     # ext. law (in  $LEPHAREDIR/ext/*)\n')
        f.write('#EXTINC_LAW	NONE	     # ext. law (in  $LEPHAREDIR/ext/*)\n')

        base_Av = [0.,0.025,0.05,0.075,0.1,0.125,0.15,0.175,0.2,0.225,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1.0]
        base_Av = [str(x) for x in base_Av]
        Av_string = ','.join(base_Av)

        if run_dusty:
            extra_Av = [1.2,1.4,1.6,1.8,2.0,2.2,2.4,2.6,2.8,3.0,3.4,3.6,3.8,4.0]
            extra_Av = [str(x) for x in extra_Av]
            Av_string += ',' + ','.join(extra_Av)

        Av_string = Av_string + '# E(B-V) (<50 values)'
        f.write(f'EB_V            {Av_string}\n')
        f.write('#EB_V		0.\n')
        f.write('EM_LINES        NO\n')
        f.write('# Z_FORM 	8,7,6,5,4,3 	     # Zformation for each SED in GAL_LIB_IN\n')
        f.write('#\n')
        f.write('#-------   ASCII OUTPUT FILES OPTION\n')
        f.write('LIB_ASCII       NO		     # Writes output in ASCII\n')
        f.write('                                     # in working directory\n')
        f.write('#\n')
        f.write('############################################################################\n')
        f.write('#              PHOTOMETRIC REDSHIFTS                                       #\n')
        f.write('# $LEPHAREDIR/source/zphota -c $LEPHAREDIR/config/rohan_test.para          #\n')
        f.write('# help: $LEPHAREDIR/source/zphot -h (or -help)                             #\n')
        f.write('############################################################################\n')
        f.write('#\n')
        f.write('#-------    Input Catalog Information \n')

        det_string = 'det_' + '_'.join(det_filters)
        if run_brown_dwarfs:
            det_string += '_bd.in'
        if run_lya:
            det_string += '_lya.in'
        if run_dusty:
            det_string += '_dusty.in'
        if not run_brown_dwarfs and not run_lya and not run_dusty:
            det_string += '.in'
        print('####################')
        print(det_string)

            
        cat_in = '/mnt/hoy/temporaryFilesROHAN/lephare/inputs/euclid/'  + det_string 

        f.write(f'CAT_IN		{cat_in}		# Input Catalog (full path)\n')
        f.write('INP_TYPE     F		          # Input type      (F:Flux or M:MAG)\n')
        f.write('CAT_MAG      AB                   # Input Magnitude (AB or VEGA)\n')
        f.write('CAT_FMT      MEME		  # MEME: (Mag,Err)i , MMEE: (Mag)i,(Err)i\n')
        f.write('CAT_LINES    1,300000             #  MIN and MAX RANGE of ROWS used in input cat [def:-99,-99]\n')
        f.write('CAT_TYPE     SHORT	          # Input Format    (LONG,SHORT-def)\n')

        # Use the input name to determine the .out name (replace .in with .out)
        cat_out = det_string.replace('.in', '.out')

        f.write(f'CAT_OUT	     $LEPHAREDIR/test/{cat_out}	# Output catalog (full path)\n')
        f.write('PARA_OUT     $LEPHAREDIR/config/zphot_output.para  # Ouput parameter (full path)\n')
        f.write('BD_SCALE     0		          # Bands used for scaling, (Sum 2^n; n=0->nbd-1, 0[-def]:all bands) \n')
        f.write('GLB_CONTEXT  -1		          # Overwrite Context (Sum 2^n; n=0->nbd-1, 0 : all bands used, -1[-def]: used context per object) \n')
        f.write('# FORB_CONTEXT -1                   # context for forbitten bands \n')
        f.write('# ERR_SCALE  0.03,0.02,0.02,0.02,0.04,0.04,0.04  # errors per band added in quadrature\n')
        f.write('ERR_FACTOR  1.0                    # error scaling factor 1.0 [-def] \n')
        f.write('#\n')
        f.write('#-------    Theoretical libraries  \n')
        f.write('ZPHOTLIB    GAL_EUC,STAR_EUC\n')
        f.write('ADD_EMLINES NO \n')
        f.write('#\n')
        f.write('########    PHOTOMETRIC REDSHIFTS OPTIONS      ###########\n')
        f.write('# FIR LIBRARY \n')
        f.write('FIR_LIB     NONE\n')
        f.write('FIR_LMIN         7.0           # Lambda Min (micron) for FIR analysis \nFIR_CONT        -1 \nFIR_SCALE       -1 \nFIR_FREESCALE    YES             # ALLOW FOR FREE SCALING \nFIR_SUBSTELLAR   NO\n')        
        f.write('# PHYSICAL LIBRARY with Stochastic models from  BC07 \nPHYS_LIB      NONE \nPHYS_CONT    -1 \nPHYS_SCALE   -1 \nPHYS_NMAX     100000 \n')
        f.write('#\n')

        #! Priors
        f.write('#-------     Priors \n')
        f.write('# MASS_SCALE	6.,16.		 # Lg(Scaling) min,max [0,0-def]\n')
        f.write('MAG_ABS 	-10.,-30.	 # Mabs_min , Mabs_max [0,0-def]\n')
        f.write('MAG_REF 	4		 # Reference number for band used by Mag_abs\n')
        f.write('# ZFORM_MIN	5,5,5,5,5,5,3,1	 # Min. Zformation per SED -> Age constraint\n')
        f.write('Z_RANGE        0.,99.99          # Z min-max used for the Galaxy library \n')
        f.write('EBV_RANGE      0,9               # E(B-V) MIN-MAX RANGE of E(B-V) used \n')
        f.write('# NZ_PRIOR      4,2,4                # I Band for prior on N(z)\n')
        f.write('# \n#-------     Fixed Z   (need format LONG for input Cat) \nZFIX		NO		 # fixed z and search best model [YES,NO-def] \n# \n#-------     Parabolic interpolation for Zbest \nZ_INTERP	YES		 # redshift interpolation [YES,NO-def] \n#')
        f.write('#-------  Analysis of normalized ML(exp-(0.5*Chi^2)) curve \n#-------  Secondary peak analysis \nDZ_WIN          0.25              # Window search for 2nd peaks [0->5;0.25-def] \nMIN_THRES       0.0              # Lower threshold for 2nd peaks[0->1; 0.1-def] \n#')
        f.write('#-------  Probability (in %) per redshift intervals\n')
        f.write('# PROB_INTZ     0,0.5,0.5,1.,1.,1.5     # even number \n')
        f.write('#\n')

        f.write('#########    ABSOLUTE MAGNITUDES COMPUTATION   ###########\n')
        f.write('MABS_METHOD	3		 # 0[-def] : obs->Ref, 1 : best  obs->Ref, 2 : fixed obs->Ref, 3 : mag from best SED, 4 : Zbin\n')
        f.write('MABS_CONTEXT    -1               # CONTEXT for Band used for MABS \n')
        f.write('MABS_REF	4		 # 0[-def]: filter obs chosen for Mabs : ONLY USED IF MABS_METHOD=2\n')
        f.write('MABS_FILT       1,2,3,4   	 # Chosen filters per redshift bin (MABS_ZBIN) : ONLY USED IF MABS_METHOD=4 \n')
        f.write('MABS_ZBIN       0,0.5,1,1.5,2,3,3.5,4 # Redshift bins (even number) : ONLY USED IF MABS_METHOD=4\n')

        #! Output spectra
        f.write('#########   OUTPUT SPECTRA                     ###########\n')
        
        #! SPEC_OUT
        f.write('SPEC_OUT	YES	 	 # spectrum for each object?  [YES,NO-def]\n')

        #! CHI2_OUT
        f.write('CHI2_OUT        NO               # output file with all values : z,mod,chi2,E(B-V),... takes up a lot of space! \n')

        #! PDZ
        f.write('#########  OUTPUT PDZ ANALYSIS  ###########\n')
        f.write('PDZ_OUT         NONE             # pdz output file name [def-NONE]\n')
        f.write('                                 # add automatically PDZ_OUT[.pdz/.mabsx/.mod/.zph] \n')
        f.write('PDZ_MABS_FILT   2,10,14           # MABS for REF FILTERS to be extracted \n')
        f.write('#\n')
        
        f.write('#########   FAST MODE : color-space reduction        #####\n')
        f.write('FAST_MODE	NO 		 # Fast computation                [NO-def] \n')
        f.write('COL_NUM		3 		 # Number of colors used           [3-def]\n')
        f.write('COL_SIGMA	3		 # Enlarge of the obs. color-errors[3-def]\n')
        f.write('COL_SEL		AND		 # Combination between used colors [AND/OR-def]\n')
        f.write('#\n')
        f.write('#########   MAGNITUDE SHIFTS applied to libraries   ######\n')
        f.write('# APPLY_SYSSHIFT  0.  # Apply systematic shifts in each band. Used only if number of shifts matches with number of filters in the library.\n')
        f.write('#########   ADAPTIVE METHOD using Z spectro sample     ###\n')
        f.write('AUTO_ADAPT	NO		 # Adapting method with spectro [NO-def]\n')
        f.write('ADAPT_BAND 	4,2,4		 # Reference band, band1, band2 for color \n')
        f.write('ADAPT_LIM       18,22.0		 # Mag limits for spectro in Ref band [18,21.5-def]\n')
        f.write('ADAPT_POLY	1		 # Number of coef in  polynom (max=4) [1-def]\n')
        f.write('ADAPT_METH      1		 # Fit as a function of 1 : Color Model  [1-def], 2 : Redshift, 3 : Models\n')
        f.write('ADAPT_CONTEXT  -1                # Context for bands used for training. -1[-def] used context per object \n')
        f.write('ADAPT_ZBIN     0.01,6            # Redshift interval used for trainin g[0.001,6-Def]\n')
        f.write('ADAPT_MODBIN   1,1000            # Model interval    used for training [1,1000-Def]\n')
        f.write('ERROR_ADAPT     NO               # [YES,NO-def] Add error in quadrature according to the difference between observed and predicted apparent magnitudes \n')
        f.write('#\n')

        print(f'Generated LePhare config file: {file_name}')

        return None
