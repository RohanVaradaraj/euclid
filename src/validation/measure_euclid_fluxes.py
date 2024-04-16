"""
measure_euclid_fluxes.py

Get fluxes of existing objects in the Euclid filters, in a consistent manner with ground-based photometry.

Created: Tuesday 16th April 2024.
"""

from pathlib import Path
import numpy as np
import sep
import matplotlib.pyplot as plt
from astropy.io import fits, ascii
from astropy.table import Table
from astropy.wcs import WCS
import sys
sys.path.append(str(Path.cwd().parent))

from cutouts.cutout_codes import isCoordInSurveyFootprints


def grid_depths(gridTable: dict, x: np.ndarray, y: np.ndarray, faster: bool = True, verbose: bool = False, nearby: bool = False) -> np.ndarray:
    ''' 
    Code to find the closest depth measurement from my previous analysis. Faster than truly local depths

    Parameters:
    gridTable (dict): Dictionary containing grid information.
    x (np.ndarray): Array of x coordinates.
    y (np.ndarray): Array of y coordinates.
    faster (bool, optional): If True, use a faster method. Defaults to True.
    verbose (bool, optional): If True, enable verbose output. Defaults to False.
    nearby (bool, optional): If True, consider nearby pixels. Defaults to False.

    Returns:
    np.ndarray: Array of depth values.
    '''

    xgrid = gridTable['x']
    ygrid = gridTable['y']
    keys = gridTable.colnames

    depthsOverField = gridTable['depths']

    ## Make an output array
    depthArray = np.zeros(x.size)
    depthArray[:] = -99.0

    if faster:
        if verbose:
            print("Using faster method.")
            print("Input array size is ", x.size)
        deltay = np.min(ygrid)
        deltax = np.min(xgrid)

    ## loop through the grid instead of each object
        for xi in range(xgrid.size):
            
            xmin = xgrid[xi] - deltax
            xmax = xgrid[xi] + deltax
            ymin = ygrid[xi] - deltay
            ymax = ygrid[xi] + deltay

            ii = (x > xmin) & (x <= xmax) & (y > ymin) & (y <= ymax)
            
            depthArray[ii] = depthsOverField[xi]

    else:
        
    ## Find the closest point to the objects x and y positions
    ## Loop!
        for xi in range(x.size):
            
        ## make a radius array
            deltax = (xgrid - x[xi])
            deltay = (ygrid - y[xi])
            radius = np.sqrt(deltax*deltax + deltay*deltay)
            mini = np.argmin(radius)
            
        ## try using argpartition
            numpoints = 10
            idx = np.argpartition(radius, numpoints)
            
            if nearby:
                
                mini = idx[0:numpoints]
                print("The nearby depths are = ", depthsOverField[mini])
                print("Before = ", depthsOverField[mini][0])
            

            depthArray[xi] = depthsOverField[mini][0]

    return depthArray


def measure_euclid_fluxes(ra: np.ndarray, dec: np.ndarray, filter_names: list, aperture_diameter: float = 1.0) -> tuple[np.ndarray, np.ndarray]:
    """
    Measure the fluxes of objects in the Euclid filters.

    Parameters
    ----------
    ra : np.ndarray
        Array of RA values of objects.
    dec : np.ndarray
        Array of Dec values of objects.
    filters : list
        List of filters we want to measure fluxes in.
    aperture_diameter : float, optional
        Diameter of the aperture, in arcseconds. Default is 1.0.

    Returns
    -------
    np.ndarray
        Array of fluxes.
    np.ndarray
        Array of errors.
    """

    euclid_dir = Path().home() / 'euclid' / 'COSMOS'
    depth_dir = Path.cwd().parent.parent / 'data' / 'depths' / 'COSMOS' / 'phot'
    psf_dir = Path.cwd().parent.parent / 'data' / 'psf' / 'COSMOS' / 'enclosedflux'

    # Initialize a dictionary to hold the data
    data_dict = {'ra': ra, 'dec': dec}

    # Open the Euclid images.lis

    image_list = euclid_dir / 'images.lis'
    info = ascii.read(image_list)

    # Loop over the filters
    for filter_name in filter_names:

        print(f'Calculating fluxes in {filter_name}...')

        # Get desired info
        this_info = info[info['Name'] == filter_name]

        image_name = this_info['Image'][0]
        weight_name = this_info['Weight'][0]
        zeropoint = this_info['zeropoint'][0]

        # Open the image
        with fits.open(euclid_dir / image_name, memmap=True) as hdu:
            image = hdu[0].data.byteswap().newbyteorder() # Change byte order to native
            header = hdu[0].header

        # Open the depth table
        depth_table = Table.read(depth_dir / f'{filter_name}_{aperture_diameter}as_gridDepths_300_200.fits')

        # Loop through the RA and DEC
        assert len(ra) == len(dec), "RA and DEC arrays must have the same length."

        # Convert to pixel coordinates
        w = WCS(header)
        x, y = w.all_world2pix(ra, dec, 0)

        # Find the depth at this x and y
        depths_here = grid_depths(depth_table, x, y)

        # If the object is in the footprint, it may not have a depth (e.g. edge of a tile).
        # But if it does not have a depth, it is definitely not in the footprint.
        #? Add some way of removing -99s

        # x, y must be array-like
        x = np.array(x)
        y = np.array(y)

        # Use sep to measure flux of object
        flux_counts, _, _ = sep.sum_circle(image, x, y, aperture_diameter/2.)

        # Convert counts to flux
        value = -(48.6 + zeropoint)/2.5
        fluxFinal = (10**value)*flux_counts

        # Calculate error from depths
        value = -(48.6 + depths_here)/2.5
        errors = 0.2*(10**value)

        # Aperture correction
        aperture_table = Table.read(psf_dir / f'{filter_name}_peak.txt', format='ascii')
        aperture_table = aperture_table[aperture_table['apD'] == aperture_diameter]

        enclosed_flux = aperture_table['ef'][0]
        print(enclosed_flux)

        fluxFinal = np.array([f*100/enclosed_flux for f in fluxFinal])

        # Wherever the depth is -99, set the error and flux to -99
        errors[depths_here == -99] = -99
        fluxFinal[depths_here == -99] = -99

        # Store fluxes and errors in the data dictionary
        data_dict[f'flux_{filter_name}'] = fluxFinal
        data_dict[f'err_{filter_name}'] = errors

    # Convert the data dictionary to an Astropy table
    flux_table = Table(data_dict)

    # Return the table
    return flux_table

def measure_ground_fluxes(ra: np.ndarray, dec: np.ndarray, filter_names: list, aperture_diameter: float = 1.8) -> tuple[np.ndarray, np.ndarray]:
    """
    Measure the fluxes of objects in the Euclid filters.

    Parameters
    ----------
    ra : np.ndarray
        Array of RA values of objects.
    dec : np.ndarray
        Array of Dec values of objects.
    filters : list
        List of filters we want to measure fluxes in.
    aperture_diameter : float, optional
        Diameter of the aperture, in arcseconds. Default is 1.8.

    Returns
    -------
    np.ndarray
        Array of fluxes.
    np.ndarray
        Array of errors.
    """

    ground_dir = Path().home().parent.parent / 'vardy' / 'vardygroupshare' / 'data' / 'COSMOS'
    depth_dir = Path().home().parent.parent / 'vardy' / 'vardygroupshare' / 'data' / 'depths' / 'COSMOS' / 'phot'
    psf_dir = Path().home().parent.parent / 'vardy' / 'vardygroupshare' / 'data' / 'psf' / 'COSMOS' / 'enclosedflux'

    # Initialize a dictionary to hold the data
    data_dict = {'ra': ra, 'dec': dec}

    # Open the Euclid images.lis
    image_list = ground_dir / 'images.lis'
    info = ascii.read(image_list)

    # Loop over the filters
    for filter_name in filter_names:

        print(f'Calculating fluxes in {filter_name}...')

        # Get desired info
        this_info = info[info['Name'] == filter_name]

        image_name = this_info['Image'][0]
        weight_name = this_info['Weight'][0]
        zeropoint = this_info['zeropoint'][0]
        directory = this_info['directory'][0]

        # Open the image
        if directory == 'here':
            with fits.open(ground_dir / image_name, memmap=True) as hdu:
                image = hdu[0].data.byteswap().newbyteorder() # Change byte order to native
                header = hdu[0].header
        else:
            with fits.open(directory / image_name, memmap=True) as hdu:
                image = hdu[0].data.byteswap().newbyteorder()
                header = hdu[0].header

        # Open the depth table
        depth_table = Table.read(depth_dir / f'{filter_name}_{aperture_diameter}as_gridDepths_300_200.fits')

        # Loop through the RA and DEC
        assert len(ra) == len(dec), "RA and DEC arrays must have the same length."

        # Convert to pixel coordinates
        w = WCS(header)
        x, y = w.all_world2pix(ra, dec, 0)

        # Find the depth at this x and y
        depths_here = grid_depths(depth_table, x, y)

        # If the object is in the footprint, it may not have a depth (e.g. edge of a tile).
        # But if it does not have a depth, it is definitely not in the footprint.
        #? Add some way of removing -99s

        # x, y must be array-like
        x = np.array(x)
        y = np.array(y)

        # Use sep to measure flux of object
        flux_counts, _, _ = sep.sum_circle(image, x, y, aperture_diameter/2.)

        # Convert counts to flux
        value = -(48.6 + zeropoint)/2.5
        fluxFinal = (10**value)*flux_counts

        # Calculate error from depths
        value = -(48.6 + depths_here)/2.5
        errors = 0.2*(10**value)

        # Aperture correction
        aperture_table = Table.read(psf_dir / f'{filter_name}_peak.txt', format='ascii')
        aperture_table = aperture_table[aperture_table['apD'] == aperture_diameter]

        enclosed_flux = aperture_table['ef'][0]
        print(enclosed_flux)

        fluxFinal = np.array([f*100/enclosed_flux for f in fluxFinal])

        # Wherever the depth is -99, set the error and flux to -99
        errors[depths_here == -99] = -99
        fluxFinal[depths_here == -99] = -99

        # Store fluxes and errors in the data dictionary
        data_dict[f'flux_{filter_name}'] = fluxFinal
        data_dict[f'err_{filter_name}'] = errors

    # Convert the data dictionary to an Astropy table
    flux_table = Table(data_dict)

    # Return the table
    return flux_table

    

#! REBELS sources
t = ascii.read(Path.cwd().parent.parent / 'data' / 'mosaic' / 'REBELS.csv', format='csv')
t = t[t['RA'] > 140]
ra = t['RA']
dec = t['Dec']
z = t['Redshift (z)']
ID = t['Object Name']

table_euclid = measure_euclid_fluxes(ra, dec, ['VIS', 'Y', 'J', 'H'])
table_ground = measure_ground_fluxes(ra, dec, ['Y', 'J', 'H', 'K'])

# Remove -99. from these astropy tables
for filter_name in ['VIS', 'Y', 'J', 'H']:
    table_euclid = table_euclid[table_euclid[f'flux_{filter_name}'] != -99]
    table_euclid = table_euclid[table_euclid[f'err_{filter_name}'] != -99]
for filter_name in ['Y', 'J', 'H', 'K']:
    table_ground = table_ground[table_ground[f'flux_{filter_name}'] != -99]
    table_ground = table_ground[table_ground[f'err_{filter_name}'] != -99]

# Loop through euclid fluxes and plot with errors
for filter_name in ['VIS', 'Y', 'J', 'H']:
    plt.errorbar(table_euclid[f'flux_{filter_name}'], table_euclid[f'err_{filter_name}'], fmt='o', label=filter_name)

# And ground data
for filter_name in ['Y', 'J', 'H', 'K']:
    plt.errorbar(table_ground[f'flux_{filter_name}'], table_ground[f'err_{filter_name}'], fmt='o', label=filter_name)




