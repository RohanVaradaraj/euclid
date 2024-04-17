#!/usr/bin/env python3

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
from astropy.table import Table, vstack
from astropy.wcs import WCS
import sys
sys.path.append(str(Path.cwd().parent))
from astropy.nddata import Cutout2D
from matplotlib.colors import LinearSegmentedColormap

from cutouts.cutout_codes import isCoordInSurveyFootprints

plt.rcParams['axes.linewidth'] = 2.5
plt.rcParams.update({'font.size': 15})
plt.rcParams['figure.dpi'] = 100

plot_dir = Path.cwd().parent.parent / 'plots'


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



def filterCentreAndWidth(filter_name: str, instrument: str) -> tuple[float, float]:
    """
    Get the central wavelength and width of a filter.

    Parameters
    ----------
    filter_names : list
        List of filter names.
    instrument : str
        Name of the filter instrument.

    Returns
    -------

    tuple[float, float]
        Tuple of central wavelength and width.
    """

    if instrument.lower() == 'euclid':

        t = ascii.read(Path.cwd() / 'euclid_filters.txt')

        centre = t[t['filter'] == filter_name]['centre'][0]
        upper_edge = t[t['filter'] == filter_name]['upper_edge'][0]
        lower_edge = t[t['filter'] == filter_name]['lower_edge'][0]
        width = upper_edge - lower_edge

    if instrument.lower() == 'vista':
            
        t = ascii.read(Path.cwd() / 'vista_filters.txt')

        centre = t[t['filter'] == filter_name]['centre'][0]
        width = t[t['filter'] == filter_name]['width'][0]

    if instrument.lower() == 'hsc':
                
            t = ascii.read(Path.cwd() / 'hsc_filters.txt')
    
            centre = t[t['filter'] == filter_name]['centre'][0]
            width = t[t['filter'] == filter_name]['width'][0]

    return centre, width



def measure_euclid_fluxes(ra: np.ndarray, dec: np.ndarray, filter_names: list, aperture_diameter: float = 1.2) -> tuple[np.ndarray, np.ndarray]:
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
        Diameter of the aperture, in arcseconds. Default is 1.2.

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

            pix_scale = np.abs(header['CD1_1']) * 3600
            print(pix_scale)

        # Open the depth table
        depth_table = Table.read(depth_dir / f'{filter_name}_{aperture_diameter}as_gridDepths_300_200.fits')

        # Loop through the RA and DEC
        assert len(ra) == len(dec), "RA and DEC arrays must have the same length."

        # Convert to pixel coordinates
        w = WCS(header)
        x, y = w.all_world2pix(ra, dec, 0)

        # Find the depth at this x and y
        depths_here = grid_depths(depth_table, x, y)
        print(depths_here)

        # x, y must be array-like
        x = np.array(x)
        y = np.array(y)

        # Use sep to measure flux of object
        flux_counts, _, _ = sep.sum_circle(image, x, y, (aperture_diameter/2.) / pix_scale)
        print((aperture_diameter/2.) * pix_scale)
        print(flux_counts)

        # # Get a quick cutout and plot
        # if filter_name == 'J':
        #     cutout = Cutout2D(image, (x[0], y[0]), (50, 50), wcs=w)
        #     plt.imshow(cutout.data, origin='lower')
        #     # Add a circle corresponding to aperture
        #     circle = plt.Circle((25, 25), (aperture_diameter/2.) / pix_scale, color='red', fill=False)
        #     plt.gca().add_artist(circle)
        #     plt.show()

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

        fluxFinal = np.array([f/enclosed_flux for f in fluxFinal])

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

def measure_ground_fluxes(ra: np.ndarray, dec: np.ndarray, filter_names: list, aperture_diameter: float = 2.0) -> tuple[np.ndarray, np.ndarray]:
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
                pix_scale = np.abs(header['CD1_1']) * 3600

        else:
            with fits.open(directory / image_name, memmap=True) as hdu:
                image = hdu[0].data.byteswap().newbyteorder()
                header = hdu[0].header
                pix_scale = np.abs(header['CD1_1']) * 3600

        # Open the depth table
        depth_table = Table.read(depth_dir / f'{filter_name}_{aperture_diameter}as_gridDepths_300_200.fits')

        # Loop through the RA and DEC
        assert len(ra) == len(dec), "RA and DEC arrays must have the same length."

        # Convert to pixel coordinates
        w = WCS(header)
        x, y = w.all_world2pix(ra, dec, 0)

        # Find the depth at this x and y
        depths_here = grid_depths(depth_table, x, y)

        # x, y must be array-like
        x = np.array(x)
        y = np.array(y)

        # Use sep to measure flux of object
        flux_counts, _, _ = sep.sum_circle(image, x, y, (aperture_diameter/2.) / pix_scale)

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

        fluxFinal = np.array([f/enclosed_flux for f in fluxFinal])

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

    


if __name__ == '__main__':

    #! REBELS sources
    t = ascii.read(Path.cwd().parent.parent / 'data' / 'mosaic' / 'REBELS.csv', format='csv')
    t = t[t['RA'] > 148]
    ra = t['RA']
    dec = t['Dec']
    ID = t['Object Name']
    ID = [name.split('>')[1].split('<')[0] for name in ID]

    #! Stars
    # stars_dir = Path.cwd().parent.parent / 'data' / 'ref_catalogues' / 'stars'
    # stars = ascii.read(stars_dir / 'Y_vista_euclid_coords.ascii')
    # ra = stars['RA_euclid']
    # dec = stars['DEC_euclid']

    #print(isCoordInSurveyFootprints(ra, dec))
    #exit()

    z = t['Redshift (z)']
    ID = t['Object Name']

    # #! ############# Test aperture sizes #################
    # apsizes = [0.8, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6]
    # tables = []
    # for apsize in [0.8, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6]:
    #     table_euclid = measure_euclid_fluxes(ra, dec, ['VIS', 'Y', 'J', 'H'], aperture_diameter=apsize)
    #     table_euclid['ap_size'] = apsize
    #     tables.append(table_euclid)

    # # Stack tables
    # table_euclid = vstack(tables)
    # #!##################################################
    table_euclid = measure_euclid_fluxes(ra, dec, ['VIS', 'Y', 'J', 'H'], aperture_diameter=1.6)
    table_ground = measure_ground_fluxes(ra, dec, ['HSC-G_DR3', 'HSC-R_DR3', 'HSC-I_DR3', 'HSC-Z_DR3', 'HSC-Y_DR3', 'Y', 'J', 'H', 'Ks'], aperture_diameter=1.8)

    # Save these tables
    table_euclid.write(Path.cwd().parent.parent / 'data' / 'ref_catalogues' / 'fluxes' / 'REBELS_euclid_fluxes.fits', format='fits', overwrite=True)
    table_ground.write(Path.cwd().parent.parent / 'data' / 'ref_catalogues' / 'fluxes' / 'REBELS_ground_fluxes.fits', format='fits', overwrite=True)

    # Open these tables
    #table_euclid = Table.read(Path.cwd().parent.parent / 'data' / 'ref_catalogues' / 'fluxes' / 'star_fluxes.fits')
    #table_ground = Table.read(Path.cwd().parent.parent / 'data' / 'ref_catalogues' / 'fluxes' / 'star_ground_fluxes.fits')
    #table_euclid = Table.read(Path.cwd().parent.parent / 'data' / 'ref_catalogues' / 'fluxes' / 'REBELS_euclid_fluxes.fits')
    #table_ground = Table.read(Path.cwd().parent.parent / 'data' / 'ref_catalogues' / 'fluxes' / 'REBELS_ground_fluxes.fits')
    #table_euclid = Table.read(Path.cwd().parent.parent / 'data' / 'ref_catalogues' / 'fluxes' / 'apSize_test_euclid.fits')
    #table_ground = Table.read(Path.cwd().parent.parent / 'data' / 'ref_catalogues' / 'fluxes' / 'apSize_test_ground.fits')

    # cmap = plt.get_cmap('viridis')  # Choose your desired colormap
    # norm = plt.Normalize(0, len(table_euclid['flux_H']) - 1)  # Normalize values of i
    
    # Remove -99. from these astropy tables
    for filter_name in ['VIS', 'Y', 'J', 'H']:
        table_euclid = table_euclid[table_euclid[f'flux_{filter_name}'] != -99]
        table_euclid = table_euclid[table_euclid[f'err_{filter_name}'] != -99]
    for filter_name in ['HSC-G_DR3', 'HSC-R_DR3', 'HSC-I_DR3', 'HSC-Z_DR3', 'HSC-Y_DR3', 'Y', 'J', 'H', 'Ks']:
        table_ground = table_ground[table_ground[f'flux_{filter_name}'] != -99]
        table_ground = table_ground[table_ground[f'err_{filter_name}'] != -99]

    plt.figure(figsize=(10, 6))
    # Loop through objects
    for i in range(len(table_euclid['flux_H'])):

        print(ID[i])

        # Loop through euclid fluxes and plot with errors
        for filter_name in ['VIS', 'Y', 'J', 'H']:

            centre, width = filterCentreAndWidth(filter_name, 'euclid')

            # color = cmap(norm(i))
            plt.errorbar(centre+i*0.01, table_euclid[f'flux_{filter_name}'][i], yerr=table_euclid[f'err_{filter_name}'][i], xerr=width/2, 
                         fmt='o', color='black', alpha=0.7)
        #plt.errorbar([], [], yerr=0, xerr=0, fmt='o', label=f'apsize={apsizes[i]} as', alpha=0.7, color=color)

        # And ground data
        for filter_name in ['HSC-G_DR3', 'HSC-R_DR3', 'HSC-I_DR3', 'HSC-Z_DR3', 'HSC-Y_DR3', 'Y', 'J', 'H', 'Ks']:

            if filter_name[0:3] == 'HSC':
                centre, width = filterCentreAndWidth(filter_name, 'HSC')
            else:
                centre, width = filterCentreAndWidth(filter_name, 'VISTA')
            plt.errorbar(centre, table_ground[f'flux_{filter_name}'][0], yerr=table_ground[f'err_{filter_name}'][0], xerr=width/2, 
                        fmt='o', color='red', alpha=0.7)
        
        plt.xlabel(r'$\lambda \ (\mu \mathrm{m})$')
        plt.ylabel(r'flux (erg s$^{-1}$ cm$^{-2}$ Hz$^{-1}$)')
        plt.legend()

        #plt.yscale('log')
        plt.show()
        #plt.savefig(Path.cwd().parent.parent / 'plots' / 'seds' / 'REBELS_fluxes.png')


