"""
cutout_codes.py

Functions and classes used to make cutouts

Created: Thursday 21st March 2024.
"""

from shapely.geometry import Point, Polygon
import numpy as np
from pathlib import Path
from astropy.io import ascii
from typing import Optional
from astropy.io import fits
from astropy.nddata.utils import Cutout2D
from astropy.wcs import WCS
import glob
from astropy.coordinates import SkyCoord
import matplotlib.pyplot as plt
from scipy.ndimage import rotate
from astropy.table import Table
import astropy.units as u

import warnings
from astropy.utils.exceptions import AstropyWarning
warnings.simplefilter('ignore', category=AstropyWarning)

# plt.rcParams['axes.linewidth'] = 2.5
# plt.rcParams.update({'font.size': 15})
# plt.rcParams['figure.dpi'] = 100

ground_dir = Path.home().parent.parent / 'vardy' / 'vardygroupshare' / 'data'
euclid_dir = Path.home() / 'euclid'
cweb_dir = Path.home().parent.parent / 'vardy' / 'vardygroupshare' / 'data' / 'CWEB'
primer_dir = Path.home() / 'JWST'
plot_dir = Path.cwd().parent.parent / 'plots' / 'cutouts'
refcat_dir = Path.cwd().parents[1] / 'data' / 'ref_catalogues'


def readAllCOSMOSGalaxies() -> Table:

    """
    Read in all crossmatched objects found from NED from the text file to an astropy table.
    """

    # Read in the file
    cosmos_file = refcat_dir / 'all_COSMOS_highz.txt'
    cosmos = np.loadtxt(cosmos_file, dtype=str, delimiter='|', skiprows=21)

    # Convert to astropy table
    cosmos_colnames = cosmos[0]
    cosmos_data = cosmos[1:]
    cosmos_table = Table(cosmos_data, names=cosmos_colnames)

    cosmos_table.write(refcat_dir / 'all_COSMOS_highz.fits', format='fits', overwrite=True)




def findPlotLimits(data: np.ndarray) -> tuple:

    mean = np.mean(data)
    std_dev = np.std(data)

    # Sigma clip the data
    data = data[(data < mean + 3 * std_dev)]

    # Recalculate mean and std_dev
    mean = np.mean(data)
    std_dev = np.std(data)

    lower = mean - 2 * std_dev
    upper = mean + 5 * std_dev

    return lower, upper



def isCoordInSurveyFootprints(ra: np.ndarray, dec: np.ndarray) -> np.ndarray:
    """
    Check if a coordinate lies in the footprints of:
        1) Euclid On-the-fly COSMOS release
        2) COSMOS-Web
        3 )PRIMER

    Parameters
    ----------
    ra : np.array
        Right ascension of the coordinate in degrees.
    dec : np.array
        Declination of the coordinate in degrees.


    Returns
    -------
    np.ndarray
        Array with rows corresponding to the input coordinates and columns corresponding to the footprints.
        The value is '0' if it is not in the footprint.
        If it is in the footprint then: 
            for Euclid it takes a numeric string value of the tile label I have assigned (1-22).
            for CWEB it takes a numeric string value of the tile label Nathan has assigned (4A, 5A, 5B, 6A, etc.)
            for PRIMER it takes the value '1'.

    """

    # Check if ra and dec are arrays
    ra = np.array(ra)
    dec = np.array(dec)

    # Reshape ra and dec to make sure they are at least one-dimensional arrays
    ra = np.atleast_1d(ra)
    dec = np.atleast_1d(dec)

    assert ra.shape == dec.shape, "ra and dec arrays must have the same shape"

    points = [Point(x, y) for x, y in zip(ra, dec)]

    # Load the footprints
    euclid_footprint = np.load(Path.cwd().parent.parent / 'data' / 'mosaic' / 'euclid_footprint.npy')
    primer_footprint = np.load(Path.cwd().parent.parent / 'data' / 'mosaic' / 'primer_footprint.npy')
    cweb_footprint = np.load(Path.cwd().parent.parent / 'data' / 'mosaic' / 'cweb_footprint.npy')
    hubble_footprint = np.load(Path.cwd().parent.parent / 'data' / 'mosaic' / 'hubble_footprint.npy')

    # Load the labels for CWEB and Euclid
    euclid_labels = np.load(Path.cwd().parent.parent / 'data' / 'mosaic' / 'euclid_labels.npy')
    cweb_labels = np.load(Path.cwd().parent.parent / 'data' / 'mosaic' / 'cweb_labels.npy')

    # Create shapely polygons
    euclid_polygons = [Polygon(footprint) for footprint in euclid_footprint]
    cweb_polygons = [Polygon(footprint) for footprint in cweb_footprint]
    primer_polygons = Polygon(primer_footprint)
    hubble_polygons = Polygon(hubble_footprint)

    points_in_euclid = []
    points_in_cweb = []
    points_in_primer = []
    points_in_hubble = []

    for point in points:

        # 0 if not in footprint, else the tile label (1 for HST and PRIMER since it is one tile each)
        in_euclid = '0'
        in_cweb = '0'
        in_primer = '0'
        in_hubble = '0'

        #! Check if the point lies in any of the Euclid polygons
        for i, polygon in enumerate(euclid_polygons):
            if polygon.contains(point):
                in_euclid = euclid_labels[i]
                break

        #! Check if the point lies in any of the COSMOS-Web polygons
        for i, polygon in enumerate(cweb_polygons):
            if polygon.contains(point):
                in_cweb = cweb_labels[i]
                break

        #! Check if point lies in PRIMER footprint
        if primer_polygons.contains(point):
            in_primer = '1'

        #! Check if point lies in Hubble footprint
        if hubble_polygons.contains(point):
            in_hubble = '1'

        # Append results of search to lists
        points_in_euclid.append(in_euclid)
        points_in_cweb.append(in_cweb)
        points_in_primer.append(in_primer)
        points_in_hubble.append(in_hubble)
    
    # Convert lists to numpy arrays
    euclid_mask = np.array(points_in_euclid)
    cweb_mask = np.array(points_in_cweb)
    primer_mask = np.array(points_in_primer)
    hubble_mask = np.array(points_in_hubble)

    return np.column_stack((euclid_mask, cweb_mask, primer_mask, hubble_mask))


def isCoordInCWEB(ra: np.ndarray, dec: np.ndarray) -> np.ndarray:
    """
    Check if a coordinate lies in the footprints of CWEB.

    If it doesn't, return zero for the coord, else return the mosaic string.

    Parameters
    ----------
    ra : np.array
        Right ascension of the coordinate in degrees.
    dec : np.array
        Declination of the coordinate in degrees.
    
    Returns
    -------
    np.ndarray
        Array with rows corresponding to the input coordinates and columns corresponding to the footprints.
        The value is '0' if it is not in the footprint.
        If it is in the footprint then it takes a numeric string value of the tile label Nathan has assigned (4A, 5A, 5B, 6A, etc.)

    """

    # Check if ra and dec are arrays
    ra = np.array(ra)
    dec = np.array(dec)

    # Reshape ra and dec to make sure they are at least one-dimensional arrays
    ra = np.atleast_1d(ra)
    dec = np.atleast_1d(dec)

    assert ra.shape == dec.shape, "ra and dec arrays must have the same shape"

    points = [Point(x, y) for x, y in zip(ra, dec)]

    mosaics = ['0A', '0B', '1A', '1B', '2A', '2B', '3A', '3B', '4A', '4B', '5A', '5B', '6A', '6B', '7A', '7B']

    polygons = []

    for mosaic in mosaics:
        footprint = np.load(Path.cwd().parent.parent / 'data' / 'mosaic' / f'CWEB_footprint_{mosaic}.npy')

        # Create shapely polygons
        cweb_polygon = Polygon(footprint)
        polygons.append(cweb_polygon)

    points_in_cweb = []
    for point in points:

        # 0 if not in footprint, else the tile label
        in_cweb = '0'

        #! Check if the point lies in any of the COSMOS-Web polygons
        #! If so, find all of them
        for i, polygon in enumerate(polygons):
            if polygon.contains(point):
                in_cweb = mosaics[i]
                break

        # Append results of search to lists
        points_in_cweb.append(in_cweb)

    # Convert lists to numpy arrays
    points_in_cweb = np.array(points_in_cweb)

    return points_in_cweb




def Cutout(ra: float, dec:float, contained_in: Optional[np.array] = None, size: float = 10.0, 
           save_cutout: bool = True, save_dir: Path = Path.cwd().parent.parent / 'data' / 'cutouts',
           plot_title: Optional[str] = None,
           add_centre_lines: Optional[bool] = False) -> None:

    """
    Create cutouts from Euclid, ground-based and JWST imaging.

    Parameters
    ----------

    ra : float
        Right ascension of the cutout center in degrees.
    dec : float
        Declination of the cutout center in degrees.
    contained_in : np.array (str, str, str)
        Array with three strings. The strings are the values returned by isCoordInSurveyFootprints.
        The first string is the Euclid tile label (0 if absent), the second is the CWEB tile label (0 if absent) and the third is '1' if in PRIMER.
    size : float, optional
        Size of the cutout in arcseconds. The default is 10.0
    save_cutout : bool, optional
        Whether to save the cutout. The default is False.
    save_dir : Path, optional
        Directory to save the cutout. The default is Path.cwd().parent.parent / 'data' / 'cutouts'.
    plot_title: str, optional
        Can pass a title to the stamp plot. The default is None.
    add_centre_lines: bool, optional
        Whether to add lines to the plot showing the centre of the cutout. The default is False.

    
    Returns
    -------
    None

    """

    # If contained_in is None, do the check for coordinates in the footprints
    if contained_in is None:
        contained_in = isCoordInSurveyFootprints(ra, dec)
        print('contained_in array: ', contained_in)

    # Check if all values of contained_in are '0'.
    if np.all(contained_in[0] == '0'):
        print("Cutout not in any survey footprint")
        #return None
    
    if contained_in[0][0] == '0':
        print('Cutout not in Euclid footprint')
        #return None

    # Convert RA, DEC to skycoords.
    c = SkyCoord(ra, dec, unit='deg')
    
    #! First get the ground-based cutouts (I, Y, J, H)
    vista_dir = ground_dir / 'COSMOS'

    #### I ####
    with fits.open(vista_dir / 'HSC-I_DR3.fits') as hdu_I:
            
            data_I = hdu_I[0].data
            hdr_I = hdu_I[0].header
            pix_scale = np.abs(hdr_I['CD1_1']) * 3600
            wcs_I = WCS(hdr_I)
    
            cutout_grI = Cutout2D(data_I, c, size=size/pix_scale, wcs=wcs_I)

    #### Y ####
    with fits.open(vista_dir / 'UVISTA_Y_DR6.fits') as hdu_Y:
    #with fits.open(vista_dir / 'UVISTA_Y_dr5_rc1.fits') as hdu_Y:

        data_Y = hdu_Y[0].data
        hdr_Y = hdu_Y[0].header
        pix_scale = np.abs(hdr_Y['CD1_1']) * 3600
        wcs_Y = WCS(hdr_Y)

        cutout_grY = Cutout2D(data_Y, c, size=size/pix_scale, wcs=wcs_Y)

    #### J ####
    with fits.open(vista_dir / 'UVISTA_J_DR6.fits') as hdu_J:
    #with fits.open(vista_dir / 'UVISTA_J_dr5_rc1.fits') as hdu_J:

        data_J = hdu_J[0].data
        hdr_J = hdu_J[0].header
        pix_scale = np.abs(hdr_J['CD1_1']) * 3600
        wcs_J = WCS(hdr_J)

        cutout_grJ = Cutout2D(data_J, c, size=size/pix_scale, wcs=wcs_J)

    #### H ####
    with fits.open(vista_dir / 'UVISTA_H_DR6.fits') as hdu_H:
    #with fits.open(vista_dir / 'UVISTA_H_dr5_rc1.fits') as hdu_H:
            
        data_H = hdu_H[0].data
        hdr_H = hdu_H[0].header
        pix_scale = np.abs(hdr_H['CD1_1']) * 3600
        wcs_H = WCS(hdr_H)

        cutout_grH = Cutout2D(data_H, c, size=size/pix_scale, wcs=wcs_H)

    #### Ks ####
    # with fits.open(vista_dir / 'UVISTA_K_DR6_cropped.fits') as hdu_K:
            
    #     data_K = hdu_K[0].data
    #     hdr_K = hdu_K[0].header
    #     pix_scale = np.abs(hdr_K['CD1_1']) * 3600
    #     wcs_K = WCS(hdr_K)

    #     cutout_grK = Cutout2D(data_K, c, size=size/pix_scale, wcs=wcs_K)

    
    #! Next get the Euclid cutouts
    euclid_tile = contained_in[0][0]

    if euclid_tile != '0':

        #### VIS ####
        tile_VIS = glob.glob(str(euclid_dir / 'VIS' / 'COSMOS' / f'*BGSUB*_{euclid_tile}.fits*'))[0]
        with fits.open(tile_VIS) as hdu_VIS:

            data_VIS = hdu_VIS[0].data
            hdr_VIS = hdu_VIS[0].header
            pix_scale = np.abs(hdr_VIS['CD1_1']) * 3600
            wcs_VIS = WCS(hdr_VIS)

            cutout_euVIS = Cutout2D(data_VIS, c, size=size/pix_scale, wcs=wcs_VIS)

        #### Y ####
        tile_Y = glob.glob(str(euclid_dir / 'Y' / 'COSMOS' / f'*BGSUB*_{euclid_tile}.fits*'))[0]
        with fits.open(tile_Y) as hdu_Y:

            data_Y = hdu_Y[0].data
            hdr_Y = hdu_Y[0].header
            pix_scale = np.abs(hdr_Y['CD1_1']) * 3600
            wcs_Y = WCS(hdr_Y)

            cutout_euY = Cutout2D(data_Y, c, size=size/pix_scale, wcs=wcs_Y)

        #### J ####
        tile_J = glob.glob(str(euclid_dir / 'J' / 'COSMOS' / f'*BGSUB*_{euclid_tile}.fits*'))[0]
        with fits.open(tile_J) as hdu_J:

            data_J = hdu_J[0].data
            hdr_J = hdu_J[0].header
            pix_scale = np.abs(hdr_J['CD1_1']) * 3600
            wcs_J = WCS(hdr_J)

            cutout_euJ = Cutout2D(data_J, c, size=size/pix_scale, wcs=wcs_J)

        #### H ####
        tile_H = glob.glob(str(euclid_dir / 'H' / 'COSMOS' / f'*BGSUB*_{euclid_tile}.fits*'))[0]
        with fits.open(tile_H) as hdu_H:

            data_H = hdu_H[0].data
            hdr_H = hdu_H[0].header
            pix_scale = np.abs(hdr_H['CD1_1']) * 3600
            wcs_H = WCS(hdr_H)

            cutout_euH = Cutout2D(data_H, c, size=size/pix_scale, wcs=wcs_H)

        # If all of the Euclid cutouts are zero, then the cutout is empty so return none
        # if np.all(cutout_euVIS.data == 0.) and np.all(cutout_euY.data == 0.) and np.all(cutout_euJ.data == 0.) and np.all(cutout_euH.data == 0.):
        #     print("Euclid cutouts are empty")
        #     return None#

    # If all the Euclid cutouts are zero, return an empty cutout
    if euclid_tile == '0':

        cutout_euVIS = Cutout2D(np.zeros((size, size)), c, size=size, wcs=WCS(naxis=2))
        cutout_euY = Cutout2D(np.zeros((size, size)), c, size=size, wcs=WCS(naxis=2))
        cutout_euJ = Cutout2D(np.zeros((size, size)), c, size=size, wcs=WCS(naxis=2))
        cutout_euH = Cutout2D(np.zeros((size, size)), c, size=size, wcs=WCS(naxis=2))


    #! Next check for and get the Hubble cutouts
    hubble_tile = contained_in[0][3]

    if hubble_tile == '1':

        dash = Path.home().parent.parent / 'hoy' / 'DASH' / 'hlsp_3d-dash_hst_wfc3_combined-cosmos_f160w_v1.0_drz-sci.fits'

        with fits.open(dash) as hdu_dash:

            data_hubble = hdu_dash[0].data
            hdr_hubble = hdu_dash[0].header
            pix_scale = np.abs(hdr_hubble['CDELT1']) * 3600
            wcs_hubble = WCS(hdr_hubble)

            cutout_hubble = Cutout2D(data_hubble, c, size=size/pix_scale, wcs=wcs_hubble)


    #! Next check for and get the JWST cutouts
    
    #### CWEB ####
    #cweb_tile = contained_in[0][1]
    cweb_tile = isCoordInCWEB(ra, dec)[0] # Updated version for full CWEB

    if cweb_tile != '0':

        tile_f115w = glob.glob(str(cweb_dir / f'mosaic_{cweb_tile}' / f'CWEB-F115W-{cweb_tile}_i2dnobg_small.fits'))[0]
        tile_f150w = glob.glob(str(cweb_dir / f'mosaic_{cweb_tile}' / f'CWEB-F150W-{cweb_tile}_i2dnobg_small.fits'))[0]
        tile_f277w = glob.glob(str(cweb_dir / f'mosaic_{cweb_tile}' / f'CWEB-F277W-{cweb_tile}_i2dnobg_small.fits'))[0]
        tile_f444w = glob.glob(str(cweb_dir / f'mosaic_{cweb_tile}' / f'CWEB-F444W-{cweb_tile}_i2dnobg_small.fits'))[0]

        #### F115W ####
        with fits.open(tile_f115w) as hdu_CWEB:

            data_CWEB = hdu_CWEB[1].data
            hdr_CWEB = hdu_CWEB[1].header
            pix_scale = np.abs(hdr_CWEB['CDELT1']) * 3600
            wcs_CWEB = WCS(hdr_CWEB)

            cutout_f115w_cweb = Cutout2D(data_CWEB, c, size=size/pix_scale, wcs=wcs_CWEB)

        #### F150W ####
        with fits.open(tile_f150w) as hdu_CWEB:

            data_CWEB = hdu_CWEB[1].data
            hdr_CWEB = hdu_CWEB[1].header
            pix_scale = np.abs(hdr_CWEB['CDELT1']) * 3600
            wcs_CWEB = WCS(hdr_CWEB)

            cutout_f150w_cweb = Cutout2D(data_CWEB, c, size=size/pix_scale, wcs=wcs_CWEB)

        #### F277W ####
        with fits.open(tile_f277w) as hdu_CWEB:

            data_CWEB = hdu_CWEB[1].data
            hdr_CWEB = hdu_CWEB[1].header
            pix_scale = np.abs(hdr_CWEB['CDELT1']) * 3600
            wcs_CWEB = WCS(hdr_CWEB)

            cutout_f277w_cweb = Cutout2D(data_CWEB, c, size=size/pix_scale, wcs=wcs_CWEB)
        
        #### F444W ####
        with fits.open(tile_f444w) as hdu_CWEB:
            
            data_CWEB = hdu_CWEB[1].data
            hdr_CWEB = hdu_CWEB[1].header
            pix_scale = np.abs(hdr_CWEB['CDELT1']) * 3600
            wcs_CWEB = WCS(hdr_CWEB)
            pa = hdr_CWEB['PA_APER']

            cutout_f444w_cweb = Cutout2D(data_CWEB, c, size=size/pix_scale, wcs=wcs_CWEB)


    #### PRIMER ####
    primer_tile = contained_in[0][2]

    if primer_tile == '1':

        #### F277W ####
        with fits.open(primer_dir / 'primer_cosmos_nircam_v0.5_f277w_30mas_sci.fits') as hdu_PRIMER:

            data_PRIMER = hdu_PRIMER[0].data
            hdr_PRIMER = hdu_PRIMER[0].header
            pix_scale = np.abs(hdr_PRIMER['CD1_1']) * 3600
            wcs_PRIMER = WCS(hdr_PRIMER)

            cutout_f277w_prim = Cutout2D(data_PRIMER, c, size=size/pix_scale, wcs=wcs_PRIMER)


        #### F444W ####
        with fits.open(primer_dir / 'primer_cosmos_nircam_v0.5_f444w_30mas_sci.fits') as hdu_PRIMER:

            data_PRIMER = hdu_PRIMER[0].data
            hdr_PRIMER = hdu_PRIMER[0].header
            pix_scale = np.abs(hdr_PRIMER['CD1_1']) * 3600
            wcs_PRIMER = WCS(hdr_PRIMER)

            cutout_f444w_prim = Cutout2D(data_PRIMER, c, size=size/pix_scale, wcs=wcs_PRIMER)


    # Convert footprint array into bools: if it is '0' then False, else True
    footprint_bools = [0 if i == '0' else 1 for i in contained_in[0]]
    footprint_bools = np.array(footprint_bools, dtype=bool)

    #! Plots based on different footprint bools

    #* Define a lambda function to plot each cutout, instead of repeating this line every time!
    plot_cutout = lambda ax, data, lims, title: ax.imshow(data, origin='lower', cmap='gist_yarg', vmin=lims[0], vmax=lims[1]) and ax.set_title(title)

    #? There must be a better way of doing the if statements below, but I can't think of it right now.

    # Only in Euclid
    if footprint_bools[0] and not footprint_bools[1] and not footprint_bools[2] and not footprint_bools[3]:

        fig, ax = plt.subplots(2, 4, figsize=(15, 10))

        # Turn off axis labels and ticks
        for axis in ax.flat:
            axis.set_axis_off()
            axis.tick_params(axis='both', which='both', bottom=False, top=False, left=False, right=False, labelbottom=False, labelleft=False)

        # Top row: ground-based cutouts
        for i, (data, title) in enumerate(zip([cutout_grI.data, cutout_grY.data, cutout_grJ.data, cutout_grH.data], ['HSC-I', 'VISTA-Y', 'VISTA-J', 'VISTA-H'])):
            lims = findPlotLimits(data)
            plot_cutout(ax[0, i], data, lims, title)

        # Bottom row: Euclid cutouts
        for i, (data, title) in enumerate(zip([cutout_euVIS.data, cutout_euY.data, cutout_euJ.data, cutout_euH.data], [r'$I_{E}$', r'$Y_{E}$', r'$J_{E}$', r'$H_{E}$'])):
            lims = findPlotLimits(data)
            plot_cutout(ax[1, i], data, lims, title)

    # In Euclid and CWEB
    elif footprint_bools[0] and cweb_tile != '0' and not footprint_bools[2] and not footprint_bools[3]: 

        fig, ax = plt.subplots(2, 6, figsize=(15, 10))

        # Turn off axis labels and ticks
        for axis in ax.flat:
            axis.set_axis_off()
            axis.tick_params(axis='both', which='both', bottom=False, top=False, left=False, right=False, labelbottom=False, labelleft=False)

        # Top row: ground-based cutouts
        for i, (data, title) in enumerate(zip([cutout_grI.data, cutout_grY.data, cutout_grJ.data, cutout_grH.data], ['HSC-I', 'VISTA-Y', 'VISTA-J', 'VISTA-H'])):
            lims = findPlotLimits(data)
            plot_cutout(ax[0, i], data, lims, title)


        # Bottom row: Euclid cutouts
        for i, (data, title) in enumerate(zip([cutout_euVIS.data, cutout_euY.data, cutout_euJ.data, cutout_euH.data], [r'$I_{E}$', r'$Y_{E}$', r'$J_{E}$', r'$H_{E}$'])):
            lims = findPlotLimits(data)
            plot_cutout(ax[1, i], data, lims, title)


        # CWEB cutouts in final column. Reproject to Euclid WCS
        for i, (data, title) in enumerate(zip([cutout_f115w_cweb.data, cutout_f150w_cweb.data, cutout_f277w_cweb.data, cutout_f444w_cweb.data], ['F115W', 'F150W', 'F277W', 'F444W'])):
            data = rotate(data, -pa, reshape=False)
            lims = findPlotLimits(data)
            if i == 0:
                plot_cutout(ax[0, 4], data, lims, title)
            if i == 1:
                plot_cutout(ax[0, 5], data, lims, title)
            if i == 2:
                plot_cutout(ax[1, 4], data, lims, title)
            if i == 3:
                plot_cutout(ax[1, 5], data, lims, title)

    # In Euclid and CWEB and PRIMER
    elif footprint_bools[0] and not footprint_bools[1] and footprint_bools[2] and not footprint_bools[3]:

        fig, ax = plt.subplots(2, 5, figsize=(15, 10))

        # Turn off axis labels and ticks
        for axis in ax.flat:
            axis.set_axis_off()
            axis.tick_params(axis='both', which='both', bottom=False, top=False, left=False, right=False, labelbottom=False, labelleft=False)

        # Top row: ground-based cutouts
        for i, (data, title) in enumerate(zip([cutout_grI.data, cutout_grY.data, cutout_grJ.data, cutout_grH.data], ['HSC-I', 'VISTA-Y', 'VISTA-J', 'VISTA-H'])):
            lims = findPlotLimits(data)
            plot_cutout(ax[0, i], data, lims, title)


        # Bottom row: Euclid cutouts
        for i, (data, title) in enumerate(zip([cutout_euVIS.data, cutout_euY.data, cutout_euJ.data, cutout_euH.data], [r'$I_{E}$', r'$Y_{E}$', r'$J_{E}$', r'$H_{E}$'])):
            lims = findPlotLimits(data)
            plot_cutout(ax[1, i], data, lims, title)

        # PRIMER cutouts in final column
        for i, (data, title) in enumerate(zip([cutout_f277w_prim.data, cutout_f444w_prim.data], ['F277W', 'F444W'])):
            lims = findPlotLimits(data)
            plot_cutout(ax[0, 4], data, lims, title) if i == 0 else plot_cutout(ax[1, 4], data, lims, title)

    # In Euclid and Hubble, not JWST
    elif footprint_bools[0] and footprint_bools[3] and not footprint_bools[1] and not footprint_bools[2]:

        fig, ax = plt.subplots(2, 5, figsize=(15, 10))

        # Turn off axis labels and ticks
        for axis in ax.flat:
            axis.set_axis_off()
            axis.tick_params(axis='both', which='both', bottom=False, top=False, left=False, right=False, labelbottom=False, labelleft=False)

        # Top row: ground-based cutouts
        for i, (data, title) in enumerate(zip([cutout_grI.data, cutout_grY.data, cutout_grJ.data, cutout_grH.data], ['HSC-I', 'VISTA-Y', 'VISTA-J', 'VISTA-H'])):
            lims = findPlotLimits(data)
            plot_cutout(ax[0, i], data, lims, title)

        # Bottom row: Euclid cutouts
        for i, (data, title) in enumerate(zip([cutout_euVIS.data, cutout_euY.data, cutout_euJ.data, cutout_euH.data], [r'$I_{E}$', r'$Y_{E}$', r'$J_{E}$', r'$H_{E}$'])):
            lims = findPlotLimits(data)
            plot_cutout(ax[1, i], data, lims, title)

        # Plot Hubble data in final column
        lims = findPlotLimits(cutout_hubble.data)
        plot_cutout(ax[1, 4], cutout_hubble.data, lims, 'F160W')

        # Could also plot Ks above F160W
        #lims = findPlotLimits(cutout_grK.data)
        #plot_cutout(ax[0, 4], cutout_grK.data, lims, 'VISTA-Ks')

    # In Euclid and Hubble and JWST
    elif footprint_bools[0] and footprint_bools[3] and (footprint_bools[1] or footprint_bools[2]):

        fig, ax = plt.subplots(2, 7, figsize=(15, 10))

        # Turn off axis labels and ticks
        for axis in ax.flat:
            axis.set_axis_off()
            axis.tick_params(axis='both', which='both', bottom=False, top=False, left=False, right=False, labelbottom=False, labelleft=False)

        # Top row: ground-based cutouts
        for i, (data, title) in enumerate(zip([cutout_grI.data, cutout_grY.data, cutout_grJ.data, cutout_grH.data], ['HSC-I', 'VISTA-Y', 'VISTA-J', 'VISTA-H'])):
            lims = findPlotLimits(data)
            plot_cutout(ax[0, i], data, lims, title)

        # Bottom row: Euclid cutouts
        for i, (data, title) in enumerate(zip([cutout_euVIS.data, cutout_euY.data, cutout_euJ.data, cutout_euH.data], [r'$I_{E}$', r'$Y_{E}$', r'$J_{E}$', r'$H_{E}$'])):
            lims = findPlotLimits(data)
            plot_cutout(ax[1, i], data, lims, title)

        # Plot Hubble data in penultimate column
        lims = findPlotLimits(cutout_hubble.data)
        plot_cutout(ax[1, 4], cutout_hubble.data, lims, 'F160W')

        # CWEB cutouts in final column. Reproject to Euclid WCS
        if cweb_tile != '0':
            for i, (data, title) in enumerate(zip([cutout_f115w_cweb.data, cutout_f150w_cweb.data, cutout_f277w_cweb.data, cutout_f444w_cweb.data], ['F115W', 'F150W', 'F277W', 'F444W'])):
                data = rotate(data, -pa, reshape=False)
                lims = findPlotLimits(data)
                if i == 0:
                    plot_cutout(ax[0, 5], data, lims, title)
                if i == 1:
                    plot_cutout(ax[0, 6], data, lims, title)
                if i == 2:
                    plot_cutout(ax[1, 5], data, lims, title)
                if i == 3:
                    plot_cutout(ax[1, 6], data, lims, title)

        if footprint_bools[2]:
            for i, (data, title) in enumerate(zip([cutout_f277w_prim.data, cutout_f444w_prim.data], ['F277W', 'F444W'])):
                lims = findPlotLimits(data)
                #plot_cutout(ax[0, 5], data, lims, title) if i == 0 else plot_cutout(ax[1, 5], data, lims, title)
                plot_cutout(ax[1, 5], data, lims, title) if i == 0 else plot_cutout(ax[1, 6], data, lims, title)

    if add_centre_lines:
        # Draw a plus symbol on each subplot to show the centre of the cutout
        for axis in ax[0, :4]:
            center_x, center_y = cutout_grY.wcs.world_to_pixel(SkyCoord(c.ra.deg, c.dec.deg, unit='deg'))
            axis.plot(center_x, center_y, 'r+', markersize=20)            

        for axis in ax[1, :4]:
            center_x, center_y = cutout_euY.wcs.world_to_pixel(SkyCoord(c.ra.deg, c.dec.deg, unit='deg'))
            axis.plot(center_x, center_y, 'r+', markersize=20)  

        center_x, center_y = cutout_hubble.wcs.world_to_pixel(SkyCoord(c.ra.deg, c.dec.deg, unit='deg'))
        ax[1,4].plot(center_x, center_y, 'r+', markersize=20)

        if len(ax[0, :]) > 5:
            # if in cweb
            if footprint_bools[1] != 0:
                center_x, center_y = cutout_f277w_cweb.wcs.world_to_pixel(SkyCoord(c.ra.deg, c.dec.deg, unit='deg'))
                ax[0,5].plot(center_x, center_y, 'r+', markersize=20)
                ax[1,5].plot(center_x, center_y, 'r+', markersize=20)
            else:
                # Primer
                center_x, center_y = cutout_f277w_prim.wcs.world_to_pixel(SkyCoord(c.ra.deg, c.dec.deg, unit='deg'))
                ax[0,5].plot(center_x, center_y, 'r+', markersize=20)
                ax[1,5].plot(center_x, center_y, 'r+', markersize=20)


    # Set title if we pass a name
    if plot_title is not None:
        plt.suptitle(plot_title)

    plt.tight_layout()
    #plt.savefig(plot_dir / f'LAE_stamps_Euclid_CWEB_6arcsec.pdf')
    #plt.show()

    if save_cutout:
        plot_title = str(plot_title)
        save_name = plot_title.split(',')[0]
        plt.savefig(plot_dir / f'{save_name}.png')

    return fig, ax
    #plt.show()
    #plt.close()
    #return None



if __name__ == '__main__':
    
    #! PRIMER STARS
    t = Table.read(Path.cwd().parents[3] / 'data' / 'psf' / 'COSMOS' / 'catalogues' / 'f410m_stars.fits')

    ra = t['ALPHA_J2000']
    dec = t['DELTA_J2000']
    ID = t['NUMBER']

    #! REBELS sources
    # t = ascii.read(Path.cwd().parent.parent / 'data' / 'mosaic' / 'REBELS.csv', format='csv')
    # t = t[t['RA'] > 148]

    # # # Skip first few
    # # #t = t[:]

    # ra = t['RA']
    # dec = t['Dec']
    # z = t['Redshift (z)']
    # ID = t['Object Name']
    # ID = [name.split('>')[1].split('<')[0] for name in ID]

    #! Strong lens
    #ra = [150.00280406167596]
    #dec = [2.2002751856804608]

    #! Nathan's z=3 sources
    # nathan_dir = Path.home().parent.parent / 'vardy' / 'vardygroupshare' / 'HSC_SSP_DR3' / 'ref_catalogues' / 'nathan'
    # cat = Table.read(nathan_dir / 'Z3_FinalSample.fits')
    # cat.sort('MUV', reverse=False)
    # cat = cat[cat['RA'] > 148]
    # ra = cat['RA']
    # dec = cat['DEC']
    # muv = cat['MUV']
    # z = cat['ZPhot']
    # ID = cat['UID']

    #! DEVILS sources
    # devils_dir = Path.home() / 'DEVILS' / 'dr1cats' / 'data' / 'catalogues'
    # cat = Table.read(devils_dir / 'D10VisualMorphology.csv', format='csv')
    # print(cat.colnames)
    # cat = cat[(cat['RAcen'] > 148) & (cat['DECcen'] > 1.8)]
    # # Drop 'NA' strings from zBest column
    # cat = cat[cat['zBest'] != 'NA']
    # # Convert zBest to float
    # cat['zBest'] = cat['zBest'].astype(float)
    # # Remove zBest values of -99
    # cat = cat[cat['zBest'] > 0.]

    # cat = cat[cat['FIRST_CLASS'] != 'NA']
    # cat.sort('zBest')
    # ra = cat['RAcen']
    # dec = cat['DECcen']

    #! Harikane z=12-16 sources
    # hd1 = '10:01:51.31 02:32:50.0'
    # hd2 = '02:18:52.44 -05:08:36.1'

    # # Use skycoord to convert these coordinates to degrees
    # c1 = SkyCoord(hd1, unit=(u.hourangle, u.deg))
    # c2 = SkyCoord(hd2, unit=(u.hourangle, u.deg))

    # ra = [c1.ra.deg, c2.ra.deg]
    # dec = [c1.dec.deg, c2.dec.deg]

    # ID = ['HD1', 'HD2']

    #! Rebecca's z>8.5 sources
    # uvista_1212 = '10:02:31.81 02:31:17.10'
    # uvista_237 = '10:00:31.88 01:57:50.04'

    # c1 = SkyCoord(uvista_1212, unit=(u.hourangle, u.deg))
    # c2 = SkyCoord(uvista_237, unit=(u.hourangle, u.deg))

    # ra = [c1.ra.deg, c2.ra.deg]
    # dec = [c1.dec.deg, c2.dec.deg]

    # ID = ['uvista_1212', 'uvista_237']

    #! Big three dragons
    # b1 = '10:01:40.69 01:54:52.42'
    # b1 = SkyCoord(b1, unit=(u.hourangle, u.deg))
    # ra = [b1.ra.deg]
    # dec = [b1.dec.deg]

    #! Bad astrometry sources
    # stars_dir = Path.cwd().parent.parent / 'data' / 'ref_catalogues' / 'stars'
    # stars = ascii.read(stars_dir / f'H_outside_pixscale_vista_euclid_coords.ascii')
    # ra = stars['RA_euclid']
    # dec = stars['DEC_euclid']

    #! Normal stars
    # stars_dir = Path.cwd().parent.parent / 'data' / 'ref_catalogues' / 'stars'
    # stars = ascii.read(stars_dir / f'Y_vista_euclid_coords.ascii')
    # ra = stars['RA_euclid']
    # dec = stars['DEC_euclid']

    #! Couple bad JWST astrometry sources
    # stars_dir = Path.cwd().parent.parent / 'data' / 'ref_catalogues' / 'stars'
    # stars = ascii.read(stars_dir / f'H_outside_pixscale_jwst_euclid_coords.ascii')
    # ra = stars['RA_euclid']
    # dec = stars['DEC_euclid'] 

    #! Spectra from Vincent!
    #ra = [150.21552279971348, 149.98832200622317, 150.26771198215147, 150.24988047861532]
    #dec = [2.02579740958649, 2.2070238197331493, 2.6213922569606107, 2.4663986583043647]

    #ID = [1502155319020257969, 1499883241022070278, 1502677115026213801, 1502498770024664099]

    #! Casey CWEB z>10 sources

    # Super bright 10<z<12 sources
    # ra = ['10:01:26.00', '09:58:55.21', '09:59:59.91', '09:59:49.04']
    # dec = ['01:55:59.70', '02:07:16.77', '02:06:59.90', '01:53:26.19']
    # z = [10.27, 12.54, 11.92, 12.03]
    # Muv = [-21.53, -22.19, -21.89, -21.58]
    # ID = ['COS-z10-1', 'COS-z12-1', 'COS-z12-2', 'COS-z12-3']

    # Bright 10<z<12 sources
    # ra = ['09:59:51.77', '09:59:57.50', '10:00:37.96', '09:59:52.53', '10:01:34.80']
    # dec = ['2:07:15.02', '02:06:20.06', '01:49:32.43', '02:00:23.53', '02:05:41.48']
    # z = [10.06, 10.17, 11.07, 11.70, 11.50]
    # Muv = [-20.62, -20.97, -20.85, -21.13, -20.77]
    # ID = ['COS-z10-2', 'COS-z10-3', 'COS-z11-1', 'COS-z11-2', 'COS-z11-3']

    # # z>13
    # ra = ['09:59:05.75', '10:00:04.24', '10:01:31.17']
    # dec = ['02:04:04.39', '02:02:11.19', '01:58:45.00']
    # z = [13.10, 12.50, 14.7]
    # Muv = [-21.27, -21.03, -20.75]
    # ID = ['COS-z13-1', 'COS-z13-2', 'COS-z14-1']

    # Probable contaminants
    # ra = ['09:59:30.49', '09:59:31.30', '10:00:20.38']
    # dec = ['02:14:44.10', '02:08:33.85', '01:49:58.33']
    # z = [12.63, 13.8, 14.7]
    # Muv = [-21.90, -20.97, -21.32]
    # ID = ['COS-z12-4', 'COS-z13-3', 'COS-z14-2']

    # Convert ra, dec to degrees
    # ra = [SkyCoord(r, d, unit=(u.hourangle, u.deg)).ra.deg for r, d in zip(ra, dec)]
    # dec = [SkyCoord(r, d, unit=(u.hourangle, u.deg)).dec.deg for r, d in zip(ra, dec)]

    #! CR7
    # cr7 = '10:00:58.005 01:48:15.251'
    # cr7 = SkyCoord(cr7, unit=(u.hourangle, u.deg))

    # ra = [cr7.ra.deg]
    # dec = [cr7.dec.deg]
    # ID = ['CR7']

    #! All COSMOS crossmatched sources
    # t = Table.read(refcat_dir / 'all_COSMOS_highz.fits')

    # # Sort by decreasing redshift
    # t.sort('Redshift')

    # # Flip
    # t = t[::-1]


    # ra = t['RA']
    # dec = t['DEC']

    # z = t['Redshift']
    # ID = t['Object Name']

    #! Initial det_YJH candidates
    # t = Table.read(Path.cwd().parents[1] / 'data' / 'catalogues' / 'XMATCH_COSMOS_5sig_Ye_2sig_VISTA_Y_nonDet_HSC_G_nonDet_HSC_R_nonDet_HSC_I.fits')
    # ra = t['RA_1']
    # dec = t['DEC_1']
    # ID = t['Object Name']
    # z = t['Redshift']

    # # Sort by flux in Y, reversed
    # t.sort('flux_Ye', reverse=True)
    # print(t['flux_Ye'])

    #! My LAE candidate
    ra = [150.11833152095758]
    dec = [2.2522416552619746]
    ID = [178396]

    #! Stars
    # t = Table.read(Path.cwd().parents[1] / 'data' / 'depths' / 'COSMOS' / 'catalogues' / 'df444w_locus_stars.fits')

    # # Sort by class_star
    # t.sort('CLASS_STAR', reverse=True)
    # print(len(t))

    # t = t[t['CLASS_STAR'] < 0.99]

    # STARS = (t['FLAGS'] < 2) & (t['ELONGATION'] < 1.5) & (t['FWHM_IMAGE'] < 6) & (t['FWHM_IMAGE'] > 5) & (t['CLASS_STAR'] < 0.99) & (t['CLASS_STAR'] > 0.8)
    # t = t[STARS]

    # print(len(t))

    # #t = t[(t['FWHM_IMAGE'] > 5.5) & (t['FWHM_IMAGE'] < 7)]

    # ra = t['RA']
    # dec = t['DEC']
    # class_star = t['CLASS_STAR']
    # flag = t['FLAGS']
    # elong = t['ELONGATION']
    # fwhm = t['FWHM_IMAGE']


    ############! GET CUTOUTS ############
    for i in range(len(ra)):

        print(f'Object number {i+1} of {len(ra)}')

    #     print(f'Object number {i+1} of {len(ra)}')

    #     if isCoordInCWEB(ra[i], dec[i])[0] == '0':
    #         print('Not in CWEB')
    #         continue

        #print(cat[i]['zBest'])
        #print(cat[i]['FIRST_CLASS'])
        #print(ID[i])
        #print(muv[i])

        #Cutout(ra[i], dec[i], size=10., plot_title=str(ID[i]) + ', z=' + str(z[i]), save_cutout=False)
        #print(class_star[i], flag[i], elong[i], fwhm[i])
        #Cutout(ra[i], dec[i], size=6., save_cutout=False)
        #Cutout(ra[i], dec[i], size=6., add_centre_lines=True)
        Cutout(ra[i], dec[i], size=6., plot_title=ID[i])
        #Cutout(ra[i], dec[i], size=10., plot_title='Big Three Dragons')   #
        #Cutout(ra[i], dec[i], size=4., plot_title=ID[i] + ', z=' + str(z[i]) + ', Muv=' + str(Muv[i]))


    
