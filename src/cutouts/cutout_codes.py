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
from astropy.visualization import simple_norm
from reproject import reproject_interp
from scipy.ndimage import rotate
from astropy.table import Table
import astropy.units as u

import warnings
from astropy.utils.exceptions import AstropyWarning
warnings.simplefilter('ignore', category=AstropyWarning)

ground_dir = Path.home().parent.parent / 'vardy' / 'vardygroupshare' / 'data'
euclid_dir = Path.home() / 'euclid'
cweb_dir = Path.home().parent.parent / 'extraspace' / 'varadaraj' / 'CWEB'
primer_dir = Path.home() / 'JWST'
plot_dir = Path.cwd().parent.parent / 'plots' / 'cutouts'


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



def Cutout(ra: float, dec:float, contained_in: Optional[np.array] = None, size: float = 10.0, 
           save_cutout: bool = False, save_dir: Path = Path.cwd().parent.parent / 'data' / 'cutouts',
           plot_title: Optional[str] = None) -> None:

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

    
    Returns
    -------
    None

    """

    # If contained_in is None, do the check for coordinates in the footprints
    if contained_in is None:
        contained_in = isCoordInSurveyFootprints(ra, dec)
        print(contained_in)

    # Check if all values of contained_in are '0'.
    if np.all(contained_in[0] == '0'):
        print("Cutout not in any survey footprint")
        return None
    
    if contained_in[0][0] == '0':
        print('Cutout not in Euclid footprint')
        return None

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
    with fits.open(vista_dir / 'UVISTA_Y_DR6_cropped.fits') as hdu_Y:

        data_Y = hdu_Y[0].data
        hdr_Y = hdu_Y[0].header
        pix_scale = np.abs(hdr_Y['CD1_1']) * 3600
        wcs_Y = WCS(hdr_Y)

        cutout_grY = Cutout2D(data_Y, c, size=size/pix_scale, wcs=wcs_Y)

    #### J ####
    with fits.open(vista_dir / 'UVISTA_J_DR6_cropped.fits') as hdu_J:

        data_J = hdu_J[0].data
        hdr_J = hdu_J[0].header
        pix_scale = np.abs(hdr_J['CD1_1']) * 3600
        wcs_J = WCS(hdr_J)

        cutout_grJ = Cutout2D(data_J, c, size=size/pix_scale, wcs=wcs_J)

    #### H ####
    with fits.open(vista_dir / 'UVISTA_H_DR6_cropped.fits') as hdu_H:
            
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
    cweb_tile = contained_in[0][1]

    if cweb_tile != '0':

        tile_f277w = glob.glob(str(cweb_dir / f'*F277W*{cweb_tile}*'))[0]
        tile_f444w = glob.glob(str(cweb_dir / f'*F444W*{cweb_tile}*'))[0]

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

    #* Define a lambda function to plot each cutout
    plot_cutout = lambda ax, data, lims, title: ax.imshow(data, origin='lower', cmap='gist_yarg', vmin=lims[0], vmax=lims[1]) and ax.set_title(title)

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
        for i, (data, title) in enumerate(zip([cutout_euVIS.data, cutout_euY.data, cutout_euJ.data, cutout_euH.data], ['VIS', 'NISP-Y', 'NISP-J', 'NISP-H'])):
            lims = findPlotLimits(data)
            plot_cutout(ax[1, i], data, lims, title)

    # In Euclid and CWEB
    elif footprint_bools[0] and footprint_bools[1] and not footprint_bools[2] and not footprint_bools[3]: 

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
        for i, (data, title) in enumerate(zip([cutout_euVIS.data, cutout_euY.data, cutout_euJ.data, cutout_euH.data], ['VIS', 'NISP-Y', 'NISP-J', 'NISP-H'])):
            lims = findPlotLimits(data)
            plot_cutout(ax[1, i], data, lims, title)


        # CWEB cutouts in final column. Reproject to Euclid WCS
        for i, (data, title) in enumerate(zip([cutout_f277w_cweb.data, cutout_f444w_cweb.data], ['F277W', 'F444W'])):
            data = rotate(data, -pa, reshape=False)
            lims = findPlotLimits(data)
            plot_cutout(ax[0, 4], data, lims, title) if i == 0 else plot_cutout(ax[1, 4], data, lims, title)


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
        for i, (data, title) in enumerate(zip([cutout_euVIS.data, cutout_euY.data, cutout_euJ.data, cutout_euH.data], ['VIS', 'NISP-Y', 'NISP-J', 'NISP-H'])):
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
        for i, (data, title) in enumerate(zip([cutout_euVIS.data, cutout_euY.data, cutout_euJ.data, cutout_euH.data], ['VIS', 'NISP-Y', 'NISP-J', 'NISP-H'])):
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
        for i, (data, title) in enumerate(zip([cutout_euVIS.data, cutout_euY.data, cutout_euJ.data, cutout_euH.data], ['VIS', 'NISP-Y', 'NISP-J', 'NISP-H'])):
            lims = findPlotLimits(data)
            plot_cutout(ax[1, i], data, lims, title)

        # Plot Hubble data in penultimate column
        lims = findPlotLimits(cutout_hubble.data)
        plot_cutout(ax[1, 4], cutout_hubble.data, lims, 'F160W')

        # CWEB cutouts in final column. Reproject to Euclid WCS
        if footprint_bools[1]:
            for i, (data, title) in enumerate(zip([cutout_f277w_cweb.data, cutout_f444w_cweb.data], ['F277W', 'F444W'])):
                data = rotate(data, -pa, reshape=False)
                lims = findPlotLimits(data)
                plot_cutout(ax[0, 5], data, lims, title) if i == 0 else plot_cutout(ax[1, 5], data, lims, title)
        if footprint_bools[2]:
            for i, (data, title) in enumerate(zip([cutout_f277w_prim.data, cutout_f444w_prim.data], ['F277W', 'F444W'])):
                lims = findPlotLimits(data)
                plot_cutout(ax[0, 5], data, lims, title) if i == 0 else plot_cutout(ax[1, 5], data, lims, title)


    #plt.savefig(plot_dir / f'cutout_test.png')
            
    if plot_title is not None:
        plt.suptitle(plot_title)
    plt.tight_layout()
    plt.show()

    return None


if __name__ == '__main__':
    #! REBELS sources
    # t = ascii.read(Path.cwd().parent.parent / 'data' / 'mosaic' / 'REBELS.csv', format='csv')
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
    # cat.sort('MUV', reverse=True)
    # cat = cat[cat['RA'] > 148]
    # ra = cat['RA']
    # dec = cat['DEC']

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
    hd1 = '10:01:51.31 02:32:50.0'
    hd2 = '02:18:52.44 -05:08:36.1'

    # Use skycoord to convert these coordinates to degrees
    c1 = SkyCoord(hd1, unit=(u.hourangle, u.deg))
    c2 = SkyCoord(hd2, unit=(u.hourangle, u.deg))

    ra = [c1.ra.deg, c2.ra.deg]
    dec = [c1.dec.deg, c2.dec.deg]


    for i in range(len(ra)):

        #print(cat[i]['zBest'])
        #print(cat[i]['FIRST_CLASS'])

        #Cutout(ra[i], dec[i], size=6., plot_title=ID[i] + ', z=' + str(z[i]))
        Cutout(ra[i], dec[i], size=10.)


    
