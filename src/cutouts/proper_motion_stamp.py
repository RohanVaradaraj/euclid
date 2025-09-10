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

plt.rcParams['axes.linewidth'] = 2.5
plt.rcParams.update({'font.size': 15})
plt.rcParams['figure.dpi'] = 100

ground_dir = Path.home().parent.parent / 'vardy' / 'vardygroupshare' / 'data'
euclid_dir = Path.home() / 'euclid'
cweb_dir = Path.home().parent.parent / 'vardy' / 'vardygroupshare' / 'data' / 'CWEB'
primer_dir = Path.home().parent.parent / 'vardy' / 'vardygroupshare' / 'data' / 'COSMOS'
plot_dir = Path.cwd().parent.parent / 'plots' / 'cutouts'
refcat_dir = Path.cwd().parents[1] / 'data' / 'ref_catalogues'

def add_crosshairs(x, y, ax, data, length_frac=0.05, gap_frac=0.02,
                   color="red", lw=1.2):
    """
    Add short crosshairs with a central gap to the centre of the image.
    
    Parameters
    ----------
    ax : matplotlib axis
    data : 2D numpy array
    length_frac : float
        Fraction of the image size for each arm of the crosshair.
    gap_frac : float
        Fraction of the image size for the central gap.
    color : str
        Crosshair color.
    lw : float
        Line width.
    """
    ny, nx = data.shape
    #cx, cy = nx // 2, ny // 2  # image centre in pixels
    cx, cy, = x, y
    
    Lx = int(nx * length_frac)
    Ly = int(ny * length_frac)
    Gx = int(nx * gap_frac)
    Gy = int(ny * gap_frac)

    # Horizontal left
    ax.plot([cx - Lx, cx - Gx], [cy, cy], color=color, lw=lw)
    # Horizontal right
    ax.plot([cx + Gx, cx + Lx], [cy, cy], color=color, lw=lw)
    # Vertical bottom
    ax.plot([cx, cx], [cy - Ly, cy - Gy], color=color, lw=lw)
    # Vertical top
    ax.plot([cx, cx], [cy + Gy, cy + Ly], color=color, lw=lw)



def findPlotLimits(data: np.ndarray) -> tuple:

    mean = np.mean(data)
    std_dev = np.std(data)

    # Sigma clip the data
    data = data[(data < mean + 1 * std_dev)]

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
    
    vista_dir = ground_dir / 'COSMOS'

    #### J ####
    with fits.open(vista_dir / 'UVISTA_J_DR6.fits') as hdu_J:
    #with fits.open(vista_dir / 'UVISTA_J_dr5_rc1.fits') as hdu_J:

        data_J = hdu_J[0].data
        hdr_J = hdu_J[0].header
        pix_scale = np.abs(hdr_J['CD1_1']) * 3600
        wcs_J = WCS(hdr_J)

        cutout_grJ = Cutout2D(data_J, c, size=size/pix_scale, wcs=wcs_J)

        # Also get x,y of the source
        center_x, center_y = cutout_grJ.wcs.world_to_pixel(SkyCoord(c.ra.deg, c.dec.deg, unit='deg'))

    #! Next get the Euclid cutouts
    euclid_tile = contained_in[0][0]

    if euclid_tile != '0':


        #### J ####
        tile_J = glob.glob(str(euclid_dir / 'J' / 'COSMOS' / f'*BGSUB*_{euclid_tile}.fits*'))[0]
        with fits.open(tile_J) as hdu_J:

            data_J = hdu_J[0].data
            hdr_J = hdu_J[0].header
            pix_scale = np.abs(hdr_J['CD1_1']) * 3600
            wcs_J = WCS(hdr_J)

            cutout_euJ = Cutout2D(data_J, c, size=size/pix_scale, wcs=wcs_J)




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
            pa = hdr_CWEB['PA_APER']

            cutout_f150w_cweb = Cutout2D(data_CWEB, c, size=size/pix_scale, wcs=wcs_CWEB)

            #cutout_f150w_cweb.data = rotate(cutout_f150w_cweb.data, angle=-pa, reshape=False)

  
    #! Plots based on different footprint bools

    plot_cutout = lambda ax, data, lims, title: ax.imshow(data, origin='lower', cmap='gist_yarg', vmin=lims[0], vmax=lims[1]) and ax.set_title(title, fontsize=30) 

    # fIG WITH 4 subplots for F160W, VISTA-J, Euclid J, and F150W, in a row
    fig, ax = plt.subplots(1, 3, figsize=(15, 5)) #, sharex=True, sharey=True)

    plt.subplots_adjust(wspace=0.1, hspace=0.1)

    for axis in ax:
            #axis.set_axis_off()
            axis.tick_params(axis='both', which='both', bottom=False, top=False, left=False, right=False, labelbottom=False, labelleft=False)

    # titles = ['F160W', r'$J$', r'$J_{\rm E}$', 'F150W']
    # cutouts = [cutout_hubble, cutout_grJ, cutout_euJ, cutout_f150w_cweb]
    titles = [r'$J$', r'$J_{\rm E}$', 'F115W']
    cutouts = [cutout_grJ, cutout_euJ, cutout_f115w_cweb]
    for i, (cutout, title) in enumerate(zip(cutouts, titles)):

        data = cutout.data
        plot_lims = findPlotLimits(data)

        if i == 2:
            # For F150W, we need to rotate the cutout
            data = rotate(data, -pa, reshape=False)

        plot_cutout(ax[i], data, plot_lims, title)

        if add_centre_lines:
            center_x, center_y = cutout.wcs.world_to_pixel(SkyCoord(c.ra.deg, c.dec.deg, unit='deg'))

            ax[i].plot(center_x, center_y, '+', markersize=20, markeredgewidth=3, color='red', alpha=1)

            # ax[i].plot(center_x, center_y, marker=1, markersize=15, markeredgewidth=2.5, color='tab:red')  
            # ax[i].plot(center_x, center_y, marker=2, markersize=15, markeredgewidth=2.5, color='tab:red')  

        # ADD THE frame around the cutout

    # Save fig as pdf
    if save_cutout:
        plot_dir = Path.cwd().parents[1] / 'plots' / 'cutouts'
        plt.savefig(plot_dir / f'proper_motion.pdf', bbox_inches='tight')
    #plt.show()
    return fig, ax

    #plt.close()
    #return None



if __name__ == '__main__':
    
    #! Brown dwarfs - any proper motion?
    t = Table.read('/mnt/vardy/vardygroupshare/rohan/euclid/data/catalogues/candidates/COSMOS_5sig_Y_J_nonDet_HSC_G_nonDet_HSC_R_nonDet_HSC_I_really_good_BDs_INTERLOPERS_2025_08_22_with_euclid.fits')
    t = t[t['ID'] == 829044]
    ra = t['RA']
    dec = t['DEC']
    ID = t['ID']


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
        Cutout(ra[i], dec[i], size=5., add_centre_lines=True, save_cutout=True)
        #Cutout(ra[i], dec[i], size=6., plot_title=ID[i])
        #Cutout(ra[i], dec[i], size=10., plot_title='Big Three Dragons')   #
        #Cutout(ra[i], dec[i], size=4., plot_title=ID[i] + ', z=' + str(z[i]) + ', Muv=' + str(Muv[i]))


    
