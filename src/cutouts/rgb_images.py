"""
Make RGB images of candidate galaxies which have Euclid imaging.

Created: Monday 3rd February
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
#from astropy.visualization import make_lupton_rgb, LogStretch, MinMaxInterval, ManualInterval# make_rgb
from astropy.stats import sigma_clip
import matplotlib.gridspec as gridspec
from astropy.visualization.wcsaxes import WCSAxes

plt.rcParams['axes.linewidth'] = 2.5
plt.rcParams.update({'font.size': 15})
plt.rcParams['figure.dpi'] = 100

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


# Sigma-clipping function
def sigma_clip_band(band, sigma=3):
    clipped = sigma_clip(band, sigma=sigma, maxiters=5, masked=False)  # Masked=False keeps it as an array
    vmin, vmax = np.nanmin(clipped), np.nanmax(clipped)  # Use min/max of clipped data
    return np.clip(band, vmin, vmax)  # Clip values to the new range


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
    upper = mean + 3.5 * std_dev

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
           save_dir: Path = Path.cwd().parent.parent / 'data' / 'cutouts',
           plot_title: Optional[str] = None,
           add_centre_lines: Optional[bool] = False,
           save_cutout=False,
           ID=None) -> None:

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
        return None

    # Convert RA, DEC to skycoords.
    c = SkyCoord(ra, dec, unit='deg')
    
    #! Get the Euclid cutouts
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

    #! Get VISTA cutout
    vista_dir = ground_dir / 'COSMOS'

    with fits.open(vista_dir / 'UVISTA_YJHK_DR6.fits') as hdu_Y:

        data_Y = hdu_Y[0].data
        hdr_Y = hdu_Y[0].header
        pix_scale = np.abs(hdr_Y['CD1_1']) * 3600
        wcs_Y = WCS(hdr_Y)

        cutout_gr = Cutout2D(data_Y, c, size=size/pix_scale, wcs=wcs_Y)

    cutouts = [cutout_euVIS, cutout_euY, cutout_euJ, cutout_euH, cutout_gr]

    # Save the cutouts as fits file
    cutout_dir = Path.cwd().parents[1] / 'data' / 'stamps'
    if save_cutout:

        IDs_to_ignore = [] # [11029, 26088, 30883, 210827, 216577, 563004, 603546, 765174, 984767]
        if cutouts != None and ID not in IDs_to_ignore:

            # Check if any cutout has nonzero pixels
            if any(np.sum(cutout.data) > 0.01 for cutout in cutouts[:-1]):

                hdu_Y = fits.PrimaryHDU(cutout_euY.data, header=cutout_euY.wcs.to_header())
                hdu_Y.writeto(cutout_dir / f'Y_{ID}.fits', overwrite=True)    

                hdu_J = fits.PrimaryHDU(cutout_euJ.data, header=cutout_euJ.wcs.to_header())
                hdu_J.writeto(cutout_dir / f'J_{ID}.fits', overwrite=True)

                hdu_H = fits.PrimaryHDU(cutout_euH.data, header=cutout_euH.wcs.to_header())
                hdu_H.writeto(cutout_dir / f'H_{ID}.fits', overwrite=True)



    return cutouts



    # if [np.sum(cutout.data) == 0. for cutout in cutouts]:
    #     print('No data')
    #     return None,None

    # Plot the 4 cutouts as subplots in a single figure
    # fig, axs = plt.subplots(1, 5, figsize=(20, 5))  # Create 1 row, 4 columns of subplots

    # for i, cutout in enumerate(cutouts):
    #     ax = axs[i]
    #     vmin, vmax = findPlotLimits(cutout.data)
    #     ax.imshow(cutout.data, origin='lower', cmap='gist_yarg', vmin=vmin, vmax=vmax)
    #     #ax.set_title(f'{filter_name} Band')
    #     ax.axis('off')  # Turn off axis labels

    # plt.subplots_adjust(wspace=0.05)
    # plt.tight_layout()
    # # Add a main title to the entire subplot
    # # if plot_title:
    # #     fig.suptitle(plot_title, fontsize=16)

    # # Save the figure if required
    # # if save_cutout:
    # #     save_path = save_dir / f"cutout_{ra}_{dec}.png"
    # #     plt.savefig(save_path, bbox_inches='tight')
    # #     print(f"Saved cutout to {save_path}")

    # # Show the plot
    # #plt.show()
    # return fig, axs



if __name__ == '__main__':
    
    #! U+E CANDIDATES
    t = Table.read(Path.cwd().parents[1] / 'data' / 'catalogues' / 'candidates' / 'COSMOS_5sig_Y_J_nonDet_HSC_G_nonDet_HSC_R_nonDet_HSC_I_candidates_2025_02_14_with_euclid.fits')

    #t_xmatch = Table.read(Path.cwd().parents[1] / 'data' / 'catalogues' / 'candidates' / 'COSMOS_5sig_Y_J_nonDet_HSC_G_nonDet_HSC_R_nonDet_HSC_I_candidates_2025_01_31_XMATCH_WITH_LITERATURE.fits')
    # print(t_xmatch.colnames)
    # exit()
    t.sort('Muv')
    #t = t[t['ID'] == 910976]
    #t = t[50:60]
    t = t[0:30]

    #! BROWN DWARFS
    # t = Table.read(Path.cwd().parents[1] / 'data' / 'catalogues' / 'COSMOS_5sig_Y_J_nonDet_HSC_G_nonDet_HSC_R_nonDet_HSC_I.fits')

    # IDs_wanted = [887134, 329431]
    # t = t[np.isin(t['ID'], IDs_wanted)]
    # print(t)

    ra = t['RA']
    dec = t['DEC']
    ID = t['ID']
    # zphot = t['Zphot']
    #Muv = t['Muv']

    #! Stavrox z=6.8
    # ra = [150.1200752011961]
    # dec = [2.0898737362171156]


    # Empty arrays for collecting figures
    all_cutouts = []
    IDs = []
    zs = []
    Muvs = []

    IDs_to_ignore = [] #[11029, 26088, 30883, 210827, 216577, 563004, 603546, 765174, 984767]

    ############! GET CUTOUTS ############
    for i in range(len(ra)):

        print(f'Object number {i+1} of {len(ra)}')

        cutouts = Cutout(ra[i], dec[i], size=6., plot_title=ID[i], save_cutout=True, ID=ID[i])

        if cutouts != None and ID[i] not in IDs_to_ignore:

            # Check if any cutout has nonzero pixels
            if any(np.sum(cutout.data) > 0.01 for cutout in cutouts[:-1]):
                all_cutouts.append(cutouts)
                IDs.append(ID[i])  # Keep track of valid IDs
                #zs.append(zphot[i])
                #Muvs.append(Muv[i])
            else:
                print(f"Skipping object {ID[i]} (empty cutout)")

    ############! PLOT CUTOUTS #############

    # Cutout titles
    titles = [r'$I_{\rm{E}}$', r'$Y_{\rm{E}}$', r'$J_{\rm{E}}$', r'$H_{\rm{E}}$', r'$YJHK_s$']

    # Save each object (with 5 cutouts) as a single image
    for i, cutouts in enumerate(all_cutouts):
        # Create a figure for the object
        fig, axs = plt.subplots(1, 5, figsize=(15, 3))  # 1 row, 5 columns for each object

        # Loop through each cutout and plot it
        for j, cutout in enumerate(cutouts):
            vmin, vmax = findPlotLimits(cutout.data)
            axs[j].imshow(cutout.data, cmap='gist_yarg', origin='lower', vmin=vmin, vmax=vmax)
            #axs[j].coords.grid(color='white', ls='dotted')  # Add grid to check alignment
            #axs[j].axis('off')  # Remove axis for a cleaner look

            # Remove axis tick marks and labels
            axs[j].get_xaxis().set_ticks([])
            axs[j].get_yaxis().set_ticks([])

            axs[j].set_title(titles[j], fontsize=35)

            # See if this object is in the xmatch cat
            # if IDs[i] in t_xmatch['ID']:
            #     xmatch_id = t_xmatch[t_xmatch['ID'] == IDs[i]]['Object Name'][0]
            #     # if len(xmatch_id) > 0:
            #     #     xmatch_id = xmatch_id[0]
            #     if j == 2:
            #         axs[j].set_title(f'ID {IDs[i]}, ' + r'$z=$'+f'{zs[i]:.2f}, ' + r'$M_{\rm{UV}}=$'+f'{Muvs[i]:.2f}, ' + f'{xmatch_id}', fontsize=30)
            # # Add title only to the first subplot (j == 0) for each object
            # else:
            #     if j == 2:
            #         axs[j].set_title(f'ID {IDs[i]}, ' + r'$z=$'+f'{zs[i]:.2f}, ' + r'$M_{\rm{UV}}=$'+f'{Muvs[i]:.2f}', fontsize=30)

        # Save the entire figure (with 5 cutouts) as a single image
        plt.subplots_adjust(wspace=0.1, hspace=0.1)
        cutout_filename = f"stamps/{i}_ID_{IDs[i]}.pdf"
        #plt.tight_layout()
        plt.savefig(plot_dir / 'bright_candidates' / cutout_filename, bbox_inches='tight')
        #plt.savefig(plot_dir / 'UCDs' / cutout_filename, bbox_inches='tight')
        plt.close(fig)  # Close the figure to free up memory

        print(f"Saved {cutout_filename}")

    





    
